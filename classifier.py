#!/usr/bin/env python3
"""

Takes a preprocessed .h5ad (log-normalized, with cell type labels in obs) and
trains a simple feedforward net on the HVG expression matrix. Outputs metrics,
a confusion matrix, and gradient-based gene importance rankings.

Usage:
    python scrna_celltype_classifier.py --input pbmc_processed.h5ad --label_key cell_type
"""

import argparse
import warnings
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import torch
import torch.nn as nn
from torch.utils.data import DataLoader, TensorDataset
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import (
    accuracy_score,
    classification_report,
    confusion_matrix,
    f1_score,
)

warnings.filterwarnings("ignore")

SEED = 13
DEVICE = "cuda" if torch.cuda.is_available() else "cpu"


def load_and_prepare(h5ad_path, label_key="cell_type", n_hvg=2000):
    """Load h5ad, subset to HVGs, encode labels, return train/val/test splits."""
    print(f"Loading {h5ad_path} ...")
    adata = sc.read_h5ad(h5ad_path)

    if label_key not in adata.obs.columns:
        available = ", ".join(adata.obs.columns.tolist()[:20])
        raise KeyError(f"'{label_key}' not in adata.obs. Available: {available}")

    # drop unlabeled cells
    mask = adata.obs[label_key].notna()
    if (~mask).sum() > 0:
        print(f"  Dropping {(~mask).sum()} unlabeled cells")
        adata = adata[mask].copy()

    if "highly_variable" not in adata.var.columns:
        print(f"  Computing {n_hvg} HVGs ...")
        sc.pp.highly_variable_genes(adata, n_top_genes=n_hvg, flavor="seurat_v3")
    adata_hvg = adata[:, adata.var["highly_variable"]].copy()

    X = adata_hvg.X
    if hasattr(X, "toarray"):
        X = X.toarray()
    X = X.astype(np.float32)
    gene_names = adata_hvg.var_names.tolist()
    print(f"  {X.shape[0]} cells x {X.shape[1]} genes")

    le = LabelEncoder()
    y = le.fit_transform(adata.obs[label_key].values)
    n_classes = len(le.classes_)
    print(f"  {n_classes} classes: {le.classes_.tolist()}")

    # stratified splits: 70/10/20
    X_tmp, X_test, y_tmp, y_test = train_test_split(
        X, y, test_size=0.2, stratify=y, random_state=SEED
    )
    X_train, X_val, y_train, y_val = train_test_split(
        X_tmp, y_tmp, test_size=0.125, stratify=y_tmp, random_state=SEED
    )
    print(f"  Split sizes: train={len(y_train)} val={len(y_val)} test={len(y_test)}")

    tensors = {
        "X_train": torch.tensor(X_train),
        "X_val": torch.tensor(X_val),
        "X_test": torch.tensor(X_test),
        "y_train": torch.tensor(y_train, dtype=torch.long),
        "y_val": torch.tensor(y_val, dtype=torch.long),
        "y_test": torch.tensor(y_test, dtype=torch.long),
    }
    return tensors, le, gene_names, n_classes


class CellTypeClassifier(nn.Module):
    """Feedforward net: input -> 256 -> 128 -> 64 -> n_classes."""

    def __init__(self, n_genes, n_classes, dropout=0.3):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(n_genes, 256),
            nn.BatchNorm1d(256),
            nn.ReLU(),
            nn.Dropout(dropout),

            nn.Linear(256, 128),
            nn.BatchNorm1d(128),
            nn.ReLU(),
            nn.Dropout(dropout),

            nn.Linear(128, 64),
            nn.BatchNorm1d(64),
            nn.ReLU(),
            nn.Dropout(dropout),

            nn.Linear(64, n_classes),
        )

    def forward(self, x):
        return self.net(x)


def train_model(model, tensors, n_epochs=50, batch_size=256, lr=1e-3, patience=10):
    """Train with early stopping on val loss. Returns loss histories."""
    print(f"\nTraining on {DEVICE}, up to {n_epochs} epochs ...")
    model = model.to(DEVICE)

    train_dl = DataLoader(
        TensorDataset(tensors["X_train"], tensors["y_train"]),
        batch_size=batch_size, shuffle=True,
    )
    val_dl = DataLoader(
        TensorDataset(tensors["X_val"], tensors["y_val"]),
        batch_size=batch_size,
    )

    criterion = nn.CrossEntropyLoss()
    opt = torch.optim.Adam(model.parameters(), lr=lr)
    sched = torch.optim.lr_scheduler.ReduceLROnPlateau(opt, factor=0.5, patience=5)

    best_val_loss = float("inf")
    best_state = None
    stale = 0
    train_hist, val_hist = [], []

    for ep in range(1, n_epochs + 1):
        model.train()
        running = 0.0
        for xb, yb in train_dl:
            xb, yb = xb.to(DEVICE), yb.to(DEVICE)
            opt.zero_grad()
            loss = criterion(model(xb), yb)
            loss.backward()
            opt.step()
            running += loss.item() * len(xb)
        t_loss = running / len(tensors["X_train"])

        model.eval()
        running = 0.0
        with torch.no_grad():
            for xb, yb in val_dl:
                xb, yb = xb.to(DEVICE), yb.to(DEVICE)
                running += criterion(model(xb), yb).item() * len(xb)
        v_loss = running / len(tensors["X_val"])

        train_hist.append(t_loss)
        val_hist.append(v_loss)
        sched.step(v_loss)

        if ep % 10 == 0 or ep == 1:
            print(f"  [{ep:>3d}/{n_epochs}] train={t_loss:.4f} val={v_loss:.4f}")

        if v_loss < best_val_loss:
            best_val_loss = v_loss
            best_state = {k: v.cpu().clone() for k, v in model.state_dict().items()}
            stale = 0
        else:
            stale += 1
            if stale >= patience:
                print(f"  Early stop at epoch {ep}")
                break

    model.load_state_dict(best_state)
    model = model.to(DEVICE)
    return model, train_hist, val_hist


def evaluate(model, tensors, label_encoder):
    """Evaluate on test set. Returns predictions and confusion matrix."""
    model.eval()
    X_test = tensors["X_test"].to(DEVICE)
    y_true = tensors["y_test"].numpy()

    with torch.no_grad():
        y_pred = model(X_test).argmax(dim=1).cpu().numpy()

    acc = accuracy_score(y_true, y_pred)
    f1 = f1_score(y_true, y_pred, average="weighted")
    names = label_encoder.classes_

    print(f"\nTest accuracy: {acc:.4f}  |  Weighted F1: {f1:.4f}")
    print(classification_report(y_true, y_pred, target_names=names))

    cm = confusion_matrix(y_true, y_pred)
    return y_pred, cm, acc, f1


def get_gene_importance(model, tensors, gene_names, top_n=20):
    """
    Gradient-based feature attribution: for each test cell, compute
    |d(predicted_logit)/d(input)| and average across cells.
    """
    model.eval()
    X = tensors["X_test"].clone().to(DEVICE).requires_grad_(True)

    logits = model(X)
    pred_logits = logits.gather(1, logits.argmax(dim=1, keepdim=True)).squeeze()
    pred_logits.sum().backward()

    imp = X.grad.abs().mean(dim=0).cpu().numpy()
    ranked = np.argsort(imp)[::-1][:top_n]

    gene_arr = np.array(gene_names)
    top_genes = gene_arr[ranked].tolist()
    top_scores = imp[ranked].tolist()

    print(f"\nTop {top_n} genes by gradient importance:")
    for g, s in zip(top_genes, top_scores):
        print(f"  {g:<18s} {s:.6f}")

    return imp, top_genes, top_scores


def save_plots(train_hist, val_hist, cm, class_names, top_genes, top_scores, outdir):
    """Loss curves, confusion matrix, gene importance barplot."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    axes[0].plot(train_hist, label="Train")
    axes[0].plot(val_hist, label="Val")
    axes[0].set(xlabel="Epoch", ylabel="Loss", title="Loss curves")
    axes[0].legend()
    axes[0].grid(alpha=0.3)

    cm_norm = cm.astype(float) / cm.sum(axis=1, keepdims=True)
    sns.heatmap(cm_norm, annot=True, fmt=".2f", cmap="Blues",
                xticklabels=class_names, yticklabels=class_names, ax=axes[1])
    axes[1].set(xlabel="Predicted", ylabel="True", title="Confusion matrix (normalized)")

    axes[2].barh(range(len(top_genes)), top_scores[::-1])
    axes[2].set_yticks(range(len(top_genes)))
    axes[2].set_yticklabels(top_genes[::-1])
    axes[2].set(xlabel="Mean |gradient|", title="Top genes by importance")
    axes[2].grid(axis="x", alpha=0.3)

    plt.tight_layout()
    outpath = outdir / "classifier_results.png"
    fig.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"\nFigure saved: {outpath}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, help="Preprocessed .h5ad")
    ap.add_argument("--label_key", default="cell_type")
    ap.add_argument("--n_hvg", type=int, default=2000)
    ap.add_argument("--epochs", type=int, default=50)
    ap.add_argument("--batch_size", type=int, default=256)
    ap.add_argument("--lr", type=float, default=1e-3)
    ap.add_argument("--output_dir", default=".")
    args = ap.parse_args()

    tensors, le, gene_names, n_classes = load_and_prepare(
        args.input, args.label_key, args.n_hvg
    )

    n_genes = tensors["X_train"].shape[1]
    model = CellTypeClassifier(n_genes, n_classes)
    print(f"\n{model}")

    model, train_hist, val_hist = train_model(
        model, tensors, n_epochs=args.epochs,
        batch_size=args.batch_size, lr=args.lr,
    )

    y_pred, cm, acc, f1 = evaluate(model, tensors, le)

    imp, top_genes, top_scores = get_gene_importance(model, tensors, gene_names)

    save_plots(train_hist, val_hist, cm, le.classes_.tolist(),
               top_genes, top_scores, args.output_dir)

    ckpt_path = Path(args.output_dir) / "celltype_classifier.pt"
    torch.save({
        "state_dict": model.cpu().state_dict(),
        "classes": le.classes_.tolist(),
        "genes": gene_names,
        "n_genes": n_genes,
        "n_classes": n_classes,
        "test_acc": acc,
        "test_f1": f1,
    }, ckpt_path)
    print(f"Checkpoint saved: {ckpt_path}")


if __name__ == "__main__":
    main()
