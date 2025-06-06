from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional

import matplotlib.pyplot as plt
import mokapot
import numpy as np
import polars as pl
import torch
import torch.nn as nn
import torch.nn.functional as F
from rich.pretty import pprint
from torch.utils.data import DataLoader, TensorDataset
from uniplot import histogram

from ..datamodels import Report
from ..folding import to_folds
from ..plotting import plot_importances

# FFN_ACTIVATION = nn.SELU
FFN_ACTIVATION = nn.ReLU


class AsymmetricMarginBCELoss(nn.Module):
    def __init__(self, *, margin_0=0.1, margin_1=0.4):
        """
        margin_0: margin for negative class (0s)
        margin_1: margin for positive class (1s)

        This loss pushes negative predictions further from decision boundary

        Larger margin_0 makes the model more conservative about predicting 1s
            (reduces false positives)
        Smaller margin_1 means we're more lenient about false negatives
        Both margins are clamped to keep predictions in [0,1] range
        """
        super().__init__()
        self.margin_0 = margin_0
        self.margin_1 = margin_1

    def forward(self, predictions, targets):
        # Add margins to predictions based on true class
        adjusted_preds = torch.where(
            targets == 1, predictions + self.margin_1, predictions - self.margin_0
        )

        # Clamp to valid probability range
        adjusted_preds = torch.clamp(adjusted_preds, 0, 1)

        # Compute BCE loss
        loss = F.binary_cross_entropy(adjusted_preds, targets, reduction="mean")

        return loss


class WeightedBCELoss(nn.Module):
    def __init__(self, pos_weight=0.2, neg_weight=1.0):
        """
        pos_weight: weight for positive class (1s)
        neg_weight: weight for negative class (0s)
        """
        super().__init__()
        self.pos_weight = pos_weight
        self.neg_weight = neg_weight

    def forward(self, predictions, targets):
        # Create weight tensor based on targets
        weights = torch.where(
            targets == 1,
            torch.tensor(self.pos_weight, device=targets.device),
            torch.tensor(self.neg_weight, device=targets.device),
        )

        # Standard BCE loss
        bce_loss = F.binary_cross_entropy(predictions, targets, reduction="none")

        # Apply weights
        weighted_loss = weights * bce_loss

        return weighted_loss.mean()


class FocalLoss3(nn.Module):
    def __init__(self, alpha=0.25, gamma=2.0, reduce=True):
        """
        Focal Loss: (1 - p)^gamma * log(p) for positive class
                   p^gamma * log(1-p) for negative class

        alpha: weighing factor for positive class
        gamma: focusing parameter that reduces the loss contribution from easy examples
        """
        super().__init__()
        if alpha < 0 or alpha > 1:
            raise ValueError("Alpha must be in [0, 1]")
        self.alpha = alpha
        self.gamma = gamma
        self.reduce = reduce

    def forward(self, predictions, targets):
        # BCE loss
        bce_loss = F.binary_cross_entropy(predictions, targets, reduction="none")

        # Focal term
        pt = torch.where(targets == 1, predictions, 1 - predictions)
        focal_term = (1 - pt) ** self.gamma

        # Alpha weighing
        alpha_weight = torch.where(
            targets == 1,
            torch.tensor(self.alpha, device=targets.device),
            torch.tensor(1 - self.alpha, device=targets.device),
        )

        loss = alpha_weight * focal_term * bce_loss
        if self.reduce:
            return loss.mean()
        else:
            return loss


class ResidualMLPBlock(nn.Module):
    def __init__(self, input_dim: int, dropout: float = 0.0):
        super().__init__()
        self.input_layer = nn.Linear(input_dim, input_dim)
        # self.activation = nn.SELU()
        self.activation = FFN_ACTIVATION()
        self.dropout = nn.Dropout(dropout)
        self.batchnorm = nn.BatchNorm1d(input_dim)

    def forward(self, x):
        residual = x
        x = self.input_layer(x)
        x = self.batchnorm(x)
        x = self.activation(x)
        x = self.dropout(x)
        return x + residual


class ResidualMLPBlockI(nn.Module):
    def __init__(self, input_dim: int, output_dim: int, dropout: float = 0.0):
        super().__init__()
        self.input_layer = nn.Linear(input_dim, output_dim)
        # self.activation = nn.SELU()
        self.activation = FFN_ACTIVATION()
        self.dropout = nn.Dropout(dropout)
        self.batchnorm = nn.BatchNorm1d(output_dim)
        self.input_projection = nn.Linear(input_dim, output_dim)

    def forward(self, x):
        residual = x
        x = self.input_layer(x)
        x = self.activation(x)
        x = self.dropout(x)
        x = self.batchnorm(x)
        return x + self.input_projection(residual)


class BinaryClassifier(nn.Module):
    def __init__(
        self,
        input_dim: int,
        nhidden_layers: int = 2,
        hidden_dims: int = 48,
        dropout: float = 0.00,
    ):
        super().__init__()

        # ACTIVATION = nn.ReLU
        # Selu seems to be critical to train deeper networks.
        # ACTIVATION = nn.SELU

        layers = []
        layers.append(nn.BatchNorm1d(input_dim))
        layers.append(
            ResidualMLPBlockI(
                input_dim=input_dim, output_dim=hidden_dims, dropout=dropout
            )
        )
        self.input_layer = nn.Sequential(*layers)

        layers = []
        for _ in range(nhidden_layers):
            layers.append(ResidualMLPBlock(input_dim=hidden_dims, dropout=dropout))

        self.hidden_layers = nn.ModuleList(layers)

        layers = []
        layers.append(nn.Linear(hidden_dims, 1))
        layers.append(nn.Sigmoid())
        self.output_layer = nn.Sequential(*layers)

        pprint(self)
        nparams = sum(p.numel() for p in self.parameters())
        nparams_trainable = sum(p.numel() for p in self.parameters() if p.requires_grad)
        pprint(f"Number of parameters: {nparams}")
        pprint(f"Number of trainable parameters: {nparams_trainable}")

    def forward(self, x):
        x = self.input_layer(x)
        for layer in self.hidden_layers:
            # x = x + layer(layer_norm(x))
            x = layer(x)

        return self.output_layer(x)

    def train_epoch(self, dataloader, device, optimizer, criterion) -> float:
        self.train()
        train_losses = []
        for x_batch, y_batch in dataloader:
            x_batch, y_batch = x_batch.to(device), y_batch.to(device)

            optimizer.zero_grad()
            y_pred = self(x_batch)
            loss = criterion(y_pred, y_batch.view(-1, 1))
            loss.backward()
            optimizer.step()

            train_losses.append(loss.item())

        self.eval()
        return np.mean(train_losses).item()

    def val_epoch(self, dataloader, device, criterion) -> float:
        val_losses = []
        with torch.no_grad():
            for x_batch, y_batch in dataloader:
                x_batch, y_batch = (
                    x_batch.to(device),
                    y_batch.to(device),
                )
                y_pred = self(x_batch)
                val_loss = criterion(y_pred, y_batch.view(-1, 1))
                val_losses.append(val_loss.item())

        avg_val_loss = np.mean(val_losses)
        return avg_val_loss


@dataclass
class MLPKFoldModel:
    """
    PyTorch implementation of K-Fold cross validation.
    Each fold contains:
    1. One fold for training
    2. One fold for validation (early stopping)
    3. Remaining folds for inference
    """

    folds: List[tuple[torch.Tensor, torch.Tensor]]  # List of (features, targets) tuples
    models: List[Optional[BinaryClassifier]]
    scores: List[Optional[torch.Tensor]]
    device: torch.device

    def determine_device(self):
        # Determine the device
        if torch.cuda.is_available():
            self.device = torch.device("cuda")
            print("Using CUDA (GPU) for training.")
        elif hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
            self.device = torch.device("mps")
            print("Using MPS (Apple Silicon GPU) for training.")
        else:
            self.device = torch.device("cpu")
            print("Using CPU for training.")

    @staticmethod
    def from_folds(
        folds: List[tuple[torch.Tensor, torch.Tensor]], device: torch.device
    ):
        return MLPKFoldModel(
            folds=folds,
            models=[None] * len(folds),
            scores=[None] * len(folds),
            device=device,
        )

    def train(
        self,
        batch_size: int = 250,
        epochs: int = 20,
        learning_rate: float = 1e-4,
        pos_weight: float = 0.2,
        **kwargs,
    ):
        for i in range(len(self.folds)):
            print(f"Training model {i}/{len(self.folds)}")

            # Prepare data
            train_data = self.folds[i]
            val_data = self.folds[(i + 1) % len(self.folds)]

            train_dataset = TensorDataset(train_data[0], train_data[1])
            val_dataset = TensorDataset(val_data[0], val_data[1])

            train_loader = DataLoader(
                train_dataset,
                batch_size=batch_size,
                shuffle=True,
                num_workers=0,
            )
            val_loader = DataLoader(
                val_dataset,
                batch_size=batch_size,
                num_workers=0,
            )

            # Initialize model
            model = BinaryClassifier(input_dim=train_data[0].shape[1], **kwargs).to(
                self.device
            )
            optimizer = torch.optim.AdamW(model.parameters(), lr=learning_rate)
            # Choose one of the following loss functions based on your needs:
            # criterion = WeightedBCELoss(pos_weight=pos_weight, neg_weight=1.0)
            # criterion = AsymmetricMarginBCELoss(margin_0=0.5, margin_1=0.2)
            # criterion = WeightedBCELoss2(fneg_weight=pos_weight)

            # Usually gamma is positive bc the desire is to emphasize well classified
            # BUT ... since we know we have a lot of false positives, we want to give
            # "false targets" a lower weight.
            criterion = FocalLoss3(alpha=pos_weight, gamma=0.5)

            best_val_loss = float("inf")
            patience = 5
            patience_counter = 0
            best_model = None

            for epoch in range(epochs):
                train_loss_mean = model.train_epoch(
                    dataloader=train_loader,
                    device=self.device,
                    optimizer=optimizer,
                    criterion=criterion,
                )
                val_loss_mean = model.val_epoch(
                    dataloader=val_loader,
                    device=self.device,
                    criterion=criterion,
                )

                # Early stopping
                print(
                    f"Epoch {epoch}: train_loss = {train_loss_mean:.4f}, val_loss = {val_loss_mean:.4f}"
                )
                if val_loss_mean < best_val_loss:
                    best_val_loss = val_loss_mean
                    best_model = model.state_dict()
                    patience_counter = 0
                else:
                    patience_counter += 1
                    if patience_counter >= patience:
                        print(
                            f"Early stopping at epoch {epoch}."
                            f" Best validation loss: {best_val_loss:.4f}"
                        )
                        break

            # Save best model
            model.load_state_dict(best_model)
            self.models[i] = model

    def score(self, batch_size: int = 32):
        for i in range(len(self.folds)):
            print(f"Scoring fold {i}/{len(self.folds)}")
            fold_data = self.folds[i]
            dataset = TensorDataset(fold_data[0], fold_data[1])
            loader = DataLoader(dataset, batch_size=batch_size)

            scores_list = []

            for j, model in enumerate(self.models):
                if j == i or j == (i + 1) % len(self.folds):
                    continue

                model.eval()
                fold_scores = []

                with torch.no_grad():
                    for x_batch, _ in loader:
                        x_batch = x_batch.to(self.device)
                        predictions = model(x_batch)
                        fold_scores.append(predictions.cpu())

                scores_list.append(torch.cat(fold_scores))

            self.scores[i] = torch.stack(scores_list).mean(dim=0)

    def get_importances2(self, feat_names: list[str] | None = None):
        """
        Calculate feature importance using gradients after BatchNorm layer.
        Uses hooks to capture intermediate gradients.
        """
        importances = []

        for model in self.models:
            importance_dict = {}
            post_bn_gradients = []

            # Register hook to capture gradients after BatchNorm
            def hook_fn(module, grad_input, grad_output):
                # post_bn_gradients.append(grad_input[0].detach().cpu())
                post_bn_gradients.append(grad_output[0].detach().cpu())

            # Get the BatchNorm layer
            bn_layer = model.input_layer[0]  # First layer is BatchNorm
            hook = bn_layer.register_backward_hook(hook_fn)

            # Compute gradients for each fold
            for i, (features, targets) in enumerate(self.folds):
                features = features.to(self.device)
                features.requires_grad_(True)
                post_bn_gradients = []  # Reset for each batch

                # Forward and backward pass
                output = model(features)
                output.sum().backward()

                # Average gradients across samples in the fold
                fold_grads = post_bn_gradients[0].abs().mean(dim=0)

                # Update importance dictionary
                if not importance_dict:
                    if feat_names is None:
                        feat_names = [f"feature_{j}" for j in range(features.shape[1])]
                    importance_dict = {
                        k: fold_grads[j].item() for j, k in enumerate(feat_names)
                    }
                else:
                    for j, k in enumerate(feat_names):
                        importance_dict[k] += fold_grads[j].item()

                features.requires_grad_(False)
                model.zero_grad()

            for key in importance_dict:
                importance_dict[key] /= len(self.folds)

            importances.append(importance_dict)

            hook.remove()

        imps_order = sorted(importances[0].items(), key=lambda x: x[1], reverse=True)
        out = {k[0]: [w.get(k[0], 0) for w in importances] for k in imps_order}
        return out

    def concat_scores(self):
        if self.scores[0] is None:
            raise ValueError("Scores not computed")
        return torch.cat(self.scores)

    def concat_targets(self):
        return torch.cat([fold[1] for fold in self.folds])


def to_torch_folds(shuffled_df: pl.LazyFrame, num_folds: int, cols):
    tmp = to_folds(shuffled_df=shuffled_df, cols=cols, num_folds=num_folds)
    out = [
        (torch.from_numpy(x[0]).float(), torch.from_numpy(x[1]).float()) for x in tmp
    ]
    return out


def mlp_stuff(
    df_use: pl.LazyFrame, cols, output_dir: Path, pos_weight: float = 0.05, **kwargs
):
    pprint("Shuffling")
    df_use = df_use.sample(frac=1).reset_index(drop=True, inplace=False)
    folds = to_torch_folds(shuffled_df=df_use, num_folds=5, cols=cols)
    fold_model = MLPKFoldModel.from_folds(folds, device="cpu")
    # Moving the data back and forth is slower ...
    # fold_model.determine_device()
    fold_model.train(pos_weight=pos_weight, **kwargs)
    fold_model.score()
    importances = fold_model.get_importances2(feat_names=cols.feature_columns)
    pprint(importances)
    fig = plot_importances(importances)
    outfile = output_dir / "importances_nn.png"
    fig.savefig(outfile)
    pprint(f"Wrote {outfile}")
    plt.close()

    ctargs = fold_model.concat_targets().numpy().flatten()
    cscores = fold_model.concat_scores().numpy().flatten()
    df_use["rescore_score"] = cscores
    df_use["qvalue"] = mokapot.qvalues.qvalues_from_scores(cscores, ctargs == 1)
    df_use = df_use.sort_values("rescore_score", ascending=False)
    outfile = output_dir / "rescored_values_nn.parquet"
    df_use.to_parquet(outfile, index=False)
    pprint(f"Wrote {outfile}")
    order = np.argsort(-cscores)
    ctargs = ctargs[order]
    cscores = cscores[order]
    qvals = mokapot.qvalues.qvalues_from_scores(cscores, ctargs == 1)
    for ct in [0.01, 0.05, 0.1, 0.5, 1.0]:
        num_at_thresh = int(np.sum(ctargs[qvals < ct]))
        ssc = cscores[qvals < ct]
        if len(ssc) == 0:
            pprint(f"No scores at {ct}")
            continue
        score_at_thresh = np.min(ssc)
        pprint(f"Score at {ct}: {score_at_thresh}")
        pprint(f"Number of targets at {ct}: {num_at_thresh}")

    target_preds = cscores[ctargs == 1]
    decoy_preds = cscores[ctargs == 0]
    histogram(target_preds, title="Target scores")
    histogram(decoy_preds, title="Decoy scores")

    report = Report(
        targets_at_1=np.sum(ctargs[qvals < 0.01]).item(),
        targets_at_5=np.sum(ctargs[qvals < 0.05]).item(),
        targets_at_10=np.sum(ctargs[qvals < 0.1]).item(),
    )
    pprint(report)
    outfile = output_dir / "report_nn.toml"
    report.save_to_toml(outfile)
    pprint(f"Wrote {outfile}")
    return np.sum(ctargs[qvals < 0.01]).item()
