from dataclasses import dataclass, field
import numpy as np
from uniplot import plot
import math

from .adam import AdamOptimizer


@dataclass
class LinearRegression:
    weights: np.ndarray | None = None
    bias: float | None = None

    def _compute_loss(self, y_true, y_pred):
        return np.mean((y_true - y_pred) ** 2)

    def fit(
        self,
        X,
        y,
        *,
        learning_rate: float = 0.01,
        n_iterations: int = 100_000,
        min_diff: float = 1e-6,
        patience: int = 100,
        plateau: int = 10,
        verbose: bool = True,
    ):
        n_samples, n_features = X.shape
        optimizer = AdamOptimizer(learning_rate=learning_rate)
        optimizer.init_params(n_features)
        loss_history = []
        self.weights = np.zeros(n_features)
        self.bias = 0.0
        best_loss = np.inf
        best_iter = 0

        for i in range(n_iterations):
            y_predicted = np.dot(X, self.weights) + self.bias
            current_loss = self._compute_loss(y, y_predicted)
            loss_history.append(current_loss)

            if i % 200 == 0 and verbose:
                mae = np.mean(np.abs(y_predicted - y))
                print(f"Iteration {i}, Loss: {current_loss:.4f} MAE: {mae:.4f}")

            dw = (1 / math.sqrt(n_samples)) * np.dot(X.T, y_predicted - y)
            db = (1 / math.sqrt(n_samples)) * np.sum(y_predicted - y)

            w_update, b_update = optimizer.update(dw, db)
            self.weights -= w_update
            self.bias -= b_update

            if current_loss < (best_loss - min_diff):
                best_loss = current_loss
                best_iter = i

            if i - best_iter > plateau:
                if verbose:
                    print(
                        f"Plateau reached at iteration {i} with loss {current_loss:.6f}"
                    )
                learning_rate /= 2
                optimizer = AdamOptimizer(learning_rate=learning_rate)
                optimizer.init_params(n_features)

            patience_count = i - best_iter
            if patience_count > patience:
                if verbose:
                    print(
                        f"Early stopping at iteration {i} with loss {current_loss:.6f}"
                    )
                break

        if verbose:
            plot(np.log1p(loss_history), title="Loss History")
            plot(y_predicted, y, title="Prediction vs Actual")

        return loss_history

    def predict(self, X):
        if self.weights is None or self.bias is None:
            raise ValueError("Model not fitted")
        return np.dot(X, self.weights) + self.bias
