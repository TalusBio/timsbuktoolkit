from dataclasses import dataclass

import numpy as np
from uniplot import plot

from .linear_regression import LinearRegression


@dataclass
class MultiPivotRegression:
    offsets: np.ndarray | None = None
    regressors: list[LinearRegression] | None = None

    def predict(self, X):
        predictions = self._predict(X)
        return predictions.sum(axis=0)

    def _predict(self, X):
        predictions = np.array([
            reg.predict(X - np.delete(self.offsets.T, i))
            for i, reg in enumerate(self.regressors)
        ])
        return predictions

    def fit(
        self,
        X,
        y,
        x_range: tuple[float, float],
        n_kernels: int = 10,
        learning_rate: float = 0.001,
        n_iterations: int = 100,
    ):
        n_samples, n_features = X.shape
        if self.offsets is None:
            self.offsets = np.linspace(x_range[0], x_range[1], n_kernels)

        if self.regressors is None:
            self.regressors = [LinearRegression() for _ in range(n_kernels)]

        for i, reg in enumerate(self.regressors):
            reg.fit(X - np.delete(self.offsets.T, i), y, verbose=False)

        for i in range(n_iterations):
            cp = self._predict(X)
            for j, reg in enumerate(self.regressors):
                rest = y - np.sum(np.delete(cp, j, axis=0), axis=0)
                reg.fit(
                    X - np.delete(self.offsets.T, j),
                    rest,
                    verbose=False,
                    learning_rate=1e-4,
                    n_iterations=200,
                )

            pred = self.predict(X)
            if i % 5 == 0:
                print(f"Iteration {i}, Loss: {np.mean(np.abs(pred - y)):.4f}")

        plot(X.flatten(), y.flatten(), title="Loess")
        plot(pred, y.flatten(), title="Loess")
