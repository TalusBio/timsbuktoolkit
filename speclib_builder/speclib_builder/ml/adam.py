from dataclasses import dataclass, field
import numpy as np


@dataclass
class AdamOptimizer:
    learning_rate: float = 0.001
    beta1: float = 0.9
    beta2: float = 0.99
    epsilon: float = 1e-8
    m_w: np.ndarray = field(init=False)
    v_w: np.ndarray = field(init=False)
    m_b: float = field(init=False)
    v_b: float = field(init=False)
    t: int = field(default=0, init=False)
    n_features: int = field(default=0, init=False)

    def init_params(self, n_features):
        self.m_w = np.zeros(n_features)
        self.v_w = np.zeros(n_features)
        self.m_b = 0.0
        self.v_b = 0.0
        self.t = 0
        self.n_features = n_features

    def update(self, dw, db):
        assert dw.shape == (self.n_features,), (
            f"Expected dw to have shape ({self.n_features},), got {dw.shape}"
        )
        assert db.shape == ()
        self.t += 1

        # Update moments
        self.m_w = self.beta1 * self.m_w + (1 - self.beta1) * dw
        self.v_w = self.beta2 * self.v_w + (1 - self.beta2) * np.square(dw)
        self.m_b = self.beta1 * self.m_b + (1 - self.beta1) * db
        self.v_b = self.beta2 * self.v_b + (1 - self.beta2) * (db**2)

        # Bias correction
        m_w_hat = self.m_w / (1 - self.beta1**self.t)
        v_w_hat = self.v_w / (1 - self.beta2**self.t)
        m_b_hat = self.m_b / (1 - self.beta1**self.t)
        v_b_hat = self.v_b / (1 - self.beta2**self.t)

        assert m_w_hat.shape == (self.n_features,)
        assert v_w_hat.shape == (self.n_features,)
        assert m_b_hat.shape == ()
        assert v_b_hat.shape == ()

        # Compute updates
        w_update = self.learning_rate * m_w_hat / (np.sqrt(v_w_hat) + self.epsilon)
        b_update = self.learning_rate * m_b_hat / (np.sqrt(v_b_hat) + self.epsilon)

        assert w_update.shape == (self.n_features,)

        return w_update, b_update
