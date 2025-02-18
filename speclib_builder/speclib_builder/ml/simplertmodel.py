import numpy as np
from rustyms import LinearPeptide

from .linear_regression import LinearRegression
from .ohe_peptide_embed import FEATURES, PeptideOHEEmbedder


class SimpleRTModel(LinearRegression):
    def __init__(self, weights: np.ndarray | None = None, bias: float | None = None):
        super().__init__(weights, bias)
        self.embedder = PeptideOHEEmbedder()

    def predict_peptide(self, peptide: LinearPeptide) -> float:
        return self.predict(self.embedder.embed(peptide)).item()

    def preict_peptides(self, peptides: list[LinearPeptide]) -> np.ndarray:
        embeddings = np.stack([self.embedder.embed(peptide) for peptide in peptides], axis=0)
        assert embeddings.shape[0] == len(peptides)
        assert embeddings.shape[1] == FEATURES
        assert len(embeddings.shape) == 2
        return self.predict(embeddings)

    def fit_peptides(self, peptides: list[LinearPeptide], rts: np.ndarray):
        embeddings = np.stack([self.embedder.embed(peptide) for peptide in peptides], axis=0)
        assert embeddings.shape[0] == len(peptides)
        assert embeddings.shape[1] == FEATURES
        assert len(embeddings.shape) == 2
        self.fit(embeddings, rts)

    def fit_stripped_sequence(
        self, peptides: list[str], masses: list[float], rts: list[float], **kwargs
    ):
        embeddings = np.stack(
            [
                self.embedder.embed_stripped_sequence(peptide, mass)
                for peptide, mass in zip(peptides, masses, strict=True)
            ],
            axis=0,
        )
        assert embeddings.shape[0] == len(peptides)
        assert embeddings.shape[1] == FEATURES
        assert len(embeddings.shape) == 2
        self.fit(embeddings, np.array(rts), **kwargs)


def default_rt_model() -> SimpleRTModel:
    model = SimpleRTModel(
        weights=np.array([
            -26.36735781,
            0.0,
            -53.56934945,
            -55.21925322,
            178.27425793,
            -66.75431489,
            -209.06799874,
            125.83324634,
            -150.87179673,
            148.44146624,
            63.25018127,
            -91.4846459,
            -17.6613223,
            -83.61554405,
            -196.30278774,
            -68.4594332,
            -50.41868969,
            50.06798268,
            187.24671262,
            20.80340359,
            -486.74342897,
            0.0,
            -433.18876836,
            -468.34916512,
            -527.09505792,
            -448.14074773,
            -467.64202839,
            -531.7258505,
            -499.37551192,
            -538.98567687,
            -493.35713573,
            -440.11155882,
            -466.97560832,
            -473.17200096,
            -462.76910048,
            -454.20171251,
            -472.42193567,
            -516.75912556,
            -521.42637442,
            -488.67922573,
            2.08151486,
            0.0,
            -3.84035216,
            0.75315676,
            9.55433749,
            20.95398651,
            40.11952704,
            -10.65939963,
            -18.59012557,
            -1.93303538,
            22.55617111,
            14.19808588,
            12.09946731,
            13.13722242,
            29.0289534,
            12.82120924,
            10.56736659,
            -1.14792772,
            -7.66936585,
            37.4589648,
            328.94046503,
            -567.69946649,
        ]),
        bias=np.float64(-567.6994664944469),
    )

    return model
