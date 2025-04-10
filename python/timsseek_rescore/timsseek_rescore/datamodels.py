from dataclasses import dataclass
from pathlib import Path


@dataclass
class Report:
    targets_at_1: int
    targets_at_5: int
    targets_at_10: int

    def save_to_toml(self, path: Path):
        with open(path, "w") as f:
            f.write("[report]\n")
            f.write(f"targets_at_1 = {self.targets_at_1}\n")
            f.write(f"targets_at_5 = {self.targets_at_5}\n")
            f.write(f"targets_at_10 = {self.targets_at_10}\n")
