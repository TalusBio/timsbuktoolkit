"""Pydantic schema and loader for fixture TOML files."""

from __future__ import annotations

import tomllib
from pathlib import Path
from typing import Any, Literal

from pydantic import (
    BaseModel,
    ConfigDict,
    Field,
    ValidationInfo,
    field_validator,
    model_validator,
)


def _require_uri(value: str, field: str) -> str:
    """Accept s3:// URIs or absolute local paths. Expands `~` to home.

    Relative paths and other URI schemes (file://, http://, etc.) are
    rejected. Existence is NOT checked here — the runner will surface a
    clearer error at use time.
    """
    if value.startswith("s3://"):
        return value
    expanded = str(Path(value).expanduser())
    if Path(expanded).is_absolute():
        return expanded
    raise ValueError(
        f"{field} must be an s3:// URI or an absolute local path, got {value!r}"
    )


class FixtureInputs(BaseModel):
    model_config = ConfigDict(extra="forbid")

    target_peptides: str | None = None
    speclib: str
    raw: str
    entrapment_peptides: str | None = None
    entrapment_ratio: float | None = None
    entrapment_mode: Literal["foreign", "shuffled"] | None = None
    pairing: str | None = None
    calibration_speclib: str | None = None

    @field_validator("speclib", "raw")
    @classmethod
    def _required_uri(cls, v: str, info: ValidationInfo) -> str:
        return _require_uri(v, info.field_name or "")

    @field_validator(
        "target_peptides", "entrapment_peptides", "calibration_speclib", "pairing"
    )
    @classmethod
    def _optional_uri(cls, v: str | None, info: ValidationInfo) -> str | None:
        if v is None:
            return v
        return _require_uri(v, info.field_name or "")

    @model_validator(mode="after")
    def _entrap_consistency(self) -> "FixtureInputs":
        has_pep = self.entrapment_peptides is not None
        has_ratio = self.entrapment_ratio is not None
        has_mode = self.entrapment_mode is not None
        if has_pep and self.target_peptides is None:
            raise ValueError(
                "entrapment_peptides requires target_peptides to also be set"
            )
        if has_pep and not (has_ratio and has_mode):
            raise ValueError(
                "entrapment_peptides requires both entrapment_ratio and entrapment_mode"
            )
        if (has_ratio or has_mode) and not has_pep:
            raise ValueError(
                "entrapment_ratio / entrapment_mode set without entrapment_peptides"
            )
        if self.pairing is not None and self.entrapment_mode != "shuffled":
            raise ValueError(
                "pairing field only valid when entrapment_mode == 'shuffled'"
            )
        if self.entrapment_ratio is not None and self.entrapment_ratio < 1.0:
            raise ValueError(
                f"entrapment_ratio must be >= 1.0, got {self.entrapment_ratio}"
            )
        return self


class Fixture(BaseModel):
    model_config = ConfigDict(extra="forbid")

    name: str
    description: str = ""
    tags: list[str] = Field(default_factory=list)
    inputs: FixtureInputs
    config: dict[str, Any] = Field(default_factory=dict)

    def has_entrapment(self) -> bool:
        return self.inputs.entrapment_peptides is not None

    def has_pairing(self) -> bool:
        return self.inputs.pairing is not None

    def has_calibration_speclib(self) -> bool:
        return self.inputs.calibration_speclib is not None


def load_fixture(path: str | Path) -> Fixture:
    """Load a fixture TOML and validate it against the schema."""
    p = Path(path)
    with p.open("rb") as f:
        raw = tomllib.load(f)
    return Fixture.model_validate(raw)
