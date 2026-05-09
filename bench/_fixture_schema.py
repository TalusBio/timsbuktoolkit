"""Pydantic schema and loader for fixture TOML files."""

from __future__ import annotations

import tomllib
from pathlib import Path
from typing import Any

from pydantic import BaseModel, ConfigDict, Field, ValidationInfo, field_validator


def _require_s3(value: str, field: str) -> str:
    if not value.startswith("s3://"):
        raise ValueError(f"{field} must be an s3:// URI, got {value!r}")
    return value


class FixtureInputs(BaseModel):
    model_config = ConfigDict(extra="forbid")

    fasta: str
    speclib: str
    raw: str
    entrapment_fasta: str | None = None
    calibration_speclib: str | None = None

    @field_validator("fasta", "speclib", "raw")
    @classmethod
    def _required_s3(cls, v: str, info: ValidationInfo) -> str:
        return _require_s3(v, info.field_name)

    @field_validator("entrapment_fasta", "calibration_speclib")
    @classmethod
    def _optional_s3(cls, v: str | None, info: ValidationInfo) -> str | None:
        if v is None:
            return v
        return _require_s3(v, info.field_name)


class Fixture(BaseModel):
    model_config = ConfigDict(extra="forbid")

    name: str
    description: str = ""
    inputs: FixtureInputs
    config: dict[str, Any] = Field(default_factory=dict)

    def has_entrapment(self) -> bool:
        return self.inputs.entrapment_fasta is not None

    def has_calibration_speclib(self) -> bool:
        return self.inputs.calibration_speclib is not None


def load_fixture(path: str | Path) -> Fixture:
    """Load a fixture TOML and validate it against the schema."""
    p = Path(path)
    with p.open("rb") as f:
        raw = tomllib.load(f)
    return Fixture.model_validate(raw)
