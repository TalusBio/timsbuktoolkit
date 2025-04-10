import io
import json
from dataclasses import dataclass
from pathlib import Path

from speclib_builder.base import SpeclibElement


@dataclass
class SpeclibWriter:
    path: Path
    _handle: io.TextIOWrapper | None = None
    _opened: bool = False
    _num_written = 0

    def __post_init__(self):
        if isinstance(self.path, str):
            self.path = Path(self.path)

        if not isinstance(self.path, Path):
            raise TypeError(
                f"Expected path to be a string or Path object, got {type(self.path)}"
            )

    def append(self, data: SpeclibElement):
        if not self._opened:
            raise RuntimeError("Writer is not open")
        if self._handle is None:
            raise RuntimeError(
                "Unregistered handle. Make sure you are using this as a context manager"
            )
        json.dump(data.as_dict(), self._handle)
        self._handle.write("\n")
        self._num_written += 1

    def __enter__(self):
        self.path.parent.mkdir(parents=True, exist_ok=True)
        print("Writing to", self.path)
        self._handle = open(self.path, "w")
        self._opened = True
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self._opened = False
        self._handle.flush()
        self._handle.close()
        print(f"Wrote {self._num_written} entries to {self.path}")
