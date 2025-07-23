import io
import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Dict, Any
import time

from loguru import logger
import msgpack
import zstandard as zstd

from speclib_builder.base import SpeclibElement


@dataclass
class FormatWriter:
    """Individual format writer"""

    path: Path
    format_name: str
    _handle: Any = None
    _num_written: int = 0
    _write_time: float = 0.0

    def open(self):
        self.path.parent.mkdir(parents=True, exist_ok=True)

        # TODO: This is super ugly ... make this a protocol...
        if self.format_name == "ndjson":
            warn_extension(self.path, "ndjson")
            self._handle = open(self.path, "w")
        elif self.format_name == "ndjson_zstd":
            warn_extension(self.path, "ndjson.zst")
            self._handle = zstd.open(self.path, "wt", encoding="utf-8")
        elif self.format_name == "msgpack":
            warn_extension(self.path, "msgpack")
            self._handle = open(self.path, "wb")
        elif self.format_name == "msgpack_zstd":
            warn_extension(self.path, "msgpack.zst")
            self._handle = zstd.open(self.path, "wb")
        else:
            raise ValueError(f"Unknown format: {self.format_name}")

    def append(self, data: SpeclibElement):
        start_time = time.time()
        data_dict = data.model_dump()

        if self.format_name == "ndjson":
            json.dump(data_dict, self._handle)
            self._handle.write("\n")

        elif self.format_name == "ndjson_zstd":
            json.dump(data_dict, self._handle)
            self._handle.write("\n")

        elif self.format_name == "msgpack":
            packed = msgpack.packb(data_dict, use_bin_type=True)
            self._handle.write(packed)

        elif self.format_name == "msgpack_zstd":
            packed = msgpack.packb(data_dict, use_bin_type=True)
            self._handle.write(packed)

        self._write_time += time.time() - start_time
        self._num_written += 1

    def close(self):
        self._handle.close()

    def get_stats(self) -> Dict:
        file_size = self.path.stat().st_size if self.path.exists() else 0
        return {
            "format": self.format_name,
            "path": str(self.path),
            "num_written": self._num_written,
            "write_time": self._write_time,
            "file_size": file_size,
            "avg_time_per_record": self._write_time / max(1, self._num_written),
        }


@dataclass
class SpeclibWriter:
    path: Path
    # format: str = "ndjson"
    file_format: str = "msgpack_zstd"
    _opened: bool = False
    _num_written: int = 0
    _writers: List[FormatWriter] = field(default_factory=list)

    def __post_init__(self):
        if isinstance(self.path, str):
            self.path = Path(self.path)

        if not isinstance(self.path, Path):
            raise TypeError(
                f"Expected path to be a string or Path object, got {type(self.path)}"
            )

        self._writer = FormatWriter(
            path=self.path,
            format_name=self.file_format,
        )

    def append(self, data: SpeclibElement):
        if not self._opened:
            raise RuntimeError("Writer is not open")

        if self._writer is None:
            raise RuntimeError(
                "Unregistered handle. Make sure you are using this as a context manager"
            )

        self._num_written += 1
        self._writer.append(data)

    def __enter__(self):
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self._writer.open()
        self._opened = True
        self._num_written = 0

        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self._opened = False
        self._writer.close()

        stats = self._writer.get_stats()
        logger.debug(
            f"Finished writing {stats['num_written']} records to {stats['path']} "
            f"in {stats['write_time']:.3f}s "
            f"({stats['avg_time_per_record'] * 1000:.2f}ms per record)"
        )


def warn_extension(path: Path, ext: str):
    """Warn if the file extension is not as expected."""
    if not path.name.endswith(f".{ext}"):
        logger.warning(
            f"Warning: File {path} does not have the expected .{ext} extension. "
            "This may cause issues with some tools."
        )
        time.sleep(5)  # Small delay to ensure the warning is visible in logs
