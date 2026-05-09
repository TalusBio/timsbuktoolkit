"""Tiny shellout wrapper around `aws s3 cp` / `aws s3 ls` / `aws s3 sync`."""

from __future__ import annotations

import subprocess
from pathlib import Path

from loguru import logger


def _aws_cp(src: str, dst: str, recursive: bool = False) -> None:
    cmd = ["aws", "s3", "cp", src, dst]
    if recursive:
        cmd.append("--recursive")
    logger.info("$ {}", " ".join(cmd))
    subprocess.run(cmd, check=True)


def _aws_sync(src: str, dst: str) -> None:
    cmd = ["aws", "s3", "sync", src, dst]
    logger.info("$ {}", " ".join(cmd))
    subprocess.run(cmd, check=True)


def s3_object_exists(s3_uri: str) -> bool:
    """Return True if the S3 object at `s3_uri` already exists."""
    if not s3_uri.startswith("s3://"):
        raise ValueError(f"not an s3:// URI: {s3_uri}")
    cmd = ["aws", "s3", "ls", s3_uri]
    logger.info("$ {}", " ".join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True)
    # `aws s3 ls` exits 0 with empty output when the object isn't there
    # (some configurations) or non-zero in others. Accept either signal:
    # object exists iff exit code 0 AND output non-empty.
    return result.returncode == 0 and bool(result.stdout.strip())


def s3_download_file(s3_uri: str, local_path: str) -> None:
    if not s3_uri.startswith("s3://"):
        raise ValueError(f"not an s3:// URI: {s3_uri}")
    _aws_cp(s3_uri, local_path)


def s3_upload_file(local_path: str, s3_uri: str, skip_if_exists: bool = False) -> None:
    if not s3_uri.startswith("s3://"):
        raise ValueError(f"not an s3:// URI: {s3_uri}")
    if skip_if_exists and s3_object_exists(s3_uri):
        logger.info("s3 upload: skipping {} (already exists)", s3_uri)
        return
    _aws_cp(local_path, s3_uri)


def s3_upload_dir(local_dir: str, s3_uri: str, idempotent: bool = False) -> None:
    """Upload a directory tree to S3.

    `idempotent=True` uses `aws s3 sync` which only transfers files whose
    size or modified-time differ. `idempotent=False` (default) does a fresh
    recursive copy, overwriting whatever is at the destination.
    """
    if not s3_uri.startswith("s3://"):
        raise ValueError(f"not an s3:// URI: {s3_uri}")
    if not Path(local_dir).is_dir():
        raise ValueError(f"not a directory: {local_dir}")
    if idempotent:
        _aws_sync(local_dir, s3_uri)
    else:
        _aws_cp(local_dir, s3_uri, recursive=True)
