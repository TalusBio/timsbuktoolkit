"""Tiny shellout wrapper around `aws s3 cp`. One module so all subprocess
invocations of the AWS CLI live in one place (easier to mock in tests)."""

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


def s3_download_file(s3_uri: str, local_path: str) -> None:
    if not s3_uri.startswith("s3://"):
        raise ValueError(f"not an s3:// URI: {s3_uri}")
    _aws_cp(s3_uri, local_path)


def s3_upload_file(local_path: str, s3_uri: str) -> None:
    if not s3_uri.startswith("s3://"):
        raise ValueError(f"not an s3:// URI: {s3_uri}")
    _aws_cp(local_path, s3_uri)


def s3_upload_dir(local_dir: str, s3_uri: str) -> None:
    if not s3_uri.startswith("s3://"):
        raise ValueError(f"not an s3:// URI: {s3_uri}")
    if not Path(local_dir).is_dir():
        raise ValueError(f"not a directory: {local_dir}")
    _aws_cp(local_dir, s3_uri, recursive=True)
