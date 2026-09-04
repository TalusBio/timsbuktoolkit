# Agent notes

- Format with `task fmt`. It runs nightly rustfmt plus ruff. Plain `cargo fmt`
  on stable silently drops the nightly-only options in `rustfmt.toml` and
  reports hundreds of spurious diffs. CI checks with nightly.
- Everything else: `docs/development.md`.
