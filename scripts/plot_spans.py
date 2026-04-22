# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "matplotlib",
# ]
# ///
"""Gantt plot of tracing spans from a timsseek `.spans.jsonl` file.

Produced by the JSON `fmt` layer in `init_tracing` with
`FmtSpan::NEW | FmtSpan::CLOSE`. Each CLOSE event carries `span.name`,
fields, and `time.busy` / `time.idle`. We match NEW↔CLOSE by
(thread, span name, nearest unfinished open) and draw one bar per span,
colored by span name, grouped by thread.

Usage:
  uv run scripts/plot_spans.py <spans.jsonl> [-o out.png] [--wandb run-name]
"""

from __future__ import annotations

import argparse
import json
import re
from collections import defaultdict
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt


@dataclass
class Span:
    name: str
    label: str
    thread: str
    start: float  # seconds since run start
    end: float
    busy_s: float
    idle_s: float
    depth: int


_DUR_RE = re.compile(r"^\s*([\d.]+)\s*(ns|us|µs|ms|s|m)\s*$")


def parse_duration(s: str) -> float:
    """Parse tracing duration strings like '1.23ms', '45s', '2m' → seconds."""
    m = _DUR_RE.match(s)
    if not m:
        return 0.0
    val, unit = float(m.group(1)), m.group(2)
    return {
        "ns": 1e-9,
        "us": 1e-6,
        "µs": 1e-6,
        "ms": 1e-3,
        "s": 1.0,
        "m": 60.0,
    }[unit] * val


def parse_ts(s: str) -> float:
    # tracing JSON layer emits RFC3339 UTC with `Z`
    return datetime.fromisoformat(s.replace("Z", "+00:00")).timestamp()


def load_spans(path: Path) -> list[Span]:
    opens: dict[tuple[str, str], list[tuple[float, int]]] = defaultdict(list)
    depth_by_thread: dict[str, int] = defaultdict(int)
    closed: list[Span] = []
    t0: float | None = None

    for line in path.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        try:
            ev = json.loads(line)
        except json.JSONDecodeError:
            continue

        msg = ev.get("fields", {}).get("message", "")
        span = ev.get("span") or {}
        name = span.get("name") or ev.get("target", "")
        if not name:
            continue
        thread = ev.get("threadName") or ev.get("thread_name") or "main"
        ts = parse_ts(ev["timestamp"])
        if t0 is None:
            t0 = ts

        key = (thread, name)
        if msg == "new":
            depth = depth_by_thread[thread]
            depth_by_thread[thread] += 1
            opens[key].append((ts, depth))
        elif msg == "close":
            if not opens[key]:
                continue
            start, depth = opens[key].pop()
            depth_by_thread[thread] = max(0, depth_by_thread[thread] - 1)
            fields = ev.get("fields", {})
            busy = parse_duration(str(fields.get("time.busy", "0ns")))
            idle = parse_duration(str(fields.get("time.idle", "0ns")))
            label_fields = {k: v for k, v in span.items() if k != "name"}
            # `step{label=...}` and `file{name=...}` are the dominant spans
            # — show only the field value, not the wrapper `step[label=x]`.
            if name == "step" and "label" in label_fields:
                label = str(label_fields["label"])
            elif name == "file" and "name" in label_fields:
                label = str(label_fields["name"])
            elif label_fields:
                label = (
                    f"{name}["
                    + ",".join(f"{k}={v}" for k, v in label_fields.items())
                    + "]"
                )
            else:
                label = name
            closed.append(
                Span(
                    name=name,
                    label=label,
                    thread=thread,
                    start=start - t0,
                    end=ts - t0,
                    busy_s=busy,
                    idle_s=idle,
                    depth=depth,
                )
            )
    return closed


def plot(spans: list[Span], out: Path, min_ms: float) -> None:
    if not spans:
        raise SystemExit("no spans parsed")

    total_end = max(s.end for s in spans)
    kept = [s for s in spans if (s.end - s.start) * 1000.0 >= min_ms]
    if not kept:
        raise SystemExit(f"no spans >= {min_ms}ms; lower --min-ms")

    # One lane per (thread, depth) — avoids vertical squish when depth is high.
    lane_keys = sorted({(s.thread, s.depth) for s in kept})
    lane_idx = {k: i for i, k in enumerate(lane_keys)}

    # Color by first-seen order for stable legend.
    names: list[str] = []
    seen: set[str] = set()
    for s in sorted(kept, key=lambda x: x.start):
        if s.name not in seen:
            names.append(s.name)
            seen.add(s.name)
    cmap = plt.get_cmap("tab20")
    color_for = {n: cmap(i % 20) for i, n in enumerate(names)}

    # Height scales with lanes so no squish. Extra headroom for title + stats.
    fig_h = max(4.5, 0.45 * len(lane_keys) + 2.5)
    fig, ax = plt.subplots(figsize=(16, fig_h))

    bar_h = 0.78
    for s in kept:
        y = lane_idx[(s.thread, s.depth)]
        width = max(s.end - s.start, 1e-4)
        ax.barh(
            y,
            width,
            left=s.start,
            height=bar_h,
            color=color_for[s.name],
            edgecolor="black",
            linewidth=0.4,
            alpha=0.85,
        )
        if width / total_end > 0.025:
            busy_pct = 100.0 * s.busy_s / max(width, 1e-9)
            ax.text(
                s.start + width / 2,
                y,
                f"{s.label}\n{width:.2f}s · busy {busy_pct:.0f}%",
                ha="center",
                va="center",
                fontsize=7,
                clip_on=True,
            )

    ax.set_yticks(range(len(lane_keys)))
    ax.set_yticklabels([f"{t} · d{d}" for t, d in lane_keys])
    ax.set_xlabel("time since run start (s)")
    ax.set_ylabel("thread · depth")
    ax.set_title(
        f"tracing spans · {len(kept)}/{len(spans)} shown "
        f"(>= {min_ms:g}ms) · {total_end:.1f}s total"
    )
    ax.set_xlim(left=-0.02 * total_end, right=total_end * 1.02)
    ax.invert_yaxis()
    ax.grid(axis="x", linestyle=":", alpha=0.4)

    # Per-name totals panel as stats footer.
    name_totals: dict[str, float] = defaultdict(float)
    name_counts: dict[str, int] = defaultdict(int)
    for s in kept:
        name_totals[s.name] += s.end - s.start
        name_counts[s.name] += 1
    top = sorted(name_totals.items(), key=lambda kv: -kv[1])[:8]
    stats = " · ".join(f"{n}:{t:.1f}s (x{name_counts[n]})" for n, t in top)
    fig.text(0.01, 0.01, stats, fontsize=7, color="#444")

    handles = [plt.Rectangle((0, 0), 1, 1, color=color_for[n], label=n) for n in names]
    ax.legend(
        handles=handles,
        loc="upper left",
        bbox_to_anchor=(1.01, 1.0),
        fontsize=8,
        frameon=False,
    )

    fig.tight_layout()
    fig.savefig(out, dpi=140, bbox_inches="tight")
    print(f"wrote {out}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("spans_jsonl", type=Path)
    ap.add_argument("-o", "--out", type=Path, default=None)
    ap.add_argument(
        "--min-ms",
        type=float,
        default=50.0,
        help="hide spans shorter than this (default 50ms) — filters out bucket-sort noise",
    )
    ap.add_argument("--wandb", metavar="RUN_NAME", default=None)
    ap.add_argument("--wandb-project", default="timsseek-spans")
    args = ap.parse_args()

    out = args.out or args.spans_jsonl.with_suffix(".gantt.png")
    spans = load_spans(args.spans_jsonl)
    plot(spans, out, min_ms=args.min_ms)

    if args.wandb:
        import wandb  # noqa: PLC0415

        run = wandb.init(project=args.wandb_project, name=args.wandb, reinit=True)  # ty: ignore[unresolved-attribute]
        run.log({"spans_gantt": wandb.Image(str(out))})  # ty: ignore[unresolved-attribute]
        run.finish()


if __name__ == "__main__":
    main()
