# /// script
# dependencies = [
#   "altair",
#   "polars",
#   "vl-convert-python",
# ]
# ///

from pathlib import Path
import polars as pl
import altair as alt
import argparse


parser = argparse.ArgumentParser()
parser.add_argument("benchmark_file")
args = parser.parse_args()

bench_file_path = Path(args.benchmark_file)
target_dir = bench_file_path.parent
file_stem = bench_file_path.stem

data = (
    pl.read_json(str(bench_file_path))
    .explode("results")
    .with_columns(
        bench=pl.col("results").struct.field("name"),
        context=pl.col("results").struct.field("context"),
        iterations=pl.col("results").struct.field("iterations"),
    )
    .explode("iterations")
    .with_columns(time_seconds=pl.col("iterations").struct.field("duration_seconds"))
)
print(data)

for context, sdf in data.group_by("context"):
    context = context[0]
    print(context)
    # Sample to max 100 points on each bench
    sdf = sdf.filter(pl.int_range(pl.len()).shuffle().over("bench") < 100)
    alt.Chart(sdf).mark_point().encode(
        x="bench",
        y=alt.Y("time_seconds", scale=alt.Scale(type="log")),
        color="bench",
    ).properties(
        width=480,
        height=360,
        title=f"{context}",
    ).save(target_dir / f"{file_stem}_{context}.png")
