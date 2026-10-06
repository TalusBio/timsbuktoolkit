"""One-fold peptide entrapment experiment for timsseek.

FDRBench's peptide pair TSV is the interchange format. The Java program is
optional and is used only to check the FDP calculations.
"""

from __future__ import annotations

import argparse
import csv
import math
import random
import re
import subprocess
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")
PAIR_COLUMNS = ("sequence", "decoy", "proteins", "peptide_type", "peptide_pair_index")
FDP_COLUMNS = ("peptide", "mod_peptide", "charge", "q_value", "score")
PLOT_COLUMNS = ("q_value", "lower_bound_fdp", "combined_fdp", "paired_fdp")
MODIFICATION = re.compile(r"\[[^\[\]]+\]|\([^()]+\)")


def stripped_sequence(sequence: str) -> str:
    stripped = MODIFICATION.sub("", sequence).replace("_", "").replace("-", "")
    stripped = stripped.upper().replace("I", "L")
    if not stripped or not set(stripped) <= AMINO_ACIDS:
        raise ValueError(f"cannot extract peptide residues from {sequence!r}")
    return stripped


def fasta_sequences(path: Path):
    header = None
    sequence = []
    with path.open() as stream:
        for line in stream:
            line = line.strip()
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(sequence)
                header = line[1:].split()[0]
                sequence = []
            elif line:
                if header is None:
                    raise ValueError(f"{path}: sequence before FASTA header")
                sequence.append(line.upper().replace("I", "L"))
    if header is not None:
        yield header, "".join(sequence)


def tryptic_peptides(protein: str):
    start = 0
    for i, residue in enumerate(protein):
        if residue in "KR" and (i + 1 == len(protein) or protein[i + 1] != "P"):
            yield protein[start : i + 1]
            start = i + 1
    if start < len(protein):
        yield protein[start:]


def generate(
    parent: Path, out: Path, seed: int = 1, min_length: int = 7, max_length: int = 35
):
    """Make one composition-matched entrapment per unique tryptic peptide."""
    if min_length < 3 or max_length < min_length:
        raise ValueError("length bounds must satisfy 3 <= min_length <= max_length")
    out.mkdir(parents=True, exist_ok=True)
    pair_path = out / "peptide_pairs.tsv"
    fasta_path = out / "entrapment.fasta"
    for path in (pair_path, fasta_path):
        if path.exists():
            raise FileExistsError(path)

    sources = defaultdict(set)
    for protein_id, protein in fasta_sequences(parent):
        for peptide in tryptic_peptides(protein):
            if (
                min_length <= len(peptide) <= max_length
                and set(peptide) <= AMINO_ACIDS
                and not any(residue in "KR" for residue in peptide[:-1])
            ):
                sources[peptide].add(protein_id)
    if not sources:
        raise ValueError("no eligible peptides in FASTA")

    rng = random.Random(seed)
    occupied = set(sources)
    pairs = []
    for peptide in sorted(sources):
        middle = list(peptide[1:-1])
        entrapment = None
        for _ in range(100):
            rng.shuffle(middle)
            candidate = peptide[0] + "".join(middle) + peptide[-1]
            if candidate not in occupied:
                entrapment = candidate
                break
        if entrapment is not None:
            occupied.add(entrapment)
            pairs.append((peptide, entrapment))
    if not pairs:
        raise ValueError("no distinct entrapment peptides could be generated")

    with (
        pair_path.open("w", newline="") as pair_stream,
        fasta_path.open("w") as fasta_stream,
    ):
        writer = csv.writer(pair_stream, delimiter="\t", lineterminator="\n")
        writer.writerow(PAIR_COLUMNS)
        for index, (target, entrapment) in enumerate(pairs):
            for kind, sequence in (("target", target), ("p_target", entrapment)):
                protein_id = f"entrap_{index}_{kind}"
                writer.writerow((sequence, "No", protein_id, kind, index))
                fasta_stream.write(f">{protein_id}\n{sequence}\n")
    return {
        "pairs": len(pairs),
        "skipped": len(sources) - len(pairs),
        "fasta": fasta_path,
        "pair_file": pair_path,
    }


@dataclass(frozen=True)
class Hit:
    peptide: str
    mod_peptide: str
    charge: int
    q_value: float
    score: float


def load_pairs(path: Path):
    labels = {}
    groups = defaultdict(dict)
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        if not set(PAIR_COLUMNS) <= set(reader.fieldnames or ()):
            raise ValueError(f"{path}: expected columns {PAIR_COLUMNS}")
        for row in reader:
            kind = row["peptide_type"]
            if kind not in ("target", "p_target"):
                continue
            peptide = row["sequence"]
            if "I" in peptide:
                raise ValueError(f"{path}: pair file must convert I to L: {peptide}")
            if peptide in labels and labels[peptide] != kind:
                raise ValueError(f"{path}: peptide has conflicting labels: {peptide}")
            labels[peptide] = kind
            group = groups[row["peptide_pair_index"]]
            if kind in group and group[kind] != peptide:
                raise ValueError(f"{path}: multiple {kind} peptides in pair")
            group[kind] = peptide
    pairing = {}
    paired_targets = set()
    for index, group in groups.items():
        if set(group) != {"target", "p_target"}:
            raise ValueError(f"{path}: incomplete pair {index}")
        if group["target"] in paired_targets:
            raise ValueError(f"{path}: target belongs to multiple pairs")
        paired_targets.add(group["target"])
        pairing[group["p_target"]] = group["target"]
    if not pairing:
        raise ValueError(f"{path}: no target-entrapment pairs")
    return labels, pairing


def read_results(path: Path, labels: dict[str, str], allow_unpaired: bool = False):
    import polars as pl

    required = {
        "sequence",
        "is_target",
        "precursor_charge",
        "qvalue",
        "discriminant_score",
    }
    missing = required - set(pl.read_parquet_schema(path))
    if missing:
        raise ValueError(
            f"{path}: missing columns {sorted(missing)}; cannot evaluate raw scores"
        )
    frame = pl.read_parquet(path, columns=sorted(required))
    hits = {}
    unknown = set()
    for row in frame.select(sorted(required)).iter_rows(named=True):
        if row["is_target"] is not True:
            continue
        mod_peptide = row["sequence"]
        if mod_peptide is None:
            raise ValueError(f"{path}: target result has no sequence")
        mod_peptide = mod_peptide.replace("I", "L")
        peptide = stripped_sequence(mod_peptide)
        if peptide not in labels:
            unknown.add(peptide)
            continue
        q_value = float(row["qvalue"])
        score = float(row["discriminant_score"])
        if not 0 <= q_value <= 1 or not math.isfinite(score):
            raise ValueError(
                f"{path}: non-finite score or invalid q-value for {mod_peptide}"
            )
        hit = Hit(peptide, mod_peptide, int(row["precursor_charge"]), q_value, score)
        key = (peptide, hit.charge)
        old = hits.get(key)
        if old is None or (hit.q_value, -hit.score) < (old.q_value, -old.score):
            hits[key] = hit
    if unknown and not allow_unpaired:
        raise ValueError(
            f"{path}: {len(unknown)} peptides absent from pair file: "
            f"{sorted(unknown)[:5]}"
        )
    if unknown:
        sys.stderr.write(f"Dropped {len(unknown)} peptides absent from the pair file\n")
    if not hits:
        raise ValueError(f"{path}: no target results")
    return list(hits.values())


def estimate(hits: list[Hit], labels: dict[str, str], pairing: dict[str, str]):
    """FDRBench one-fold lower, combined, and paired FDP at each q-value cutoff."""
    ordered = sorted(hits, key=lambda h: (h.q_value, -h.score, h.peptide, h.charge))
    target_rank = {
        (hit.peptide, hit.charge): rank
        for rank, hit in enumerate(ordered)
        if labels[hit.peptide] == "target"
    }
    target_to_entrapment = {
        target: entrapment for entrapment, target in pairing.items()
    }
    seen_entrapments = set()
    targets = entrapments = orphaned = outranked = 0
    values = []
    for rank, hit in enumerate(ordered):
        if labels[hit.peptide] == "target":
            targets += 1
            # An earlier entrapment moves from p>s>t to p>t>s.
            partner = target_to_entrapment.get(hit.peptide)
            if (partner, hit.charge) in seen_entrapments:
                orphaned -= 1
                outranked += 1
        else:
            entrapments += 1
            seen_entrapments.add((hit.peptide, hit.charge))
            t_rank = target_rank.get((pairing[hit.peptide], hit.charge))
            if t_rank is None or rank < t_rank:
                # Count as p>s>t until the target enters the cutoff.
                orphaned += 1
        total = targets + entrapments
        values.append((
            entrapments / total,
            2 * entrapments / total,
            (entrapments + orphaned + 2 * outranked) / total,
            targets,
            entrapments,
        ))
    # FDRBench reports one value per q-value threshold, including tied rows.
    start = 0
    while start < len(ordered):
        end = start + 1
        while end < len(ordered) and ordered[end].q_value == ordered[start].q_value:
            end += 1
        values[start:end] = [values[end - 1]] * (end - start)
        start = end
    return [
        {
            **dict(
                zip(
                    FDP_COLUMNS,
                    (h.peptide, h.mod_peptide, h.charge, h.q_value, h.score),
                )
            ),
            "lower_bound_fdp": v[0],
            "combined_fdp": v[1],
            "paired_fdp": v[2],
            "n_targets": v[3],
            "n_entrapments": v[4],
        }
        for h, v in zip(ordered, values, strict=True)
    ]


def write_tsv(path: Path, rows, columns):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(
            stream, fieldnames=columns, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def plot_fdp(rows, out: Path):
    """Plot q-value calibration without materializing the full report again."""
    import matplotlib

    matplotlib.use("Agg")
    from matplotlib import pyplot as plt

    # Keep the last reported cutoff in each q-value bin. The plot needs at most
    # 1,200 points per panel even when the report has millions of precursors.
    windows = [(0.1, {}), (1.0, {})]
    last_q = -1.0
    final_counts = None
    for row in rows:
        values = tuple(float(row[column]) for column in PLOT_COLUMNS)
        q_value = values[0]
        if not all(math.isfinite(value) for value in values) or not 0 <= q_value <= 1:
            raise ValueError(f"invalid q-value or FDP in plot input: {values}")
        if q_value < last_q:
            raise ValueError("entrapment report must be sorted by q-value")
        last_q = q_value
        final_counts = (int(row["n_targets"]), int(row["n_entrapments"]))
        for limit, buckets in windows:
            if q_value <= limit:
                bucket = min(1199, int(q_value / limit * 1200))
                buckets[bucket] = values
    if final_counts is None:
        raise ValueError("entrapment report has no rows to plot")

    fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
    colors = ("#4e79a7", "#e15759", "#59a14f")
    labels = ("Lower bound", "Combined", "Paired")
    for ax, (limit, buckets), title in zip(
        axes, windows, ("Discovery range", "Full range"), strict=True
    ):
        points = [buckets[index] for index in sorted(buckets)]
        ax.plot(
            [0, limit], [0, limit], color="#404040", linestyle="--", label="q = FDP"
        )
        for offset, (color, label) in enumerate(
            zip(colors, labels, strict=True), start=1
        ):
            ax.plot(
                [point[0] for point in points],
                [point[offset] for point in points],
                color=color,
                linewidth=1.8,
                label=label,
            )
        ax.set_xlim(0, limit)
        ax.set_ylim(
            0, max(limit, *(point[i] for point in points for i in (1, 2, 3))) * 1.05
        )
        ax.set_title(f"{title} (q ≤ {limit:g})")
        ax.set_xlabel("Reported q-value")
        ax.set_ylabel("Entrapment FDP estimate")
        ax.grid(alpha=0.2)
    axes[0].legend(loc="upper left", frameon=False)
    fig.suptitle(
        "Entrapment FDP versus reported q-value\n"
        f"{final_counts[0]:,} target and {final_counts[1]:,} entrapment precursors"
    )
    out.mkdir(parents=True, exist_ok=True)
    paths = (out / "entrapment_fdp.png", out / "entrapment_fdp.pdf")
    for path in paths:
        fig.savefig(path, dpi=200)
    plt.close(fig)
    return paths


def plot_report(report: Path, out: Path | None = None):
    """Plot an existing `entrapment_fdp.tsv` without rerunning the search."""
    with report.open(newline="") as stream:
        rows = csv.DictReader(stream, delimiter="\t")
        required = set(PLOT_COLUMNS) | {"n_targets", "n_entrapments"}
        if not required <= set(rows.fieldnames or ()):
            raise ValueError(f"{report}: expected columns {sorted(required)}")
        return plot_fdp(rows, out or report.parent)


def analyze(
    results: Path,
    pairs: Path,
    out: Path,
    reference_jar: Path | None = None,
    allow_unpaired: bool = False,
):
    labels, pairing = load_pairs(pairs)
    hits = read_results(results, labels, allow_unpaired)
    report = estimate(hits, labels, pairing)
    out.mkdir(parents=True, exist_ok=True)
    input_file = out / "fdrbench_input.tsv"
    report_file = out / "entrapment_fdp.tsv"
    write_tsv(
        input_file,
        [
            dict(
                zip(
                    FDP_COLUMNS,
                    (h.peptide, h.mod_peptide, h.charge, h.q_value, h.score),
                )
            )
            for h in hits
        ],
        FDP_COLUMNS,
    )
    write_tsv(report_file, report, tuple(report[0]))
    plot_fdp(report, out)
    if reference_jar is not None:
        reference_file = out / "fdrbench_reference.csv"
        subprocess.run(
            [
                "java",
                "-jar",
                str(reference_jar),
                "-i",
                str(input_file),
                "-level",
                "precursor",
                "-pep",
                str(pairs),
                "-fold",
                "1",
                "-score",
                "score:1",
                "-o",
                str(reference_file),
            ],
            check=True,
        )
        compare_reference(report, reference_file, pairing)
    return report_file


def compare_reference(
    report: list[dict],
    reference_file: Path,
    pairing: dict[str, str],
    tolerance: float = 1e-6,
):
    ours = {(r["peptide"], r["charge"]): r for r in report}
    tied_q_values = []
    for row in report:
        target = pairing.get(row["peptide"])
        if target is None:
            continue
        mate = ours.get((target, row["charge"]))
        if mate and (row["q_value"], row["score"]) == (
            mate["q_value"],
            mate["score"],
        ):
            tied_q_values.append(row["q_value"])
    if tied_q_values:
        sys.stderr.write(
            f"Reference comparison allows sort-order differences for "
            f"{len(tied_q_values)} exactly tied pairs\n"
        )
    with reference_file.open(newline="") as stream:
        reference = list(csv.DictReader(stream))
    if len(reference) != len(ours):
        raise ValueError("FDRBench returned a different number of precursors")
    for row in reference:
        key = (row["peptide"].replace("I", "L"), int(row["charge"]))
        if key not in ours:
            raise ValueError(f"FDRBench returned unknown precursor {key}")
        for field in ("lower_bound_fdp", "combined_fdp", "paired_fdp"):
            limit = tolerance
            if field == "paired_fdp":
                item = ours[key]
                tied = sum(q <= item["q_value"] for q in tied_q_values)
                limit += 2 * tied / (item["n_targets"] + item["n_entrapments"])
            if abs(float(row[field]) - ours[key][field]) > limit:
                raise ValueError(
                    f"FDRBench {field} differs for {key}: "
                    f"{row[field]} vs {ours[key][field]}"
                )


def run_search(
    timsseek: Path,
    fasta: Path,
    raw: list[Path],
    out: Path,
    library: Path,
    config: Path | None,
    reuse_library: bool = False,
):
    if not str(library).endswith((".mzspeclib.txt", ".mzspeclib.txt.gz")):
        raise ValueError("library path must end in .mzspeclib.txt[.gz]")
    if library.exists() and not reuse_library:
        raise FileExistsError(f"{library} exists; pass --reuse-library to search it")
    if not library.exists() and reuse_library:
        raise FileNotFoundError(library)
    if not reuse_library:
        lengths = [len(sequence) for _, sequence in fasta_sequences(fasta)]
        if not lengths:
            raise ValueError(f"{fasta}: no FASTA sequences")
        command = [
            str(timsseek),
            "build-library",
            "--fasta",
            str(fasta),
            "--out",
            str(library),
            "--missed-cleavages",
            "0",
            "--min-length",
            str(min(lengths)),
            "--max-length",
            str(max(lengths)),
            "--decoys",
        ]
        if config is not None:
            command.extend(("--config", str(config)))
        subprocess.run(command, check=True)
    command = [
        str(timsseek),
        "--speclib-uri",
        str(library),
        "--output-uri",
        str(out),
        "--max-qvalue",
        "1",
        "--decoy-strategy",
        "never",
    ]
    for path in raw:
        command.extend(("--raw-inputs", str(path)))
    if config is not None:
        command.extend(("--config", str(config)))
    subprocess.run(command, check=True)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    gen = commands.add_parser(
        "generate", help="Make peptide pairs and a FASTA from a parent proteome"
    )
    gen.add_argument("--fasta", type=Path, required=True)
    gen.add_argument("--out", type=Path, required=True)
    gen.add_argument("--seed", type=int, default=1)
    gen.add_argument("--min-length", type=int, default=7)
    gen.add_argument("--max-length", type=int, default=35)
    run = commands.add_parser("run", help="Build a library and search raw files")
    run.add_argument("--fasta", type=Path, required=True)
    run.add_argument("--raw", type=Path, action="append", required=True)
    run.add_argument("--out", type=Path, required=True)
    run.add_argument("--library", type=Path, required=True)
    run.add_argument("--timsseek", type=Path, default=Path("target/release/timsseek"))
    run.add_argument("--config", type=Path)
    run.add_argument("--reuse-library", action="store_true")
    calc = commands.add_parser("analyze", help="Calculate FDP from one results.parquet")
    calc.add_argument("--results", type=Path, required=True)
    calc.add_argument("--pairs", type=Path, required=True)
    calc.add_argument("--out", type=Path, required=True)
    calc.add_argument("--reference-jar", type=Path)
    calc.add_argument("--allow-unpaired", action="store_true")
    plot = commands.add_parser("plot", help="Plot an existing entrapment FDP report")
    plot.add_argument("--report", type=Path, required=True)
    plot.add_argument(
        "--out", type=Path, help="Output directory (default: report directory)"
    )
    args = parser.parse_args(argv)
    if args.command == "generate":
        sys.stdout.write(
            "Entrapment generation: trypsin (K/R unless followed by P), "
            f"missed cleavages=0, length={args.min_length}..{args.max_length}, "
            f"I→L, modifications=none, seed={args.seed}\n"
        )
        sys.stdout.write(
            "Library prediction defaults for `run`: fixed C[UNIMOD:4], "
            "variable M[UNIMOD:35], max variable mods=1 "
            "(--config can override these)\n"
        )
        sys.stdout.flush()
        result = generate(
            args.fasta, args.out, args.seed, args.min_length, args.max_length
        )
        sys.stdout.write(f"{result}\n")
    elif args.command == "run":
        run_search(
            args.timsseek,
            args.fasta,
            args.raw,
            args.out,
            args.library,
            args.config,
            args.reuse_library,
        )
    elif args.command == "analyze":
        result = analyze(
            args.results,
            args.pairs,
            args.out,
            args.reference_jar,
            args.allow_unpaired,
        )
        sys.stdout.write(f"{result}\n")
    else:
        png, pdf = plot_report(args.report, args.out)
        sys.stdout.write(f"{png}\n{pdf}\n")


if __name__ == "__main__":
    main()
