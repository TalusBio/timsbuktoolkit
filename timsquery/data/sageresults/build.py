# /// script
# dependencies = [
#   "polars",
# ]
# ///

import polars as pl
import json

PROTON_MASS = 1.007276
NEUTRON_MASS = 1.008664

ub_peptides_df = (
    pl.scan_csv("./results.sage.tsv", separator="\t")
    .filter(pl.col("proteins").str.contains("UBB_"), pl.col("peptide_q") < 0.01)
    .group_by("peptide")
    .agg(pl.all().sort_by("sage_discriminant_score", descending=True).head(1))
    .explode(pl.all().exclude("peptide"))
    .select(
        [
            "peptide",
            "charge",
            "psm_id",
            "rt",
            "ion_mobility",
            "calcmass",
        ]
    )
    .with_columns(
        mz=(pl.col("calcmass") + (pl.col("charge") * PROTON_MASS)) / pl.col("charge")
    )
    .collect()
)

psm_ids = ub_peptides_df["psm_id"].unique()

fragments = (
    pl.scan_csv("matched_fragments.sage.tsv", separator="\t")
    .filter(pl.col("psm_id").is_in(psm_ids))
    .select(
        [
            "psm_id",
            "fragment_charge",
            "fragment_mz_calculated",
            "fragment_intensity",
        ]
    )
    .group_by(["psm_id"])
    .agg(pl.all())
    .collect()
)

df = ub_peptides_df.join(fragments, on="psm_id", how="inner")

# Convert to json
out = []
for x in df.iter_rows(named=True):
    mzs = [x["mz"] + (z * (NEUTRON_MASS / x["charge"])) for z in range(4)]
    expected_ms1_intensities = (0.491, 0.334, 0.129, 0.036)
    expected_intensities = {str(i): y for i, y in enumerate(x["fragment_intensity"])}

    out.append(
        {
            "id": x["psm_id"],
            "human_id": x["peptide"] + "_" + str(x["charge"]),
            "mobility": x["ion_mobility"],
            # Note: the rt in sage is minutes ...
            "rt_seconds": x["rt"] * 60,
            "precursor_mzs": mzs,
            "precursor_charge": x["charge"],
            "fragment_mzs": {
                str(i): y for i, y in enumerate(x["fragment_mz_calculated"])
            },
            "expected_fragment_intensity": expected_intensities,
            "expected_precursor_intensity": expected_ms1_intensities,
        }
    )


out.sort(key=lambda x: x["rt_seconds"])
print(out)

with open("ubb_elution_groups.json", "w") as f:
    json.dump(out, f, indent=2)
