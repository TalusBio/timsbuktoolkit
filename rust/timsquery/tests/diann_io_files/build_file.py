import polars as pl
from pathlib import Path
from polars import Schema, String, Int64, Float32
from typing import Any


schema_use = Schema({
    "Precursor.Id": String,
    "Modified.Sequence": String,
    "Stripped.Sequence": String,
    "Precursor.Charge": Int64,
    "Proteotypic": Int64,
    "Decoy": Int64,
    "N.Term": Int64,
    "C.Term": Int64,
    "RT": Float32,
    "IM": Float32,
    "Q.Value": Float32,
    "Peptidoform.Q.Value": Float32,
    "PTM.Site.Confidence": Float32,
    "PG.Q.Value": Float32,
    "Precursor.Mz": Float32,
    "Product.Mz": Float32,
    "Relative.Intensity": Float32,
    "Fragment.Type": String,
    "Fragment.Charge": Int64,
    "Fragment.Series.Number": Int64,
    "Fragment.Loss.Type": String,
    "Exclude.From.Quant": Int64,
    "Protein.Ids": String,
    "Protein.Group": String,
    "Protein.Names": String,
    "Genes": String,
    "Flags": Int64,
})


def decompress_rle(data: dict[str, list[tuple[int, Any]]]):
    decompressed = {}
    for key, rle_list in data.items():
        decompressed_list = []
        for value, count in rle_list:
            decompressed_list.extend([value] * count)
        decompressed[key] = decompressed_list
    return decompressed


def compress_rle(data: dict[str, list[Any]]) -> dict[str, list[tuple[Any, int]]]:
    compressed = {}
    for key, data_list in data.items():
        if not data_list:
            compressed[key] = []
            continue
        rle_list = []
        prev_value = data_list[0]
        count = 1
        for value in data_list[1:]:
            if value == prev_value:
                count += 1
            else:
                rle_list.append((prev_value, count))
                prev_value = value
                count = 1
        rle_list.append((prev_value, count))
        compressed[key] = rle_list
    return compressed


NAME_USE = "sample_pq_speclib.parquet"
SAVE_LOC = Path(__file__).parent / NAME_USE
COMPRESSED = {
    "Precursor.Id": [
        ("GREEWESAALQNANTK3", 4),
        ("SGGLGGSHALLLLR2", 4),
        ("ILQGAPEILDR2", 4),
    ],
    "Modified.Sequence": [
        ("GREEWESAALQNANTK", 4),
        ("SGGLGGSHALLLLR", 4),
        ("ILQGAPEILDR", 4),
    ],
    "Stripped.Sequence": [
        ("GREEWESAALQNANTK", 4),
        ("SGGLGGSHALLLLR", 4),
        ("ILQGAPEILDR", 4),
    ],
    "Precursor.Charge": [(3, 4), (2, 8)],
    "Proteotypic": [(1, 12)],
    "Decoy": [(0, 12)],
    "N.Term": [(0, 12)],
    "C.Term": [(0, 12)],
    "RT": [(44.62127685546875, 4), (58.92487716674805, 4), (56.13072967529297, 4)],
    "IM": [(0.9339450001716614, 4), (1.066414475440979, 4), (0.9522402286529541, 4)],
    "Q.Value": [
        (0.0011079016840085387, 4),
        (0.001989169977605343, 4),
        (0.0008191685192286968, 4),
    ],
    "Peptidoform.Q.Value": [
        (0.0011079016840085387, 4),
        (0.001989169977605343, 4),
        (0.0008191685192286968, 4),
    ],
    "PTM.Site.Confidence": [(1.0, 12)],
    "PG.Q.Value": [(0.00027824152493849397, 12)],
    "Precursor.Mz": [
        (601.9588623046875, 4),
        (675.896240234375, 4),
        (612.8509521484375, 4),
    ],
    "Product.Mz": [
        (675.342041015625, 1),
        (788.4260864257812, 1),
        (547.283447265625, 1),
        (658.2943725585938, 1),
        (698.4923095703125, 1),
        (1036.626220703125, 1),
        (627.4552001953125, 1),
        (835.5512084960938, 1),
        (742.4093627929688, 1),
        (870.4679565429688, 1),
        (998.5265502929688, 1),
        (499.76690673828125, 1),
    ],
    "Relative.Intensity": [
        (1.0, 1),
        (0.7401171326637268, 1),
        (0.6205617189407349, 1),
        (0.5351390838623047, 1),
        (1.0, 1),
        (0.9322429895401001, 1),
        (0.579563558101654, 1),
        (0.5688385367393494, 1),
        (1.0, 1),
        (0.7854613065719604, 1),
        (0.32938897609710693, 1),
        (0.12993720173835754, 1),
    ],
    "Fragment.Type": [("y", 3), ("b", 1), ("y", 8)],
    "Fragment.Charge": [(1, 11), (2, 1)],
    "Fragment.Series.Number": [
        (6, 1),
        (7, 1),
        (5, 2),
        (6, 1),
        (10, 1),
        (5, 1),
        (7, 1),
        (6, 1),
        (8, 1),
        (9, 2),
    ],
    "Fragment.Loss.Type": [("unknown", 12)],
    "Exclude.From.Quant": [(0, 12)],
    "Protein.Ids": [("Q5T4S7-3", 4), ("P42704", 4), ("Q86TP1", 4)],
    "Protein.Group": [("Q5T4S7-3", 4), ("P42704", 4), ("Q86TP1", 4)],
    "Protein.Names": [("", 12)],
    "Genes": [("", 12)],
    # WTF does flags mean?
    "Flags": [(17, 1), (1, 3), (17, 1), (1, 3), (17, 1), (1, 3)],
}

if __name__ == "__main__":
    decompressed = decompress_rle(COMPRESSED)
    recompressed = compress_rle(decompressed)
    assert COMPRESSED == recompressed

    df = pl.DataFrame(decompressed, schema=schema_use)
    df.write_parquet(SAVE_LOC)
