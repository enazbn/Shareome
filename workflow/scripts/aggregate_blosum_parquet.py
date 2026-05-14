#!/usr/bin/env python3

"""
aggregate_blosum_parquet.py

Aggregate BLOSUM-based k-mer match results into a peptide-level core table.

Input:
    Parquet file produced by match_kmers_blosum.py / csv_to_parquet.py

Expected input columns:
    virus_kmer
    human_kmer
    blosum_raw_score
    blosum_self_virus
    blosum_self_human
    blosum_similarity_fraction
    blosum_similarity_percent
    mismatches
    identity_fraction
    identity_percent
    human_accession
    human_position
    virus_accession
    virus_position

Output:
    Peptide-level parquet file, grouped by virus_kmer.

Usage:
    python workflow/scripts/aggregate_blosum_parquet.py \
        results/Rota_seq_nr_clean_blosum_matched.parquet \
        results/Rota_seq_nr_clean_blosum_peptide_core.parquet
"""

import sys
from pathlib import Path

import pandas as pd


REQUIRED_COLUMNS = {
    "virus_kmer",
    "human_kmer",
    "blosum_raw_score",
    "blosum_self_virus",
    "blosum_self_human",
    "blosum_similarity_fraction",
    "blosum_similarity_percent",
    "mismatches",
    "identity_fraction",
    "identity_percent",
    "human_accession",
    "human_position",
    "virus_accession",
    "virus_position",
}


def collapse_locations(accessions: pd.Series, positions: pd.Series) -> str:
    """
    Collapse accession and position columns into a semicolon-separated string.

    Example:
        XP_001:45;XP_002:102
    """
    pairs = (
        pd.DataFrame(
            {
                "accession": accessions.astype(str),
                "position": positions.astype(str),
            }
        )
        .drop_duplicates()
        .sort_values(["accession", "position"])
    )

    return ";".join(
        f"{row.accession}:{row.position}"
        for row in pairs.itertuples(index=False)
    )


def collapse_unique_values(values: pd.Series) -> str:
    """
    Collapse unique non-null values into a semicolon-separated string.
    """
    unique_values = (
        values.dropna()
        .astype(str)
        .drop_duplicates()
        .sort_values()
        .tolist()
    )

    return ";".join(unique_values)


def validate_columns(df: pd.DataFrame) -> None:
    """
    Check whether all required columns are present.
    """
    missing = REQUIRED_COLUMNS - set(df.columns)

    if missing:
        missing_str = ", ".join(sorted(missing))
        raise ValueError(
            f"Input file is missing required columns: {missing_str}"
        )


def aggregate_blosum_matches(df: pd.DataFrame) -> pd.DataFrame:
    """
    Aggregate BLOSUM match results at the viral k-mer level.
    """

    validate_columns(df)

    # Make sure numeric columns are actually numeric.
    numeric_cols = [
        "blosum_raw_score",
        "blosum_self_virus",
        "blosum_self_human",
        "blosum_similarity_fraction",
        "blosum_similarity_percent",
        "mismatches",
        "identity_fraction",
        "identity_percent",
        "human_position",
        "virus_position",
    ]

    for col in numeric_cols:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    grouped_rows = []

    for virus_kmer, g in df.groupby("virus_kmer", sort=False):
        human_locs = collapse_locations(
            g["human_accession"],
            g["human_position"],
        )

        virus_locs = collapse_locations(
            g["virus_accession"],
            g["virus_position"],
        )

        row = {
            "virus_kmer": virus_kmer,

            # Counts
            "total_matches": len(g),
            "human_hits": len(g[["human_accession", "human_position"]].drop_duplicates()),
            "virus_hits": len(g[["virus_accession", "virus_position"]].drop_duplicates()),
            "human_proteins": g["human_accession"].nunique(dropna=True),
            "virus_proteins": g["virus_accession"].nunique(dropna=True),
            "matched_human_kmers": g["human_kmer"].nunique(dropna=True),

            # BLOSUM score summaries
            "max_blosum_raw_score": g["blosum_raw_score"].max(),
            "mean_blosum_raw_score": g["blosum_raw_score"].mean(),
            "median_blosum_raw_score": g["blosum_raw_score"].median(),

            "max_blosum_similarity_percent": g["blosum_similarity_percent"].max(),
            "mean_blosum_similarity_percent": g["blosum_similarity_percent"].mean(),
            "median_blosum_similarity_percent": g["blosum_similarity_percent"].median(),

            "max_blosum_similarity_fraction": g["blosum_similarity_fraction"].max(),
            "mean_blosum_similarity_fraction": g["blosum_similarity_fraction"].mean(),

            # Identity / mismatch summaries retained for interpretation
            "min_mismatches": g["mismatches"].min(),
            "mean_mismatches": g["mismatches"].mean(),
            "max_identity_percent": g["identity_percent"].max(),
            "mean_identity_percent": g["identity_percent"].mean(),

            # Matched human sequence diversity
            "human_kmers": collapse_unique_values(g["human_kmer"]),

            # Location strings
            "human_locs": human_locs,
            "virus_locs": virus_locs,
        }

        grouped_rows.append(row)

    out = pd.DataFrame(grouped_rows)

    # A useful default sorting:
    # strongest BLOSUM similarity first, then more viral recurrence, then more human overlap.
    if not out.empty:
        out = out.sort_values(
            [
                "max_blosum_similarity_percent",
                "virus_proteins",
                "human_proteins",
                "total_matches",
            ],
            ascending=[False, False, False, False],
        ).reset_index(drop=True)

    return out


def main() -> None:
    if len(sys.argv) != 3:
        sys.stderr.write(
            "Usage:\n"
            "  python aggregate_blosum_parquet.py "
            "<input_blosum_matched.parquet> <output_blosum_peptide_core.parquet>\n"
        )
        sys.exit(1)

    input_path = Path(sys.argv[1])
    output_path = Path(sys.argv[2])

    if not input_path.exists():
        raise FileNotFoundError(f"Input file does not exist: {input_path}")

    output_path.parent.mkdir(parents=True, exist_ok=True)

    print(f"[INFO] Reading input parquet: {input_path}")
    df = pd.read_parquet(input_path)

    print(f"[INFO] Input rows: {len(df):,}")
    print(f"[INFO] Input columns: {list(df.columns)}")

    if df.empty:
        print("[WARNING] Input file is empty. Writing empty output table.")
        empty_out = pd.DataFrame(
            columns=[
                "virus_kmer",
                "total_matches",
                "human_hits",
                "virus_hits",
                "human_proteins",
                "virus_proteins",
                "matched_human_kmers",
                "max_blosum_raw_score",
                "mean_blosum_raw_score",
                "median_blosum_raw_score",
                "max_blosum_similarity_percent",
                "mean_blosum_similarity_percent",
                "median_blosum_similarity_percent",
                "max_blosum_similarity_fraction",
                "mean_blosum_similarity_fraction",
                "min_mismatches",
                "mean_mismatches",
                "max_identity_percent",
                "mean_identity_percent",
                "human_kmers",
                "human_locs",
                "virus_locs",
            ]
        )
        empty_out.to_parquet(output_path, index=False)
        print(f"[INFO] Wrote empty output parquet: {output_path}")
        return

    print("[INFO] Aggregating BLOSUM matches by virus_kmer...")
    out = aggregate_blosum_matches(df)

    print(f"[INFO] Output rows: {len(out):,}")
    print(f"[INFO] Writing output parquet: {output_path}")
    out.to_parquet(output_path, index=False)

    print("[INFO] Done.")


if __name__ == "__main__":
    main()