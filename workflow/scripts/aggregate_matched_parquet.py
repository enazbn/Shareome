import sys
import pandas as pd


def aggregate_matched_parquet(input_parquet, output_parquet):
    df = pd.read_parquet(input_parquet)

    required = [
        "virus_kmer",
        "human_kmer",
        "mismatches",
        "identity_percent",
        "human_accession",
        "human_position",
        "virus_accession",
        "virus_position",
    ]

    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    results = []

    group_cols = ["virus_kmer", "mismatches"]

    for (virus_peptide, mismatches), sub in df.groupby(group_cols, sort=False):
        human_sites = sub[["human_kmer", "human_accession", "human_position"]].drop_duplicates()
        virus_sites = sub[["virus_accession", "virus_position"]].drop_duplicates()

        human_hits = len(human_sites)
        virus_hits = len(virus_sites)

        human_proteins = human_sites["human_accession"].nunique()
        virus_proteins = virus_sites["virus_accession"].nunique()

        human_kmers = ";".join(sorted(sub["human_kmer"].dropna().astype(str).unique()))

        human_locs = ";".join(
            human_sites["human_kmer"].astype(str)
            + "|"
            + human_sites["human_accession"].astype(str)
            + ":"
            + human_sites["human_position"].astype(str)
        )

        virus_locs = ";".join(
            virus_sites["virus_accession"].astype(str)
            + ":"
            + virus_sites["virus_position"].astype(str)
        )

        results.append(
            {
                "virus_peptide": virus_peptide,
                "mismatches": mismatches,
                "identity_percent": sub["identity_percent"].max(),
                "matched_human_kmers": human_kmers,
                "n_matched_human_kmers": sub["human_kmer"].nunique(),
                "human_hits": human_hits,
                "virus_hits": virus_hits,
                "human_proteins": human_proteins,
                "virus_proteins": virus_proteins,
                "human_locs": human_locs,
                "virus_locs": virus_locs,
            }
        )

    out_df = pd.DataFrame(results)
    out_df.to_parquet(output_parquet, index=False)

    print(f"Wrote {len(out_df)} aggregated virus peptide groups to {output_parquet}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(
            "Usage: python aggregate_matched_parquet.py <input.parquet> <output.parquet>"
        )
        sys.exit(1)

    aggregate_matched_parquet(sys.argv[1], sys.argv[2])