import sys
import pandas as pd


def make_loc_strings(df, acc_col, pos_col):
    tmp = df[[acc_col, pos_col]].drop_duplicates().copy()
    tmp["loc"] = tmp[acc_col].astype(str) + ":" + tmp[pos_col].astype(str)
    return tmp["loc"].tolist()


def aggregate_matched_parquet(input_parquet, output_parquet):
    df = pd.read_parquet(input_parquet)

    required = [
        "kmer",
        "human_accession",
        "human_position",
        "virus_accession",
        "virus_position",
    ]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    results = []

    for peptide, sub in df.groupby("kmer", sort=False):
        human_sites = sub[["human_accession", "human_position"]].drop_duplicates()
        virus_sites = sub[["virus_accession", "virus_position"]].drop_duplicates()

        human_hits = len(human_sites)
        virus_hits = len(virus_sites)

        human_proteins = human_sites["human_accession"].nunique()
        virus_proteins = virus_sites["virus_accession"].nunique()

        human_locs = ";".join(
            human_sites["human_accession"].astype(str)
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
                "peptide": peptide,
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


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(
            "Usage: python aggregate_matched_parquet.py <input.parquet> <output.parquet>"
        )
        sys.exit(1)

    input_parquet = sys.argv[1]
    output_parquet = sys.argv[2]

    aggregate_matched_parquet(input_parquet, output_parquet)
