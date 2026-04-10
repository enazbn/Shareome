import sys
import csv
from collections import defaultdict


def load_human_kmers(human_csv):
    human_dict = defaultdict(list)

    with open(human_csv, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            kmer = row["kmer"]
            accession = row["accession"]
            position = row["position"]
            human_dict[kmer].append((accession, position))

    return human_dict


def match_kmers(human_csv, virus_csv, output_csv):
    human_dict = load_human_kmers(human_csv)

    with open(virus_csv, newline="") as vf, open(output_csv, "w", newline="") as out:
        reader = csv.DictReader(vf)
        writer = csv.writer(out)

        writer.writerow([
            "kmer",
            "human_accession",
            "human_position",
            "virus_accession",
            "virus_position"
        ])

        for row in reader:
            kmer = row["kmer"]
            virus_accession = row["accession"]
            virus_position = row["position"]

            if kmer in human_dict:
                for human_accession, human_position in human_dict[kmer]:
                    writer.writerow([
                        kmer,
                        human_accession,
                        human_position,
                        virus_accession,
                        virus_position
                    ])


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python match_kmers.py <human_kmers.csv> <virus_kmers.csv> <output.csv>")
        sys.exit(1)

    human_csv = sys.argv[1]
    virus_csv = sys.argv[2]
    output_csv = sys.argv[3]

    match_kmers(human_csv, virus_csv, output_csv)
