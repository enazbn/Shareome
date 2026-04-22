import sys
import os
import csv
from collections import defaultdict
from itertools import combinations


def generate_output_filename(human_csv, virus_csv, max_mismatches):
    human_name = os.path.splitext(os.path.basename(human_csv))[0]
    virus_name = os.path.splitext(os.path.basename(virus_csv))[0]
    return f"{human_name}_vs_{virus_name}_mm{max_mismatches}.csv"


def generate_masked_patterns(kmer, num_masks):
    """
    Generate all masked versions of a kmer by replacing num_masks positions with '*'.
    """
    kmer_list = list(kmer)
    for positions in combinations(range(len(kmer)), num_masks):
        masked = kmer_list.copy()
        for pos in positions:
            masked[pos] = "*"
        yield "".join(masked)


def count_mismatches(seq1, seq2):
    """
    Count position-wise mismatches between two equal-length strings.
    """
    return sum(a != b for a, b in zip(seq1, seq2))


def load_human_kmers_index(human_csv, max_mismatches):
    """
    Build a masked-pattern index for human kmers.

    Returns
    -------
    index[length][pattern] -> list of (kmer, accession, position)
    """
    index = defaultdict(lambda: defaultdict(list))

    with open(human_csv, newline="") as f:
        reader = csv.DictReader(f)
        required = {"kmer", "accession", "position"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Missing required columns in {human_csv}: {sorted(missing)}")

        for row in reader:
            kmer = row["kmer"].strip()
            accession = row["accession"].strip()
            position = row["position"].strip()
            k = len(kmer)

            if max_mismatches == 0:
                index[k][kmer].append((kmer, accession, position))
            else:
                for pattern in generate_masked_patterns(kmer, max_mismatches):
                    index[k][pattern].append((kmer, accession, position))

    return index


def match_kmers_threshold(human_csv, virus_csv, output_csv, max_mismatches):
    """
    Match virus kmers to human kmers allowing up to max_mismatches.
    """
    human_index = load_human_kmers_index(human_csv, max_mismatches)

    total_rows = 0
    written_rows = 0

    with open(virus_csv, newline="") as vf, open(output_csv, "w", newline="") as out:
        reader = csv.DictReader(vf)
        required = {"kmer", "accession", "position"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Missing required columns in {virus_csv}: {sorted(missing)}")

        writer = csv.writer(out)
        writer.writerow([
            "virus_kmer",
            "human_kmer",
            "matches",
            "mismatches",
            "identity_fraction",
            "identity_percent",
            "human_accession",
            "human_position",
            "virus_accession",
            "virus_position"
        ])

        for row in reader:
            total_rows += 1

            virus_kmer = row["kmer"].strip()
            virus_accession = row["accession"].strip()
            virus_position = row["position"].strip()
            k = len(virus_kmer)

            candidate_matches = set()

            if max_mismatches == 0:
                for candidate in human_index[k].get(virus_kmer, []):
                    candidate_matches.add(candidate)
            else:
                for pattern in generate_masked_patterns(virus_kmer, max_mismatches):
                    for candidate in human_index[k].get(pattern, []):
                        candidate_matches.add(candidate)

            for human_kmer, human_accession, human_position in candidate_matches:
                mismatches = count_mismatches(virus_kmer, human_kmer)

                if mismatches <= max_mismatches:
                    matches = k - mismatches
                    identity_fraction = matches / k
                    identity_percent = identity_fraction * 100

                    writer.writerow([
                        virus_kmer,
                        human_kmer,
                        matches,
                        mismatches,
                        round(identity_fraction, 4),
                        round(identity_percent, 2),
                        human_accession,
                        human_position,
                        virus_accession,
                        virus_position
                    ])
                    written_rows += 1

            if total_rows % 100000 == 0:
                print(f"[INFO] Processed {total_rows} virus kmers; wrote {written_rows} matches")

    print(f"[INFO] Finished. Processed {total_rows} virus kmers; wrote {written_rows} matches")
    print(f"[INFO] Output written to: {output_csv}")


if __name__ == "__main__":
    if len(sys.argv) not in [4, 5]:
        print(
            "Usage: python match_kmers_threshold.py "
            "<human_kmers.csv> <virus_kmers.csv> <max_mismatches> [output.csv]"
        )
        sys.exit(1)

    human_csv = sys.argv[1]
    virus_csv = sys.argv[2]

    try:
        max_mismatches = int(sys.argv[3])
    except ValueError:
        print("Error: <max_mismatches> must be an integer")
        sys.exit(1)

    if max_mismatches < 0:
        print("Error: <max_mismatches> must be >= 0")
        sys.exit(1)

    if len(sys.argv) == 5:
        output_csv = sys.argv[4]
    else:
        output_csv = generate_output_filename(human_csv, virus_csv, max_mismatches)

    match_kmers_threshold(human_csv, virus_csv, output_csv, max_mismatches)