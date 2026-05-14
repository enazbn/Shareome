#!/usr/bin/env python3

"""
match_kmers_blosum.py

BLOSUM-based similarity matching between viral and human k-mers.

This is intended as a secondary Shareome matching mode, separate from exact/mismatch
matching. It uses the same upstream k-mer CSV files produced by kmerslicer.

Input CSV columns:
    kmer
    accession
    position

Output CSV columns:
    virus_kmer
    human_kmer
    matches
    mismatches
    identity_fraction
    identity_percent
    blosum_raw_score
    blosum_self_virus
    blosum_self_human
    blosum_similarity_fraction
    blosum_similarity_percent
    human_accession
    human_position
    virus_accession
    virus_position

Example:
    python workflow/scripts/match_kmers_blosum.py \
        results/human_refseq_nr_cleaned_kmers.csv \
        results/Rota_seq_nr_clean_kmers.csv \
        results/Rota_seq_nr_clean_blosum_matched.csv \
        --matrix BLOSUM62 \
        --min-similarity-percent 70 \
        --max-candidate-mismatches 3
"""

import argparse
import csv
import os
import sys
from collections import defaultdict
from itertools import combinations


VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")


def generate_masked_patterns(kmer, num_masks):
    """
    Generate masked versions of a k-mer by replacing num_masks positions with '*'.

    Example:
        kmer = ABCD, num_masks = 2
        patterns include: **CD, *B*D, *BC*, A**D, A*C*, AB**
    """
    kmer_list = list(kmer)

    for positions in combinations(range(len(kmer)), num_masks):
        masked = kmer_list.copy()

        for pos in positions:
            masked[pos] = "*"

        yield "".join(masked)


def count_mismatches(seq1, seq2):
    """
    Count position-wise mismatches between two equal-length sequences.
    """
    return sum(a != b for a, b in zip(seq1, seq2))


def is_valid_kmer(kmer):
    """
    Keep only canonical amino-acid k-mers.
    """
    return bool(kmer) and set(kmer).issubset(VALID_AA)


def load_blosum_matrix(matrix_name):
    """
    Load a substitution matrix from Biopython.

    Requires:
        biopython

    Example:
        conda install -c conda-forge biopython
    """
    try:
        from Bio.Align import substitution_matrices
    except ImportError as exc:
        raise ImportError(
            "Biopython is required for BLOSUM mode. Install it with:\n"
            "  conda install -c conda-forge biopython\n"
            "or:\n"
            "  pip install biopython"
        ) from exc

    try:
        return substitution_matrices.load(matrix_name)
    except Exception as exc:
        raise ValueError(
            f"Could not load substitution matrix '{matrix_name}'. "
            "Common choices include BLOSUM62, BLOSUM80, and BLOSUM45."
        ) from exc


def substitution_score(matrix, aa1, aa2):
    """
    Get substitution score for an amino-acid pair.

    Biopython matrices usually support tuple access, but this helper also tries
    the reverse pair for safety.
    """
    try:
        return matrix[aa1, aa2]
    except Exception:
        return matrix[aa2, aa1]


def calculate_blosum_scores(virus_kmer, human_kmer, matrix):
    """
    Calculate raw and normalized BLOSUM similarity.

    Normalization used here:
        blosum_similarity_fraction = raw_score / virus_self_score

    This asks:
        How similar is the human k-mer to the viral k-mer relative to the
        viral k-mer's own maximum self-score?
    """
    raw_score = 0.0
    virus_self_score = 0.0
    human_self_score = 0.0

    for v_aa, h_aa in zip(virus_kmer, human_kmer):
        raw_score += substitution_score(matrix, v_aa, h_aa)
        virus_self_score += substitution_score(matrix, v_aa, v_aa)
        human_self_score += substitution_score(matrix, h_aa, h_aa)

    if virus_self_score == 0:
        similarity_fraction = 0.0
    else:
        similarity_fraction = raw_score / virus_self_score

    similarity_percent = similarity_fraction * 100

    return (
        raw_score,
        virus_self_score,
        human_self_score,
        similarity_fraction,
        similarity_percent,
    )


def load_human_kmers_index(human_csv, max_candidate_mismatches):
    """
    Build a masked-pattern index for human k-mers.

    This is only a candidate retrieval step. The final filtering is done using
    the BLOSUM similarity score.

    Returns:
        index[k][pattern] -> list of (human_kmer, accession, position)
    """
    index = defaultdict(lambda: defaultdict(list))

    total_rows = 0
    kept_rows = 0

    with open(human_csv, newline="") as f:
        reader = csv.DictReader(f)

        required = {"kmer", "accession", "position"}
        missing = required - set(reader.fieldnames or [])

        if missing:
            raise ValueError(
                f"Missing required columns in {human_csv}: {sorted(missing)}"
            )

        for row in reader:
            total_rows += 1

            kmer = row["kmer"].strip().upper()
            accession = row["accession"].strip()
            position = row["position"].strip()

            if not is_valid_kmer(kmer):
                continue

            k = len(kmer)

            if max_candidate_mismatches == 0:
                index[k][kmer].append((kmer, accession, position))
            else:
                for pattern in generate_masked_patterns(kmer, max_candidate_mismatches):
                    index[k][pattern].append((kmer, accession, position))

            kept_rows += 1

            if total_rows % 1_000_000 == 0:
                print(
                    f"[INFO] Indexed {total_rows:,} human rows; "
                    f"kept {kept_rows:,}",
                    flush=True,
                )

    print(
        f"[INFO] Finished indexing human k-mers. "
        f"Read {total_rows:,}; kept {kept_rows:,}.",
        flush=True,
    )

    return index


def match_kmers_blosum(
    human_csv,
    virus_csv,
    output_csv,
    matrix_name="BLOSUM62",
    min_similarity_percent=70.0,
    max_candidate_mismatches=3,
):
    """
    Match virus k-mers to human k-mers using BLOSUM similarity.

    Important:
        max_candidate_mismatches is not the biological cutoff.
        It only controls which human k-mers are considered as candidates.

        Final filtering is:
            blosum_similarity_percent >= min_similarity_percent
    """
    print(f"[INFO] Loading matrix: {matrix_name}", flush=True)
    matrix = load_blosum_matrix(matrix_name)

    print(
        f"[INFO] Building human candidate index using "
        f"max_candidate_mismatches={max_candidate_mismatches}",
        flush=True,
    )
    human_index = load_human_kmers_index(human_csv, max_candidate_mismatches)

    total_virus_rows = 0
    total_candidates = 0
    written_rows = 0

    with open(virus_csv, newline="") as vf, open(output_csv, "w", newline="") as out:
        reader = csv.DictReader(vf)

        required = {"kmer", "accession", "position"}
        missing = required - set(reader.fieldnames or [])

        if missing:
            raise ValueError(
                f"Missing required columns in {virus_csv}: {sorted(missing)}"
            )

        writer = csv.writer(out)

        writer.writerow(
            [
                "virus_kmer",
                "human_kmer",
                "matches",
                "mismatches",
                "identity_fraction",
                "identity_percent",
                "blosum_raw_score",
                "blosum_self_virus",
                "blosum_self_human",
                "blosum_similarity_fraction",
                "blosum_similarity_percent",
                "human_accession",
                "human_position",
                "virus_accession",
                "virus_position",
            ]
        )

        for row in reader:
            total_virus_rows += 1

            virus_kmer = row["kmer"].strip().upper()
            virus_accession = row["accession"].strip()
            virus_position = row["position"].strip()

            if not is_valid_kmer(virus_kmer):
                continue

            k = len(virus_kmer)
            candidate_matches = set()

            if max_candidate_mismatches == 0:
                for candidate in human_index[k].get(virus_kmer, []):
                    candidate_matches.add(candidate)
            else:
                for pattern in generate_masked_patterns(
                    virus_kmer,
                    max_candidate_mismatches,
                ):
                    for candidate in human_index[k].get(pattern, []):
                        candidate_matches.add(candidate)

            total_candidates += len(candidate_matches)

            for human_kmer, human_accession, human_position in candidate_matches:
                if len(human_kmer) != k:
                    continue

                mismatches = count_mismatches(virus_kmer, human_kmer)
                matches = k - mismatches
                identity_fraction = matches / k
                identity_percent = identity_fraction * 100

                (
                    blosum_raw_score,
                    blosum_self_virus,
                    blosum_self_human,
                    blosum_similarity_fraction,
                    blosum_similarity_percent,
                ) = calculate_blosum_scores(
                    virus_kmer,
                    human_kmer,
                    matrix,
                )

                if blosum_similarity_percent >= min_similarity_percent:
                    writer.writerow(
                        [
                            virus_kmer,
                            human_kmer,
                            matches,
                            mismatches,
                            round(identity_fraction, 4),
                            round(identity_percent, 2),
                            round(blosum_raw_score, 4),
                            round(blosum_self_virus, 4),
                            round(blosum_self_human, 4),
                            round(blosum_similarity_fraction, 4),
                            round(blosum_similarity_percent, 2),
                            human_accession,
                            human_position,
                            virus_accession,
                            virus_position,
                        ]
                    )

                    written_rows += 1

            if total_virus_rows % 100_000 == 0:
                print(
                    f"[INFO] Processed {total_virus_rows:,} virus k-mers; "
                    f"tested {total_candidates:,} candidates; "
                    f"wrote {written_rows:,} BLOSUM matches",
                    flush=True,
                )

    print(
        f"[INFO] Finished. Processed {total_virus_rows:,} virus k-mers; "
        f"tested {total_candidates:,} candidates; "
        f"wrote {written_rows:,} BLOSUM matches.",
        flush=True,
    )
    print(f"[INFO] Output written to: {output_csv}", flush=True)


def parse_args():
    parser = argparse.ArgumentParser(
        description="BLOSUM-based k-mer similarity matcher."
    )

    parser.add_argument(
        "human_csv",
        help="Human k-mer CSV file.",
    )

    parser.add_argument(
        "virus_csv",
        help="Virus k-mer CSV file.",
    )

    parser.add_argument(
        "output_csv",
        help="Output matched CSV file.",
    )

    parser.add_argument(
        "--matrix",
        default="BLOSUM62",
        help="Substitution matrix name. Default: BLOSUM62.",
    )

    parser.add_argument(
        "--min-similarity-percent",
        type=float,
        default=70.0,
        help="Minimum normalized BLOSUM similarity percent to keep. Default: 70.",
    )

    parser.add_argument(
        "--max-candidate-mismatches",
        type=int,
        default=3,
        help=(
            "Masked-index candidate retrieval depth. "
            "This is not the final cutoff. Default: 3."
        ),
    )

    return parser.parse_args()


def main():
    args = parse_args()

    if not os.path.exists(args.human_csv):
        raise FileNotFoundError(f"Human CSV does not exist: {args.human_csv}")

    if not os.path.exists(args.virus_csv):
        raise FileNotFoundError(f"Virus CSV does not exist: {args.virus_csv}")

    if args.max_candidate_mismatches < 0:
        raise ValueError("--max-candidate-mismatches must be >= 0")

    if args.min_similarity_percent < 0:
        raise ValueError("--min-similarity-percent must be >= 0")

    output_dir = os.path.dirname(args.output_csv)

    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    match_kmers_blosum(
        human_csv=args.human_csv,
        virus_csv=args.virus_csv,
        output_csv=args.output_csv,
        matrix_name=args.matrix,
        min_similarity_percent=args.min_similarity_percent,
        max_candidate_mismatches=args.max_candidate_mismatches,
    )


if __name__ == "__main__":
    main()