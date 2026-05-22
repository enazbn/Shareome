#!/usr/bin/env python3

"""
match_kmers_blosum.py

Optimized BLOSUM-based similarity matching between viral and human k-mers.

This version is designed for large Shareome runs.

Main optimization:
    - Human k-mer locations are collapsed by unique k-mer.
    - The masked candidate index stores unique human k-mers only.
    - Virus k-mer locations are also collapsed by unique k-mer.
    - BLOSUM scores are calculated once per unique virus_kmer/human_kmer pair.
    - Location expansion happens only after a pair passes the BLOSUM threshold.

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
"""

import argparse
import csv
import os
from collections import defaultdict
from itertools import combinations


VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")


# ------------------------------------------------------------
# Basic helpers
# ------------------------------------------------------------

def is_valid_kmer(kmer):
    return bool(kmer) and set(kmer).issubset(VALID_AA)


def count_mismatches(seq1, seq2):
    return sum(a != b for a, b in zip(seq1, seq2))


def precompute_mask_positions(k, num_masks):
    """
    Precompute mask-combination positions once per k.
    """
    if num_masks == 0:
        return [()]

    return list(combinations(range(k), num_masks))


def generate_masked_patterns_from_positions(kmer, mask_positions):
    """
    Generate masked patterns using precomputed mask positions.
    """
    if mask_positions == [()]:
        yield kmer
        return

    kmer_list = list(kmer)

    for positions in mask_positions:
        masked = kmer_list.copy()

        for pos in positions:
            masked[pos] = "*"

        yield "".join(masked)


# ------------------------------------------------------------
# BLOSUM helpers
# ------------------------------------------------------------

def load_blosum_matrix(matrix_name):
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


def matrix_to_fast_lookup(matrix):
    """
    Convert Biopython substitution matrix to plain dictionaries for faster scoring.
    """
    pair_score = {}
    self_score = {}

    aas = sorted(VALID_AA)

    for aa1 in aas:
        for aa2 in aas:
            try:
                score = matrix[aa1, aa2]
            except Exception:
                score = matrix[aa2, aa1]

            pair_score[(aa1, aa2)] = float(score)

    for aa in aas:
        self_score[aa] = pair_score[(aa, aa)]

    return pair_score, self_score


def calculate_blosum_scores(virus_kmer, human_kmer, pair_score, self_score):
    raw_score = 0.0
    virus_self_score = 0.0
    human_self_score = 0.0

    for v_aa, h_aa in zip(virus_kmer, human_kmer):
        raw_score += pair_score[(v_aa, h_aa)]
        virus_self_score += self_score[v_aa]
        human_self_score += self_score[h_aa]

    if virus_self_score == 0:
        similarity_fraction = 0.0
    else:
        similarity_fraction = raw_score / virus_self_score

    similarity_percent = similarity_fraction * 100.0

    return (
        raw_score,
        virus_self_score,
        human_self_score,
        similarity_fraction,
        similarity_percent,
    )


# ------------------------------------------------------------
# Loading and indexing
# ------------------------------------------------------------

def load_human_unique_index(human_csv, max_candidate_mismatches):
    """
    Load human k-mers and build a masked-pattern index using unique human k-mers.

    Returns:
        human_locs:
            human_locs[human_kmer] = [(accession, position), ...]

        human_index:
            human_index[k][pattern] = [human_kmer1, human_kmer2, ...]

        mask_cache:
            mask_cache[k] = precomputed mask position combinations
    """
    human_locs = defaultdict(list)
    human_index = defaultdict(lambda: defaultdict(list))
    mask_cache = {}

    total_rows = 0
    valid_rows = 0
    unique_kmers = 0

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

            if not is_valid_kmer(kmer):
                continue

            accession = row["accession"].strip()
            position = row["position"].strip()

            is_new_kmer = kmer not in human_locs

            human_locs[kmer].append((accession, position))
            valid_rows += 1

            if is_new_kmer:
                unique_kmers += 1
                k = len(kmer)

                if k not in mask_cache:
                    mask_cache[k] = precompute_mask_positions(
                        k,
                        max_candidate_mismatches,
                    )

                for pattern in generate_masked_patterns_from_positions(
                    kmer,
                    mask_cache[k],
                ):
                    human_index[k][pattern].append(kmer)

            if total_rows % 1_000_000 == 0:
                print(
                    f"[INFO] Human rows read: {total_rows:,}; "
                    f"valid rows: {valid_rows:,}; "
                    f"unique human k-mers: {unique_kmers:,}",
                    flush=True,
                )

    print(
        f"[INFO] Finished human loading/indexing. "
        f"Rows read: {total_rows:,}; "
        f"valid rows: {valid_rows:,}; "
        f"unique human k-mers: {unique_kmers:,}",
        flush=True,
    )

    return human_locs, human_index, mask_cache


def load_virus_unique_locs(virus_csv):
    """
    Collapse virus rows by unique virus k-mer.

    Returns:
        virus_locs[virus_kmer] = [(accession, position), ...]
    """
    virus_locs = defaultdict(list)

    total_rows = 0
    valid_rows = 0

    with open(virus_csv, newline="") as f:
        reader = csv.DictReader(f)

        required = {"kmer", "accession", "position"}
        missing = required - set(reader.fieldnames or [])

        if missing:
            raise ValueError(
                f"Missing required columns in {virus_csv}: {sorted(missing)}"
            )

        for row in reader:
            total_rows += 1

            kmer = row["kmer"].strip().upper()

            if not is_valid_kmer(kmer):
                continue

            accession = row["accession"].strip()
            position = row["position"].strip()

            virus_locs[kmer].append((accession, position))
            valid_rows += 1

            if total_rows % 1_000_000 == 0:
                print(
                    f"[INFO] Virus rows read: {total_rows:,}; "
                    f"valid rows: {valid_rows:,}; "
                    f"unique virus k-mers: {len(virus_locs):,}",
                    flush=True,
                )

    print(
        f"[INFO] Finished virus loading. "
        f"Rows read: {total_rows:,}; "
        f"valid rows: {valid_rows:,}; "
        f"unique virus k-mers: {len(virus_locs):,}",
        flush=True,
    )

    return virus_locs


# ------------------------------------------------------------
# Main matching
# ------------------------------------------------------------

def match_kmers_blosum(
    human_csv,
    virus_csv,
    output_csv,
    matrix_name="BLOSUM62",
    min_similarity_percent=80.0,
    max_candidate_mismatches=1,
):
    print(f"[INFO] Loading matrix: {matrix_name}", flush=True)
    matrix = load_blosum_matrix(matrix_name)
    pair_score, self_score = matrix_to_fast_lookup(matrix)

    print(
        f"[INFO] Loading human k-mers using unique-kmer index. "
        f"max_candidate_mismatches={max_candidate_mismatches}",
        flush=True,
    )

    human_locs, human_index, mask_cache = load_human_unique_index(
        human_csv,
        max_candidate_mismatches,
    )

    print("[INFO] Loading virus k-mers as unique k-mer groups...", flush=True)
    virus_locs = load_virus_unique_locs(virus_csv)

    total_unique_virus = 0
    total_candidates = 0
    retained_unique_pairs = 0
    written_rows = 0

    with open(output_csv, "w", newline="") as out:
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

        for virus_kmer, v_locs in virus_locs.items():
            total_unique_virus += 1

            k = len(virus_kmer)

            if k not in mask_cache:
                mask_cache[k] = precompute_mask_positions(
                    k,
                    max_candidate_mismatches,
                )

            candidate_human_kmers = set()

            for pattern in generate_masked_patterns_from_positions(
                virus_kmer,
                mask_cache[k],
            ):
                candidate_human_kmers.update(
                    human_index[k].get(pattern, [])
                )

            total_candidates += len(candidate_human_kmers)

            for human_kmer in candidate_human_kmers:
                if len(human_kmer) != k:
                    continue

                mismatches = count_mismatches(virus_kmer, human_kmer)
                matches = k - mismatches
                identity_fraction = matches / k
                identity_percent = identity_fraction * 100.0

                (
                    blosum_raw_score,
                    blosum_self_virus,
                    blosum_self_human,
                    blosum_similarity_fraction,
                    blosum_similarity_percent,
                ) = calculate_blosum_scores(
                    virus_kmer,
                    human_kmer,
                    pair_score,
                    self_score,
                )

                if blosum_similarity_percent < min_similarity_percent:
                    continue

                retained_unique_pairs += 1
                h_locs = human_locs[human_kmer]

                for virus_accession, virus_position in v_locs:
                    for human_accession, human_position in h_locs:
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

            if total_unique_virus % 100_000 == 0:
                print(
                    f"[INFO] Unique virus k-mers processed: {total_unique_virus:,}; "
                    f"unique candidates tested: {total_candidates:,}; "
                    f"retained unique pairs: {retained_unique_pairs:,}; "
                    f"expanded rows written: {written_rows:,}",
                    flush=True,
                )

    print(
        f"[INFO] Finished BLOSUM matching. "
        f"Unique virus k-mers processed: {total_unique_virus:,}; "
        f"unique candidates tested: {total_candidates:,}; "
        f"retained unique pairs: {retained_unique_pairs:,}; "
        f"expanded rows written: {written_rows:,}",
        flush=True,
    )

    print(f"[INFO] Output written to: {output_csv}", flush=True)


# ------------------------------------------------------------
# CLI
# ------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description="Optimized BLOSUM-based k-mer similarity matcher."
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
        default=80.0,
        help="Minimum normalized BLOSUM similarity percent to keep. Default: 80.",
    )

    parser.add_argument(
        "--max-candidate-mismatches",
        type=int,
        default=1,
        help=(
            "Masked-index candidate retrieval depth. "
            "This is not the final biological cutoff. Default: 1."
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