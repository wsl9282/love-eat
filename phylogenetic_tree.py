#!/usr/bin/env python3
import argparse
from typing import Dict, List, Tuple


def read_fasta(path: str) -> Dict[str, str]:
    sequences: Dict[str, str] = {}
    current_name = None
    current_seq: List[str] = []

    with open(path, "r", encoding="utf-8") as f:
        for raw_line in f:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_name is not None:
                    sequences[current_name] = "".join(current_seq).upper()
                current_name = line[1:].strip()
                if not current_name:
                    raise ValueError("FASTA header cannot be empty.")
                if current_name in sequences:
                    raise ValueError(f"Duplicated FASTA header: {current_name}")
                current_seq = []
            else:
                if current_name is None:
                    raise ValueError("Invalid FASTA format: sequence found before header.")
                current_seq.append(line)

    if current_name is not None:
        sequences[current_name] = "".join(current_seq).upper()

    if len(sequences) < 2:
        raise ValueError("At least 2 sequences are required to build a phylogenetic tree.")

    lengths = {len(seq) for seq in sequences.values()}
    if len(lengths) != 1:
        raise ValueError("All sequences must have the same length (aligned FASTA required).")

    return sequences


def p_distance(seq1: str, seq2: str) -> float:
    differences = sum(1 for a, b in zip(seq1, seq2) if a != b)
    return differences / len(seq1)


def pair_key(a: int, b: int) -> Tuple[int, int]:
    return (a, b) if a < b else (b, a)


def build_upgma_newick(sequences: Dict[str, str]) -> str:
    names = list(sequences.keys())
    n = len(names)

    clusters = {
        i: {"newick": names[i], "height": 0.0, "size": 1}
        for i in range(n)
    }
    active = set(range(n))
    distances: Dict[Tuple[int, int], float] = {}

    for i in range(n):
        for j in range(i + 1, n):
            distances[(i, j)] = p_distance(sequences[names[i]], sequences[names[j]])

    next_id = n
    while len(active) > 1:
        best_pair = None
        best_dist = float("inf")
        ordered_active = sorted(active)
        for idx, i in enumerate(ordered_active):
            for j in ordered_active[idx + 1:]:
                d = distances[pair_key(i, j)]
                if d < best_dist:
                    best_dist = d
                    best_pair = (i, j)

        if best_pair is None:
            raise RuntimeError("Failed to find a pair to merge.")

        i, j = best_pair
        ci, cj = clusters[i], clusters[j]
        new_height = best_dist / 2.0
        branch_i = max(0.0, new_height - ci["height"])
        branch_j = max(0.0, new_height - cj["height"])

        newick = f"({ci['newick']}:{branch_i:.6f},{cj['newick']}:{branch_j:.6f})"
        new_cluster = {
            "newick": newick,
            "height": new_height,
            "size": ci["size"] + cj["size"],
        }

        clusters[next_id] = new_cluster

        for k in list(active):
            if k in (i, j):
                continue
            dik = distances[pair_key(i, k)]
            djk = distances[pair_key(j, k)]
            weighted = (dik * ci["size"] + djk * cj["size"]) / new_cluster["size"]
            distances[pair_key(next_id, k)] = weighted

        active.remove(i)
        active.remove(j)
        active.add(next_id)
        next_id += 1

    root_id = next(iter(active))
    return clusters[root_id]["newick"] + ";"


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build a phylogenetic tree (UPGMA) from aligned FASTA sequences and output Newick."
    )
    parser.add_argument("-i", "--input", required=True, help="Input aligned FASTA file.")
    parser.add_argument("-o", "--output", help="Output Newick file path (optional).")
    args = parser.parse_args()

    sequences = read_fasta(args.input)
    newick = build_upgma_newick(sequences)

    if args.output:
        with open(args.output, "w", encoding="utf-8") as f:
            f.write(newick + "\n")
    else:
        print(newick)


if __name__ == "__main__":
    main()
