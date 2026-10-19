from __future__ import annotations

from functools import lru_cache
from Bio.SeqUtils import MeltingTemp as mt


DNA_ALPHABET = frozenset("ACGT")


def normalize_dna(sequence: str) -> str:
    seq = "".join(sequence.split()).upper()
    invalid = sorted(set(seq) - DNA_ALPHABET)
    if invalid:
        raise ValueError(f"Invalid DNA symbols: {invalid}")
    return seq


_COMPLEMENT = str.maketrans("ACGT", "TGCA")


@lru_cache(maxsize=100_000)
def reverse_complement(sequence: str) -> str:
    return normalize_dna(sequence).translate(_COMPLEMENT)[::-1]


def snapgene_tm(sequence: str) -> float:
    """Match the Tm scale stored in the existing SnapGene primer packets."""
    seq = normalize_dna(sequence)
    if len(seq) < 2:
        raise ValueError("Tm requires at least 2 nt")
    return float(
        mt.Tm_NN(
            seq,
            nn_table=mt.DNA_NN4,
            dnac1=250,
            dnac2=0,
            Na=50,
            Mg=0,
            dNTPs=0,
            saltcorr=5,
        )
    )


def gc_fraction(sequence: str) -> float:
    seq = normalize_dna(sequence)
    return (seq.count("G") + seq.count("C")) / len(seq)


def longest_homopolymer(sequence: str) -> int:
    seq = normalize_dna(sequence)
    longest = current = 1 if seq else 0
    for left, right in zip(seq, seq[1:]):
        current = current + 1 if left == right else 1
        longest = max(longest, current)
    return longest


def max_contiguous_complementarity(left: str, right: str) -> int:
    """Maximum contiguous match between left and reverse-complement(right)."""
    a = normalize_dna(left)
    b = reverse_complement(right)
    previous = [0] * (len(b) + 1)
    best = 0
    for char_a in a:
        current = [0]
        for index, char_b in enumerate(b, start=1):
            value = previous[index - 1] + 1 if char_a == char_b else 0
            current.append(value)
            best = max(best, value)
        previous = current
    return best


def three_prime_complementarity(left: str, right: str, tail: int = 12) -> int:
    a = normalize_dna(left)
    b = normalize_dna(right)
    reverse_b = reverse_complement(b)
    reverse_a = reverse_complement(a)
    best = 0
    for length in range(1, min(tail, len(a), len(b)) + 1):
        if a[-length:] in reverse_b:
            best = length
        if b[-length:] in reverse_a:
            best = length
    return best
