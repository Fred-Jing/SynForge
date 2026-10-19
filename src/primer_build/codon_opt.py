from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Iterable

from .thermo import gc_fraction, longest_homopolymer


# Relative synonymous-codon frequencies for a practical E. coli starting model.
# Values are normalized within each amino acid. They are intentionally used as
# soft preferences rather than hard claims about expression yield.
ECOLI_CODON_USAGE: dict[str, tuple[tuple[str, float], ...]] = {
    "A": (("GCG", 0.36), ("GCC", 0.27), ("GCA", 0.21), ("GCT", 0.16)),
    "R": (
        ("CGT", 0.36),
        ("CGC", 0.36),
        ("CGG", 0.11),
        ("CGA", 0.07),
        ("AGA", 0.07),
        ("AGG", 0.03),
    ),
    "N": (("AAC", 0.55), ("AAT", 0.45)),
    "D": (("GAT", 0.63), ("GAC", 0.37)),
    "C": (("TGC", 0.55), ("TGT", 0.45)),
    "Q": (("CAG", 0.65), ("CAA", 0.35)),
    "E": (("GAA", 0.68), ("GAG", 0.32)),
    "G": (("GGC", 0.37), ("GGT", 0.35), ("GGG", 0.15), ("GGA", 0.13)),
    "H": (("CAT", 0.57), ("CAC", 0.43)),
    "I": (("ATT", 0.49), ("ATC", 0.39), ("ATA", 0.12)),
    "L": (
        ("CTG", 0.50),
        ("TTG", 0.13),
        ("TTA", 0.13),
        ("CTT", 0.12),
        ("CTC", 0.10),
        ("CTA", 0.02),
    ),
    "K": (("AAA", 0.76), ("AAG", 0.24)),
    "M": (("ATG", 1.0),),
    "F": (("TTC", 0.57), ("TTT", 0.43)),
    "P": (("CCG", 0.53), ("CCA", 0.19), ("CCT", 0.16), ("CCC", 0.12)),
    "S": (
        ("AGC", 0.28),
        ("TCT", 0.17),
        ("TCC", 0.15),
        ("TCG", 0.15),
        ("AGT", 0.14),
        ("TCA", 0.11),
    ),
    "T": (("ACC", 0.40), ("ACG", 0.24), ("ACT", 0.19), ("ACA", 0.17)),
    "W": (("TGG", 1.0),),
    "Y": (("TAT", 0.57), ("TAC", 0.43)),
    "V": (("GTG", 0.37), ("GTT", 0.28), ("GTC", 0.20), ("GTA", 0.15)),
}


@dataclass(frozen=True)
class CodonCandidate:
    cds: str
    intrinsic_score: float


def _tail_gc_penalty(dna: str, rules: dict) -> float:
    window = int(rules["gc_window_nt"])
    if len(dna) < window:
        return 0.0
    tail = dna[-window:]
    observed = gc_fraction(tail)
    low = float(rules["gc_window_min"])
    high = float(rules["gc_window_max"])
    target = float(rules["gc_target"])
    weight = float(rules["gc_penalty_weight"])
    penalty = 0.25 * weight * (observed - target) ** 2
    if observed < low:
        penalty += weight * (low - observed) ** 2 * 10
    elif observed > high:
        penalty += weight * (observed - high) ** 2 * 10
    return penalty


def _repeat_penalty(dna: str, rules: dict) -> float:
    k = int(rules["repeat_kmer"])
    if len(dna) < 2 * k:
        return 0.0
    current = dna[-k:]
    previous = dna[:-k]
    return float(rules["repeat_penalty"]) if current in previous else 0.0


def _homopolymer_penalty(dna: str, rules: dict) -> float:
    longest = longest_homopolymer(dna[-12:])
    if longest <= 4:
        return 0.0
    return float(rules["homopolymer_penalty"]) * (longest - 4) ** 2


def _sd_like_search_penalty(
    dna: str,
    motifs: Iterable[str],
    penalty: float,
) -> float:
    # A newly appended codon contributes three bases, so a motif may end at
    # any of those three positions rather than exactly at the codon boundary.
    return (
        penalty
        if any(motif in dna[-(len(motif) + 2) :] for motif in motifs)
        else 0.0
    )


def _codons_for_residue(residue: str, ncaa_codon: str) -> Iterable[tuple[str, float]]:
    if residue == "X":
        return ((ncaa_codon, 1.0),)
    try:
        return ECOLI_CODON_USAGE[residue]
    except KeyError as error:
        raise ValueError(f"Unsupported amino-acid symbol: {residue!r}") from error


def _score_complete_cds(
    protein: str,
    cds: str,
    rules: dict,
    ncaa_codon: str,
) -> float:
    if len(cds) != 3 * len(protein):
        raise ValueError("CDS length does not match protein length")
    usage_weight = float(rules["codon_usage_weight"])
    frequency_floor = float(rules["rare_codon_frequency_floor"])
    score = 0.0
    partial = ""
    for index, residue in enumerate(protein):
        codon = cds[index * 3 : index * 3 + 3]
        frequencies = dict(_codons_for_residue(residue, ncaa_codon))
        if codon not in frequencies:
            raise ValueError(
                f"Codon {codon!r} does not encode residue {residue!r} "
                f"at position {index + 1}"
            )
        partial += codon
        score -= usage_weight * math.log(
            max(float(frequencies[codon]), frequency_floor)
        )
        score += _tail_gc_penalty(partial, rules)
        score += _repeat_penalty(partial, rules)
        score += _homopolymer_penalty(partial, rules)
    return score


def _select_diverse_beam(
    expanded: list[tuple[float, str]],
    beam_width: int,
    five_prime_nt: int,
) -> list[tuple[float, str]]:
    expanded.sort(key=lambda item: (item[0], item[1]))
    if not expanded:
        return []

    selected: list[tuple[float, str]] = []
    selected_dna: set[str] = set()

    def add(item: tuple[float, str]) -> bool:
        if item[1] in selected_dna:
            return False
        selected.append(item)
        selected_dna.add(item[1])
        return True

    # Preserve the globally best half of the beam.
    best_quota = max(1, beam_width // 2)
    for item in expanded:
        add(item)
        if len(selected) >= best_quota:
            break

    # Once the translation-initiation window has been synthesized, reserve a
    # quarter of the beam for distinct synonymous 5' sequences. Without this,
    # a conventional score-only beam converges to one identical start and
    # leaves no sequence space for downstream RNA-structure optimization.
    if len(expanded[0][1]) >= five_prime_nt:
        prefix_quota = max(1, beam_width // 4)
        prefixes = {item[1][:five_prime_nt] for item in selected}
        for item in expanded:
            prefix = item[1][:five_prime_nt]
            if prefix in prefixes:
                continue
            if add(item):
                prefixes.add(prefix)
            if len(prefixes) >= prefix_quota:
                break

    # Fill the remainder with diverse recent sequence contexts, retaining the
    # original assembly-aware behavior of the optimizer.
    seen_suffixes = {item[1][-18:] for item in selected}
    for item in expanded:
        suffix = item[1][-18:]
        if suffix in seen_suffixes:
            continue
        if add(item):
            seen_suffixes.add(suffix)
        if len(selected) >= beam_width:
            return selected

    for item in expanded:
        add(item)
        if len(selected) >= beam_width:
            break
    return selected


def _select_diverse_candidates(
    ordered: list[tuple[str, float]],
    candidate_count: int,
    five_prime_nt: int,
) -> list[tuple[str, float]]:
    if len(ordered) <= candidate_count:
        return ordered
    selected: list[tuple[str, float]] = []
    selected_dna: set[str] = set()
    prefixes: set[str] = set()

    best_quota = max(1, candidate_count // 3)
    for item in ordered[:best_quota]:
        selected.append(item)
        selected_dna.add(item[0])
        prefixes.add(item[0][:five_prime_nt])

    for item in ordered:
        prefix = item[0][:five_prime_nt]
        if item[0] in selected_dna or prefix in prefixes:
            continue
        selected.append(item)
        selected_dna.add(item[0])
        prefixes.add(prefix)
        if len(selected) >= candidate_count:
            break

    if len(selected) >= candidate_count:
        return sorted(selected, key=lambda item: (item[1], item[0]))

    for item in ordered:
        if item[0] in selected_dna:
            continue
        selected.append(item)
        selected_dna.add(item[0])
        if len(selected) >= candidate_count:
            break
    return sorted(selected, key=lambda item: (item[1], item[0]))


def optimize_codons(
    protein_sequence: str,
    rules: dict,
    ncaa_codon: str = "TAG",
    sd_like_motifs: Iterable[str] = (),
    sd_like_search_penalty: float = 0.0,
) -> list[CodonCandidate]:
    """Return a diverse beam of expression/manufacturability-aware CDS candidates."""
    protein = protein_sequence.strip().upper()
    beam_width = int(rules["beam_width"])
    candidate_count = int(rules["candidate_count"])
    usage_weight = float(rules["codon_usage_weight"])
    frequency_floor = float(rules["rare_codon_frequency_floor"])
    five_prime_nt = int(rules["five_prime_diversity_nt"])
    normalized_sd_motifs = tuple(
        motif.upper().replace("U", "T") for motif in sd_like_motifs
    )
    states: list[tuple[float, str]] = [(0.0, "")]

    for residue in protein:
        expanded: list[tuple[float, str]] = []
        for score, dna in states:
            for codon, frequency in _codons_for_residue(residue, ncaa_codon):
                next_dna = dna + codon
                use_frequency = max(float(frequency), frequency_floor)
                next_score = score - usage_weight * math.log(use_frequency)
                next_score += _tail_gc_penalty(next_dna, rules)
                next_score += _repeat_penalty(next_dna, rules)
                next_score += _homopolymer_penalty(next_dna, rules)
                next_score += _sd_like_search_penalty(
                    next_dna,
                    normalized_sd_motifs,
                    sd_like_search_penalty,
                )
                expanded.append((next_score, next_dna))

        states = _select_diverse_beam(expanded, beam_width, five_prime_nt)

    unique: dict[str, float] = {}
    for _, cds in states:
        intrinsic_score = _score_complete_cds(
            protein,
            cds,
            rules,
            ncaa_codon,
        )
        unique[cds] = min(
            intrinsic_score,
            unique.get(cds, float("inf")),
        )

    # Recombine retained synonymous 5' prefixes with the best-scoring tail.
    # This isolates translation-initiation diversity from unrelated downstream
    # codon changes, avoiding an artificial score surcharge on every open-start
    # candidate.
    if states and len(states[0][1]) > five_prime_nt:
        best_tail_source = min(states, key=lambda item: (item[0], item[1]))[1]
        for _, cds in states:
            hybrid = cds[:five_prime_nt] + best_tail_source[five_prime_nt:]
            if hybrid in unique:
                continue
            unique[hybrid] = _score_complete_cds(
                protein,
                hybrid,
                rules,
                ncaa_codon,
            )

    ordered = sorted(
        unique.items(),
        key=lambda item: (
            item[1]
            + (
                sd_like_search_penalty
                if any(motif in item[0] for motif in normalized_sd_motifs)
                else 0.0
            ),
            item[1],
            item[0],
        ),
    )
    selected = _select_diverse_candidates(
        ordered,
        candidate_count,
        five_prime_nt,
    )
    return [
        CodonCandidate(cds=cds, intrinsic_score=score)
        for cds, score in selected
    ]
