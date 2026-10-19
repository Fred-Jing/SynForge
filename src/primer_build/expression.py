from __future__ import annotations

from dataclasses import asdict, dataclass
from functools import lru_cache
from typing import Iterable

try:
    import RNA
except ImportError as error:  # pragma: no cover - exercised by deployment checks
    raise RuntimeError(
        "ViennaRNA is required for 5' translation-initiation analysis. "
        "Install the pinned project dependencies with: "
        "uv sync --locked"
    ) from error

from .codon_opt import CodonCandidate


@dataclass(frozen=True)
class SdLikeSite:
    start_nt: int
    end_nt: int
    motif: str


@dataclass(frozen=True)
class ExpressionMetrics:
    five_prime_mfe_kcal_mol: float
    five_prime_structure: str
    five_prime_upstream_nt: int
    five_prime_cds_nt: int
    initiation_unpaired_probability: float
    local_mfe_min_kcal_mol: float
    local_mfe_min_start_nt: int
    local_window_nt: int
    sd_like_sites: tuple[SdLikeSite, ...]
    expression_penalty: float
    risk_flags: tuple[str, ...]

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class RankedCodonCandidate:
    cds: str
    intrinsic_score: float
    combined_score: float
    metrics: ExpressionMetrics


@lru_cache(maxsize=16_384)
def _rna_mfe(sequence: str, temperature_c: float) -> tuple[str, float]:
    rna = sequence.upper().replace("T", "U")
    model = RNA.md()
    model.temperature = temperature_c
    fold_compound = RNA.fold_compound(rna, model)
    structure, mfe = fold_compound.mfe()
    return structure, round(float(mfe), 3)


@lru_cache(maxsize=4_096)
def _rna_mfe_and_unpaired(
    sequence: str,
    temperature_c: float,
) -> tuple[str, float, tuple[float, ...]]:
    rna = sequence.upper().replace("T", "U")
    model = RNA.md()
    model.temperature = temperature_c
    fold_compound = RNA.fold_compound(rna, model)
    structure, mfe = fold_compound.mfe()
    fold_compound.pf()
    pair_probabilities = fold_compound.bpp()
    unpaired: list[float] = []
    for position in range(1, len(rna) + 1):
        paired_probability = sum(
            pair_probabilities[min(position, other)][max(position, other)]
            for other in range(1, len(rna) + 1)
            if other != position
        )
        unpaired.append(max(0.0, min(1.0, 1.0 - paired_probability)))
    return (
        structure,
        round(float(mfe), 3),
        tuple(unpaired),
    )


def _local_mfe_min(
    cds: str,
    window_nt: int,
    step_nt: int,
    temperature_c: float,
) -> tuple[float, int, int]:
    effective_window = min(window_nt, len(cds))
    last_start = len(cds) - effective_window
    starts = list(range(0, last_start + 1, step_nt))
    if starts[-1] != last_start:
        starts.append(last_start)
    values = [
        (_rna_mfe(cds[start : start + effective_window], temperature_c)[1], start)
        for start in starts
    ]
    mfe, start = min(values, key=lambda item: (item[0], item[1]))
    return mfe, start + 1, effective_window


def _find_sd_like_sites(
    cds: str,
    motifs: Iterable[str],
    scan_start_nt: int,
) -> tuple[SdLikeSite, ...]:
    sequence = cds.upper().replace("U", "T")
    configured_motifs = list(
        dict.fromkeys(motif.upper().replace("U", "T") for motif in motifs)
    )
    ordered_motifs = sorted(
        configured_motifs,
        key=lambda motif: (-len(motif), configured_motifs.index(motif)),
    )
    motif_priority = {
        motif: index for index, motif in enumerate(configured_motifs)
    }
    raw_hits: list[SdLikeSite] = []
    for start in range(max(0, scan_start_nt), len(sequence)):
        for motif in ordered_motifs:
            if sequence.startswith(motif, start):
                raw_hits.append(
                    SdLikeSite(
                        start_nt=start + 1,
                        end_nt=start + len(motif),
                        motif=motif,
                    )
                )
                break

    # A six-base SD-like sequence also contains overlapping five-base motifs.
    # Keep one representative per overlapping cluster so the penalty is not
    # counted twice.
    collapsed: list[SdLikeSite] = []
    for hit in raw_hits:
        if collapsed and hit.start_nt <= collapsed[-1].end_nt:
            previous = collapsed[-1]
            if (
                len(hit.motif) > len(previous.motif)
                or (
                    len(hit.motif) == len(previous.motif)
                    and motif_priority[hit.motif]
                    < motif_priority[previous.motif]
                )
            ):
                collapsed[-1] = hit
            continue
        collapsed.append(hit)
    return tuple(collapsed)


def evaluate_expression_metrics(
    cds: str,
    upstream_sequence: str,
    rules: dict,
) -> ExpressionMetrics:
    temperature_c = float(rules["temperature_c"])
    upstream_nt = min(
        int(rules["five_prime_upstream_nt"]),
        len(upstream_sequence),
    )
    cds_nt = min(int(rules["five_prime_cds_nt"]), len(cds))
    upstream_context = upstream_sequence[-upstream_nt:]
    five_prime_sequence = upstream_context + cds[:cds_nt]
    five_prime_structure, five_prime_mfe, unpaired = _rna_mfe_and_unpaired(
        five_prime_sequence,
        temperature_c,
    )
    rbs_motif = str(rules["rbs_motif"]).upper().replace("U", "T")
    rbs_start = upstream_context.rfind(rbs_motif)
    if rbs_start < 0:
        raise ValueError(
            f"Configured RBS motif {rbs_motif!r} was not found within the "
            f"last {upstream_nt} nt of the fixed vector context"
        )
    initiation_positions = [
        *range(rbs_start, rbs_start + len(rbs_motif)),
        *range(upstream_nt, min(upstream_nt + 3, len(five_prime_sequence))),
    ]
    initiation_unpaired_probability = round(
        sum(unpaired[position] for position in initiation_positions)
        / len(initiation_positions),
        4,
    )
    local_mfe, local_start, local_window = _local_mfe_min(
        cds,
        int(rules["local_window_nt"]),
        int(rules["local_step_nt"]),
        temperature_c,
    )
    sd_like_sites = _find_sd_like_sites(
        cds,
        rules["sd_like_motifs"],
        int(rules["sd_scan_start_nt"]),
    )

    five_prime_floor = float(rules["five_prime_mfe_soft_floor_kcal_mol"])
    local_floor = float(rules["local_mfe_soft_floor_kcal_mol"])
    accessibility_floor = float(
        rules["initiation_unpaired_probability_soft_floor"]
    )
    five_prime_penalty = (
        max(0.0, -five_prime_mfe)
        * float(rules["five_prime_mfe_base_weight"])
        + max(0.0, five_prime_floor - five_prime_mfe)
        * float(rules["five_prime_mfe_below_floor_weight"])
    )
    local_penalty = (
        max(0.0, local_floor - local_mfe)
        * float(rules["local_mfe_below_floor_weight"])
    )
    accessibility_penalty = (
        max(0.0, accessibility_floor - initiation_unpaired_probability)
        * float(rules["initiation_accessibility_penalty_weight"])
    )
    sd_penalty = len(sd_like_sites) * float(rules["sd_like_penalty"])
    expression_penalty = round(
        five_prime_penalty
        + local_penalty
        + accessibility_penalty
        + sd_penalty,
        4,
    )

    risk_flags: list[str] = []
    if five_prime_mfe < five_prime_floor:
        risk_flags.append("stable_5prime_structure")
    if initiation_unpaired_probability < accessibility_floor:
        risk_flags.append("low_initiation_accessibility")
    if local_mfe < local_floor:
        risk_flags.append("stable_internal_structure")
    if sd_like_sites:
        risk_flags.append("internal_sd_like")
    return ExpressionMetrics(
        five_prime_mfe_kcal_mol=five_prime_mfe,
        five_prime_structure=five_prime_structure,
        five_prime_upstream_nt=upstream_nt,
        five_prime_cds_nt=cds_nt,
        initiation_unpaired_probability=initiation_unpaired_probability,
        local_mfe_min_kcal_mol=local_mfe,
        local_mfe_min_start_nt=local_start,
        local_window_nt=local_window,
        sd_like_sites=sd_like_sites,
        expression_penalty=expression_penalty,
        risk_flags=tuple(risk_flags),
    )


def rank_expression_candidates(
    candidates: Iterable[CodonCandidate],
    upstream_sequence: str,
    rules: dict,
) -> list[RankedCodonCandidate]:
    ranked: list[RankedCodonCandidate] = []
    for candidate in candidates:
        metrics = evaluate_expression_metrics(
            candidate.cds,
            upstream_sequence,
            rules,
        )
        ranked.append(
            RankedCodonCandidate(
                cds=candidate.cds,
                intrinsic_score=candidate.intrinsic_score,
                combined_score=candidate.intrinsic_score
                + metrics.expression_penalty,
                metrics=metrics,
            )
        )
    return sorted(
        ranked,
        key=lambda candidate: (
            candidate.combined_score,
            len(candidate.metrics.risk_flags),
            -candidate.metrics.initiation_unpaired_probability,
            -candidate.metrics.five_prime_mfe_kcal_mol,
            -candidate.metrics.local_mfe_min_kcal_mol,
            candidate.intrinsic_score,
            candidate.cds,
        ),
    )
