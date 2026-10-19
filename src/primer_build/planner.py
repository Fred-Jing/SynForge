from __future__ import annotations

from dataclasses import replace
from functools import lru_cache
from itertools import product
from typing import Iterable

from .models import AssemblyPlan, Junction, Oligo, SyntheticFragment
from .thermo import (
    gc_fraction,
    longest_homopolymer,
    reverse_complement,
    snapgene_tm,
    three_prime_complementarity,
)


@lru_cache(maxsize=200_000)
def _cached_tm(sequence: str) -> float:
    return snapgene_tm(sequence)


def _circular_occurrences(sequence: str, motif: str) -> int:
    motif = motif.upper()
    extended = sequence.upper() + sequence.upper()[: len(motif) - 1]
    # str.find runs the search in C; advance by one to retain overlapping hits.
    count = 0
    position = extended.find(motif)
    while 0 <= position < len(sequence):
        count += 1
        position = extended.find(motif, position + 1)
    return count



def _oligo_risks(sequence: str, rules: dict, compare: str | None = None) -> tuple[str, ...]:
    return _cached_oligo_risks(sequence, int(rules["oligos"]["homopolymer_max"]), compare)


@lru_cache(maxsize=100_000)
def _cached_oligo_risks(sequence: str, max_homopolymer: int, compare: str | None) -> tuple[str, ...]:
    flags: list[str] = []
    if longest_homopolymer(sequence) > max_homopolymer:
        flags.append("homopolymer")
    if three_prime_complementarity(sequence, sequence) >= 5:
        flags.append("3prime_self_complementarity")
    if compare and three_prime_complementarity(sequence, compare) >= 6:
        flags.append("3prime_pair_complementarity")
    return tuple(flags)


def _tier(length: int, rules: dict) -> str:
    if length <= int(rules["oligos"]["standard_tier_max_nt"]):
        return "standard"
    if length <= int(rules["oligos"]["long_tier_max_nt"]):
        return "long"
    return "unsupported"


def _tm_window_score(tm_c: float, rules: dict) -> float:
    thermo = rules["thermodynamics"]
    target = float(thermo["annealing_tm_target_c"])
    low = float(thermo["annealing_tm_min_c"])
    high = float(thermo["annealing_tm_max_c"])
    score = (tm_c - target) ** 2
    if tm_c < low:
        score += 30.0 * (low - tm_c) ** 2
    elif tm_c > high:
        score += 18.0 * (tm_c - high) ** 2
    return score


def _length_score(length: int, rules: dict) -> float:
    thermo = rules["thermodynamics"]
    low = int(thermo["annealing_length_preferred_min"])
    high = int(thermo["annealing_length_preferred_max"])
    if low <= length <= high:
        return 0.015 * (length - 26) ** 2
    distance = low - length if length < low else length - high
    return 5.0 * distance**2


def _valid_annealing_segments(sequence: str, rules: dict) -> tuple[tuple[int, float], ...]:
    thermo = rules["thermodynamics"]
    return _annealing_segments(sequence,
        int(thermo["annealing_length_hard_min"]), int(thermo["annealing_length_hard_max"]),
        float(thermo["annealing_tm_min_c"]), float(thermo["annealing_tm_max_c"]))


@lru_cache(maxsize=100_000)
def _annealing_segments(sequence: str, hard_min: int, hard_max: int,
                        low: float, high: float) -> tuple[tuple[int, float], ...]:
    values = []
    for length in range(hard_min, min(hard_max, len(sequence)) + 1):
        tm_c = _cached_tm(sequence[:length])
        if low <= tm_c <= high:
            values.append((length, tm_c))
    return tuple(values)


def _best_fragment_builder(
    insert_sequence: str,
    start: int,
    end: int,
    rules: dict,
) -> SyntheticFragment | None:
    fragment = insert_sequence[start:end]
    if len(fragment) < 60:
        return None
    oligo_rules = rules["oligos"]
    standard_max = int(oligo_rules["standard_tier_max_nt"])
    long_max = int(oligo_rules["long_tier_max_nt"])
    core_min = int(oligo_rules["preferred_core_min_nt"])
    best: tuple[float, SyntheticFragment] | None = None

    core_length_groups = [list(range(core_min, standard_max + 1))]
    if bool(oligo_rules["allow_long_core_pair"]):
        core_length_groups.append(list(range(standard_max + 1, long_max + 1)))

    for core_lengths in core_length_groups:
        for core_length in core_lengths:
            if core_length >= len(fragment):
                continue
            # Each outer primer consists of an optional extension plus a
            # 17–40 nt annealing segment and must remain <=59 nt.
            maximum_tail = standard_max - int(
                rules["thermodynamics"]["annealing_length_hard_min"]
            )
            core_start_min = max(0, len(fragment) - core_length - maximum_tail)
            core_start_max = min(maximum_tail, len(fragment) - core_length)
            for core_start_local in range(core_start_min, core_start_max + 1):
                core_end_local = core_start_local + core_length
                core_seq = fragment[core_start_local:core_end_local]
                try:
                    core_tm = _cached_tm(core_seq)
                except ValueError:
                    continue
                if core_tm < 62.0:
                    continue

                forward_segments = _valid_annealing_segments(core_seq, rules)
                reverse_segments = _valid_annealing_segments(reverse_complement(core_seq), rules)
                forward_choices = [
                    (annealing_length, tm_c)
                    for annealing_length, tm_c in forward_segments
                    if core_start_local + annealing_length <= standard_max
                ]
                reverse_choices = [
                    (annealing_length, tm_c)
                    for annealing_length, tm_c in reverse_segments
                    if (
                        len(fragment) - core_end_local + annealing_length
                        <= standard_max
                    )
                ]
                if not forward_choices or not reverse_choices:
                    continue
                forward_length, forward_tm = min(
                    forward_choices,
                    key=lambda item: _tm_window_score(item[1], rules)
                    + _length_score(item[0], rules),
                )
                reverse_length, reverse_tm = min(
                    reverse_choices,
                    key=lambda item: _tm_window_score(item[1], rules)
                    + _length_score(item[0], rules),
                )
                outer_forward_length = core_start_local + forward_length
                outer_reverse_length = (
                    len(fragment) - core_end_local + reverse_length
                )
                outer_forward_seq = fragment[:outer_forward_length]
                outer_reverse_seq = reverse_complement(
                    fragment[core_end_local - reverse_length :]
                )
                core_forward_seq = core_seq
                core_reverse_seq = reverse_complement(core_seq)
                core_tier = _tier(core_length, rules)
                if core_tier == "unsupported":
                    continue

                score = _tm_window_score(forward_tm, rules)
                score += _tm_window_score(reverse_tm, rules)
                score += _length_score(forward_length, rules)
                score += _length_score(reverse_length, rules)
                if core_tier == "long":
                    score += float(oligo_rules["long_core_pair_penalty"])
                    if oligo_rules.get("prefer_short_long_oligos", False):
                        # Synthesis turnaround and error risk matter more than
                        # filling a supplier's long-oligo price tier.
                        score += float(oligo_rules["short_long_core_penalty_weight"]) * (
                            core_length - standard_max
                        ) ** 2
                    else:
                        score += 0.02 * (core_length - long_max) ** 2
                else:
                    score += 0.02 * (core_length - 55) ** 2
                risks = set()
                for sequence in (
                    outer_forward_seq,
                    outer_reverse_seq,
                    core_forward_seq,
                    core_reverse_seq,
                ):
                    risks.update(_oligo_risks(sequence, rules))
                score += (100.0 if oligo_rules.get("prefer_short_long_oligos", False) else 8.0) * len(risks)

                global_core_start = start + core_start_local
                global_core_end = start + core_end_local
                candidate = SyntheticFragment(
                    name="fragment",
                    start=start,
                    end=end,
                    core_start=global_core_start,
                    core_end=global_core_end,
                    core_forward=Oligo(
                        name="core-FW",
                        sequence=core_forward_seq,
                        purpose="core_template",
                        orientation="FW",
                        start=global_core_start,
                        end=global_core_end,
                        effective_sequence=core_forward_seq,
                        effective_tm_c=core_tm,
                        full_tm_c=core_tm,
                        tier=core_tier,
                        risk_flags=_oligo_risks(core_forward_seq, rules),
                    ),
                    core_reverse=Oligo(
                        name="core-RV",
                        sequence=core_reverse_seq,
                        purpose="core_template",
                        orientation="RV",
                        start=global_core_start,
                        end=global_core_end,
                        effective_sequence=core_reverse_seq,
                        effective_tm_c=core_tm,
                        full_tm_c=core_tm,
                        tier=core_tier,
                        risk_flags=_oligo_risks(core_reverse_seq, rules),
                    ),
                    outer_forward=Oligo(
                        name="outer-FW",
                        sequence=outer_forward_seq,
                        purpose="fragment_extension",
                        orientation="FW",
                        start=start,
                        end=start + outer_forward_length,
                        effective_sequence=core_seq[:forward_length],
                        effective_tm_c=forward_tm,
                        full_tm_c=_cached_tm(outer_forward_seq),
                        tier="standard",
                        risk_flags=_oligo_risks(outer_forward_seq, rules),
                    ),
                    outer_reverse=Oligo(
                        name="outer-RV",
                        sequence=outer_reverse_seq,
                        purpose="fragment_extension",
                        orientation="RV",
                        start=end - outer_reverse_length,
                        end=end,
                        effective_sequence=reverse_complement(
                            core_seq[-reverse_length:]
                        ),
                        effective_tm_c=reverse_tm,
                        full_tm_c=_cached_tm(outer_reverse_seq),
                        tier="standard",
                        risk_flags=_oligo_risks(outer_reverse_seq, rules),
                    ),
                    score=score,
                    risk_flags=tuple(sorted(risks)),
                )
                if best is None or score < best[0]:
                    best = (score, candidate)
        if best is not None:
            # Do not pay for a 60–89 nt core pair when a standard-tier core
            # yields a valid fragment.
            break
    return best[1] if best else None


def _fragment_overlap_candidates_near(
    insert_sequence: str,
    target_center: float,
    rules: dict,
    radius: int = 24,
    keep: int = 12,
) -> list[tuple[Junction, float]]:
    thermo = rules["thermodynamics"]
    candidates: list[tuple[Junction, float]] = []
    lower = max(0, int(target_center) - radius)
    upper = min(len(insert_sequence), int(target_center) + radius)
    for start in range(lower, upper + 1):
        for length in range(
            int(thermo["annealing_length_preferred_min"]),
            int(thermo["annealing_length_preferred_max"]) + 1,
        ):
            end = start + length
            if end > len(insert_sequence):
                continue
            center = (start + end) / 2
            if abs(center - target_center) > radius:
                continue
            sequence = insert_sequence[start:end]
            tm_c = _cached_tm(sequence)
            if not (
                float(thermo["annealing_tm_min_c"])
                <= tm_c
                <= float(thermo["annealing_tm_max_c"])
            ):
                continue
            score = _tm_window_score(tm_c, rules) + _length_score(length, rules)
            score += 0.03 * (center - target_center) ** 2
            candidates.append(
                (
                    Junction(
                        name="fragment-overlap",
                        start=start,
                        end=end,
                        sequence=sequence,
                        tm_c=tm_c,
                        purpose="overlap_extension",
                    ),
                    score,
                )
            )
    candidates.sort(key=lambda item: item[1])
    # Preserve positional diversity instead of returning many nearly
    # identical windows from one local optimum.
    selected: list[tuple[Junction, float]] = []
    seen_center_bins: set[int] = set()
    for candidate in candidates:
        junction = candidate[0]
        center_bin = round(((junction.start + junction.end) / 2) / 2)
        if center_bin in seen_center_bins:
            continue
        seen_center_bins.add(center_bin)
        selected.append(candidate)
        if len(selected) >= keep:
            return selected
    for candidate in candidates:
        if candidate not in selected:
            selected.append(candidate)
            if len(selected) >= keep:
                break
    return selected


def plan_synthetic_fragments(
    insert_sequence: str,
    design_id: str,
    rules: dict,
    primer_design_id: str | None = None,
) -> tuple[list[SyntheticFragment], list[Junction], float] | None:
    """Find a 2–4 fragment plan using duplex cores plus extension primers.

    Boundary windows are searched around equal-length targets.  This retains
    the useful sequence-level optimization while avoiding an exhaustive
    O(n^4) enumeration of nearly equivalent fragment boundaries.
    """
    length = len(insert_sequence)
    minimum_length = 75
    maximum_length = 145
    max_count = int(rules["assembly"]["max_synthetic_fragments"])
    min_count = int(rules["assembly"]["min_synthetic_fragments"])
    preferred_count = int(rules["assembly"]["preferred_synthetic_fragments"])

    @lru_cache(maxsize=None)
    def build(start: int, end: int):
        return _best_fragment_builder(insert_sequence, start, end, rules)

    completed: list[tuple[float, list[SyntheticFragment], list[Junction]]] = []
    count_order = sorted(
        range(min_count, max_count + 1),
        key=lambda count: (abs(count - preferred_count), count),
    )
    for count in count_order:
        count_completed: list[
            tuple[float, list[SyntheticFragment], list[Junction]]
        ] = []
        effective_mean = length / count
        if not (
            minimum_length - 35 <= effective_mean <= maximum_length
        ):
            continue
        boundary_options = [
            _fragment_overlap_candidates_near(
                insert_sequence,
                target_center=length * boundary_index / count,
                rules=rules,
                radius=60 if count == preferred_count else 28,
                keep=16 if count == preferred_count else 3,
            )
            for boundary_index in range(1, count)
        ]
        if any(not options for options in boundary_options):
            continue
        for selected in product(*boundary_options):
            junctions = [item[0] for item in selected]
            if any(
                left.end >= right.start
                for left, right in zip(junctions, junctions[1:])
            ):
                continue
            fragment_ranges: list[tuple[int, int]] = [(0, junctions[0].end)]
            fragment_ranges.extend(
                (left.start, right.end)
                for left, right in zip(junctions, junctions[1:])
            )
            fragment_ranges.append((junctions[-1].start, length))
            if any(
                end - start < minimum_length or end - start > maximum_length
                for start, end in fragment_ranges
            ):
                continue
            fragments: list[SyntheticFragment] = []
            for start, end in fragment_ranges:
                fragment = build(start, end)
                if fragment is None:
                    break
                fragments.append(fragment)
            if len(fragments) != count:
                continue
            junction_score = sum(item[1] for item in selected)
            fragment_score = sum(fragment.score for fragment in fragments)
            count_completed.append(
                (
                    junction_score + fragment_score,
                    fragments,
                    junctions,
                )
            )
        if count_completed:
            # Fragment count is a rule-level preference. Scores only compare
            # plans within the same fragment-count class.
            completed = count_completed
            break

    if not completed:
        return None
    completed.sort(key=lambda item: item[0])
    score, fragments, junctions = completed[0]
    renamed: list[SyntheticFragment] = []
    primer_root = primer_design_id or design_id
    for index, fragment in enumerate(fragments, start=1):
        fragment_name = f"{primer_root}-{index}"
        renamed.append(
            replace(
                fragment,
                name=fragment_name,
                core_forward=replace(
                    fragment.core_forward, name=f"{fragment_name}-FW"
                ),
                core_reverse=replace(
                    fragment.core_reverse, name=f"{fragment_name}-RV"
                ),
                outer_forward=replace(
                    fragment.outer_forward, name=f"{fragment_name}-1-FW"
                ),
                outer_reverse=replace(
                    fragment.outer_reverse, name=f"{fragment_name}-1-RV"
                ),
            )
        )
    renamed_junctions = [
        replace(junction, name=f"{design_id}-OE{i}-{i+1}")
        for i, junction in enumerate(junctions, start=1)
    ]
    return renamed, renamed_junctions, score


def _homology_windows(
    sequence: str,
    start_min: int,
    start_max: int,
    length_min: int,
    length_max: int,
    tm_min: float,
) -> list[Junction]:
    values: list[Junction] = []
    for start in range(start_min, start_max + 1):
        for length in range(length_min, length_max + 1):
            end = start + length
            if end > len(sequence):
                continue
            motif = sequence[start:end]
            tm_c = _cached_tm(motif)
            if tm_c < tm_min:
                continue
            if _circular_occurrences(sequence, motif) != 1:
                continue
            values.append(
                Junction(
                    name="homology",
                    start=start,
                    end=end,
                    sequence=motif,
                    tm_c=tm_c,
                    purpose="gibson",
                )
            )
    return values


def _choose_vector_split(
    plasmid: str,
    insert_start: int,
    insert_end: int,
    rules: dict,
) -> Junction | None:
    thermo = rules["thermodynamics"]
    candidates = _homology_windows(
        plasmid,
        200,
        max(201, insert_start - int(rules["assembly"]["vector_split_min_distance_from_insert"])),
        int(thermo["gibson_length_min"]),
        int(thermo["gibson_length_max"]),
        float(thermo["gibson_tm_min_c"]),
    )
    if not candidates:
        return None
    target_tm = float(thermo["annealing_tm_target_c"])
    target_position = insert_start // 2
    return min(
        candidates,
        key=lambda j: (
            (j.tm_c - target_tm) ** 2
            + 0.000001 * (j.start - target_position) ** 2
            + 0.01 * (j.length - 26) ** 2
            + 8.0
            * len(_oligo_risks(j.sequence, rules))
        ),
    )


def _junction_pair_candidates(
    plasmid: str,
    poi_start: int,
    poi_end: int,
    rules: dict,
) -> list[tuple[Junction, Junction, float]]:
    thermo = rules["thermodynamics"]
    assembly = rules["assembly"]
    length_min = int(thermo["gibson_length_min"])
    length_max = int(thermo["gibson_length_max"])
    tm_min = float(thermo["gibson_tm_min_c"])
    upstream: list[Junction] = []
    for end in range(
        int(assembly["insert_upstream_scan_start"]) + length_min,
        min(poi_start, int(assembly["insert_upstream_scan_end"])) + 1,
    ):
        for length in range(length_min, length_max + 1):
            start = end - length
            if start < int(assembly["insert_upstream_scan_start"]):
                continue
            motif = plasmid[start:end]
            tm_c = _cached_tm(motif)
            if tm_c >= tm_min and _circular_occurrences(plasmid, motif) == 1:
                upstream.append(
                    Junction(
                        name="vector1-insert",
                        start=start,
                        end=end,
                        sequence=motif,
                        tm_c=tm_c,
                        purpose="gibson",
                    )
                )

    downstream: list[Junction] = []
    downstream_start_limit = poi_end + int(assembly["insert_downstream_extra_min"])
    downstream_end_limit = min(
        len(plasmid), poi_end + int(assembly["insert_downstream_extra_max"])
    )
    for start in range(poi_end, downstream_end_limit - length_min + 1):
        for length in range(length_min, length_max + 1):
            end = start + length
            if end < downstream_start_limit or end > downstream_end_limit:
                continue
            motif = plasmid[start:end]
            tm_c = _cached_tm(motif)
            if tm_c >= tm_min and _circular_occurrences(plasmid, motif) == 1:
                downstream.append(
                    Junction(
                        name="insert-vector2",
                        start=start,
                        end=end,
                        sequence=motif,
                        tm_c=tm_c,
                        purpose="gibson",
                    )
                )

    pairs: list[tuple[Junction, Junction, float]] = []
    target = float(thermo["annealing_tm_target_c"])
    preferred_right_end = poi_end + int(assembly["insert_downstream_extra_min"])
    for left in upstream:
        if left.end > poi_start:
            continue
        for right in downstream:
            if right.start < poi_end:
                continue
            insert_length = right.end - left.start
            score = (left.tm_c - target) ** 2 + (right.tm_c - target) ** 2
            score += 0.02 * (left.length - 26) ** 2
            score += 0.02 * (right.length - 26) ** 2
            # Prefer compact inserts close to the expression cassette while
            # still allowing the sequence to move a boundary when that
            # materially improves Tm or oligo manufacturability.
            score += 0.002 * (poi_start - left.end) ** 2
            score += 0.002 * (right.end - preferred_right_end) ** 2
            score += 0.002 * insert_length
            score += 8.0 * len(
                _oligo_risks(reverse_complement(left.sequence), rules)
            )
            score += 8.0 * len(_oligo_risks(right.sequence, rules))
            pairs.append((left, right, score))
    pairs.sort(key=lambda item: item[2])
    return pairs[:6]


def plan_assembly(
    design_id: str,
    protein_sequence: str,
    cds: str,
    final_plasmid: str,
    poi_start: int,
    intrinsic_score: float,
    rules: dict,
    primer_design_id: str | None = None,
) -> AssemblyPlan | None:
    poi_end = poi_start + len(cds)
    best: AssemblyPlan | None = None
    for left, right, junction_score in _junction_pair_candidates(
        final_plasmid, poi_start, poi_end, rules
    ):
        insert_sequence = final_plasmid[left.start : right.end]
        planned = plan_synthetic_fragments(
            insert_sequence,
            design_id,
            rules,
            primer_design_id=primer_design_id,
        )
        if not planned:
            continue
        fragments, fragment_junctions, fragment_score = planned
        split = _choose_vector_split(
            final_plasmid, left.start, right.end, rules
        )
        if split is None:
            continue
        vector_oligos = [
            Oligo(
                name="vec-1-FW",
                sequence=split.sequence,
                purpose="vector_fragment_1",
                orientation="FW",
                start=split.start,
                end=split.end,
                effective_sequence=split.sequence,
                effective_tm_c=split.tm_c,
                full_tm_c=split.tm_c,
                risk_flags=_oligo_risks(split.sequence, rules),
            ),
            Oligo(
                name="vec-1-RV",
                sequence=reverse_complement(left.sequence),
                purpose="vector_fragment_1",
                orientation="RV",
                start=left.start,
                end=left.end,
                effective_sequence=reverse_complement(left.sequence),
                effective_tm_c=left.tm_c,
                full_tm_c=left.tm_c,
                risk_flags=_oligo_risks(reverse_complement(left.sequence), rules),
            ),
            Oligo(
                name="vec-2-FW",
                sequence=right.sequence,
                purpose="vector_fragment_2",
                orientation="FW",
                start=right.start,
                end=right.end,
                effective_sequence=right.sequence,
                effective_tm_c=right.tm_c,
                full_tm_c=right.tm_c,
                risk_flags=_oligo_risks(right.sequence, rules),
            ),
            Oligo(
                name="vec-2-RV",
                sequence=reverse_complement(split.sequence),
                purpose="vector_fragment_2",
                orientation="RV",
                start=split.start,
                end=split.end,
                effective_sequence=reverse_complement(split.sequence),
                effective_tm_c=split.tm_c,
                full_tm_c=split.tm_c,
                risk_flags=_oligo_risks(reverse_complement(split.sequence), rules),
            ),
        ]
        gibson = [
            left,
            right,
            replace(split, name="vector2-vector1"),
        ]
        total = intrinsic_score + junction_score + fragment_score
        total += 3.0 * sum(bool(oligo.risk_flags) for oligo in vector_oligos)
        plan = AssemblyPlan(
            design_id=design_id,
            protein_sequence=protein_sequence,
            cds=cds,
            final_plasmid_sequence=final_plasmid,
            insert_start=left.start,
            insert_end=right.end,
            insert_sequence=insert_sequence,
            synthetic_fragments=fragments,
            fragment_junctions=fragment_junctions,
            gibson_junctions=gibson,
            vector_oligos=vector_oligos,
            score=total,
        )
        if best is None or plan.score < best.score:
            best = plan
        # Junction pairs are pre-ranked by Tm, length, and compactness.  The
        # first feasible pair avoids spending minutes optimizing negligible
        # boundary differences across a batch.
        break
    return best
