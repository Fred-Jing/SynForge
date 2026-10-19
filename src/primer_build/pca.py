"""Bounded, whole-target PCA layout search and explicit thermodynamic screening.

Coordinates are zero-based half-open. Screening conditions are a reference
model, not a claim about the composition of any proprietary PCR master mix.
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
import heapq
import math
import random
import primer3
from Bio.Seq import Seq

from .thermo import reverse_complement, longest_homopolymer
from .codon_opt import ECOLI_CODON_USAGE


@dataclass(frozen=True)
class Conditions:
    mv_conc: float = 50.0
    dv_conc: float = 1.5
    dntp_conc: float = 0.6
    dna_conc: float = 50.0

    def kwargs(self):
        return vars(self)


@dataclass(frozen=True)
class Rules:
    conditions: Conditions = Conditions()
    min_overlap: int = 20
    max_overlap: int = 40
    min_tm: float = 60.0
    max_tm: float = 64.0
    max_oligo: int = 59
    min_oligo: int = 0
    max_hairpin: float = 50.0
    max_dimer: float = 45.0
    max_end: float = 40.0
    beam_per_state: int = 12
    max_final_candidates: int = 400


@lru_cache(maxsize=100_000)
def tm(seq: str, c: Conditions = Conditions()) -> float:
    return primer3.calc_tm(seq, **c.kwargs())


@lru_cache(maxsize=100_000)
def single_qc(seq: str, c: Conditions = Conditions()) -> tuple[float, float, float]:
    return (primer3.calc_hairpin(seq, **c.kwargs()).tm,
            primer3.calc_homodimer(seq, **c.kwargs()).tm,
            primer3.calc_end_stability(seq, seq, **c.kwargs()).tm)


@lru_cache(maxsize=100_000)
def pair_qc(a: str, b: str, c: Conditions = Conditions()) -> tuple[float, float, float]:
    return (primer3.calc_heterodimer(a, b, **c.kwargs()).tm,
            primer3.calc_end_stability(a, b, **c.kwargs()).tm,
            primer3.calc_end_stability(b, a, **c.kwargs()).tm)


def single_ok(seq: str, rules: Rules) -> bool:
    if len(seq) > rules.max_oligo or len(seq) < rules.min_oligo or longest_homopolymer(seq) > 5:
        return False
    hp, hd, end = single_qc(seq, rules.conditions)
    return hp <= rules.max_hairpin and hd <= rules.max_dimer and end <= rules.max_end


def nonadjacent_qc(oligos: list[str], rules: Rules) -> tuple[bool, list[dict]]:
    results = []
    for i, a in enumerate(oligos):
        for j in range(i + 2, len(oligos)):
            h, e1, e2 = pair_qc(a, oligos[j], rules.conditions)
            results.append(dict(i=i, j=j, heterodimer_tm=h, end_i_tm=e1, end_j_tm=e2))
            if h > rules.max_dimer or max(e1, e2) > rules.max_end:
                return False, results
    return True, results


def design_layout(target: str, rules: Rules = Rules()) -> dict:
    """Search a junction DAG, retaining bounded alternate paths per state.

    Junctions cannot overlap each other. Thus oligos on the same strand have
    no overlapping target coordinates. Inter-junction gaps are allowed.
    Search is deterministic, heuristic (bounded beam), not proof of global
    optimality or of infeasibility when no path is found.
    """
    n = len(target)
    nodes = []
    for start in range(n - rules.min_overlap + 1):
        for end in range(start + rules.min_overlap, min(n, start + rules.max_overlap) + 1):
            temp = tm(target[start:end], rules.conditions)
            if rules.min_tm <= temp <= rules.max_tm:
                nodes.append((start, end, temp))
    nodes.sort()
    # A state records last junction and count parity; each path stores all oligos.
    states: dict[tuple[int, int], list[tuple]] = {}
    finals = []

    @lru_cache(maxsize=20000)
    def oligo(start, end, orientation):
        seq = target[start:end]
        if orientation == 1:
            seq = reverse_complement(seq)
        if not single_ok(seq, rules):
            return None
        # Exact 3' 12-mer uniqueness in the target, including its other strand.
        tail = seq[-12:]
        if target.count(tail) + reverse_complement(target).count(tail) != 1:
            return None
        hp, hd, et = single_qc(seq, rules.conditions)
        penalty = 30 + 0.002 * (hp ** 2 + hd ** 2) + .02 * (len(seq) - 52) ** 2
        return seq, penalty

    for idx, (s, e, temp) in enumerate(nodes):
        incoming = {0: [], 1: []}
        if e <= rules.max_oligo:
            candidate = oligo(0, e, 0)
            if candidate:
                seq, score = candidate
                incoming[1].append((score + (temp-62)**2, ((0, e),), (idx,)))
        for prev in range(idx-1, -1, -1):
            ps, pe, _ = nodes[prev]
            if e-ps > rules.max_oligo:
                break
            if pe > s:
                continue
            for parity in (0, 1):
                paths = states.get((prev, parity))
                if not paths:
                    continue
                candidate = oligo(ps, e, parity)
                if candidate is None:
                    continue
                seq, delta = candidate
                delta += (temp-62)**2
                for score, intervals, junctions in paths:
                    incoming[1-parity].append((score+delta, intervals+((ps,e),), junctions+(idx,)))
        for parity in (0, 1):
            chosen = heapq.nsmallest(rules.beam_per_state, incoming[parity])
            if chosen:
                states[idx, parity] = chosen
            # End on an R oligo, with no separate requirement to reuse termini.
            if parity == 1 and n-s <= rules.max_oligo and e < n:
                last = oligo(s, n, 1)
                if last:
                    for score, intervals, junctions in chosen:
                        finals.append((score+last[1], intervals+((s,n),), junctions))
    rejected = 0
    for score, intervals, junctions in heapq.nsmallest(rules.max_final_candidates, finals):
        seqs = [target[s:e] if i%2 == 0 else reverse_complement(target[s:e])
                for i, (s,e) in enumerate(intervals)]
        passed, pairs = nonadjacent_qc(seqs, rules)
        if not passed:
            rejected += 1
            continue
        return dict(score=score, intervals=intervals, sequences=seqs,
                    junctions=[nodes[j] for j in junctions], pairs=pairs,
                    nodes=len(nodes), final_candidates=len(finals), rejected_pair_sets=rejected)
    furthest = max((nodes[i][1] for i, p in states), default=0)
    raise ValueError(f"No accepted layout: nodes={len(nodes)}, reached={furthest}/{n}, finals={len(finals)}, rejected={rejected}")


def synonymous_candidates(original: str, coding_nt: int, seed: int, steps: int = 1600):
    """Deterministic, local synonymous search; yields intermediate best sequences.

    Freeze all original stop codons (including TAG), Met and the downstream
    cassette. Windows are a search surrogate; final oligos must pass real QC.
    """
    rng = random.Random(seed)
    codons = [original[i:i+3] for i in range(0, coding_nt, 3)]
    choices = {}
    for i, codon in enumerate(codons):
        aa = str(Seq(codon).translate())
        if aa in ('*', 'M', 'W'):
            continue
        choices[i] = [(s, f) for s, f in ECOLI_CODON_USAGE[aa] if f >= .10]
    starts = list(range(0, max(1, len(original)-59), 6))
    if len(original)-60 not in starts:
        starts.append(len(original)-60)
    @lru_cache(maxsize=30000)
    def window_score(s):
        hp = primer3.calc_hairpin(s, **Conditions().kwargs()).tm
        gc = (s.count('G')+s.count('C'))/len(s)
        return max(0, hp-45)**2 + 2000*max(0, gc-.65)**2 + 100*max(0,longest_homopolymer(s)-4)**2
    seq = original
    scores = [window_score(seq[s:s+60]) for s in starts]
    score = sum(scores)
    best_score, best = score, seq
    for step in range(steps):
        # Bias proposals towards problematic windows, but retain global exploration.
        wi = rng.choices(range(len(starts)), weights=[x+5 for x in scores])[0]
        eligible = [i for i in choices if starts[wi] <= i*3 < starts[wi]+60]
        if not eligible:
            eligible = list(choices)
        i = rng.choice(eligible)
        codon = rng.choices([x[0] for x in choices[i]], weights=[x[1] for x in choices[i]])[0]
        if seq[i*3:i*3+3] == codon:
            continue
        candidate = seq[:i*3]+codon+seq[i*3+3:]
        changed = [j for j,s in enumerate(starts) if s < i*3+3 and s+60 > i*3]
        updates = {j:window_score(candidate[starts[j]:starts[j]+60]) for j in changed}
        candidate_score = score + sum(updates[j]-scores[j] for j in changed)
        temperature = 25*(1-step/steps)+.1
        if candidate_score <= score or rng.random() < math.exp(min(0,(score-candidate_score)/temperature)):
            seq,score = candidate,candidate_score
            for j,v in updates.items(): scores[j]=v
            if score < best_score:
                best_score,best = score,seq
        if (step+1)%200 == 0:
            yield best, best_score, step+1
    yield best,best_score,steps
