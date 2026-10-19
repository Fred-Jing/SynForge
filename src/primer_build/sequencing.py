"""Local circular-consensus comparison and conservative BAM evidence checks."""
from __future__ import annotations

from collections import Counter
import statistics

from Bio import Align
from Bio.Seq import Seq
import pysam


def aligner():
    return Align.PairwiseAligner(mode="global", match_score=2, mismatch_score=-4,
                                open_gap_score=-7, extend_gap_score=-1)


def query_boundary(alignment, target_boundary: int) -> int:
    """Map a target boundary to query; deleted positions map to the gap boundary."""
    for (a, c), (b, d) in zip(alignment.coordinates.T, alignment.coordinates.T[1:]):
        if int(a) <= target_boundary < int(b):
            return int(c) + (target_boundary - int(a) if d > c else 0)
    if target_boundary == int(alignment.coordinates[0, -1]):
        return int(alignment.coordinates[1, -1])
    raise ValueError("Target boundary outside alignment")


def normalize_circular(reference: str, query: str):
    """Return query in reference orientation/origin; never use length alone as QC.

    A unique exact backbone anchor chooses orientation and rotation. When an indel
    occurs before that anchor, align the rotated molecules to recover base 1.
    Failure to find a unique anchor is explicit, rather than a guessed mapping.
    """
    reference, query = reference.upper(), query.upper()
    if not reference or not query:
        raise ValueError("Empty sequence")
    width = min(80, max(12, len(reference) // 4))
    for strand, oriented in [("+", query), ("-", str(Seq(query).reverse_complement()))]:
        for start in range(0, len(reference) - width + 1, width):
            anchor = reference[start:start + width]
            if set(anchor) - set("ACGT") or reference.count(anchor) != 1:
                continue
            doubled = oriented + oriented
            pos = doubled.find(anchor)
            if pos < 0 or pos >= len(oriented) or doubled.find(anchor, pos + 1) < len(oriented) and doubled.find(anchor, pos + 1) >= 0:
                continue
            qr = oriented[pos:] + oriented[:pos]
            rr = reference[start:] + reference[:start]
            if start:
                boundary = len(reference) - start
                q_origin = boundary if qr == rr else query_boundary(aligner().align(rr, qr)[0], boundary)
            else:
                q_origin = 0
            offset = (pos + q_origin) % len(oriented)
            return oriented[offset:] + oriented[:offset], strand, offset
    raise ValueError("No unique exact circular anchor; manual/reference review required")


def compare_sequences(reference: str, query: str):
    """Global alignment, canonical reference coordinates and minimal differences."""
    if reference == query:
        return None, [], dict(enumerate(range(len(reference))))
    al = aligner().align(reference, query)[0]
    changes = []
    for (a, c), (b, d) in zip(al.coordinates.T, al.coordinates.T[1:]):
        a, b, c, d = map(int, (a, b, c, d))
        if b - a == d - c:
            for i, (r, q) in enumerate(zip(reference[a:b], query[c:d])):
                if r != q:
                    changes.append(dict(kind="SNV", start0=a+i, end0=a+i+1,
                                        query_start0=c+i, query_end0=c+i+1, ref=r, alt=q))
        else:
            # Repeat-associated indels have equivalent positions. Left-normalize.
            r, q = reference[a:b], query[c:d]
            allele = r or q
            while a > 0 and c > 0 and reference[a-1] == query[c-1] == allele[-1]:
                allele = allele[-1] + allele[:-1]
                a, b, c, d = a-1, b-1, c-1, d-1
            changes.append(dict(kind="DEL" if r else "INS", start0=a, end0=b,
                                query_start0=c, query_end0=d,
                                ref=allele if r else "", alt=allele if q else ""))
    # Base mapping uses the original valid alignment, independently of display normalization.
    mapping = {}
    for (a, c), (b, d) in zip(al.coordinates.T, al.coordinates.T[1:]):
        if b > a and d > c:
            mapping.update((int(c+i), int(a+i)) for i in range(int(b-a)))
    return al, changes, mapping


def keep_read(read):
    return not (read.is_unmapped or read.is_secondary or read.is_supplementary
                or read.is_duplicate or read.is_qcfail) and read.mapping_quality >= 20


def original_position(normalized_position, length, strand, offset):
    oriented = (normalized_position + offset) % length
    return oriented if strand == "+" else length - 1 - oriented


def bam_qc(bam, consensus, strand, offset, core_start, core_end):
    if len(bam.references) != 1 or bam.lengths[0] != len(consensus):
        raise ValueError("BAM reference length does not match supplied consensus")
    counts = bam.count_coverage(bam.references[0], 0, len(consensus),
                               quality_threshold=15, read_callback=keep_read)
    depths = list(map(sum, zip(*counts)))
    # Supplementary records are excluded, so depth is deliberately conservative
    # at a circular alignment split. All samples in this batch still exceed 20x.
    normalized = [depths[original_position(i, len(consensus), strand, offset)]
                  for i in range(len(consensus))]
    snv_candidates = []
    for pos, base in enumerate(consensus):
        if base not in "ACGT" or depths[pos] < 20:
            continue
        alt_index = max((i for i in range(4) if "ACGT"[i] != base),
                        key=lambda i: counts[i][pos])
        alt_count = counts[alt_index][pos]
        if alt_count >= 5 and alt_count / depths[pos] >= .10:
            alt = "ACGT"[alt_index]
            fwd, rev = 0, 0
            for col in bam.pileup(bam.references[0], pos, pos+1, truncate=True,
                                  stepper="nofilter", min_base_quality=0,
                                  max_depth=100000, ignore_overlaps=False):
                for pr in col.pileups:
                    read, qp = pr.alignment, pr.query_position
                    if keep_read(read) and qp is not None and read.query_qualities[qp] >= 15 and read.query_sequence[qp] == alt:
                        if read.is_reverse: rev += 1
                        else: fwd += 1
            if min(fwd, rev) >= 3:
                snv_candidates.append(dict(consensus_pos1=pos+1, ref=base, alt=alt,
                                           depth=depths[pos], frequency=alt_count/depths[pos],
                                           alt_forward=fwd, alt_reverse=rev))
    vals = sorted(normalized)
    qc = dict(depth_min=min(vals), depth_median=statistics.median(vals),
              depth_q05=vals[int(.05*len(vals))], bases_below20=sum(d<20 for d in vals),
              core_depth_min=min(normalized[core_start:core_end]),
              primary_reads=sum(1 for r in bam.fetch() if keep_read(r)),
              independent_snv_flags=snv_candidates)
    return qc, normalized


def haplotype_support(bam, consensus, normalized_query, expected_window,
                      query_start, query_end, strand, offset):
    """Rank reads between two local alleles, allowing long-read sequencing errors.

    This is evidence for a consensus difference, not a calibrated allele fraction.
    Exact agreement is also counted. Windows crossing the artificial BAM origin
    are explicitly unavailable. Supplementary records never inflate support.
    """
    if query_start < 0 or query_end > len(normalized_query):
        return {"unavailable": "variant too close to canonical origin"}
    original = [original_position(i, len(consensus), strand, offset)
                for i in (query_start, query_end-1)]
    lo, hi = min(original), max(original)+1
    if hi-lo != query_end-query_start:
        return {"unavailable": "window crosses artificial BAM origin"}
    observed_window = normalized_query[query_start:query_end]
    if strand == "-":
        observed_window = str(Seq(observed_window).reverse_complement())
        expected_window = str(Seq(expected_window).reverse_complement())
    scores = Counter()
    identities = set()
    scorer = aligner()
    for read in bam.fetch(bam.references[0], lo, hi):
        if not keep_read(read) or read.query_name in identities:
            continue
        pairs = {r: q for q, r in read.get_aligned_pairs() if r in (lo, hi-1) and q is not None}
        if len(pairs) < 2:
            continue
        identities.add(read.query_name)
        a, b = pairs[lo], pairs[hi-1]+1
        if not read.query_qualities or min(read.query_qualities[a:b], default=0) < 5:
            continue
        seq = read.query_sequence[a:b]
        diff = scorer.score(observed_window, seq) - scorer.score(expected_window, seq)
        key = "consensus_preferred" if diff > 0 else "target_preferred" if diff < 0 else "tie"
        scores[key] += 1
        scores["consensus_exact"] += seq == observed_window
        scores["target_exact"] += seq == expected_window
    scores["informative_reads"] = scores["consensus_preferred"] + scores["target_preferred"]
    return dict(scores)
