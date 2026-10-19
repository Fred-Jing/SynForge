#!/usr/bin/env python3
"""Design a bounded, screened alternating PCA layout for one exact DNA target."""
from __future__ import annotations

import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import re
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from Bio import SeqIO
from primer_build.pca import Conditions, Rules, design_layout, tm


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--target', type=Path, required=True, help='Exactly one FASTA record, A/C/G/T only')
    parser.add_argument('--output', type=Path, required=True, help='New JSON file; existing files are never overwritten')
    parser.add_argument('--max-oligo', type=int, default=59, help='20–60 nt; Primer3 structural screens are limited to 60 nt')
    parser.add_argument('--min-overlap', type=int, default=20)
    parser.add_argument('--max-overlap', type=int, default=40)
    parser.add_argument('--min-tm', type=float, default=60)
    parser.add_argument('--max-tm', type=float, default=64)
    parser.add_argument('--monovalent-mm', type=float, default=50)
    parser.add_argument('--divalent-mm', type=float, default=1.5)
    parser.add_argument('--dntp-mm', type=float, default=0.6)
    parser.add_argument('--dna-nm', type=float, default=50)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Output already exists; choose a new file')
    if not 20 <= args.max_oligo <= 60:
        parser.error('--max-oligo must be 20–60 for the current structural screening engine')
    if not 12 <= args.min_overlap <= args.max_overlap < args.max_oligo:
        parser.error('Require 12 <= min-overlap <= max-overlap < max-oligo')
    if not 0 < args.min_tm <= args.max_tm < 100:
        parser.error('Require 0 < min-tm <= max-tm < 100')
    concentrations = (args.monovalent_mm, args.divalent_mm, args.dntp_mm, args.dna_nm)
    if any(x < 0 for x in concentrations) or args.dna_nm <= 0:
        parser.error('Concentrations must be nonnegative; DNA concentration must be positive')
    record = SeqIO.read(args.target, 'fasta')
    target = str(record.seq).upper()
    if re.fullmatch('[ACGT]+', target) is None or len(target) <= args.max_oligo:
        parser.error('Target must contain only A/C/G/T and be longer than max-oligo')
    rules = Rules(max_oligo=args.max_oligo, min_overlap=args.min_overlap,
                  max_overlap=args.max_overlap, min_tm=args.min_tm, max_tm=args.max_tm,
                  conditions=Conditions(*concentrations))
    layout = design_layout(target, rules)
    rows = [dict(name=f'{record.id}-{("F" if i % 2 == 0 else "R")}{i // 2 + 1}',
                 sequence_5to3=seq, length=len(seq), start0=start, end0=end,
                 strand='+' if i % 2 == 0 else '-')
            for i, (seq, (start, end)) in enumerate(zip(layout['sequences'], layout['intervals']))]
    result = dict(schema_version=1, target_id=record.id, target_length=len(target),
                  target_sha256=hashlib.sha256(target.encode()).hexdigest(),
                  input_sha256=hashlib.sha256(args.target.read_bytes()).hexdigest(),
                  rules=asdict(rules), oligos=rows,
                  junctions=[dict(start0=s, end0=e, tm_c=tm(target[s:e], rules.conditions))
                             for s, e, _ in layout['junctions']],
                  search=layout, scope='PCA assembly oligos only; no synonymous changes, outer amplification primers or vector assembly design',
                  code_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest()
                               for p in (Path(__file__), ROOT/'src/primer_build/pca.py', ROOT/'uv.lock')})
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x') as handle:
        json.dump(result, handle, indent=2)
        handle.write('\n')
    print(f'{len(rows)} screened assembly oligos -> {args.output}')


if __name__ == '__main__':
    main()
