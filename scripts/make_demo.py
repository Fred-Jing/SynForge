#!/usr/bin/env python3
"""Generate deterministic synthetic inputs. These are not experimental results."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import random
import struct
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from Bio import SeqIO
from Bio.Seq import Seq
from primer_build.pca import Rules, design_layout


def write_ab1(path, sequence):
    """Minimal processed ABIF 1.01 with ideal synthetic four-channel peaks."""
    positions = [8*i+4 for i in range(len(sequence))]
    size = 8*len(sequence)+8
    entries = [(b'PBAS', 2, 2, 1, len(sequence), sequence.encode()),
               (b'PCON', 2, 2, 1, len(sequence), bytes([40]*len(sequence))),
               (b'PLOC', 2, 4, 2, len(sequence), struct.pack('>'+str(len(sequence))+'h', *positions)),
               (b'FWO_', 1, 2, 1, 4, b'GATC')]
    for i, base in enumerate('GATC', 9):
        signal = [0]*size
        for pos, called in zip(positions, sequence):
            for dx, weight in ((-2,100),(-1,500),(0,1000),(1,500),(2,100)):
                signal[pos+dx] = weight if called == base else 2
        entries.append((b'DATA', i, 4, 2, size, struct.pack('>'+str(size)+'h', *signal)))
    offset = 34 + 28*len(entries)
    directories = []
    payload = bytearray()
    for tag, number, kind, element_size, count, data in entries:
        if len(data) <= 4:
            location = int.from_bytes(data.ljust(4,b'\0'),'big')
        else:
            location = offset+len(payload)
            payload.extend(data)
        directories.append(struct.pack('>4sI2H4I', tag, number, kind, element_size, count, len(data), location, 0))
    header = b'ABIF'+struct.pack('>H4sI2H3I',101,b'tdir',1,1023,28,len(entries),28*len(entries),34)+struct.pack('>I',0)
    path.write_bytes(header+b''.join(directories)+payload)
    record = SeqIO.read(path, 'abi')
    if str(record.seq) != sequence or record.letter_annotations['phred_quality'] != [40]*len(sequence):
        raise ValueError('Synthetic ABIF roundtrip failed')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True, help='New output directory')
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Output directory already exists; choose a new directory')
    args.output.mkdir(parents=True)
    rng = random.Random(20261019)
    reference = ''.join(rng.choices('ACGT', k=1000))
    edit0 = 400
    alt = next(b for b in 'ACGT' if b != reference[edit0])
    expected = reference[:edit0]+alt+reference[edit0+1:]
    (args.output/'reference.fasta').write_text('>synthetic_circular_reference\n'+reference+'\n')
    reads = args.output/'reads'
    reads.mkdir()
    # Three independent synthetic passing clones ensure there is no top-two quota.
    for colony in range(1,5):
        template = expected if colony <= 3 else reference
        for primer,seq in [('F',template[120:820]),('R',str(Seq(template[180:880]).reverse_complement()))]:
            write_ab1(reads/f'DEMO-1-{colony}_{primer}.ab1',seq)
    config = dict(schema_version=1,batch_name='Synthetic demo — not wet-lab evidence',
        input_directory='reads',sample_regex=r'(?P<sample>DEMO-(?P<construct>\d+)-(?P<colony>\d+))_(?P<primer>F|R)\.ab1',
        primers={'F':reference[100:120],'R':str(Seq(reference[880:900]).reverse_complement())},
        constructs={'1':dict(reference='reference.fasta',edits=[dict(position1=401,ref=reference[400],alt=alt)],
                            required_regions=[dict(start1=251,end1=750)])},
        expected_samples=[f'DEMO-1-{i}' for i in range(1,6)])
    (args.output/'batch.json').write_text(json.dumps(config,indent=2)+'\n')
    # Random noncoding DNA exercises the layout engine without real protein designs.
    rng = random.Random(2)
    target = ''.join(rng.choices('ACGT',k=150))
    design_layout(target,Rules())
    (args.output/'pca_target.fasta').write_text('>synthetic_PCA_target\n'+target+'\n')
    # Sequence-only mmCIF fixture: no structural prediction or coordinates are claimed.
    protein = 'MVTETATVQLPSPEVQQLLQAGVEAGDNVGLLLR'
    names = dict(zip('ACDEFGHIKLMNPQRSTVWY','ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR'.split()))
    structures = args.output/'design/structures'
    structures.mkdir(parents=True)
    cif = 'data_synthetic\nloop_\n_entity_poly_seq.entity_id\n_entity_poly_seq.num\n_entity_poly_seq.mon_id\n'
    cif += ''.join(f'1 {i} {names[aa]}\n' for i,aa in enumerate(protein,1))
    cif += '#\nloop_\n_entity_poly.entity_id\n_entity_poly.pdbx_strand_id\n1 A\n#\n'
    (structures/'synthetic.cif').write_text(cif)
    (args.output/'PROVENANCE.json').write_text(json.dumps(dict(synthetic=True,seed=20261019,
        expected_pass=['DEMO-1-1','DEMO-1-2','DEMO-1-3'],expected_reject=['DEMO-1-4'],
        expected_missing=['DEMO-1-5'],description='Ideal artificial ABIF traces; not an estimate of empirical accuracy'),indent=2)+'\n')
    print(f'Synthetic demo -> {args.output}')


if __name__ == '__main__':
    main()
