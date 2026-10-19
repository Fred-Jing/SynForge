"""Independent reference checks for optimized hot paths."""
import itertools
import random
import sys
import unittest
from pathlib import Path
from Bio.Seq import Seq

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from primer_build.planner import _circular_occurrences, _valid_annealing_segments, _oligo_risks
from primer_build.thermo import reverse_complement, snapgene_tm, longest_homopolymer, three_prime_complementarity
from primer_build.io import load_rules
from primer_build.expression import _rna_mfe, _rna_mfe_and_unpaired

class PerformanceEquivalenceTests(unittest.TestCase):
    def test_circular_search_overlaps_and_origin(self):
        sequences = ['', 'AAAA', 'ACGTAC', 'GTACG', 'aaccggtt']
        motifs = ['', 'A', 'AAA', 'ACG', 'CGTA', 'AACCGGTTA', 'TTA']
        for sequence,motif in itertools.product(sequences,motifs):
            extended=sequence.upper()+sequence.upper()[:len(motif)-1]
            expected=sum(extended.startswith(motif.upper(),i) for i in range(len(sequence)))
            self.assertEqual(_circular_occurrences(sequence,motif),expected,(sequence,motif))

    def test_reverse_complement_matches_biopython(self):
        rng=random.Random(12)
        for length in (0,1,20,35,59,89,240):
            seq=''.join(rng.choices('ACGT',k=length))
            self.assertEqual(reverse_complement(seq),str(Seq(seq).reverse_complement()))
            self.assertEqual(reverse_complement(seq.lower()+'\n'),str(Seq(seq).reverse_complement()))
        with self.assertRaises(ValueError):
            reverse_complement('ACGN')

    def test_cached_segments_and_risks_respect_rule_changes(self):
        rules=load_rules(ROOT/'config/design_rules.json')
        sequence='GCGATCGTACGATGCGATGCTAGCTAGCATCGATCGACTAGCTAGCGATCGATGAAAAAA'
        for low,high,limit in [(62,66,5),(50,60,8),(68,72,4)]:
            t=rules['thermodynamics'];t['annealing_tm_min_c']=low;t['annealing_tm_max_c']=high
            expected=[]
            for length in range(t['annealing_length_hard_min'],min(t['annealing_length_hard_max'],len(sequence))+1):
                tm=snapgene_tm(sequence[:length])
                if low<=tm<=high:expected.append((length,tm))
            self.assertEqual(_valid_annealing_segments(sequence,rules),tuple(expected))
            rules['oligos']['homopolymer_max']=limit
            flags=[]
            if longest_homopolymer(sequence)>limit:flags.append('homopolymer')
            if three_prime_complementarity(sequence,sequence)>=5:flags.append('3prime_self_complementarity')
            self.assertEqual(_oligo_risks(sequence,rules),tuple(flags))

    def test_rna_cache_preserves_result_and_temperature(self):
        seq='ACGCGCAUACGCGCAUACGCGC'
        for fn in (_rna_mfe,_rna_mfe_and_unpaired):
            fn.cache_clear()
            for temperature in (25.0,37.0):
                self.assertEqual(fn(seq,temperature),fn.__wrapped__(seq,temperature))
                self.assertEqual(fn(seq,temperature),fn.__wrapped__(seq,temperature))
            self.assertEqual(fn.cache_info().hits,2)

if __name__=='__main__':unittest.main()
