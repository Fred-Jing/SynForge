import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'src'))
from primer_build.pca import Conditions, Rules, tm, single_ok, nonadjacent_qc
from primer_build.thermo import reverse_complement


class PCAScreeningTests(unittest.TestCase):
    def test_concentration_is_part_of_cache_key(self):
        seq='ATGAGCGCCGAACTGGCAAC'
        low=tm(seq,Conditions(dna_conc=25))
        high=tm(seq,Conditions(dna_conc=100))
        self.assertGreater(high,low)
        self.assertEqual(tm(seq,Conditions(dna_conc=25)),low)

    def test_old_operator_oligo_fails_hairpin_gate(self):
        old='AGGGGAATTGTGAGCGGATAACAATTCCCCTCTAGAAATAATTTTGT'
        self.assertFalse(single_ok(old,Rules()))

    def test_intended_adjacent_pair_is_not_classified_as_off_target(self):
        seq='ATGAGCGCCGAACTGGCAAC'
        partner=reverse_complement(seq)
        self.assertTrue(nonadjacent_qc([seq,partner],Rules())[0])
        # The same complementary oligo at a nonadjacent position must fail.
        self.assertFalse(nonadjacent_qc([seq,'AATGTTAGCTATGACTTGCA',partner],Rules())[0])

    def test_overlength_oligo_fails_before_thermodynamic_call(self):
        self.assertFalse(single_ok('ACGT'*16,Rules()))


if __name__=='__main__':unittest.main()
