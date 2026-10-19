from pathlib import Path
import sys
import tempfile
import unittest

from Bio.Seq import Seq
import pysam

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from primer_build.sequencing import (normalize_circular, compare_sequences,
                                    original_position, bam_qc)


class SequencingTests(unittest.TestCase):
    REF = ("AGTCAGGACTTACCGATGCTAGATCGTTGACCTAGGTACGATCAGCTAGGTCAT"
           "CCGATTCAGATGCAATGTCCTAGGATCGATCTGAGGATCCATGCTTAGCAGTGC"
           "ATCAGATGGTCCATTGACCTCGATACGTGATCCTAGATTCGACATGCTTAGCAC")

    def test_rotation_reverse_complement_and_coordinate_mapping(self):
        for reverse in (False, True):
            q = self.REF[37:] + self.REF[:37]
            if reverse: q = str(Seq(q).reverse_complement())
            norm, strand, offset = normalize_circular(self.REF, q)
            self.assertEqual(norm, self.REF)
            for i in (0, 10, len(q)-1):
                base = q[original_position(i, len(q), strand, offset)]
                self.assertEqual(str(Seq(base).complement()) if reverse else base, norm[i])

    def test_first_base_mutation_and_indel_before_anchor(self):
        for q in ("C" + self.REF[1:], self.REF[:5] + "AAA" + self.REF[5:], self.REF[1:]):
            rotated = q[60:] + q[:60]
            norm, _, _ = normalize_circular(self.REF, rotated)
            self.assertEqual(norm, q)

    def test_snv_deletion_and_duplication(self):
        q = self.REF[:30]+"A"+self.REF[31:80]+self.REF[81:100]+self.REF[100:125]+self.REF[100:]
        _, variants, mapping = compare_sequences(self.REF, q)
        self.assertEqual([v["kind"] for v in variants], ["SNV", "DEL", "INS"])
        self.assertEqual(sum(len(v["alt"])-len(v["ref"]) for v in variants), len(q)-len(self.REF))
        self.assertEqual(len(variants[-1]["alt"]), 25)
        self.assertEqual(mapping[len(q)-1], len(self.REF)-1)

    def test_repeat_indel_is_left_normalized(self):
        _, variants, _ = compare_sequences("CTAAAAGC", "CTAAAGC")
        self.assertEqual(variants[0]["start0"], 2)

    def test_unrelated_sequence_requires_review(self):
        with self.assertRaises(ValueError): normalize_circular(self.REF, "A"*150)

    def test_ambiguous_base_is_a_difference(self):
        _, variants, _ = compare_sequences(self.REF, self.REF[:60]+"N"+self.REF[61:])
        self.assertEqual(variants[0]["alt"], "N")

    def test_bam_filters_strand_flags_and_depth(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d)/"reads.bam"
            ref = "A"*40
            header = {"HD":{"VN":"1.6","SO":"coordinate"}, "SQ":[{"SN":"sample","LN":40}]}
            with pysam.AlignmentFile(path, "wb", header=header) as f:
                for i, (flag, mapq, quality, base) in enumerate([(0,60,40,"A")]*20 + [(16,60,40,"C")]*4 + [(0,60,40,"C")]*4 + [(2048,60,40,"C"), (0,5,40,"C"), (0,60,5,"C")]):
                    r=pysam.AlignedSegment();r.query_name=f"read{i}";r.flag=flag;r.reference_id=0;r.reference_start=0;r.mapping_quality=mapq;r.cigarstring="40M";r.query_sequence=base*40;r.query_qualities=[quality]*40;f.write(r)
            pysam.index(str(path))
            with pysam.AlignmentFile(path,"rb") as f:
                qc, dp=bam_qc(f,ref,"+",0,0,40)
            self.assertEqual(min(dp),28)
            self.assertEqual(qc["primary_reads"],29)  # Low BQ read is kept at read level only.
            self.assertEqual(qc["independent_snv_flags"][0]["alt_reverse"],4)


if __name__ == "__main__": unittest.main()
