from copy import deepcopy
from pathlib import Path
import random
import sys
import tempfile
import unittest

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from primer_build.sanger import (DEFAULTS, aggregate_clone, build_reference,
    directed_window, mott_bounds, primer_site, natural_key, trace_evidence, valid_trace, align_read)


class SangerValidationTests(unittest.TestCase):
    def make_inputs(self):
        reference=dict(expected='ACGT'*10,required=set(range(10,30)),edits=[dict(position1=15,ref='A',alt='G')])
        reads=[]
        for name,strand in [('forward',1),('reverse',-1)]:
            reads.append(dict(name=name,summary=dict(trace_available=True,seq_sidecar_matches=True),gaps=[],
                observations={p:dict(base=b,quality=40,read_index0=p,strand=strand,peak=None,mixed=False)
                              for p,b in enumerate(reference['expected'])}))
        return reference,reads

    def call(self,reference,reads):
        return aggregate_clone('sample','construct','colony',reads,reference,DEFAULTS)

    def test_trim_retains_raw_trace_indices(self):
        self.assertEqual(mott_bounds([0]*5+[40]*200+[0]*5),(5,205))
        self.assertEqual(mott_bounds([0]*30),(0,0))

    def test_primer_direction_and_origin_wrap(self):
        rng=random.Random(37)
        sequence=''.join(rng.choices('ACGT',k=160))
        for strand,primer in [(1,sequence[145:]+sequence[:5]),(-1,str(Seq(sequence[5:25]).reverse_complement()))]:
            site=primer_site(sequence,primer)
            self.assertEqual(site[2],strand)
            window,positions=directed_window(sequence,site,80)
            self.assertEqual(window,''.join(sequence[p] if strand==1 else str(Seq(sequence[p]).complement()) for p in positions))
            self.assertEqual(len(set(positions)),80)
        with self.assertRaises(ValueError):primer_site('A'*100,'A'*20)

    def test_reference_edit_is_explicit_and_source_is_read_only(self):
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/'ref.fasta';path.write_text('>ref\nACGTACGT\n');original=path.read_bytes()
            spec=dict(reference='ref.fasta',edits=[dict(position1=3,ref='G',alt='C')],required_regions=[dict(start1=1,end1=8)])
            ref=build_reference(spec,Path(d))
            self.assertEqual(ref['expected'],'ACCTACGT')
            self.assertEqual(path.read_bytes(),original)
            spec['edits'][0]['ref']='T'
            with self.assertRaises(ValueError):build_reference(spec,Path(d))

    def test_clean_and_complementary_coverage_are_distinct(self):
        reference,reads=self.make_inputs()
        self.assertEqual(self.call(reference,reads)['status'],'PASS_BIDIRECTIONAL')
        del reads[0]['observations'][11]
        c=self.call(reference,reads)
        self.assertEqual(c['status'],'PASS_COMPLEMENTARY')
        self.assertEqual(c['covered_expected_bases'],20)
        self.assertEqual(c['bidirectional_clean_bases'],19)

    def test_no_reference_filling_for_missing_or_low_quality_bases(self):
        reference,reads=self.make_inputs()
        for read in reads:read['observations'][11]['quality']=19
        c=self.call(reference,reads)
        self.assertEqual(c['status'],'INCOMPLETE_COVERAGE')
        self.assertIn(12,c['missing_positions1'])

    def test_mutation_needs_two_high_quality_directions(self):
        reference,reads=self.make_inputs()
        reads[0]['observations'][14]['quality']=21
        c=self.call(reference,reads)
        self.assertEqual(c['status'],'SINGLE_STRAND_OR_LOW_Q')
        self.assertFalse(c['edit_results'][0]['bidirectional_hq'])

    def test_one_high_quality_conflict_is_not_outvoted(self):
        reference,reads=self.make_inputs()
        reads[0]['observations'][14]['base']='A'
        self.assertEqual(self.call(reference,reads)['status'],'REVIEW_TRACE_OR_DISCORDANCE')
        reads[1]['observations'][14]['base']='A'
        self.assertEqual(self.call(reference,reads)['status'],'REJECT_DIFFERENCE')

    def test_unilateral_peak_is_a_flag_not_mixture_fraction(self):
        reference,reads=self.make_inputs()
        for index in [0,1]:
            obs=reads[index]['observations'][11]
            obs.update(mixed=True,peak=dict(secondary_base='A',secondary_ratio=0.4))
            c=self.call(reference,reads)
            self.assertEqual(c['status'],'PASS_BIDIRECTIONAL' if index==0 else 'REVIEW_TRACE_OR_DISCORDANCE')
        self.assertEqual(len(c['trace_flags']),2)

    def test_high_quality_deletion_is_not_called_missing_only(self):
        reference,reads=self.make_inputs()
        for read in reads:
            del read['observations'][18]
            read['gaps'].append(dict(kind='DEL',positions1=[19],position1=19,ref='G',alt='',quality=40,
                strand=1 if read['name']=='forward' else -1,read_index0=18))
        c=self.call(reference,reads)
        self.assertEqual(c['status'],'REJECT_DIFFERENCE')
        self.assertTrue(all(g['confirmed_difference'] for g in c['issues']))

    def test_external_discrepancy_is_visible_and_not_prioritized(self):
        reference,reads=self.make_inputs()
        reads[0]['observations'][35]['base']='A'
        c=self.call(reference,reads)
        self.assertEqual(c['region_status'],'PASS_BIDIRECTIONAL')
        self.assertEqual(c['status'],'REVIEW_OUTSIDE_REGION')

    def test_missing_direction_and_sidecar_disagreement_do_not_pass(self):
        reference,reads=self.make_inputs()
        self.assertEqual(self.call(reference,reads[:1])['status'],'SINGLE_STRAND_OR_LOW_Q')
        reads[0]['summary']['seq_sidecar_matches']=False
        self.assertEqual(self.call(reference,reads)['status'],'REVIEW_INPUT_OR_TRACE')

    def test_channel_order_and_mixed_peak_values(self):
        raw={'FWO_1':b'GATC','PLOC2':[1],
             'DATA9':[0,25,0],'DATA10':[0,5,0],
             'DATA11':[0,10,0],'DATA12':[0,100,0]}
        peak=trace_evidence(raw,0,'C',DEFAULTS)
        self.assertEqual(peak['secondary_base'],'G')
        self.assertEqual(peak['secondary_ratio'],0.25)

    def test_mixed_numeric_construct_names_sort_without_type_error(self):
        names=sorted(['A10','12','A2','3'],key=natural_key)
        self.assertLess(names.index('3'),names.index('12'))
        self.assertLess(names.index('A2'),names.index('A10'))

    def test_long_five_prime_quality_trim_does_not_truncate_reference_window(self):
        rng=random.Random(517)
        sequence=''.join(rng.choices('ACGT',k=1200))
        record=SeqRecord(Seq(sequence[120:900]))
        record.letter_annotations['phred_quality']=[0]*300+[40]*480
        record.annotations['abif_raw']={}
        with tempfile.TemporaryDirectory() as d, patch('primer_build.sanger.SeqIO.read',return_value=record):
            read=align_read(Path(d)/'read.ab1',dict(expected=sequence),sequence[100:120],DEFAULTS)
        self.assertEqual(read['summary']['trim_start0'],300)
        self.assertEqual(set(read['observations']),set(range(420,900)))

    def test_high_quality_external_indel_requires_review_or_rejection(self):
        reference,reads=self.make_inputs()
        for i,read in enumerate(reads):
            del read['observations'][35]
            read['gaps'].append(dict(kind='DEL',positions1=[36],position1=36,ref='T',alt='',quality=40,
                strand=1 if read['name']=='forward' else -1,read_index0=35))
            c=self.call(reference,reads)
            self.assertEqual(c['status'],'REVIEW_OUTSIDE_REGION' if i==0 else 'REJECT_DIFFERENCE')

    def test_trace_metadata_validation_allows_duplicate_low_quality_tail_positions(self):
        raw={'FWO_1':b'GATC','PLOC2':[1,2,2]}
        raw.update({f'DATA{i}':[0]*5 for i in range(9,13)})
        self.assertTrue(valid_trace(raw,3))
        self.assertFalse(valid_trace(raw,4))
        raw['DATA12']=[0]*2
        self.assertFalse(valid_trace(raw,3))


if __name__=='__main__':unittest.main()
