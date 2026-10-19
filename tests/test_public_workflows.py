"""End-to-end synthetic examples; no laboratory files or account access required."""
from pathlib import Path
import csv
import hashlib
import json
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from Bio import SeqIO
from Bio.Seq import Seq


class PublicWorkflowTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temporary = tempfile.TemporaryDirectory()
        cls.demo = Path(cls.temporary.name)/'demo'
        cls.run_cli('make_demo.py','--output',cls.demo)

    @classmethod
    def tearDownClass(cls):
        cls.temporary.cleanup()

    @staticmethod
    def run_cli(script,*args):
        return subprocess.run([sys.executable,str(ROOT/'scripts'/script),*map(str,args)],
                              cwd=ROOT,check=True,capture_output=True,text=True)

    def test_sanger_all_passing_clones_no_quota_and_parent_rejected(self):
        originals={p:hashlib.sha256(p.read_bytes()).hexdigest() for p in self.demo.rglob('*.ab1')}
        self.run_cli('analyze_sanger_batch.py','--config',self.demo/'batch.json')
        out=self.demo/'results'
        rows=list(csv.DictReader((out/'usable_samples.csv').read_text(encoding='utf-8-sig').splitlines()))
        self.assertEqual({r['sample'] for r in rows},{f'DEMO-1-{i}' for i in (1,2,3)})
        samples=list(csv.DictReader((out/'samples.csv').read_text(encoding='utf-8-sig').splitlines()))
        self.assertEqual(next(r for r in samples if r['sample']=='DEMO-1-4')['status'],'REJECT_DIFFERENCE')
        manifest=json.loads((out/'manifest.json').read_text())
        self.assertEqual(manifest['missing_samples'],['DEMO-1-5'])
        self.assertIn('Synthetic demo',(out/'report.html').read_text())
        for p,h in originals.items():self.assertEqual(hashlib.sha256(p.read_bytes()).hexdigest(),h)
        with self.assertRaises(subprocess.CalledProcessError):
            self.run_cli('analyze_sanger_batch.py','--config',self.demo/'batch.json')

    def test_structure_design_replaces_placeholder_and_preserves_target(self):
        self.run_cli('design_project.py','--project',self.demo/'design','--chain','A',
                     '--map-start','DEMO-1','--primer-start','DEMO-1-1-1-FW',
                     '--candidates',2,'--workers',1)
        record=SeqIO.read(self.demo/'design/maps/PET-22B-DEMO-1.dna','snapgene')
        from primer_build.io import get_poi_feature
        feature=get_poi_feature(record)
        self.assertEqual(str(feature.extract(record.seq).translate()),'MVTETATVQLPSPEVQQLLQAGVEAGDNVGLLLR')
        self.assertTrue((self.demo/'design/all_primers.xlsx').is_file())
        manifest=json.loads((self.demo/'design/run_manifest.json').read_text())
        self.assertTrue(manifest)

    def test_pca_roundtrip_complete_target_and_screen_limits(self):
        out=self.demo/'pca.json'
        self.run_cli('design_pca.py','--target',self.demo/'pca_target.fasta','--output',out)
        data=json.loads(out.read_text())
        target=str(SeqIO.read(self.demo/'pca_target.fasta','fasta').seq)
        covered=set()
        for i,row in enumerate(data['oligos']):
            expected=target[row['start0']:row['end0']]
            if i%2:expected=str(Seq(expected).reverse_complement())
            self.assertEqual(row['sequence_5to3'],expected)
            self.assertLessEqual(row['length'],59)
            covered.update(range(row['start0'],row['end0']))
        self.assertEqual(covered,set(range(len(target))))
        self.assertEqual(len(data['oligos'])%2,0)
        for j in data['junctions']:
            self.assertTrue(20<=j['end0']-j['start0']<=40)
            self.assertTrue(60<=j['tm_c']<=64)
        with self.assertRaises(subprocess.CalledProcessError):
            self.run_cli('design_pca.py','--target',self.demo/'pca_target.fasta','--output',out)
        with self.assertRaises(subprocess.CalledProcessError):
            self.run_cli('design_pca.py','--target',self.demo/'pca_target.fasta','--output',self.demo/'invalid.json','--max-oligo',89)


if __name__=='__main__':unittest.main()
