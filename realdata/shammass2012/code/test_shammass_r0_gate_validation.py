import importlib.util
import unittest
import copy
from pathlib import Path
from unittest.mock import patch

SCRIPT = Path(__file__).with_name('shammass_r0_v1_0_3_reviewed.py')
spec = importlib.util.spec_from_file_location('r0_gate_validation', SCRIPT)
r0 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(r0)

class GateValidationTests(unittest.TestCase):
    def setUp(self):
        self.root = Path('.')
        self.hashes = {'exp_large':'a'*64,'exp_small':'b'*64,'gpe_large':'c'*64,'gpe_small':'d'*64,'wpd.json':'e'*64}
        cases = {f'q{q}_m{m}_s{s}' for q in r0.SHAPES for m in r0.MS for s in r0.SIGNS}
        results = {}
        for series in ('large','small'):
            deltas = {}
            for d in ('1','2','4','8'):
                deltas[d] = {
                    'null': {'n':r0.N_NULL},
                    'signals': {case:{'n':r0.N_SIGNAL} for case in cases},
                    'calibration_jitter_summary': {'n':r0.N_NULL + len(cases)*r0.N_SIGNAL},
                }
            results[series] = {'deltas':deltas}
        self.data = {
            'spec':r0.SPEC, 'stage':'synthetic-gate', 'seed':r0.SEED,
            'replicates_per_series_delta':{'null':r0.N_NULL,'signal':r0.N_SIGNAL},
            'total_replicates':176000, 'input_sha256':copy.deepcopy(self.hashes),
            'spec_sha256':'f'*64, 'code_sha256':'1'*64,
            'results':results,
            'verdicts':{'large':{'verdict':'TESTABLE'},'small':{'verdict':'SIGN_ONLY'}},
        }
        self.patches = [
            patch.object(r0,'input_hash_record',return_value=self.hashes),
            patch.object(r0,'spec_file_hash',return_value='f'*64),
            patch.object(r0,'running_code_hash',return_value='1'*64),
        ]
        for p in self.patches: p.start()

    def tearDown(self):
        for p in reversed(self.patches): p.stop()

    def test_accepts_complete_structural_fixture(self):
        self.assertEqual(r0.validate_gate_file(self.data,self.root),self.hashes)

    def test_rejects_wrong_seed(self):
        self.data['seed'] += 1
        with self.assertRaises(r0.R0Error): r0.validate_gate_file(self.data,self.root)

    def test_rejects_wrong_code_hash(self):
        self.data['code_sha256'] = '2'*64
        with self.assertRaises(r0.R0Error): r0.validate_gate_file(self.data,self.root)

    def test_rejects_missing_delta(self):
        del self.data['results']['large']['deltas']['8']
        with self.assertRaises(r0.R0Error): r0.validate_gate_file(self.data,self.root)

    def test_rejects_wrong_jitter_count(self):
        self.data['results']['small']['deltas']['2']['calibration_jitter_summary']['n'] -= 1
        with self.assertRaises(r0.R0Error): r0.validate_gate_file(self.data,self.root)

    def test_rejects_missing_signal_case(self):
        self.data['results']['large']['deltas']['1']['signals'].pop('q2_m5_s1')
        with self.assertRaises(r0.R0Error): r0.validate_gate_file(self.data,self.root)

    def test_rejects_wrong_input_hashes(self):
        self.data['input_sha256']['wpd.json'] = '0'*64
        with self.assertRaises(r0.R0Error): r0.validate_gate_file(self.data,self.root)

if __name__ == '__main__':
    unittest.main(verbosity=2)
