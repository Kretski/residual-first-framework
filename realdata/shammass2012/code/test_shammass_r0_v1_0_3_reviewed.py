import importlib.util
import tempfile
import unittest
from pathlib import Path
import numpy as np

SCRIPT = Path(__file__).with_name("shammass_r0_v1_0_3_reviewed.py")
spec = importlib.util.spec_from_file_location("r0_reviewed", SCRIPT)
r0 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(r0)

class ReviewedR0Tests(unittest.TestCase):
    def test_curve_only_csv_does_not_parse_experimental_frequency(self):
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)/"exp.csv"
            p.write_text("k_um_inv,omega_over_2pi_kHz,px,py\n0.1,NOT_A_NUMBER,also_bad,also_bad\n",encoding="utf-8")
            d=r0.load_csv(p,{"k_um_inv"})
            self.assertEqual(d["k_um_inv"].tolist(),[0.1])
            self.assertNotIn("omega_over_2pi_kHz",d)
            self.assertNotIn("px",d)

    def test_holm_adjusted_values_in_original_q_order(self):
        got=r0.holm_adjusted([0.001,0.03,0.02])
        self.assertAlmostEqual(got[1],0.003)
        self.assertAlmostEqual(got[2],0.04)
        self.assertAlmostEqual(got[3],0.04)

    def test_holm_decision_uses_adjusted_p_values(self):
        fits={1:{"p":0.001,"gamma":2.0},2:{"p":0.03,"gamma":1.0},3:{"p":0.02,"gamma":-1.0}}
        d=r0.holm_decision(fits)
        self.assertEqual(d["outcome"],"DETECTED")
        self.assertEqual(d["q_star"],1)
        self.assertAlmostEqual(fits[1]["p_holm"],0.003)
        self.assertAlmostEqual(fits[2]["p_holm"],0.04)
        self.assertAlmostEqual(fits[3]["p_holm"],0.04)

    def test_g1_rule_does_not_add_not_testable_count_condition(self):
        # §8.5 G1 is determined by the two CP bounds only.
        self.assertTrue(r0.g1_cell_pass([0.009,0.014]))
        self.assertFalse(r0.g1_cell_pass([0.011,0.014]))
        self.assertFalse(r0.g1_cell_pass([0.009,0.016]))
        self.assertFalse(r0.g1_cell_pass([None,None]))

    def test_clopper_pearson_edges(self):
        self.assertEqual(r0.clopper_pearson(0,100)[0],0.0)
        self.assertEqual(r0.clopper_pearson(100,100)[1],1.0)
        with self.assertRaises(r0.R0Error):
            r0.clopper_pearson(101,100)

    def test_seedsequence_case_count(self):
        children=np.random.SeedSequence(r0.SEED).spawn(2*4*(1+3*4*2))
        self.assertEqual(len(children),200)

    def test_observed_design_diagnostics_include_base_condition(self):
        curve={"k":np.linspace(0.1,4.0,30),"f":np.linspace(0.05,1.0,30)**1.2}
        corr,conds=r0.design_diagnostics(curve,np.linspace(0.2,3.8,12))
        self.assertEqual(set(conds),{"base","1","2","3"})
        self.assertTrue(all(np.isfinite(v) for v in conds.values()))

if __name__ == "__main__":
    unittest.main(verbosity=2)
