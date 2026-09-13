import tempfile
import unittest
from pathlib import Path

from cef_util import combine_cef_results


class MergeCompoundsTests(unittest.TestCase):
    def combine(self, records, tolerance=0.1):
        # Exercise parsing and merging together using one CEF per observation.
        with tempfile.TemporaryDirectory(dir=Path(__file__).resolve().parent) as directory:
            for i, (name, rt, area) in enumerate(records):
                Path(directory, f'{i}.cef').write_text(
                    f'<CEF><Compound><Location rt="{rt}" ri="{1000 + i}" a="{area}"/>'
                    f'<Results><Molecule name="{name}" formula="C2H6">'
                    '<Database><Accession id="74-84-0"/></Database></Molecule></Results>'
                    f'<Spectrum><MSPeaks><p x="{50 + i}" y="999"/></MSPeaks></Spectrum>'
                    '</Compound></CEF>', encoding='utf-8')
            return combine_cef_results(directory, rt_tolerance=tolerance, exclude_elements=[])

    def test_retains_strongest_observation_anywhere_in_group(self):
        for areas, winner in [([10, 20, 30], 2), ([30, 10, 20], 0),
                              ([10, 30, 20], 1), ([30, 30, 30], 0)]:
            with self.subTest(areas=areas):
                result = self.combine([('A', 1 + i * 0.01, area)
                                       for i, area in enumerate(areas)])
                self.assertEqual(len(result), 1)
                retained = result.iloc[0]
                self.assertEqual(retained['MaxArea'], 30)
                self.assertEqual(retained['MS_Peaks'], [(50 + winner, 999)])
                self.assertEqual(retained['File'], f'{winner}.cef')
                self.assertAlmostEqual(retained['RT'], 1 + winner * 0.01)
                self.assertEqual(retained['RI'], 1000 + winner)

    def test_adjacent_rt_chain_remains_one_group(self):
        result = self.combine([('A', 1, 30), ('A', 1.075, 10), ('A', 1.15, 20)])
        self.assertEqual(len(result), 1)
        self.assertEqual(result.iloc[0]['MaxArea'], 30)

    def test_distinct_names_and_rt_groups_remain_separate(self):
        result = self.combine([('A', 1, 10), ('A', 1.01, 20), ('A', 1.02, 30),
                               ('A', 2, 40), ('B', 2.01, 50)])
        self.assertEqual(result['Chemical_Name'].tolist(), ['A', 'A peak 2', 'B'])
        self.assertEqual(result['MaxArea'].tolist(), [30, 40, 50])

    def test_rt_difference_equal_to_tolerance_starts_new_group(self):
        result = self.combine([('A', 1, 10), ('A', 1.125, 20)], tolerance=0.125)
        self.assertEqual(result['Chemical_Name'].tolist(), ['A', 'A peak 2'])


if __name__ == '__main__':
    unittest.main()
