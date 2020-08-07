import unittest
import tempfile
import datetime
from cotton2k.io import parse_profile


class ProfileFileTestCase(unittest.TestCase):
    content = """\
HAKB1.PRO           HAZOR 1984 experiment, treatment KB1               
08-APR-1984    01-APR-1984    28-SEP-1984                                       
REHA84.ACT                                                                      
HAZOR.HYD           HAZOR.INT           HAKB1.AGI                               
    32.000    35.000    50.000         1
    96.520     0.000    10.000         4
        10    20-APR-1984    20-SEP-1984        10    01-JUN-1984    20-SEP-1984    
  0  0  1  1  0  1  1  1  1  1  1  1  0  0  0  0  1  0  0  0  0  0  0"""

    def test_line_1(self):
        self.assertEqual(parse_profile(self.content)['description'], 'HAZOR 1984 experiment, treatment KB1')
    
    def test_line_2(self):
        result = parse_profile(self.content)
        self.assertEqual(result.get('dateEmerge'), datetime.date(1984, 4, 8))
        self.assertEqual(result.get('dateSimStart'), datetime.date(1984, 4, 1))
        self.assertEqual(result.get('dateSimEnd'), datetime.date(1984, 9, 28))
        self.assertNotIn('datePlant', result)
        self.assertEqual(result.get('CO2EnrichmentFactor'), 0)

    def test_line_3(self):
        result = parse_profile(self.content)
        self.assertEqual(result.get('actualWeatherFileName'), 'REHA84.ACT')

    def test_line_4(self):
        result = parse_profile(self.content)
        self.assertEqual(result.get('soilHydraulicFileName'), 'HAZOR.HYD')
        self.assertEqual(result.get('soilInitFileName'), 'HAZOR.INT')
        self.assertEqual(result.get('agriculturalInputFileName'), 'HAKB1.AGI')
        self.assertEqual(result.get('plantmapFileName'), '')

    def test_line_5(self):
        result = parse_profile(self.content)
        self.assertEqual(result.get('latitude'), 32.000)
        self.assertEqual(result.get('longitude'), 35.000)
        self.assertEqual(result.get('elevation'), 50.000)
        self.assertEqual(result.get('siteNumber'), 1)

    def test_line_6(self):
        result = parse_profile(self.content)
        self.assertEqual(result['rowSpace'], 96.520)
        self.assertEqual(result['skipRowWidth'], 0.000)
        self.assertEqual(result['plantsPerMeter'], 10.000)
        self.assertEqual(result['varNumber'], 4)

    def test_line_7(self):
        result = parse_profile(self.content)
        self.assertEqual(result['soilMapFrequency'], 10)
        self.assertEqual(result['soilMapStartDate'], datetime.date(1984, 4, 20))
        self.assertEqual(result['soilMapEndDate'], datetime.date(1984, 9, 20))
        self.assertEqual(result['plantMapFrequency'], 10)
        self.assertEqual(result['plantMapStartDate'], datetime.date(1984, 6, 1))
        self.assertEqual(result['plantMapEndDate'], datetime.date(1984, 9, 20))
