import unittest
import tempfile
import datetime
from unittest.case import skip
from cotton2k.io import parse_profile, parse_profile_carbon_dioxide, parse_profile_description, parse_profile_output_flags, parse_profile_simulation_dates


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

    def test_description(self):
        line = 'HAKB1.PRO           HAZOR 1984 experiment, treatment KB1               '
        self.assertEqual(parse_profile_description(line)['description'], 'HAZOR 1984 experiment, treatment KB1')

    def test_simulation_dates(self):
        line = '08-APR-1984    01-APR-1984    28-SEP-1984                                       '
        result = parse_profile_simulation_dates(line)
        self.assertEqual(result.get('dateEmerge'), datetime.date(1984, 4, 8))
        self.assertEqual(result.get('dateSimStart'), datetime.date(1984, 4, 1))
        self.assertEqual(result.get('dateSimEnd'), datetime.date(1984, 9, 28))
        self.assertIsNone(result.get('datePlant'))
        line = '               01-APR-1984    28-SEP-1984                                       '
        with self.assertRaises(TypeError):
            parse_profile_simulation_dates(line)

    def test_carbon_dioxide(self):
        result = parse_profile_carbon_dioxide('')
        self.assertEqual(result.get('CO2EnrichmentFactor'), 0)
        result = parse_profile_carbon_dioxide('     1.000  100  101')
        self.assertEqual(result.get('CO2EnrichmentFactor'), 1.000)
        self.assertEqual(result.get('DayStartCO2'), 100)
        self.assertEqual(result.get('DayEndCO2'), 101)
        with self.assertRaises(ValueError):
            parse_profile_carbon_dioxide('     1.000  101  100')

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

    def test_output_flags(self):
        flags = '  0  0  1  1  0  1  1  1  1  1  1  1  0  0  0  0  1  0  0  0  0  0  0'
        result = parse_profile_output_flags(flags)
        self.assertFalse(result['UnitedStatesCustomarySystemOfUnitsOrInternationalSystemOfUnits'])
        self.assertFalse(result['perSquareMeterOrPerPlant'])
        self.assertTrue(result['outputDryWeight'])
