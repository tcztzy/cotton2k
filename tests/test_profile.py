import datetime
import unittest

from cotton2k.io import (
    parse_profile,
    parse_profile_carbon_dioxide,
    parse_profile_description,
    parse_profile_field,
    parse_profile_geometry,
    parse_profile_output_flags,
    parse_profile_output_options,
    parse_profile_parameter_files,
    parse_profile_simulation_dates,
    parse_profile_soil_mulch,
    parse_profile_weather,
)


class ProfileFileTestCase(unittest.TestCase):
    def test_parsing(self):
        result = parse_profile(
            """test.pro            Test profile
01-MAY-2020    20-APR-2020    15-OCT-2020                        1.000  100  101
test.act                                         1     0.000     0.000  122    0
test.hyd            test.int            test.agi
    40.548    81.296  1013.000         0
    75.000     0.000    40.000         0
        10    10-APR-2020    20-OCT-2020        10    01-JUN-2020    20-OCT-2020
  0  0  1  1  0  1  0  1  1  1  1  1  0  0  0  0  0  1  0  0  0  0  0"""
        )
        self.assertEqual(result['dayEndMulch'], 289)

    def test_description(self):
        line = "HAKB1.PRO           HAZOR 1984 experiment, treatment KB1               "
        self.assertEqual(
            parse_profile_description(line)["description"],
            "HAZOR 1984 experiment, treatment KB1",
        )

    def test_simulation_dates(self):
        line = "08-APR-1984    01-APR-1984    28-SEP-1984                                       "
        result = parse_profile_simulation_dates(line)
        self.assertEqual(result.get("dateEmerge"), datetime.date(1984, 4, 8))
        self.assertEqual(result.get("dateSimStart"), datetime.date(1984, 4, 1))
        self.assertEqual(result.get("dateSimEnd"), datetime.date(1984, 9, 28))
        self.assertIsNone(result.get("datePlant"))
        line = "               01-APR-1984    28-SEP-1984                                       "
        with self.assertRaises(TypeError):
            parse_profile_simulation_dates(line)

    def test_carbon_dioxide(self):
        result = parse_profile_carbon_dioxide("")
        self.assertEqual(result.get("CO2EnrichmentFactor"), 0)
        result = parse_profile_carbon_dioxide("     1.000  100  101")
        self.assertEqual(result.get("CO2EnrichmentFactor"), 1.000)
        self.assertEqual(result.get("DayStartCO2"), 100)
        self.assertEqual(result.get("DayEndCO2"), 101)
        with self.assertRaises(ValueError):
            parse_profile_carbon_dioxide("     1.000  101  100")

    def test_weather(self):
        result = parse_profile_weather(
            "REHA84.ACT                                                                      "
        )
        self.assertEqual(result.get("actualWeatherFileName"), "REHA84.ACT")

    def test_soil_mulch(self):
        result = parse_profile_soil_mulch("")
        self.assertEqual(result["mulchIndicator"], 0)

    def test_parameter_files(self):
        result = parse_profile_parameter_files(
            "HAZOR.HYD           HAZOR.INT           HAKB1.AGI                               "
        )
        self.assertEqual(result.get("soilHydraulicFileName"), "HAZOR.HYD")
        self.assertEqual(result.get("soilInitFileName"), "HAZOR.INT")
        self.assertEqual(result.get("agriculturalInputFileName"), "HAKB1.AGI")
        self.assertEqual(result.get("plantmapFileName"), "")

    def test_geometry(self):
        result = parse_profile_geometry("    32.000    35.000    50.000         1")
        self.assertEqual(result.get("latitude"), 32.000)
        self.assertEqual(result.get("longitude"), 35.000)
        self.assertEqual(result.get("elevation"), 50.000)
        self.assertEqual(result.get("siteNumber"), 1)

    def test_field(self):
        result = parse_profile_field("    96.520     0.000    10.000         4")
        self.assertEqual(result["rowSpace"], 96.520)
        self.assertEqual(result["skipRowWidth"], 0.000)
        self.assertEqual(result["plantsPerMeter"], 10.000)
        self.assertEqual(result["varNumber"], 4)

    def test_output_options(self):
        result = parse_profile_output_options(
            "        10    20-APR-1984    20-SEP-1984        10    01-JUN-1984    20-SEP-1984    "
        )
        self.assertEqual(result["soilMapFrequency"], 10)
        self.assertEqual(result["soilMapStartDate"], datetime.date(1984, 4, 20))
        self.assertEqual(result["soilMapEndDate"], datetime.date(1984, 9, 20))
        self.assertEqual(result["plantMapFrequency"], 10)
        self.assertEqual(result["plantMapStartDate"], datetime.date(1984, 6, 1))
        self.assertEqual(result["plantMapEndDate"], datetime.date(1984, 9, 20))

    def test_output_flags(self):
        flags = "  0  0  1  1  0  1  1  1  1  1  1  1  0  0  0  0  1  0  0  0  0  0  0"
        result = parse_profile_output_flags(flags)
        self.assertFalse(
            result["UnitedStatesCustomarySystemOfUnitsOrInternationalSystemOfUnits"]
        )
        self.assertFalse(result["perSquareMeterOrPerPlant"])
        self.assertTrue(result["outputDryWeight"])
