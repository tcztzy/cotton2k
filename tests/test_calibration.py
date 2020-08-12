import unittest

from cotton2k.io import parse_list_dat, parse_parameter


class CalibrationTestCase(unittest.TestCase):
    def test_list_dat(self):
        result = parse_list_dat(
            """   2 GC510                              VPARGC510.DAT
   3 MAXXA                              VPARMAXXA.DAT
   4 ACALA SJ2                          VPARSJ2.DAT
   5 SIVON                              VPARSIVON.DAT
  11 DELTAPINE 61                       VPARDP61.DAT
  12 DELTAPINE 77                       VPARDP77.DAT
  17 新陆早8号                              VPARxlz8.DAT
  18 China Test1                        VPARTest1_xlz8.DAT
  19 China Test2                        VPARTest2_xlz8.DAT
  20 China Test3                        VPARTest3_xlz8.DAT"""
        )
        self.assertEqual(result[2], ("GC510", "VPARGC510.DAT"))
        self.assertEqual(result[3], ("MAXXA", "VPARMAXXA.DAT"))
        self.assertEqual(result[4], ("ACALA SJ2", "VPARSJ2.DAT"))
        self.assertEqual(result[5], ("SIVON", "VPARSIVON.DAT"))
        self.assertEqual(result[11], ("DELTAPINE 61", "VPARDP61.DAT"))
        self.assertEqual(result[12], ("DELTAPINE 77", "VPARDP77.DAT"))
        self.assertEqual(result[17], ("新陆早8号", "VPARxlz8.DAT"))
        self.assertEqual(result[18], ("China Test1", "VPARTest1_xlz8.DAT"))
        self.assertEqual(result[19], ("China Test2", "VPARTest2_xlz8.DAT"))
        self.assertEqual(result[20], ("China Test3", "VPARTest3_xlz8.DAT"))
        result = parse_list_dat(
            """   0 CALIFORNIA, West SJ                USACALSJWEST.DAT   
   1 ISRAEL, Coastal                    ILCOASTAL.DAT
   2 ARIZONA, Phoenix                   USAAZCENTRAL.DAT
   3 ISRAEL, Galil                      ILGALIL.DAT 
   4 ISRAEL, Avdat                      ILAVDAT.DAT 
"""
        )
        self.assertEqual(result[0], ("CALIFORNIA, West SJ", "USACALSJWEST.DAT"))
        self.assertEqual(result[1], ("ISRAEL, Coastal", "ILCOASTAL.DAT"))
        self.assertEqual(result[2], ("ARIZONA, Phoenix", "USAAZCENTRAL.DAT"))
        self.assertEqual(result[3], ("ISRAEL, Galil", "ILGALIL.DAT"))
        self.assertEqual(result[4], ("ISRAEL, Avdat", "ILAVDAT.DAT"))

    def test_parameter(self):
        result = parse_parameter(
            """Test
   0.050       Whatever
   0.049       Yet another""",
            2,
        )
        self.assertEqual(result[0], 0.050)
        self.assertEqual(result[1], 0.049)
