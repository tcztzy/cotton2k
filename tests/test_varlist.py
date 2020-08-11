import unittest

from cotton2k.io import parse_varlist


class VarListTestCase(unittest.TestCase):
    def test_varlist(self):
        result = parse_varlist(
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
