# distutils: language=c++
# cython: language_level=3
from _cotton2k._global cimport (
    gh2oc,
    tstbd,
    impede,
    inrim,
    ncurve,
)


cdef class SoilImpedance:

    @property
    def curves(self):
        return {gh2oc[i]: {tstbd[j][i]: impede[j][i] for j in range(inrim)} for i in range(ncurve)}

    @curves.setter
    def curves(self, impedance_table):
        global ncurve, inrim
        ncurve = len(impedance_table)
        inrim = len(impedance_table[0])
        for i, row in enumerate(impedance_table):
            gh2oc[i] = row.pop("water")
            for j, pair in enumerate(sorted(row.items())):
                tstbd[j][i], impede[j][i] = pair
