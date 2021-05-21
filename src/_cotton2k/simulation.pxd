from libc.stdlib cimport malloc
from libc.math cimport exp
from libc.stdint cimport uint32_t

from _cotton2k._global cimport *
from _cotton2k._structs cimport *
from _cotton2k._io cimport *

cdef extern from "CottonPhenology.h":
    void CottonPhenology(cSimulation &, uint32_t)

cdef extern from "DailyClimate.h":
    void DayClim(cSimulation &, uint32_t u)

cdef extern from "PlantGrowth.h":
    void CheckDryMatterBal(cState &)
    void Defoliate(cSimulation &, uint32_t)
    void GetNetPhotosynthesis(cSimulation &, uint32_t, const double &)
    double PhysiologicalAge(Hour[24])
    void PlantGrowth(cSimulation &, const uint32_t &, const int &, const double &)
    void Stress(cSimulation &, unsigned int)

cdef extern from "PlantNitrogen.h":
    void PlantNitrogen(cSimulation &, uint32_t)

cdef extern from "SoilNitrogen.h":
    void SoilNitrogen(cSimulation &, unsigned int)

cdef extern from "SoilProcedures.h":
    void SoilProcedures(cSimulation &, uint32_t)
    void SoilSum(cState &, double)

cdef extern from "SoilTemperature.h":
    void ColumnShading(cState &, double[20], double, double, unsigned int)
    void SoilTemperature(cSimulation &, uint32_t, double[20])
