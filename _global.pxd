cdef extern from "global.h":
    void InitializeGlobal()
    int isw
    int Kday
    int SoilMapFreq
    int OutIndex[24]
    double PlantPopulation
    double TotalSoilNo3N
    double TotalSoilNh4N
    double TotalSoilUreaN
    double SoilNitrogenAtStart
    double PlantWeightAtStart
    double TotalRootWeight
    double TotalStemWeight
    double TotalLeafWeight
    double ReserveC
