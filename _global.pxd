cdef extern:
    double dl(unsigned int)
    double wk(unsigned int, double)


cdef extern from "global.h":
    void InitializeGlobal()
    int isw
    int Kday
    int SoilMapFreq
    int OutIndex[24]
    const int maxl
    const int maxk
    int nl
    int nk
    double VarPar[61]
    double SitePar[21]
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
    double PlantRowLocation
    double PerPlantArea
    double DensityFactor
