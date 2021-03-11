cdef extern:
    double dl(unsigned int)
    double wk(unsigned int, double)


cdef extern from "global.h":
    void InitializeGlobal()
    unsigned int ncurve
    int inrim
    int isw
    int Kday
    const int maxl
    const int maxk
    int nl
    int nk
    double gh2oc[10]
    double tstbd[10][10]
    double impede[10][10]
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
