cdef extern:
    double dl(unsigned int)
    double wk(unsigned int, double)
    double tdewest(double, double, double)
    int SlabLoc(int, double)


cdef extern from "global.h":
    ctypedef struct NitrogenFertilizer:
        int day
        int mthfrt
        int ksdr
        int lsdr
        double amtamm
        double amtnit
        double amtura
    void InitializeGlobal()
    unsigned int ncurve
    int inrim
    int isw
    int LastDayWeatherData
    const int maxl
    const int maxk
    int nl
    int nk
    double gh2oc[10]
    double tstbd[10][10]
    double impede[10][10]
    double SitePar[21]
    double TotalSoilNo3N
    double TotalSoilNh4N
    double TotalSoilUreaN
    double TotalRootWeight
    double ReserveC
    double PlantRowLocation
    double RatioImplicit
    double conmax
    double airdr[9]
    double thetas[9]
    double alpha[9]
    double vanGenuchtenBeta[9]
    double SaturatedHydCond[9]
    double BulkDensity[9]
    double DefoliantAppRate[5]
    int DefoliationDate[5]
    int DefoliationMethod[5]
    NitrogenFertilizer NFertilizer[150]
    int NumNitApps
    int NumIrrigations
