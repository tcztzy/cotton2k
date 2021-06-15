from .climate cimport ClimateStruct
from .irrigation cimport Irrigation
from .state cimport cState

cdef extern from "PlantGrowth.h":
    void LeafWaterPotential(cState &, double)

cdef extern from "Simulation.hpp":
    ctypedef struct cSimulation "Simulation":
        int year
        unsigned int day_emerge
        unsigned int day_start
        unsigned int day_finish
        unsigned int day_plant
        unsigned int day_topping
        unsigned int day_defoliate
        double latitude
        double longitude
        double elevation
        double row_space
        double plant_population
        double per_plant_area
        double density_factor
        unsigned int first_bloom
        unsigned int first_square
        unsigned int plant_row_column
        double cultivar_parameters[61]
        ClimateStruct climate[400]
        Irrigation irrigation[150]
        cState states[200]

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
    double SandVolumeFraction[40]
    double ClayVolumeFraction[40]
    double ActualStemGrowth
    double PotGroStem
    double PotGroAllRoots
    double TotalPetioleWeight
    double LwpMin
    double LwpMinX[3]
    double LwpMax
    double LwpX[3]
    double AverageLwp
    double AverageLwpMin
    double StemWeight[365]
    double NetPhotosynthesis
    double CumNetPhotosynth
    double DayTimeTemp
    double NightTimeTemp
    double PotGroAllSquares
    double PotGroAllBolls
    double PotGroAllBurrs
    double PotGroAllLeaves
    double PotGroAllPetioles
    double PotGroLeafAreaPreFru[9]
    double PotGroLeafWeightPreFru[9]
    double PotGroPetioleWeightPreFru[9]
    int DefoliationDate[5]
    int DefoliationMethod[5]
    double PercentDefoliation
    double VolWaterContent[40][20]
    NitrogenFertilizer NFertilizer[150]
    int NumNitApps
    int NumIrrigations
