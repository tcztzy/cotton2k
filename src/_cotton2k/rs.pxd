from libc.stdint cimport uint32_t, int32_t
from .fruiting_site cimport Stage

cdef extern:
    double dl(unsigned int)
    double wk(unsigned int, double)
    double tdewest(double, double, double)
    int SlabLoc(int, double)
    double daywnd(double, double, double, double, double, double)
    double AddPlantHeight(double, double, uint32_t, Stage, double, double, double, double, double, double, double,
                          double, double, double, double, double, double, double)
    double TemperatureOnFruitGrowthRate(double)
    double VaporPressure(double)
    double clearskyemiss(double, double)
    double SoilNitrateOnRootGrowth(double)
    double SoilAirOnRootGrowth(double, double, double)
    double SoilMechanicResistance(double)
    double SoilTemOnRootGrowth(double)
    double wcond(double, double, double, double, double, double)