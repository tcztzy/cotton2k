from libc.stdint cimport uint32_t, int32_t
from .fruiting_site cimport Stage

cdef extern:
    double dl(unsigned int)
    double wk(unsigned int, double)
    double tdewest(double, double, double)
    int SlabLoc(int, double)
    double dayrad(double, double, double, double)
    double dayrh(double, double)
    double daywnd(double, double, double, double, double, double)
    double PotentialStemGrowth(double, int, Stage, double, double, double, double, double, double, double, double)
    double AddPlantHeight(double, double, uint32_t, Stage, double, double, double, double, double, double, double,
                          double, double, double, double, double, double, double)
    double TemperatureOnFruitGrowthRate(double)