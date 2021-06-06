from libc.stdint cimport uint32_t, int32_t

cdef extern:
    double dl(unsigned int)
    double wk(unsigned int, double)
    double tdewest(double, double, double)
    int SlabLoc(int, double)
    double dayrad(double, double, double, double)
    double dayrh(double, double)
    double daywnd(double, double, double, double, double, double)