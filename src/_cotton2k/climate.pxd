cdef extern from "Climate.h":
    ctypedef struct ClimateStruct:
        double Rad
        double Tmax
        double Tmin
        double Rain
        double Wind
        double Tdew
