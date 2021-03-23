cdef extern from "State.h":
    ctypedef struct State:
        double plant_height
        double plant_weight
        double lint_yield

cdef extern from "Climate.h":
    ctypedef struct ClimateStruct:
        double Rad
        double Tmax
        double Tmin
        double Rain
        double Wind
        double Tdew

cdef extern from "Irrigation.h":
    ctypedef struct Irrigation:
        int day
        int method
        int LocationColumnDrip
        int LocationLayerDrip
        double amount

cdef extern from "Simulation.hpp":
    ctypedef struct Simulation:
        const char *profile_name
        int year
        unsigned int day_start
        unsigned int day_finish
        unsigned int day_emerge
        unsigned int day_plant
        double latitude
        double longitude
        double elevation
        double row_space
        unsigned int plant_row_column
        State *states
        ClimateStruct climate[400]
        Irrigation irrigation[150]
