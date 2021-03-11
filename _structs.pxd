cdef extern from "State.h":
    ctypedef struct State:
        double plant_height
        double plant_weight
        double lint_yield

cdef extern from "Climate.h":
    ctypedef struct ClimateStruct:
        pass

cdef extern from "Simulation.hpp":
    ctypedef struct Simulation:
        const char *profile_name
        int year
        unsigned int day_start
        unsigned int day_finish
        unsigned int day_emerge
        unsigned int day_plant
        unsigned int day_start_soil_maps
        unsigned int day_stop_soil_maps
        double latitude
        double longitude
        double elevation
        double row_space
        State *states
        ClimateStruct climate[400]
