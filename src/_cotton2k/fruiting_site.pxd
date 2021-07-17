cdef extern from "FruitingSite.h":
    cdef enum Stage:
        NotYetFormed
        Square
        GreenBoll
        MatureBoll
        AbscisedAsBoll
        AbscisedAsSquare
        AbscisedAsFlower
        YoungGreenBoll
    ctypedef struct Leaf:
        double age
        double potential_growth
        double area
        double weight
    struct SquareStruct:
        double potential_growth
        double weight
    ctypedef struct cBoll "Boll":
        double age
        double potential_growth
        double weight
        double cumulative_temperature
    ctypedef struct Burr:
        double potential_growth
        double weight
    ctypedef struct Petiole:
        double potential_growth
        double weight
    ctypedef struct FruitingSite:
        double age
        double fraction
        double average_temperature
        double ginning_percent
        Stage stage
        Leaf leaf
        SquareStruct square
        cBoll boll
        Burr burr
        Petiole petiole
