from libcpp cimport bool as bool_t

cdef extern from "Soil.h":

    ctypedef struct Root:
        double potential_growth
        double growth_factor
        double actual_growth
        double age
        double weight_capable_uptake
        double weight[3]

    ctypedef struct SoilCell:
        double nitrate_nitrogen_content
        double fresh_organic_matter
        Root root

    ctypedef struct SoilLayer:
        unsigned int number_of_left_columns_with_root
        unsigned int number_of_right_columns_with_root

    ctypedef struct cSoil "Soil":
        unsigned int number_of_layers_with_root
        SoilLayer layers[40]
        SoilCell cells[40][20]

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
    ctypedef struct Boll:
        double age
        double potential_growth
        double weight
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
        Stage stage
        Leaf leaf
        SquareStruct square
        Boll boll
        Burr burr
        Petiole petiole

cdef extern from "State.hpp":
    ctypedef struct Hour:
        double temperature
        double radiation
        double cloud_cov
        double cloud_cor
        double et1
        double et2
        double ref_et
        double wind_speed
        double dew_point
        double humidity
        double albedo
    ctypedef struct MainStemLeaf:
        double leaf_area
        double leaf_weight
        double petiole_weight
        double potential_growth_for_leaf_area
        double potential_growth_for_leaf_weight
        double potential_growth_for_petiole_weight
    ctypedef struct cFruitingBranch "FruitingBranch":
        unsigned int number_of_fruiting_nodes
        double delay_for_new_node
        MainStemLeaf main_stem_leaf
        FruitingSite nodes[5]
    ctypedef struct cVegetativeBranch "VegetativeBranch":
        unsigned int number_of_fruiting_branches
        cFruitingBranch fruiting_branches[30]
    ctypedef struct cState "State":
        unsigned int daynum
        double day_inc
        double lint_yield
        double bloom_weight_loss
        double abscised_fruit_sites
        double abscised_leaf_weight
        double cumulative_nitrogen_loss
        double applied_water
        double water_stress
        double water_stress_stem
        double carbon_stress
        double extra_carbon
        double day_length
        double plant_height
        double plant_weight
        double runoff
        double solar_noon
        double net_radiation
        double evapotranspiration
        double actual_transpiration
        double cumulative_transpiration
        double actual_soil_evaporation
        double cumulative_evaporation
        unsigned int number_of_vegetative_branches
        unsigned int number_of_fruiting_sites
        double number_of_squares
        double number_of_green_bolls
        double number_of_open_bolls
        double nitrogen_stress
        double total_required_nitrogen
        double leaf_area_index
        double leaf_nitrogen_concentration
        double petiole_nitrogen_concentration
        double seed_nitrogen_concentration
        double root_nitrogen_concentration
        double ginning_percent
        bool_t pollination_switch
        cVegetativeBranch vegetative_branches[3]
        Hour hours[24]
        cSoil soil

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
    ctypedef struct cSimulation "Simulation":
        int year
        unsigned int day_start
        unsigned int day_finish
        unsigned int day_emerge
        unsigned int day_plant
        unsigned int day_topping
        double latitude
        double longitude
        double elevation
        double row_space
        unsigned int plant_row_column
        cState *states
        ClimateStruct climate[400]
        Irrigation irrigation[150]
