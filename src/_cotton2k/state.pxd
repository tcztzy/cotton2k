from libcpp cimport bool as bool_t
from .fruiting_site cimport FruitingSite
from .soil cimport cSoil

cdef extern from "State.hpp":
    ctypedef struct cHour "Hour":
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

    ctypedef struct cMainStemLeaf "MainStemLeaf":
        double leaf_area
        double leaf_weight
        double petiole_weight
        double potential_growth_for_leaf_area
        double potential_growth_for_leaf_weight
        double potential_growth_for_petiole_weight

    ctypedef struct cFruitingBranch "FruitingBranch":
        unsigned int number_of_fruiting_nodes
        double delay_for_new_node
        cMainStemLeaf main_stem_leaf
        FruitingSite nodes[5]

    ctypedef struct cVegetativeBranch "VegetativeBranch":
        unsigned int number_of_fruiting_branches
        double delay_for_new_fruiting_branch
        cFruitingBranch fruiting_branches[30]

    ctypedef struct cState "State":
        unsigned int daynum
        unsigned int kday
        double day_inc
        double bloom_weight_loss
        double abscised_fruit_sites
        double abscised_leaf_weight
        double cumulative_nitrogen_loss
        double water_stress
        double water_stress_stem
        double carbon_stress
        double day_length
        double plant_height
        double stem_weight
        double root_weight
        double petiole_weight
        double square_weight
        double green_bolls_weight
        double green_bolls_burr_weight
        double open_bolls_weight
        double open_bolls_burr_weight
        double reserve_carbohydrate
        double runoff
        double solar_noon
        double net_radiation
        double evapotranspiration
        double actual_transpiration
        double potential_evaporation
        double cumulative_transpiration
        double actual_soil_evaporation
        double cumulative_evaporation
        double light_interception
        unsigned int number_of_vegetative_branches
        unsigned int number_of_fruiting_sites
        double number_of_squares
        double number_of_green_bolls
        double number_of_open_bolls
        double nitrogen_stress
        double nitrogen_stress_vegetative
        double nitrogen_stress_fruiting
        double nitrogen_stress_root
        double total_required_nitrogen
        double leaf_area_index
        double leaf_area
        double leaf_weight
        double leaf_weight_pre_fruiting[9]
        double leaf_weight_area_ratio
        double leaf_nitrogen_concentration
        double leaf_nitrogen
        double petiole_nitrogen_concentration
        double seed_nitrogen_concentration
        double seed_nitrogen
        double root_nitrogen_concentration
        double root_nitrogen
        double square_nitrogen_concentration
        double burr_nitrogen_concentration
        double burr_nitrogen
        double square_nitrogen
        double stem_nitrogen_concentration
        double stem_nitrogen
        double fruit_growth_ratio
        double ginning_percent
        double average_temperature
        double daytime_temperature
        double deep_soil_temperature
        double total_actual_leaf_growth
        double total_actual_petiole_growth
        double actual_square_growth
        double actual_stem_growth
        double actual_boll_growth
        double actual_burr_growth
        double supplied_nitrate_nitrogen
        double supplied_ammonium_nitrogen
        double petiole_nitrogen
        double petiole_nitrate_nitrogen_concentration
        bool_t pollination_switch
        double age_of_pre_fruiting_nodes[9]
        int number_of_pre_fruiting_nodes
        double leaf_area_pre_fruiting[9]
        double delay_for_new_branch[3]
        cVegetativeBranch vegetative_branches[3]
        cHour hours[24]
        cSoil soil


cdef class StateBase:
    cdef cState *_
    cdef public unsigned int seed_layer_number  # layer number where the seeds are located.
    cdef public unsigned int year
    cdef public unsigned int version
    cdef public double carbon_allocated_for_root_growth # available carbon allocated for root growth, g per plant.
    cdef public double delay_of_emergence  # effect of negative values of xt on germination rate.
    cdef public double extra_carbon  # Extra carbon, not used for plant potential growth requirements, assumed to accumulate in taproot.
    cdef public double fiber_length
    cdef public double fiber_strength
    cdef public double hypocotyl_length  # length of hypocotyl, cm.
    cdef public double lint_yield  # yield of lint, kgs per hectare.
    cdef public double net_photosynthesis  # net photosynthetic rate, g per plant per day.
    cdef public double pavail  # residual available carbon for root growth from previous day.
    cdef public double seed_moisture  # moisture content of germinating seeds, percent.
    cdef public double stem_potential_growth  # potential growth rate of stems, g plant-1 day-1.
