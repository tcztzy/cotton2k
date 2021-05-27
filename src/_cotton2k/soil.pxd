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
