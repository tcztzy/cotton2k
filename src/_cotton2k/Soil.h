#ifndef SOIL_TYPE
#define SOIL_TYPE
typedef struct RootStruct
{
    double growth_factor;    // root growth correction factor in a soil cell (0 to 1).
    double weight_capable_uptake; // root weight capable of uptake, in g per soil cell.
    double weight[3];
} Root;

typedef struct SoilCellStruct
{
    double nitrate_nitrogen_content; // volumetric nitrate nitrogen content of a soil cell, mg N cm-3.
    double fresh_organic_matter;     // fresh organic matter in the soil, mg / cm3.
    double water_content;            // volumetric water content of a soil cell, cm3 cm-3.
    Root root;
} SoilCell;

typedef struct SoilLayerStruct
{
    unsigned int number_of_left_columns_with_root;  // first column with roots in a soil layer.
    unsigned int number_of_right_columns_with_root; // last column with roots in a soil layer.
} SoilLayer;

typedef struct SoilStruct
{
    unsigned int number_of_layers_with_root;
    SoilLayer layers[40];
    SoilCell cells[40][20];
} Soil;
#endif
