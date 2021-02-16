#ifndef SOIL_TYPE
#define SOIL_TYPE
typedef struct RootStruct
{
    double potential_growth; // potential root growth in a soil cell (g per day).
    double growth_factor;    // root growth correction factor in a soil cell (0 to 1).
    double actual_growth;
    double age;
    double weight_capable_uptake; // root weight capable of uptake, in g per soil cell.
    double weight[3];
} Root;

typedef struct SoilCellStruct
{
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
