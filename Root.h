#ifndef SOIL_CELL_TYPE
#define SOIL_CELL_TYPE
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
#endif
