#ifndef STATE_TYPE
#define STATE_TYPE
#include "Root.h"
#include "FruitingSite.h"

typedef struct HourStruct
{
    double temperature; // hourly air temperatures, C.
    double radiation;   // hourly global radiation, W / m2.
    double cloud_cov;   // cloud cover ratio (0 to 1).
    double cloud_cor;   // hourly cloud type correction.
    double et1;         // part of hourly Penman evapotranspiration affected by net radiation, in mm per hour.
    double et2;         // part of hourly Penman evapotranspiration affected by wind and vapor pressure deficit, in mm per hour.
    double ref_et;      // reference evapotranspiration, mm per hour.
    double wind_speed;  // Hourly wind velocity, m per second.
    double dew_point;   // hourly dew point temperatures, C.
    double humidity;    // hourly values of relative humidity (%).
    double albedo;      // hourly albedo of a reference crop.
} Hour;
typedef struct MainStemLeafStruct
{
    double leaf_area;
    double leaf_weight;                         // mainstem leaf weight at each node, g.
    double petiole_weight;                      // weight of mainstem leaf petiole at each node, g.
    double potential_growth_for_leaf_area;      // potential growth in area of an individual main stem node leaf, dm2 day-1.
    double potential_growth_for_leaf_weight;    // potential growth in weight of an individual main stem node leaf, g day-1.
    double potential_growth_for_petiole_weight; // potential growth in weight of an individual main stem node petiole, g day-1.
} MainStemLeaf;
typedef struct FruitingBranchStruct
{
    unsigned int number_of_fruiting_nodes; // number of nodes on each fruiting branch.
    double delay_for_new_node;             // cumulative effect of stresses on delaying the formation of a new node on a fruiting branch.
    MainStemLeaf main_stem_leaf;
    FruitingSite nodes[5];
} FruitingBranch;
typedef struct VegetativeBranchStruct
{
    unsigned int number_of_fruiting_branches; // number of fruiting branches at each vegetative branch.
    FruitingBranch fruiting_branches[30];
} VegetativeBranch;
typedef struct State
{
    char date[12];
    double day_inc;                          // physiological days increment for this day. computes physiological age
    unsigned int number_of_layers_with_root; //
    double abscised_fruit_sites;             // total number of abscised fruit sites, per plant.
    double abscised_leaf_weight;             // weight of abscissed leaves, g per plant.
    double water_stress;                     // general water stress index (0 to 1).
    double carbon_stress;                    // carbohydrate stress factor.
    double day_length;                       // day length, in hours
    double plant_height;
    double runoff;
    double solar_noon;
    double net_radiation;                       // daily total net radiation, W m-2.
    double evapotranspiration;                  // daily sum of hourly reference evapotranspiration, mm per day.
    unsigned int number_of_vegetative_branches; // number of vegetative branches (including the main branch), per plant.
    unsigned int number_of_fruiting_sites;      // total number of fruiting sites per plant.
    double number_of_squares;                   // number of squares per plant.
    double number_of_green_bolls;               // average number of retained green bolls, per plant.
    VegetativeBranch vegetative_branches[3];
    Hour hours[24];
    Root root[40][20];
    FruitingSite site[3][30][5];
} State;
#endif
