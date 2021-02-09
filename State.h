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
    double wind_speed;  // Hourly wind velocity, m per second.
} Hour;
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
    Hour hours[24];
    Root root[40][20];
    FruitingSite site[3][30][5];
} State;
#endif
