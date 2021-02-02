#ifndef STATE_TYPE
#define STATE_TYPE
#include "Root.h"
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
    Root root[40][20];
    double plant_height;
    double runoff;
} State;
#endif
