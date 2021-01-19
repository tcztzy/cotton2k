#ifndef STATE_TYPE
#define STATE_TYPE
#include "Root.h"
typedef struct State
{
    uint32_t day_number = 0;   // Days from the start of the first year of simulation (day of year = DOY)
    char date[12];             // date string formatted as "dd-MMM-yyyy", for example 25-JUN-2003
    double plant_height = 4.0; // plant height, cm.
    Root root;
} State;
#endif
