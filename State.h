#ifndef STATE_TYPE
#define STATE_TYPE
#include "Root.h"
typedef struct State
{
    char date[12];
    double day_inc; // physiological days increment for this day. computes physiological age
    Root root[40][20];
} State;
#endif
