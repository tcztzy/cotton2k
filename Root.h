#ifndef ROOT_TYPE
#define ROOT_TYPE
typedef struct RootStruct
{
    double potential_growth;
    double growth_factor;
    double actual_growth;
    double age;
    double weight_capable_uptake; // root weight capable of uptake, in g per soil cell.
    double weight[3];
} Root;
#endif
