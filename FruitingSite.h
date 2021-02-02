#ifndef FRUITING_SITE_TYPE
#define FRUITING_SITE_TYPE
typedef struct LeafStruct
{
    double age; // leaf age at each fruiting site, physiological days.
} Leaf;
typedef struct FruitingSiteStruct
{
    double age; // age of each fruiting site, physiological days from square initiation.
    Leaf leaf;
} FruitingSite;
#endif