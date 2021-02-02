#ifndef FRUITING_SITE_TYPE
#define FRUITING_SITE_TYPE
typedef struct LeafStruct
{
    double age; // leaf age at each fruiting site, physiological days.
} Leaf;
typedef struct FruitingSiteStruct
{
    Leaf leaf;
} FruitingSite;
#endif