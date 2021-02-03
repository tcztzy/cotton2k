#ifndef FRUITING_SITE_TYPE
#define FRUITING_SITE_TYPE
typedef struct LeafStruct
{
    double age; // leaf age at each fruiting site, physiological days.
} Leaf;
typedef struct BollStruct
{
    double age; // age of each boll, physiological days from flowering.
} Boll;
typedef struct FruitingSiteStruct
{
    double age; // age of each fruiting site, physiological days from square initiation.
    Leaf leaf;
    Boll boll;
} FruitingSite;
#endif