#ifndef FRUITING_SITE_TYPE
#define FRUITING_SITE_TYPE
// code indicating the developmental state of each fruiting site
enum Stage
{
    NotYetFormed,
    Square,
    GreenBoll, // not susceptible to shedding
    MatureBoll,
    AbscisedAsBoll,
    AbscisedAsSquare,
    AbscisedAsFlower,
    YoungGreenBoll, // susceptible to shedding
};
typedef struct LeafStruct
{
    double age; // leaf age at each fruiting site, physiological days.
} Leaf;
typedef struct BollStruct
{
    double age;              // age of each boll, physiological days from flowering.
    double potential_growth; // age of each boll, physiological days from flowering.
    double weight;           // weight of seedcotton for each site, g per plant.
} Boll;
typedef struct BurrStruct
{
    double weight; // weight of burrs for each site, g per plant.
} Burr;
typedef struct PetioleStruct
{
    double weight; // petiole weight at each fruiting site, g.
} Petiole;
typedef struct FruitingSiteStruct
{
    double age; // age of each fruiting site, physiological days from square initiation.
    enum Stage stage;
    Leaf leaf;
    Boll boll;
    Burr burr;
    Petiole petiole;
} FruitingSite;
#endif