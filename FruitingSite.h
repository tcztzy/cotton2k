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
    double age;    // leaf age at each fruiting site, physiological days.
    double area;   // leaf area at each fruiting site, dm2.
    double weight; // leaf weight at each fruiting site, g.
} Leaf;
struct SquareStruct
{
    double potential_growth; // potential growth in weight of an individual fruiting node squares, g day-1.
    double weight;           // weight of each square, g per plant.
};
typedef struct BollStruct
{
    double age;              // age of each boll, physiological days from flowering.
    double potential_growth; // potential growth in weight of an individual fruiting node bolls, g day-1.
    double weight;           // weight of seedcotton for each site, g per plant.
} Boll;
typedef struct BurrStruct
{
    double weight; // weight of burrs for each site, g per plant.
} Burr;
typedef struct PetioleStruct
{
    double potential_growth; // potential growth in weight of an individual fruiting node petiole, g day-1.
    double weight;           // petiole weight at each fruiting site, g.
} Petiole;
typedef struct FruitingSiteStruct
{
    double age;      // age of each fruiting site, physiological days from square initiation.
    double fraction; // fraction of fruit remaining at each fruiting site (0 to 1).
    enum Stage stage;
    Leaf leaf;
    struct SquareStruct square;
    Boll boll;
    Burr burr;
    Petiole petiole;
} FruitingSite;
#endif