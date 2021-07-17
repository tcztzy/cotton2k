//  CottonPhenology.cpp
//
//   Functions in this file:
// CottonPhenology()
// SimulateFruitingSite{}
// NewBollFormation()
// BollOpening()
//
#include <algorithm>
#include <cmath>
#include "LeafAbscission.h"
#include "FruitAbscission.h"

using namespace std;

//   Declaration of file-scope variables:
double FibLength;          // fiber length
double FibStrength;        // fiber strength
//////////////////////////////////////////////////
//      The following documentation describes the implementation of
//  the simulation of the cotton plant phenology and the abscission of
//  leaves and of fruiting sites in Cotton2K.
//
//      The calling sequence is as follows:
//      LeafAbscission() calls PreFruitLeafAbscission(), MainStemLeafAbscission(),
//          FruitNodeLeafAbscission(), DefoliationLeafAbscission()
//          === see file LeafAbscission.cpp
//      FruitingSitesAbscission() calls SiteAbscissionRatio(), SquareAbscission(), BollAbscission() and ComputeSiteNumbers()
//          === see file FruitAbscission.cpp
//////////////////////////

/////////////////////////
void NewBollFormation(State &state, FruitingSite &site)
//     Function NewBollFormation() simulates the formation of a new boll at a
//   fruiting site. It is called from function SimulateFruitingSite().
//
//     The following global variable are set here:
//        BloomWeightLoss, BollWeight,
//        CumPlantNLoss, FruitFraction, FruitingCode.
//     The following arguments are used:
//        k, l, m - indices of vegetative branch, fruiting branch, and
//                  node on fruiting branch for this site.
{
    //     The following constant parameters are used:
    const double seedratio = 0.64; // ratio of seeds in seedcotton weight.
    const double vnewboll[2] = {0.31, 0.02};
    //     If bPollinSwitch is false accumulate number of blooms to be dropped,
    //  and define FruitingCode as 6.
    if (!state.pollination_switch)
    {
        site.stage = Stage::AbscisedAsFlower;
        site.fraction = 0;
        state.bloom_weight_loss += site.square.weight;
        site.square.weight = 0;
        return;
    }
    //     The initial weight of the new boll (BollWeight) and new burr (state.burr_weight)
    //  will be a fraction of the square weight, and the rest will be added
    //  to BloomWeightLoss. 80% of the initial weight will be in the burr.
    //     The nitrogen in the square is partitioned in the same proportions. The nitrogen
    //  that was in the square is transferred to the burrs. Update state.green_bolls_weight,
    //  state.green_bolls_burr_weight and state.square_weight. assign zero to SquareWeight at this site.
    double bolinit; // initial weight of boll after flowering.
    bolinit = vnewboll[0] * site.square.weight;
    site.boll.weight = 0.2 * bolinit;
    site.burr.weight = bolinit - site.boll.weight;
    state.bloom_weight_loss += site.square.weight - bolinit;
    //
    double sqr1n; // the nitrogen content of one square before flowering.
    sqr1n = state.square_nitrogen_concentration * site.square.weight;
    state.square_nitrogen -= sqr1n;
    state.cumulative_nitrogen_loss += sqr1n * (1 - vnewboll[0]);
    sqr1n = sqr1n * vnewboll[0];
    //
    double seed1n; // the nitrogen content of seeds in a new boll on flowering.
    seed1n = site.boll.weight * seedratio * vnewboll[1];
    if (seed1n > sqr1n)
        seed1n = sqr1n;
    state.seed_nitrogen += seed1n;
    state.burr_nitrogen += sqr1n - seed1n;
    //
    state.green_bolls_weight += site.boll.weight;
    state.green_bolls_burr_weight += site.burr.weight;
    state.square_weight -= site.square.weight;
    site.square.weight = 0;
}

/////////////////////////
void BollOpening(State &state, int k, int l, int m, unsigned int day_defoliate, double tmpboll, double plant_population, double var39, double var40, double var41, double var42)
//     Function BollOpening() simulates the transition of each fruiting site
//  from green to dehissed (open) boll. It is called from SimulateFruitingSite().
//     The following global variables are referenced here:
//        AgeOfBoll, BollWeight, FruitFraction,
//        LeafAreaIndex, VarPar
//     The following global variable are set here:
//        FibLength, FruitingCode, FibStrength, NumOpenBolls, LintYield.
//     The following arguments are used in this function:
//        k, l, m - indices of vegetative branch, fruiting branch, and
//                  node on fruiting branch for this site.
//        tmpboll - average temperature of this boll.
//
{
    FruitingSite &site = state.vegetative_branches[k].fruiting_branches[l].nodes[m];
    //     The following constant parameters are used:
    double ddpar1 = 1;
    double ddpar2 = 0.8; // constant parameters for computing fdhslai.
    double vboldhs[11] = {30.0, 41.189, -1.6057, 0.020743, 70.0,
                          0.994, 56.603, -2.921, 0.059, 1.219, 0.0065};
    //     Assign atn as the average boll temperature (tmpboll), and check that it
    //  is not higher than a maximum value.
    double atn; // modified average temperature of this boll.
    atn = tmpboll;
    if (atn > vboldhs[0])
        atn = vboldhs[0];
    //     Compute dehiss as a function of boll temperature.
    double dehiss; // days from flowering to boll opening.
    dehiss = var39 + atn * (vboldhs[1] + atn * (vboldhs[2] + atn * vboldhs[3]));
    dehiss = dehiss * var40;
    if (dehiss > vboldhs[4])
        dehiss = vboldhs[4];
    //     Dehiss is decreased after a defoliation.
    if (day_defoliate > 0 && state.daynum > day_defoliate)
        dehiss = dehiss * pow(vboldhs[5], (state.daynum - day_defoliate));
    //     If leaf area index is less than dpar1, decrease dehiss.
    if (state.leaf_area_index < ddpar1)
    {
        double fdhslai; // effect of small lai on dehiss
        fdhslai = ddpar2 + state.leaf_area_index * (1 - ddpar2) / ddpar1;
        if (fdhslai < 0)
            fdhslai = 0;
        if (fdhslai > 1)
            fdhslai = 1;
        dehiss = dehiss * fdhslai;
    }
    if (site.boll.age < dehiss)
        return;
    //     If green boll is old enough (AgeOfBoll greater than dehiss), make
    //  it an open boll, set FruitingCode to 3, and update state.open_bolls_weight, state.open_bolls_burr_weight,
    //  state.green_bolls_weight, state.green_bolls_burr_weight.
    site.stage = Stage::MatureBoll;
    state.open_bolls_weight += site.boll.weight;
    state.open_bolls_burr_weight += site.burr.weight;
    state.green_bolls_weight -= site.boll.weight;
    state.green_bolls_burr_weight -= site.burr.weight;
    //     Compute the ginning percentage as a function of boll temperature.
    //     Compute the average ginning percentage of all the bolls opened
    //  until now (state.ginning_percent).
    site.ginning_percent = (var41 - var42 * atn) / 100;
    state.ginning_percent = (state.ginning_percent * state.number_of_open_bolls + site.ginning_percent * site.fraction) /
             (state.number_of_open_bolls + site.fraction);
    //     Cumulative lint yield (LintYield) is computed in kg per ha.
    state.lint_yield += site.ginning_percent * site.boll.weight * plant_population * .001;
    //     Note: computation of fiber properties is as in GOSSYM, it is
    //  not used in COTTON2K, and it has not been tested. It is included here
    //  for compatibility, and it may be developed in future versions.
    //     fsx (fiber strength in g / tex at 1/8 inch) is computed, and
    //  averaged (as FibStrength) for all open bolls.
    //     flx (fiber length in inches, 2.5% span) is computed, and
    //  averaged (as FibLength) for all open bolls.
    double flx; // fiber length (inches, 2.5% span) of this boll.
    double fsx; // fiber strength (g / tex at 1/8 inch) of this boll.
    fsx = vboldhs[6] + atn * (vboldhs[7] + vboldhs[8] * atn);
    flx = vboldhs[9] - vboldhs[10] * atn;
    FibStrength = (FibStrength * state.number_of_open_bolls + fsx * site.fraction) / (state.number_of_open_bolls + site.fraction);
    FibLength = (FibLength * state.number_of_open_bolls + flx * site.fraction) / (state.number_of_open_bolls + site.fraction);
    //     Update the number of open bolls per plant (nopen).
    state.number_of_open_bolls += site.fraction;
}
