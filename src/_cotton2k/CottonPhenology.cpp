//  CottonPhenology.cpp
//
//   Functions in this file:
// CottonPhenology()
// SimulateFruitingSite{}
// NewBollFormation()
//
#include <algorithm>
#include <cmath>
#include "LeafAbscission.h"
#include "FruitAbscission.h"

using namespace std;

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
