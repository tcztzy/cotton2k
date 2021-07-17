//  PlantGrowth_3.cpp
//
//   functions in this file:
// DryMatterBalance()
// ActualFruitGrowth()
// ActualLeafGrowth()
// CheckDryMatterBal()
//
#include <cmath>
#include <string>
#include "global.h"
#include "Simulation.hpp"

using namespace std;

///////////////////////////////////////////////////////////////////////////////
void ActualFruitGrowth(State &state)
//     This function simulates the actual growth of squares and
//  bolls of cotton plants. It is called from PlantGrowth().
//
//     The following global variables are referenced here:
//        FruitingCode, NumFruitBranches, NumNodes, NumVegBranches,
//        PotGroBolls, PotGroBurrs, PotGroSquares.
//
//     The following global variables are set here:
//        BollWeight, SquareWeight.
//
{
    //     Assign zero to all the sums to be computed.
    state.square_weight = 0;
    state.green_bolls_weight = 0;
    state.green_bolls_burr_weight = 0;
    state.actual_square_growth = 0;
    state.actual_boll_growth = 0;
    state.actual_burr_growth = 0;
    //     Begin loops over all fruiting sites.
    for (int k = 0; k < state.number_of_vegetative_branches; k++) // loop of vegetative branches
        for (int l = 0; l < state.vegetative_branches[k].number_of_fruiting_branches; l++)  // loop of fruiting branches
            for (int m = 0; m < state.vegetative_branches[k].fruiting_branches[l].number_of_fruiting_nodes; m++) // loop of nodes on a fruiting branch
            {
                FruitingSite &site = state.vegetative_branches[k].fruiting_branches[l].nodes[m];
                //     If this site is a square, the actual dry weight added to it
                //  (dwsq) is proportional to its potential growth.
                //     Update the weight of this square (SquareWeight), sum of today's added dry
                //  weight to squares (state.actual_square_growth), and total weight of squares (state.square_weight).
                if (site.stage == Stage::Square)
                {
                    double dwsq = site.square.potential_growth * state.fruit_growth_ratio; // dry weight added to square.

                    site.square.weight += dwsq;
                    state.actual_square_growth += dwsq;
                    state.square_weight += site.square.weight;
                }
                //     If this site is a green boll, the actual dry weight added to seedcotton and burrs
                //  is proportional to their respective potential growth.
                if (site.stage == Stage::GreenBoll || site.stage == Stage::YoungGreenBoll)
                {
                    double dwboll; // dry weight added to seedcotton in a boll.
                    dwboll = site.boll.potential_growth * state.fruit_growth_ratio;
                    site.boll.weight += dwboll;
                    state.actual_boll_growth += dwboll;
                    state.green_bolls_weight += site.boll.weight;
                    double dwburr; // dry weight added to the burrs in a boll.
                    dwburr = site.burr.potential_growth * state.fruit_growth_ratio;
                    site.burr.weight += dwburr;
                    state.actual_burr_growth += dwburr;
                    state.green_bolls_burr_weight += site.burr.weight;
                }
            }
}

//////////////////////////
void CheckDryMatterBal(State &state)
//     This function checks the dry matter balances in the cotton model, for diagnostic
//  purposes. The units are g per plant of dry matter. It is called from SimulateThisDay().
//     The following global variables are referenced here:
//       AbscisedLeafWeight, BloomWeightLoss,
//       CumNetPhotosynth, GreenBollsLost, Kday,
//       ReserveC, RootWeightLoss.
//     The following global variable is set here:     PlantWeight.
{
    //     PlantWeight Is the total dry weight of all plant organs, including C reserves.
    state.plant_weight = state.root_weight + state.stem_weight + state.green_bolls_weight + state.green_bolls_burr_weight + state.leaf_weight + state.petiole_weight + state.square_weight + state.open_bolls_weight + state.open_bolls_burr_weight + ReserveC;
}
