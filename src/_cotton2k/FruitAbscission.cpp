//  FruitAbscission.cpp
//
//   functions in this file:
//       FruitingSitesAbscission()
//       SiteAbscissionRatio()
//       SquareAbscission()
//       BollAbscission()
//       ComputeSiteNumbers()
//
#include <algorithm>
#include <cmath>
#include <cstdint>
#include "global.h"
#include "Simulation.hpp"

using namespace std;

double SiteAbscissionRatio(State &, int, int, int, int);

void SquareAbscission(State &, FruitingSite &, double);

void BollAbscission(State &, FruitingSite &, double, double);

extern "C"
{
    void ComputeSiteNumbers(State &);
}

//////////////////////////////////////////////////
void FruitingSitesAbscission(Simulation &sim, uint32_t u)
//     This function simulates the abscission of squares and bolls.
//  It is called from function CottonPhenology().  It calls SiteAbscissionRatio(),
//	SquareAbscission(), BollAbscission() and ComputeSiteNumbers()
//
//     The following global variables are referenced here:
//  CarbonStress, DayInc, FruitingCode, ginp, Gintot, Kday, NitrogenStress,
//  NumFruitBranches, NumNodes, NumVegBranches, VarPar, WaterStress.
//
//     The following global variable are set here:
//  AbscissionLag, NumSheddingTags, ShedByCarbonStress, ShedByNitrogenStress, ShedByWaterStress.
//
{
    State &state = sim.states[u];
//     The following constant parameters are used:
    const double vabsfr[9] = {21.0, 0.42, 30.0, 0.05, 6.0, 2.25, 0.60, 5.0, 0.20};
//
//     Update tags for shedding: Increment NumSheddingTags by 1, and move the array members
// of ShedByCarbonStress, ShedByNitrogenStress, ShedByWaterStress, and AbscissionLag.
    NumSheddingTags++;
    if (NumSheddingTags > 1)
        for (int lt = NumSheddingTags - 1; lt > 0; lt--) {
            int ltm1 = lt - 1;
            ShedByCarbonStress[lt] = ShedByCarbonStress[ltm1];
            ShedByNitrogenStress[lt] = ShedByNitrogenStress[ltm1];
            ShedByWaterStress[lt] = ShedByWaterStress[ltm1];
            AbscissionLag[lt] = AbscissionLag[ltm1];
        }
//     Calculate shedding intensity: The shedding intensity due to stresses of this day is assigned
//  to the first members of the arrays ShedByCarbonStress, ShedByNitrogenStress, and ShedByWaterStress.
    if (state.carbon_stress < VarPar[43])
        ShedByCarbonStress[0] = (VarPar[43] - state.carbon_stress) / VarPar[43];
    else
        ShedByCarbonStress[0] = 0;
    if (state.nitrogen_stress < vabsfr[1])
        ShedByNitrogenStress[0] = (vabsfr[1] - state.nitrogen_stress) / vabsfr[1];
    else
        ShedByNitrogenStress[0] = 0;
    if (state.water_stress < VarPar[44])
        ShedByWaterStress[0] = (VarPar[44] - state.water_stress) / VarPar[44];
    else
        ShedByWaterStress[0] = 0;
//     Assign 0.01 to the first member of AbscissionLag.
    AbscissionLag[0] = 0.01;
//     Updating age of tags for shedding: Each member of array AbscissionLag is
//  incremented by physiological age of today. It is further increased (i.e., shedding
//  will occur sooner) when maximum temperatures are high.
    double tmax = sim.climate[u].Tmax;
    for (int lt = 0; lt < NumSheddingTags; lt++) {
        AbscissionLag[lt] += std::max(sim.states[u].day_inc, 0.40);
        if (tmax > vabsfr[2])
            AbscissionLag[lt] += (tmax - vabsfr[2]) * vabsfr[3];
    }
//     Assign zero to idecr, and start do loop over all days since
//  first relevant tagging. If AbscissionLag reaches a value of vabsfr[4],
//  calculate actual shedding of each site:
    int idecr = 0; // decrease in NumSheddingTags after shedding has been executed.
    for (int lt = 0; lt < NumSheddingTags; lt++) {
        if (AbscissionLag[lt] >= vabsfr[4] || lt >= 20) {
            double gin1;
            if (Gintot > 0)
                gin1 = Gintot;
            else
                gin1 = ginp;
//      Start loop over all possible fruiting sites. The abscission functions
//  will be called for sites that are squares or green bolls.
            for (int k = 0; k < state.number_of_vegetative_branches; k++)
                for (int l = 0; l < state.vegetative_branches[k].number_of_fruiting_branches; l++)
                    for (int m = 0; m < state.vegetative_branches[k].fruiting_branches[l].number_of_fruiting_nodes; m++) {
                        FruitingSite &site = state.vegetative_branches[k].fruiting_branches[l].nodes[m];
                        if (site.stage == Stage::Square || site.stage == Stage::YoungGreenBoll || site.stage == Stage::GreenBoll) {
                            double abscissionRatio; // ratio of abscission for a fruiting site.
                            abscissionRatio = SiteAbscissionRatio(state, k, l, m, lt);
                            if (abscissionRatio > 0) {
                                if (site.stage == Stage::Square)
                                    SquareAbscission(state, site, abscissionRatio);
                                else
                                    BollAbscission(state, site, abscissionRatio, gin1);
                            }
                        }
                    }
//      Assign zero to the array members for this day.
            ShedByCarbonStress[lt] = 0;
            ShedByNitrogenStress[lt] = 0;
            ShedByWaterStress[lt] = 0;
            AbscissionLag[lt] = 0;
            idecr++;
        } // if AbscissionLag
    } // for lt
//  Decrease NumSheddingTags. If plantmap adjustments are necessary for square
// number, or green boll number, or open boll number - call AdjustAbscission().
    NumSheddingTags = NumSheddingTags - idecr;
//
//
    ComputeSiteNumbers(state);
}

/////////////////////////
double SiteAbscissionRatio(State &state, int k, int l, int m, int lt)
//     This function computes and returns the probability of abscission of a single
//  site (k, l, m). It is called from function FruitingSitesAbscission().
//
//     The following global variables are referenced here:
//        AgeOfSite, AgeOfBoll, FruitingCode, ShedByCarbonStress,
//        ShedByNitrogenStress, ShedByWaterStress, VarPar
//
//     The following arguments are used here:
//        k, l, m - indices defining position of this site.
//        lt - lag index for this node.
//
{
    FruitingSite &site = state.vegetative_branches[k].fruiting_branches[l].nodes[m];
//     The following constant parameters are used:
    const double vabsc[5] = {21.0, 2.25, 0.60, 5.0, 0.20};
//
//     For each site, compute the probability of its abscission (pabs)
//  as afunction of site age, and the total shedding ratio (shedt) as a
//  function of plant stresses that occurred when abscission was
//  triggered.
    double pabs = 0;  // probability of abscission of a fruiting site.
    double shedt = 0; // total shedding ratio, caused by various stresses.
//     (1) Squares (FruitingCode = 1).
    if (site.stage == Stage::Square) {
        if (site.age < vabsc[3])
            pabs = 0; // No abscission of very young squares (AgeOfSite less than vabsc(3))
        else {
            double xsqage; // square age after becoming susceptible to shedding.
            xsqage = site.age - vabsc[3];
            if (xsqage >= vabsc[0])
                pabs = VarPar[46]; // Old squares have a constant probability of shedding.
            else
//     Between these limits, pabs is a function of xsqage.
                pabs = VarPar[46] + (VarPar[45] - VarPar[46])
                                    * pow(((vabsc[0] - xsqage) / vabsc[0]), vabsc[1]);
        }
//     Total shedding ratio (shedt) is a product of the effects of
//  carbohydrate stress and nitrogen stress.
        shedt = 1 - (1 - ShedByCarbonStress[lt]) * (1 - ShedByNitrogenStress[lt]);
    }
//     (2) Very young bolls (FruitingCode = 7, and AgeOfBoll less than VarPar[47]).
    else if (site.stage == Stage::YoungGreenBoll && site.boll.age <= VarPar[47]) {
//     There is a constant probability of shedding (VarPar[48]), and shedt is a product
//  of the effects carbohydrate, and nitrogen stresses. Note that nitrogen stress has only a
//  partial effect in this case, as modified by vabsc[2].
        pabs = VarPar[48];
        shedt = 1 - (1 - ShedByCarbonStress[lt]) * (1 - vabsc[2] * ShedByNitrogenStress[lt]);
    }
//     (3) Medium age bolls (AgeOfBoll between VarPar[47] and VarPar[47] + VarPar[49]).
    else if (site.boll.age > VarPar[47] && site.boll.age <= (VarPar[47] + VarPar[49])) {
//     pabs is linearly decreasing with age, and shedt is a product of the effects
//  carbohydrate, nitrogen and water stresses.  Note that nitrogen stress has only
//  a partial effect in this case, as modified by vabsc[4].
        pabs = VarPar[48] - (VarPar[48] - VarPar[50]) * (site.boll.age - VarPar[47]) / VarPar[49];
        shedt = 1 -
                (1 - ShedByCarbonStress[lt]) * (1 - vabsc[4] * ShedByNitrogenStress[lt]) * (1 - ShedByWaterStress[lt]);
    }
//     (4) Older bolls (AgeOfBoll between VarPar[47] + VarPar[49] and VarPar[47] + 2*VarPar[49]).
    else if (site.boll.age > (VarPar[47] + VarPar[49]) && site.boll.age <= (VarPar[47] + 2 * VarPar[49])) {
//     pabs is linearly decreasing with age, and shedt is affected only by water stress.
        pabs = VarPar[50] / VarPar[49] * (VarPar[47] + 2 * VarPar[49] - site.boll.age);
        shedt = ShedByWaterStress[lt];
    }
//     (5) bolls older than VarPar[47] + 2*VarPar[49]
    else if (site.boll.age > (VarPar[47] + 2 * VarPar[49]))
        pabs = 0; // no abscission
//      Actual abscission of tagged sites (abscissionRatio) is a product of pabs,
//  shedt and DayInc for this day. It can not be greater than 1.
    double abscissionRatio = pabs * shedt * state.day_inc;
    if (abscissionRatio > 1)
        abscissionRatio = 1;
    return abscissionRatio;
}

//////////////////////////////////////////////////////////////////
void SquareAbscission(State &state, FruitingSite &site, double abscissionRatio)
//     This function simulates the abscission of a single square
//  at site (k, l, m). It is called from function FruitingSitesAbscission()
//  if this site is a square.
//
//     The following global variable is referenced here:   SquareNConc
//
//     The following global variable are set here:
//   BloomWeightLoss, CumPlantNLoss, FruitingCode, FruitFraction, SquareNitrogen,
//   SquareWeight, TotalSquareWeight.
//
//     The following arguments are used in this function:
//        abscissionRatio - ratio of abscission of a fruiting site.
//        k, l, m - indices defining position of this site.
//
{
//     Compute the square weight lost by shedding (wtlos) as a proportion of SquareWeight
//  of this site. Update SquareNitrogen, CumPlantNLoss, SquareWeight[k][l][m], BloomWeightLoss,
//  TotalSquareWeight, and FruitFraction[k][l][m].
    double wtlos = site.square.weight * abscissionRatio; // weight lost by shedding at this site.
    SquareNitrogen -= wtlos * SquareNConc;
    state.cumulative_nitrogen_loss += wtlos * SquareNConc;
    site.square.weight -= wtlos;
    state.bloom_weight_loss += wtlos;
    TotalSquareWeight -= wtlos;
    site.fraction *= (1 - abscissionRatio);
//     If FruitFraction[k][l][m] is less than 0.001 make it zero, and update
//  SquareNitrogen, CumPlantNLoss, BloomWeightLoss, TotalSquareWeight, SquareWeight[k][l][m],
//  and assign 5 to FruitingCode.
    if (site.fraction <= 0.001) {
        site.fraction = 0;
        SquareNitrogen -= site.square.weight * SquareNConc;
        state.cumulative_nitrogen_loss += site.square.weight * SquareNConc;
        state.bloom_weight_loss += site.square.weight;
        TotalSquareWeight -= site.square.weight;
        site.square.weight = 0;
        site.stage = Stage::AbscisedAsSquare;
    }
}

////////////////////////
void BollAbscission(State &state, FruitingSite &site, double abscissionRatio, double gin1)
//     This function simulates the abscission of a single green boll
//  at site (k, l, m). It is called from function FruitingSitesAbscission() if this site
//  is a green boll.
//
//     The following global variables are referenced here:
//   BurrNConc, SeedNConc
//
//     The following global variable are set here:
//   BollWeight, BurrNitrogen, BurrWeight, BurrWeightGreenBolls, CottonWeightGreenBolls,
//   CumPlantNLoss, FruitFraction, FruitingCode, FruitFraction, GreenBollsLost, SeedNitrogen,
//
//     The following arguments are used in this function:
//        abscissionRatio - ratio of abscission of a fruiting site.
//        gin1 - percent of seeds in seedcotton, used to compute lost nitrogen.
//        k, l, m - location of this site on the plant
//
{
//     Update SeedNitrogen, BurrNitrogen, CumPlantNLoss, GreenBollsLost, CottonWeightGreenBolls, BurrWeightGreenBolls,
//  BollWeight[k][l][m], BurrWeight[k][l][m], and FruitFraction[k][l][m].
    SeedNitrogen -= site.boll.weight * abscissionRatio * (1 - gin1) * state.seed_nitrogen_concentration;
    BurrNitrogen -= site.burr.weight * abscissionRatio * BurrNConc;
    state.cumulative_nitrogen_loss += site.boll.weight * abscissionRatio * (1. - gin1) * state.seed_nitrogen_concentration;
    state.cumulative_nitrogen_loss += site.burr.weight * abscissionRatio * BurrNConc;
    GreenBollsLost += (site.boll.weight + site.burr.weight) * abscissionRatio;
    CottonWeightGreenBolls -= site.boll.weight * abscissionRatio;
    BurrWeightGreenBolls -= site.burr.weight * abscissionRatio;
    site.boll.weight -= site.boll.weight * abscissionRatio;
    site.burr.weight -= site.burr.weight * abscissionRatio;
    site.fraction -= site.fraction * abscissionRatio;
//
//     If FruitFraction[k][l][m] is less than 0.001 make it zero, update SeedNitrogen,
//  BurrNitrogen, CumPlantNLoss, CottonWeightGreenBolls, BurrWeightGreenBolls, GreenBollsLost,
//  BollWeight[k][l][m], BurrWeight[k][l][m], and assign 4 to FruitingCode.
//
    if (site.fraction <= 0.001) {
        site.stage = Stage::AbscisedAsBoll;
        SeedNitrogen -= site.boll.weight * (1 - gin1) * state.seed_nitrogen_concentration;
        BurrNitrogen -= site.burr.weight * BurrNConc;
        state.cumulative_nitrogen_loss += site.boll.weight * (1 - gin1) * state.seed_nitrogen_concentration;
        state.cumulative_nitrogen_loss += site.burr.weight * BurrNConc;
        site.fraction = 0;
        CottonWeightGreenBolls -= site.boll.weight;
        BurrWeightGreenBolls -= site.burr.weight;
        GreenBollsLost += site.boll.weight + site.burr.weight;
        site.boll.weight = 0;
        site.burr.weight = 0;
    }
}
