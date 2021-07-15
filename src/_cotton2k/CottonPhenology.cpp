//  CottonPhenology.cpp
//
//   Functions in this file:
// CottonPhenology()
// AddFruitingNode()
// SimulateFruitingSite{}
// NewBollFormation()
// BollOpening()
//
#include <algorithm>
#include <cmath>
#include "LeafAbscission.h"
#include "FruitAbscission.h"

using namespace std;

void NewBollFormation(State &, FruitingSite &);

void BollOpening(Simulation &, uint32_t, int, int, int, double);

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
/////////////////////////////////////////////////////////////////////////////////////////
void AddFruitingNode(State &state, int k, int l, double delayFrtByCStress, double stemNRatio, double density_factor, double VarPar[61], double PhenDelayByNStress)
//     Function AddFruitingNode() decides if a new node is to be added to a fruiting branch,
//  and forms it. It is called from function CottonPhenology().
//     The following global variables are referenced here:
//  AdjAddSitesRate, AgeOfSite, WaterStress,
//     The following global variable are set here:
//        AvrgNodeTemper, DelayNewNode, FruitFraction, FruitingCode, LeafAreaNodes, LeafWeightNodes,
//        NumNodes.
//     The following arguments are used in this function:
//        delayFrtByCStress - delay caused by carbohydrate and nitrogen stresses.
//        k, l - indices of this vegetative branch and fruiting branch.
//        stemNRatio - the ratio of n to dm in the stems.
//
{
    //     The following constant parameters are used:
    const double vfrtnod[6] = {1.32, 0.90, 33.0, 7.6725, -0.3297, 0.004657};
    //      Compute the cumulative delay for the appearance of the next
    //  node on the fruiting branch, caused by carbohydrate, nitrogen, and
    //  water stresses.
    state.vegetative_branches[k].fruiting_branches[l].delay_for_new_node += delayFrtByCStress + vfrtnod[0] * PhenDelayByNStress;
    state.vegetative_branches[k].fruiting_branches[l].delay_for_new_node += vfrtnod[1] * (1 - state.water_stress);
    //     Define nnid, and compute the average temperature of the last
    //  node of this fruiting branch, from the time it was formed.
    int nnid = state.vegetative_branches[k].fruiting_branches[l].number_of_fruiting_nodes - 1;     // the number of the last node on this fruiting branche.
    double tav = state.vegetative_branches[k].fruiting_branches[l].nodes[nnid].average_temperature; // modified daily average temperature.
    if (tav > vfrtnod[2])
        tav = vfrtnod[2];
    //     Compute TimeToNextFruNode, the time (in physiological days) needed for the
    //  formation of each successive node on the fruiting branch. This is
    //  a function of temperature, derived from data of K. R. Reddy, CSRU,
    //  adjusted for age in physiological days. It is modified for plant density.
    double TimeToNextFruNode; // time, in physiological days, for the next node on the fruiting branch to be formed
    TimeToNextFruNode = VarPar[36] + tav * (vfrtnod[3] + tav * (vfrtnod[4] + tav * vfrtnod[5]));
    TimeToNextFruNode = TimeToNextFruNode * (1 + VarPar[37] * (1 - density_factor)) + state.vegetative_branches[k].fruiting_branches[l].delay_for_new_node;
    //     Check if the the age of the last node on the fruiting branch exceeds TimeToNextFruNode.
    //  If so, form the new node:
    if (state.vegetative_branches[k].fruiting_branches[l].nodes[nnid].age < TimeToNextFruNode)
        return;
    //     Increment NumNodes, define newnod, and assign 1 to FruitFraction and FruitingCode.
    state.vegetative_branches[k].fruiting_branches[l].number_of_fruiting_nodes++;
    if (state.vegetative_branches[k].fruiting_branches[l].number_of_fruiting_nodes > 5)
    {
        state.vegetative_branches[k].fruiting_branches[l].number_of_fruiting_nodes = 5;
        return;
    }
    int newnod = nnid + 1; // the number of the new node on this fruiting branche.
    state.vegetative_branches[k].fruiting_branches[l].nodes[newnod].fraction = 1;
    state.vegetative_branches[k].fruiting_branches[l].nodes[newnod].stage = Stage::Square;
    //     Initiate a new leaf at the new node. The mass and nitrogen in
    //  the new leaf is substacted from the stem.
    state.vegetative_branches[k].fruiting_branches[l].nodes[newnod].leaf.area = VarPar[34];
    state.vegetative_branches[k].fruiting_branches[l].nodes[newnod].leaf.weight = VarPar[34] * state.leaf_weight_area_ratio;
    state.stem_weight -= state.vegetative_branches[k].fruiting_branches[l].nodes[newnod].leaf.weight;
    state.leaf_weight += state.vegetative_branches[k].fruiting_branches[l].nodes[newnod].leaf.weight;
    state.leaf_nitrogen += state.vegetative_branches[k].fruiting_branches[l].nodes[newnod].leaf.weight * stemNRatio;
    state.stem_nitrogen -= state.vegetative_branches[k].fruiting_branches[l].nodes[newnod].leaf.weight * stemNRatio;
    //     Begin computing AvrgNodeTemper of the new node, and assign zero to DelayNewNode.
    state.vegetative_branches[k].fruiting_branches[l].nodes[newnod].average_temperature = state.average_temperature;
    state.vegetative_branches[k].fruiting_branches[l].delay_for_new_node = 0;
}

//////////////////////////////////////////////////
void SimulateFruitingSite(Simulation &sim, uint32_t u, int k, int l, int m, int &NodeRecentWhiteFlower, const double &WaterStress)
//     Function SimulateFruitingSite() simulates the development of each fruiting site.
//  It is called from function CottonPhenology().
//     The following global variable are set here:
//        AgeOfSite, AgeOfBoll, AvrgNodeTemper, BollWeight,
//        FirstBloom, FruitingCode, LeafAge, NumFruitSites,
//     The following arguments are used in this function:
//        k, l, m - indices of vegetative branch, fruiting branch, and
//                  node on fruiting branch for this site.
//        NodeRecentWhiteFlower - the node of the most recent white flower is computed in this
//                  function, although it is not used in this version.
//
{
    State &state = sim.states[u];
    FruitingSite &site = state.vegetative_branches[k].fruiting_branches[l].nodes[m];
    //     The following constant parameters are used:
    const double vfrsite[15] = {0.60, 0.40, 12.25, 0.40, 33.0, 0.20, 0.04, 0.45,
                                26.10, 9.0, 0.10, 3.0, 1.129, 0.043, 0.26};
    //      FruitingCode = 0 indicates that this node has not yet been formed.
    //  In this case, assign zero to boltmp and return.
    static double boltmp[3][30][5]; // cumulative boll temperature.
    if (site.stage == Stage::NotYetFormed)
    {
        boltmp[k][l][m] = 0;
        return;
    }
    //      Assign zero to FibLength and FibStrength before any sites have been formed.
    if (state.number_of_fruiting_sites <= 0)
    {
        FibLength = 0;
        FibStrength = 0;
    }
    state.number_of_fruiting_sites++; //      Increment site number.
                     //     LeafAge(k,l,m) is the age of the leaf at this site. it is updated
                     //  by adding the physiological age of this day, the effect of water
                     //  and nitrogen stresses (agefac).
    double agefac;   // effect of water and nitrogen stresses on leaf aging.
    agefac = (1 - WaterStress) * vfrsite[0] + (1 - state.nitrogen_stress_vegetative) * vfrsite[1];
    site.leaf.age += state.day_inc + agefac;
    //  After the application of defoliation, add the effect of defoliation on leaf age.
    if (sim.day_defoliate > 0 && state.daynum > sim.day_defoliate)
        site.leaf.age += sim.cultivar_parameters[38];
    //     FruitingCode = 3, 4, 5 or 6 indicates that this node has an open boll,
    //  or has been completely abscised. Return in this case.
    if (site.stage == Stage::MatureBoll || site.stage == Stage::AbscisedAsBoll || site.stage == Stage::AbscisedAsSquare || site.stage == Stage::AbscisedAsFlower)
        return;
    //     Age of node is modified for low minimum temperatures and for high
    //  maximum temperatures.
    double tmin = sim.climate[u].Tmin;
    double tmax = sim.climate[u].Tmax;
    double ageinc = state.day_inc; // daily addition to site age.
                                   //     Adjust leaf aging for low minimum twmperatures.
    if (tmin < vfrsite[2])
        ageinc += vfrsite[3] * (vfrsite[2] - tmin);
    //     Adjust leaf aging for high maximum twmperatures.
    if (tmax > vfrsite[4])
    {
        double adjust = vfrsite[6] * (tmax - vfrsite[4]); // vfrsite [4] = 33.0  [5] 0.20  [6] 0.04
        if (adjust > vfrsite[5])
            adjust = vfrsite[5];
        ageinc -= adjust;
    }
    //
    if (ageinc < vfrsite[7])
        ageinc = vfrsite[7];
    //     Compute average temperature of this site since formation.
    site.average_temperature = (site.average_temperature * site.age + state.average_temperature * ageinc) / (site.age + ageinc);
    //     Update the age of this node, AgeOfSite(k,l,m), by adding ageinc.
    site.age += ageinc;
    //     The following is executed if this node is a square (FruitingCode =  1):
    //     If square is old enough, make it a green boll: initialize the
    //  computations of average boll temperature (boltmp) and boll age
    //  (AgeOfBoll). FruitingCode will now be 7.
    if (site.stage == Stage::Square)
    {
        if (site.age >= vfrsite[8])
        {
            boltmp[k][l][m] = state.average_temperature;
            site.boll.age = state.day_inc;
            site.stage = Stage::YoungGreenBoll;
            NewBollFormation(state, state.vegetative_branches[k].fruiting_branches[l].nodes[m]);
            //     If this is the first flower, define FirstBloom.
            if (state.green_bolls_weight > 0 && sim.first_bloom <= 1)
                sim.first_bloom = state.daynum;
            //     Determine node of most recent white flower.
            if (k == 0 && m == 0)
                NodeRecentWhiteFlower = std::max(NodeRecentWhiteFlower, l);
        }
        return;
    }
    //     If there is a boll at this site:
    //     Calculate average boll temperature (boltmp), and boll age
    //  (AgeOfBoll) which is its physiological age, modified by water stress.
    //  If leaf area index is low, dum is calculated as an intermediate
    //  variable. It is used to increase boll temperature and to accelerate
    //  boll aging when leaf cover is decreased. Boll age is also modified
    //  by nitrogen stress (state.nitrogen_stress_fruiting).
    if (site.boll.weight > 0)
    {
        double dum; // effect of leaf area index on boll temperature and age.
        if (state.leaf_area_index <= vfrsite[11] && state.kday > 100)
            dum = vfrsite[12] - vfrsite[13] * state.leaf_area_index;
        else
            dum = 1;
        double dagebol; // added physiological age of boll on this day.
        dagebol = state.day_inc * dum + vfrsite[14] * (1 - WaterStress) + vfrsite[10] * (1 - state.nitrogen_stress_fruiting);
        boltmp[k][l][m] = (boltmp[k][l][m] * site.boll.age + state.average_temperature * dagebol) / (site.boll.age + dagebol);
        site.boll.age += dagebol;
    }
    //     If this node is a young green boll (FruitingCode = 7):
    //     Check boll age and after a fixed age convert it to an "old"
    //  green boll (FruitingCode = 2).
    if (site.stage == Stage::YoungGreenBoll)
    {
        if (site.boll.age >= vfrsite[9])
            site.stage = Stage::GreenBoll;
        return;
    }
    //     If this node is an older green boll (FruitingCode = 2):
    if (site.stage == Stage::GreenBoll)
        BollOpening(sim, u, k, l, m, boltmp[k][l][m]);
}

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
void BollOpening(Simulation &sim, uint32_t u, int k, int l, int m, double tmpboll)
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
    State &state = sim.states[u];
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
    dehiss = sim.cultivar_parameters[39] + atn * (vboldhs[1] + atn * (vboldhs[2] + atn * vboldhs[3]));
    dehiss = dehiss * sim.cultivar_parameters[40];
    if (dehiss > vboldhs[4])
        dehiss = vboldhs[4];
    //     Dehiss is decreased after a defoliation.
    if (sim.day_defoliate > 0 && state.daynum > sim.day_defoliate)
        dehiss = dehiss * pow(vboldhs[5], (state.daynum - sim.day_defoliate));
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
    site.ginning_percent = (sim.cultivar_parameters[41] - sim.cultivar_parameters[42] * atn) / 100;
    state.ginning_percent = (state.ginning_percent * state.number_of_open_bolls + site.ginning_percent * site.fraction) /
             (state.number_of_open_bolls + site.fraction);
    //     Cumulative lint yield (LintYield) is computed in kg per ha.
    state.lint_yield += site.ginning_percent * site.boll.weight * sim.plant_population * .001;
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
