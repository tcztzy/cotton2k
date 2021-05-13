//  PlantGrowth_3.cpp
//
//   functions in this file:
// DryMatterBalance()
// ActualFruitGrowth()
// ActualLeafGrowth()
// CheckDryMatterBal()
// Defoliate()
//
#include <cmath>
#include <string>
#include "global.h"
#include "Simulation.hpp"

using namespace std;

double vratio; // ratio of carbohydrates supplied to leaf and petiole growth to their requirements.
//////////////////////////////////////////////////
void DryMatterBalance(State &state, double &cdstem, double &cdleaf, double &cdpet, double &cdroot, double per_plant_area)
//     This function computes the cotton plant dry matter (carbon) balance, its allocation to
//  growing plant parts, and carbon stress. It is called from PlantGrowth().
//     The following global variables are referenced here:
//        Kday, NetPhotosynthesis, NStressFruiting, NStressRoots,
//        PotGroAllBolls, PotGroAllBurrs, PotGroAllLeaves, PotGroAllPetioles,
//        PotGroAllRoots, PotGroAllSquares, PotGroStem, WaterStress.
//     The following global and file scope variables are set here:
//        ActualStemGrowth, CarbonAllocatedForRootGrowth, CarbonStress, ExtraCarbon,
//        ReserveC, TotalActualLeafGrowth, TotalActualPetioleGrowth, vratio.
{
    //     The following constant parameters are used:
    const double vchbal[15] = {6.0, 2.5, 1.0, 5.0, 0.20, 0.80, 0.48, 0.40,
                               0.2072, 0.60651, 0.0065, 1.10, 4.0, 0.25, 4.0};
    //     Assign values for carbohydrate requirements for growth of stems, roots, leaves, petioles,
    //  squares and bolls. Potential growth of all plant parts is modified by nitrogen stresses.
    double cdsqar; // carbohydrate requirement for square growth, g per plant per day.
    cdsqar = PotGroAllSquares * (NStressFruiting + vchbal[0]) / (vchbal[0] + 1);
    double cdboll; // carbohydrate requirement for boll and burr growth, g per plant per day.
    cdboll = (PotGroAllBolls + PotGroAllBurrs) * (NStressFruiting + vchbal[0]) / (vchbal[0] + 1);
    //     cdleaf is carbohydrate requirement for leaf growth, g per plant per day.
    cdleaf = PotGroAllLeaves * (state.nitrogen_stress_vegetative + vchbal[1]) / (vchbal[1] + 1);
    //     cdstem is carbohydrate requirement for stem growth, g per plant per day.
    cdstem = PotGroStem * (state.nitrogen_stress_vegetative + vchbal[2]) / (vchbal[2] + 1);
    //     cdroot is carbohydrate requirement for root growth, g per plant per day.
    cdroot = PotGroAllRoots * (NStressRoots + vchbal[3]) / (vchbal[3] + 1);
    //     cdpet is carbohydrate requirement for petiole growth, g per plant per day.
    cdpet = PotGroAllPetioles * (state.nitrogen_stress_vegetative + vchbal[14]) / (vchbal[14] + 1);
    double cdsum; // total carbohydrate requirement for plant growth, g per plant per day.
    cdsum = cdstem + cdleaf + cdpet + cdroot + cdsqar + cdboll;
    double cpool; // total available carbohydrates for growth (cpool, g per plant).
                  //     cpool is computed as: net photosynthesis plus a fraction (vchbal(13) ) of the
                  //  stored reserves (ReserveC).
    cpool = NetPhotosynthesis + ReserveC * vchbal[13];
    //     Compute CarbonStress as the ratio of available to required carbohydrates.
    if (cdsum <= 0)
    {
        state.carbon_stress = 1;
        return; // Exit function if cdsum is 0.
    }
    state.carbon_stress = cpool / cdsum;
    if (state.carbon_stress > 1)
        state.carbon_stress = 1;
    //     When carbohydrate supply is sufficient for growth requirements, CarbonStress will be
    //  assigned 1, and the carbohydrates actually supplied for plant growth (TotalActualLeafGrowth,
    //  TotalActualPetioleGrowth, ActualStemGrowth, CarbonAllocatedForRootGrowth, pdboll, pdsq)
    //  will be equal to the required amounts.
    double pdboll;         // amount of carbohydrates allocated to boll growth.
    double pdsq;           // amount of carbohydrates allocated to square growth.
    double xtrac1, xtrac2; // first and second components of ExtraCarbon.
    if (state.carbon_stress >= 1)
    {
        TotalActualLeafGrowth = cdleaf;
        TotalActualPetioleGrowth = cdpet;
        ActualStemGrowth = cdstem;
        CarbonAllocatedForRootGrowth = cdroot;
        pdboll = cdboll;
        pdsq = cdsqar;
        xtrac1 = 0;
    }
    //     When carbohydrate supply is less than the growth requirements, set priorities for
    //  allocation of carbohydrates.
    else
    {
        double cavail; // remaining available carbohydrates.
                       //     First priority is for fruit growth. Compute the ratio of available carbohydrates
                       //  to the requirements for boll and square growth (bsratio).
        if ((cdboll + cdsqar) > 0)
        {
            double bsratio; // ratio of available carbohydrates to the requirements for boll and square growth.
            bsratio = cpool / (cdboll + cdsqar);
            double ffr; // ratio of actual supply of carbohydrates to the requirement for boll and square growth.
                        //     The factor ffr is a function of bsratio and WaterStress. It is assumed that water stress
                        //  increases allocation of carbohydrates to bolls. Check that ffr is not less than zero, or
                        //  greater than 1 or than bsratio.
            ffr = (vchbal[5] + vchbal[6] * (1 - state.water_stress)) * bsratio;
            if (ffr < 0)
                ffr = 0;
            if (ffr > 1)
                ffr = 1;
            if (ffr > bsratio)
                ffr = bsratio;
            //     Now compute the actual carbohydrates used for boll and square growth, and the
            //  remaining available carbohydrates.
            pdboll = cdboll * ffr;
            pdsq = cdsqar * ffr;
            cavail = cpool - pdboll - pdsq;
        }
        else
        {
            cavail = cpool;
            pdboll = 0;
            pdsq = 0;
        } // if cdboll+cdsquar
          //     The next priority is for leaf and petiole growth. Compute the factor flf for leaf
          //  growth allocation, and check that it is not less than zero or greater than 1.
        if ((cdleaf + cdpet) > 0)
        {
            double flf; // ratio of actual supply of carbohydrates to the requirement for leaf growth.
            flf = vchbal[7] * cavail / (cdleaf + cdpet);
            if (flf < 0)
                flf = 0;
            if (flf > 1)
                flf = 1;
            //     Compute the actual carbohydrates used for leaf and petiole growth, and the
            //  remaining available carbohydrates.
            TotalActualLeafGrowth = cdleaf * flf;
            TotalActualPetioleGrowth = cdpet * flf;
            cavail -= (TotalActualLeafGrowth + TotalActualPetioleGrowth);
        }
        else
        {
            TotalActualLeafGrowth = 0;
            TotalActualPetioleGrowth = 0;
        } // if cdleaf+cdpet
          //     The next priority is for root growth.
        if (cdroot > 0)
        {
            double ratio; // ratio between carbohydrate supply to root and to stem growth.
                          //     At no water stress conditions, ratio is an exponential function of dry
                          //  weight of vegetative shoot (stem + leaves). This equation is based
                          //  on data from Avi Ben-Porath's PhD thesis.
                          //     ratio is modified (calibrated) by vchbal[11].
            ratio = vchbal[8] + vchbal[9] * exp(-vchbal[10] * (state.stem_weight + state.leaf_weight + TotalPetioleWeight) *
                                                per_plant_area);
            ratio = ratio * vchbal[11];
            //     rtmax is the proportion of remaining available carbohydrates that can be supplied to
            //  root growth. This is increased by water stress.
            double rtmax;
            rtmax = ratio / (ratio + 1);
            rtmax = rtmax * (1 + vchbal[12] * (1 - state.water_stress));
            if (rtmax > 1)
                rtmax = 1;
            //     Compute the factor frt for root growth allocation, as a function of rtmax, and check that
            //  it is not less than zero or greater than 1.
            double frt; // ratio of actual supply of carbohydrates to the requirement for root growth.
            frt = rtmax * cavail / cdroot;
            if (frt < 0)
                frt = 0;
            if (frt > 1)
                frt = 1;
            //     Compute the actual carbohydrates used for root growth, and the
            //  remaining available carbohydrates.
            CarbonAllocatedForRootGrowth = std::max((cdroot * frt), (cavail - cdstem));
            cavail -= CarbonAllocatedForRootGrowth;
        }
        else
            CarbonAllocatedForRootGrowth = 0;
        //     The remaining available carbohydrates are used for stem growth. Compute the
        //  factor fst and the actual carbohydrates used for stem growth.
        if (cdstem > 0)
        {
            double fst; // ratio of actual supply of carbohydrates to the requirement for stem growth.
            fst = cavail / cdstem;
            if (fst < 0)
                fst = 0;
            if (fst > 1)
                fst = 1;
            ActualStemGrowth = cdstem * fst;
        }
        else
            ActualStemGrowth = 0;
        //     If there are any remaining available unused carbohydrates, define them as xtrac1.
        if (cavail > ActualStemGrowth)
            xtrac1 = cavail - ActualStemGrowth;
        else
            xtrac1 = 0;
    } // if CarbonStress
      //     Check that the amounts of carbohydrates supplied to each organ
      //  will not be less than zero.
    if (ActualStemGrowth < 0)
        ActualStemGrowth = 0;
    if (TotalActualLeafGrowth < 0)
        TotalActualLeafGrowth = 0;
    if (TotalActualPetioleGrowth < 0)
        TotalActualPetioleGrowth = 0;
    if (CarbonAllocatedForRootGrowth < 0)
        CarbonAllocatedForRootGrowth = 0;
    if (pdboll < 0)
        pdboll = 0;
    if (pdsq < 0)
        pdsq = 0;
    //     Update the amount of reserve carbohydrates (ReserveC) in the leaves.
    ReserveC = ReserveC + NetPhotosynthesis - (ActualStemGrowth + TotalActualLeafGrowth + TotalActualPetioleGrowth + CarbonAllocatedForRootGrowth + pdboll + pdsq);
    double resmax; // maximum possible amount of carbohydrate reserves that can be stored in the leaves.
                   //     resmax is a fraction (vchbal[4])) of leaf weight. Excessive reserves are defined as xtrac2.
    resmax = vchbal[4] * state.leaf_weight;
    if (ReserveC > resmax)
    {
        xtrac2 = ReserveC - resmax;
        ReserveC = resmax;
    }
    else
        xtrac2 = 0;
    //     ExtraCarbon is computed as total excessive carbohydrates.
    state.extra_carbon = xtrac1 + xtrac2;
    //     Compute state.fruit_growth_ratio as the ratio of carbohydrates supplied to
    //  square and boll growth to their carbohydrate requirements.
    if ((PotGroAllSquares + PotGroAllBolls + PotGroAllBurrs) > 0)
        state.fruit_growth_ratio = (pdsq + pdboll) / (PotGroAllSquares + PotGroAllBolls + PotGroAllBurrs);
    else
        state.fruit_growth_ratio = 1;
    //     Compute vratio as the ratio of carbohydrates supplied to leaf
    //  and petiole growth to their carbohydrate requirements.
    if ((PotGroAllLeaves + PotGroAllPetioles) > 0)
        vratio = (TotalActualLeafGrowth + TotalActualPetioleGrowth) / (PotGroAllLeaves + PotGroAllPetioles);
    else
        vratio = 1;
}

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
//        ActualBollGrowth, ActualBurrGrowth, ActualSquareGrowth, BollWeight, BurrWeight,
//        BurrWeightGreenBolls, CottonWeightGreenBolls, SquareWeight, TotalSquareWeight.
//
{
    //     Assign zero to all the sums to be computed.
    TotalSquareWeight = 0;
    CottonWeightGreenBolls = 0;
    BurrWeightGreenBolls = 0;
    ActualSquareGrowth = 0;
    ActualBollGrowth = 0;
    ActualBurrGrowth = 0;
    //     Begin loops over all fruiting sites.
    for (int k = 0; k < state.number_of_vegetative_branches; k++) // loop of vegetative branches
        for (int l = 0; l < state.vegetative_branches[k].number_of_fruiting_branches; l++)  // loop of fruiting branches
            for (int m = 0; m < state.vegetative_branches[k].fruiting_branches[l].number_of_fruiting_nodes; m++) // loop of nodes on a fruiting branch
            {
                FruitingSite &site = state.vegetative_branches[k].fruiting_branches[l].nodes[m];
                //     If this site is a square, the actual dry weight added to it
                //  (dwsq) is proportional to its potential growth.
                //     Update the weight of this square (SquareWeight), sum of today's added dry
                //  weight to squares (ActualSquareGrowth), and total weight of squares (TotalSquareWeight).
                if (site.stage == Stage::Square)
                {
                    double dwsq = site.square.potential_growth * state.fruit_growth_ratio; // dry weight added to square.

                    site.square.weight += dwsq;
                    ActualSquareGrowth += dwsq;
                    TotalSquareWeight += site.square.weight;
                }
                //     If this site is a green boll, the actual dry weight added to seedcotton and burrs
                //  is proportional to their respective potential growth.
                if (site.stage == Stage::GreenBoll || site.stage == Stage::YoungGreenBoll)
                {
                    double dwboll; // dry weight added to seedcotton in a boll.
                    dwboll = site.boll.potential_growth * state.fruit_growth_ratio;
                    site.boll.weight += dwboll;
                    ActualBollGrowth += dwboll;
                    CottonWeightGreenBolls += site.boll.weight;
                    double dwburr; // dry weight added to the burrs in a boll.
                    dwburr = site.burr.potential_growth * state.fruit_growth_ratio;
                    site.burr.weight += dwburr;
                    ActualBurrGrowth += dwburr;
                    BurrWeightGreenBolls += site.burr.weight;
                }
            }
}

//////////////////////////
void ActualLeafGrowth(State &state)
//     This function simulates the actual growth of leaves of
//  cotton plants. It is called from PlantGrowth().
//
//     The following global and file scope variables are referenced here:
//       NumFruitBranches, NumNodes, NumVegBranches, PotGroLeafAreaNodes,
//       PotGroLeafAreaMainStem, PotGroLeafWeightNodes, PotGroPetioleWeightNodes,
//       PotGroLeafWeightMainStem, PotGroPetioleWeightMainStem, PotGroLeafAreaPreFru,
//       PotGroLeafWeightPreFru, PotGroPetioleWeightPreFru, vratio.
//     The following global variables are set here:
//       LeafAreaMainStem, LeafAreaNodes, LeafWeightMainStem,
//       LeafWeightNodes, PetioleWeightMainStem, PetioleWeightNodes,
//       PetioleWeightPreFru, TotalPetioleWeight, .
{
    //     Loop for all prefruiting node leaves. Added dry weight to each leaf is
    //  proportional to PotGroLeafWeightPreFru. Update leaf weight (state.leaf_weight_pre_fruiting) and
    //  leaf area (state.leaf_area_pre_fruiting) for each prefruiting node leaf. added dry weight to
    //  each petiole is proportional to PotGroPetioleWeightPreFru. update petiole weight
    //  (PetioleWeightPreFru) for each prefruiting node leaf.
    //     Compute total leaf weight (state.leaf_weight), total petiole
    //  weight (PetioleWeightNodes), and state.leaf_area.
    for (int j = 0; j < state.number_of_pre_fruiting_nodes; j++) // loop by prefruiting node.
    {
        state.leaf_weight_pre_fruiting[j] += PotGroLeafWeightPreFru[j] * vratio;
        state.leaf_weight += state.leaf_weight_pre_fruiting[j];
        PetioleWeightPreFru[j] += PotGroPetioleWeightPreFru[j] * vratio;
        TotalPetioleWeight += PetioleWeightPreFru[j];
        state.leaf_area_pre_fruiting[j] += PotGroLeafAreaPreFru[j] * vratio;
        state.leaf_area += state.leaf_area_pre_fruiting[j];
    }
    //     Loop for all fruiting branches on each vegetative branch, to
    //  compute actual growth of mainstem leaves.
    //     Added dry weight to each leaf is proportional to PotGroLeafWeightMainStem,
    //  added dry weight to each petiole is proportional to PotGroPetioleWeightMainStem,
    //  and added area to each leaf is proportional to PotGroLeafAreaMainStem.
    //     Update leaf weight (LeafWeightMainStem), petiole weight (PetioleWeightMainStem)
    //  and leaf area(LeafAreaMainStem) for each main stem node leaf.
    //     Update the total leaf weight (state.leaf_weight), total
    //  petiole weight (TotalPetioleWeight) and total area (state.leaf_area).
    for (int k = 0; k < state.number_of_vegetative_branches; k++) // loop of vegetative branches
        for (int l = 0; l < state.vegetative_branches[k].number_of_fruiting_branches; l++) // loop of fruiting branches
        {
            MainStemLeaf &main_stem_leaf = state.vegetative_branches[k].fruiting_branches[l].main_stem_leaf;
            main_stem_leaf.leaf_weight += main_stem_leaf.potential_growth_for_leaf_weight * vratio;
            state.leaf_weight += main_stem_leaf.leaf_weight;
            main_stem_leaf.petiole_weight += main_stem_leaf.potential_growth_for_petiole_weight * vratio;
            TotalPetioleWeight += main_stem_leaf.petiole_weight;
            main_stem_leaf.leaf_area += main_stem_leaf.potential_growth_for_leaf_area * vratio;
            state.leaf_area += main_stem_leaf.leaf_area;
            //     Loop for all fruiting nodes on each fruiting branch. to compute
            //  actual growth of fruiting node leaves.
            //     Added dry weight to each leaf is proportional to PotGroLeafWeightNodes,
            //  added dry weight to each petiole is proportional to PotGroPetioleWeightNodes,
            //  and added area to each leaf is proportional to PotGroLeafAreaNodes.
            //  Update leaf weight (LeafWeightNodes), petiole weight (PetioleWeightNodes) and leaf area
            //  (LeafAreaNodes) for each fruiting node leaf.
            //     Compute total leaf weight (state.leaf_weight), total petiole weight
            //  (PetioleWeightNodes) and total area (state.leaf_area).
            for (int m = 0; m < state.vegetative_branches[k].fruiting_branches[l].number_of_fruiting_nodes; m++) // loop of nodes on a fruiting branch
            {
                FruitingSite &site = state.vegetative_branches[k].fruiting_branches[l].nodes[m];
                site.leaf.weight += site.leaf.potential_growth * state.leaf_weight_area_ratio * vratio;
                state.leaf_weight += site.leaf.weight;
                site.petiole.weight += site.petiole.potential_growth * vratio;
                TotalPetioleWeight += site.petiole.weight;
                site.leaf.area += site.leaf.potential_growth * vratio;
                state.leaf_area += site.leaf.area;
            } // loop m
        }     // loopl
}

//////////////////////////
void CheckDryMatterBal(State &state)
//     This function checks the dry matter balances in the cotton model, for diagnostic
//  purposes. The units are g per plant of dry matter. It is called from SimulateThisDay().
//     The following global variables are referenced here:
//       AbscisedLeafWeight, BloomWeightLoss, BurrWeightGreenBolls, BurrWeightOpenBolls,
//       CottonWeightGreenBolls, CottonWeightOpenBolls, CumNetPhotosynth, GreenBollsLost, Kday,
//       ReserveC, RootWeightLoss, TotalPetioleWeight,
//       TotalRootWeight, TotalSquareWeight.
//     The following global variable is set here:     PlantWeight.
{
    //     PlantWeight Is the total dry weight of all plant organs, including C reserves.
    state.plant_weight = TotalRootWeight + state.stem_weight + CottonWeightGreenBolls + BurrWeightGreenBolls + state.leaf_weight + TotalPetioleWeight + TotalSquareWeight + CottonWeightOpenBolls + BurrWeightOpenBolls + ReserveC;
}

//////////////////////////
void Defoliate(Simulation &sim, uint32_t u)
//     This function simulates the effects of defoliating chemicals
//  applied on the cotton. It is called from SimulateThisDay().
//
//     The following global variables are referenced here:
//       LeafAreaIndex, LightIntercept,
//       NumGreenBolls, NumOpenBolls, LwpMin.
//
//     The following global variables are set here:
//       DefoliantAppRate, DefoliationDate, DefoliationMethod, PercentDefoliation.
//
{
    State &state = sim.states[u];
    //        constant parameters:
    const double p1 = -50.0;
    const double p2 = 0.525;
    const double p3 = 7.06;
    const double p4 = 0.85;
    const double p5 = 2.48;
    const double p6 = 0.0374;
    const double p7 = 0.0020;
    //
    static double defkgh; // amount of defoliant applied, kg per ha
    static double tdfkgh; // total cumulative amount of defoliant
    static int idsw;      // switch indicating if predicted defoliation date was defined.
                          //     If this is first day set initial values of tdfkgh, defkgh and idsw to 0.
    if (state.daynum <= sim.day_emerge)
    {
        tdfkgh = 0;
        defkgh = 0;
        idsw = 0;
    }
    //     Start a loop for five possible defoliant applications.
    for (int i = 0; i < 5; i++)
    {
        //     If there are open bolls and defoliation prediction has been set, execute the following.
        if (state.number_of_open_bolls > 0 && DefoliantAppRate[i] <= -99.9)
        {
            int OpenRatio; // percentage of open bolls in total boll number
            OpenRatio = (int)(100 * state.number_of_open_bolls / (state.number_of_open_bolls + state.number_of_green_bolls));
            if (i == 0 && idsw == 0)
            {
                //     If this is first defoliation - check the percentage of boll opening.
                //     If it is after the defined date, or the percent boll opening is greater than the
                //  defined threshold - set defoliation date as this day and set a second prediction.
                if ((state.daynum >= DefoliationDate[i] && DefoliationDate[0] > 0) ||
                    OpenRatio > DefoliationMethod[i])
                {
                    idsw = 1;
                    DefoliationDate[i] = state.daynum;
                    DefoliantAppRate[1] = -99.9;
                    if (state.daynum < sim.day_defoliate || sim.day_defoliate <= 0)
                        sim.day_defoliate = state.daynum;
                    DefoliationMethod[i] = 0;
                } // if Daynum
            }     // if i, idsw
                  //     If 10 days have passed since the last defoliation, and the
                  //  leaf area index is still greater than 0.2, set another defoliation.
            if (i >= 1)
            {
                if (state.daynum == (DefoliationDate[i - 1] + 10) && state.leaf_area_index >= 0.2)
                {
                    DefoliationDate[i] = state.daynum;
                    if (i < 4)
                        DefoliantAppRate[i + 1] = -99.9;
                    DefoliationMethod[i] = 0;
                } // if Daynum
            }     // if i>=1
        }         //  if NumOpenBolls
                  //
        if (state.daynum == DefoliationDate[i])
        {
            //     If it is a predicted defoliation, assign tdfkgh as 2.5 .
            //     Else, compute the amount intercepted by the plants in kg per ha
            //  (defkgh), and add it to tdfkgh.
            if (DefoliantAppRate[i] < -99)
                tdfkgh = 2.5;
            else
            {
                if (DefoliationMethod[i] == 0)
                    defkgh += DefoliantAppRate[i] * 0.95 * 1.12085 * 0.75;
                else
                    defkgh += DefoliantAppRate[i] * LightIntercept * 1.12085 * 0.75;
                tdfkgh += defkgh;
            }
        } // Daynum
          //     If this is after the first day of defoliant application, compute the
          //  percent of leaves to be defoliated (PercentDefoliation), as a function of
          //  average daily temperature, leaf water potential, days after first
          //  defoliation application, and tdfkgh. The regression equation is
          //  modified from the equation suggested in GOSSYM.
        if (DefoliationDate[i] > 0 && state.daynum > sim.day_defoliate)
        {
            double dum = -LwpMin * 10; // value of LwpMin in bars.
            PercentDefoliation = p1 + p2 * state.average_temperature + p3 * tdfkgh + p4 * (state.daynum - sim.day_defoliate) + p5 * dum - p6 * dum * dum + p7 * state.average_temperature * tdfkgh * (state.daynum - sim.day_defoliate) * dum;
            if (PercentDefoliation < 0)
                PercentDefoliation = 0;
            double perdmax = 40; // maximum possible percent of defoliation.
            if (PercentDefoliation > perdmax)
                PercentDefoliation = perdmax;
        } // if DefoliationDate
    }     // loop i
}
