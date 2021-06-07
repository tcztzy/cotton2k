// File  PlantGrowth_2.cpp
//
//   functions in this file:
// PotentialStemGrowth()
// PotentialLeafGrowth()
// TemperatureOnLeafGrowthRate()
// TemperatureOnFruitGrowthRate()
//
#include <cmath>
#include "global.h"
#include "State.hpp"

using namespace std;

extern "C"
{
    double TemperatureOnLeafGrowthRate(double);
}

//////////////////////////
void PotentialLeafGrowth(State &state, double density_factor, double VarPar[61])
//     This function simulates the potential growth of leaves of cotton plants.
//  It is called from PlantGrowth(). It calls function TemperatureOnLeafGrowthRate().
//
//     The following monomolecular growth function is used :
//       leaf area = smax * (1 - exp(-c * pow(t,p)))
//     where    smax = maximum leaf area.
//              t = time (leaf age).
//              c, p = constant parameters.
//     Note: p is constant for all leaves, whereas smax and c depend on leaf position.
//     The rate per day (the derivative of this function) is :
//           r = smax * c * p * exp(-c * pow(t,p)) * pow(t, (p-1))
//
//     The following global variables are referenced here:
//        LeafAreaNodes,
//        LeafAreaMainStem, NumFruitBranches, NumNodes,
//        NumVegBranches, WaterStress.
//     The following global variables are set here:
//        PotGroAllLeaves, PotGroAllPetioles,
//        PotGroLeafAreaMainStem, PotGroLeafAreaNodes,  PotGroLeafAreaPreFru,
//        PotGroLeafWeightMainStem, PotGroLeafWeightNodes, PotGroLeafWeightPreFru,
//        PotGroPetioleWeightMainStem, PotGroPetioleWeightNodes, PotGroPetioleWeightPreFru.
{
    //     The following constant parameters. are used in this function:
    const double p = 1.6; // parameter of the leaf growth rate equation.
    const double vpotlf[14] = {3.0, 0.95, 1.2, 13.5, -0.62143, 0.109365,
                               0.00137566, 0.025, 0.00005, 30., 0.02, 0.001, 2.50, 0.18};
    //     Calculate water stress reduction factor for leaf growth rate
    // (wstrlf). This has been empirically calibrated in COTTON2K.
    double wstrlf = state.water_stress * (1 + vpotlf[0] * (2 - state.water_stress)) - vpotlf[0];
    if (wstrlf < 0.05)
        wstrlf = 0.05;
    //     Calculate wtfstrs, the effect of leaf water stress on state.leaf_weight_area_ratio (the ratio
    //  of leaf dry weight to leaf area). This has also been empirically calibrated in COTTON2K.
    double wtfstrs = vpotlf[1] + vpotlf[2] * (1 - wstrlf);
    //     Compute the ratio of leaf dry weight increment to leaf area increment (g per dm2),
    //  as a function of average daily temperature and water stress. Parameters for the effect
    //  of temperature are adapted from GOSSYM.
    double tdday = state.average_temperature; // limited value of today's average temperature.
    if (tdday < vpotlf[3])
        tdday = vpotlf[3];
    state.leaf_weight_area_ratio = wtfstrs / (vpotlf[4] + tdday * (vpotlf[5] - tdday * vpotlf[6]));
    //     Assign zero to total potential growth of leaf and petiole.
    PotGroAllLeaves = 0;
    PotGroAllPetioles = 0;
    double c = 0;                            // parameter of the leaf growth rate equation.
    double smax = 0;                         // maximum possible leaf area, a parameter of the leaf growth rate equation.
    double rate;                             // growth rate of area of a leaf.
                                             //     Compute the potential growth rate of prefruiting leaves.
                                             //  smax and c are functions of prefruiting node number.
    for (int j = 0; j < state.number_of_pre_fruiting_nodes; j++) // loop by prefruiting node.
    {
        if (state.leaf_area_pre_fruiting[j] <= 0)
        {
            PotGroLeafAreaPreFru[j] = 0;
            PotGroLeafWeightPreFru[j] = 0;
            PotGroPetioleWeightPreFru[j] = 0;
        }
        else
        {
            int jp1 = j + 1;
            smax = jp1 * (VarPar[2] - VarPar[3] * jp1);
            if (smax < VarPar[4])
                smax = VarPar[4];
            c = vpotlf[7] + vpotlf[8] * jp1 * (jp1 - vpotlf[9]);
            rate = smax * c * p * exp(-c * pow(state.age_of_pre_fruiting_nodes[j], p)) * pow(state.age_of_pre_fruiting_nodes[j], (p - 1));
            //     Growth rate is modified by water stress and a function of average temperature.
            //     Compute potential growth of leaf area, leaf weight and petiole
            //  weight for leaf on node j. Add leaf weight potential growth to PotGroAllLeaves.
            //  Add potential growth of petiole weight to PotGroAllPetioles.
            if (rate >= 1e-12)
            {
                PotGroLeafAreaPreFru[j] = rate * wstrlf * TemperatureOnLeafGrowthRate(state.average_temperature);
                PotGroLeafWeightPreFru[j] = PotGroLeafAreaPreFru[j] * state.leaf_weight_area_ratio;
                PotGroPetioleWeightPreFru[j] = PotGroLeafAreaPreFru[j] * state.leaf_weight_area_ratio * vpotlf[13];
                PotGroAllLeaves += PotGroLeafWeightPreFru[j];
                PotGroAllPetioles += PotGroPetioleWeightPreFru[j];
            } // rate
        }     // state.leaf_area_pre_fruiting
    }         // state.number_of_pre_fruiting_nodes
              //     denfac is the effect of plant density on leaf growth rate.
    double denfac = 1 - vpotlf[12] * (1 - density_factor);
    for (int k = 0; k < state.number_of_vegetative_branches; k++) // loop of vegetative branches
    {
        for (int l = 0; l < state.vegetative_branches[k].number_of_fruiting_branches; l++) // loop of fruiting branches
        {
            //     smax and c are  functions of fruiting branch number.
            //     smax is modified by plant density, using the density factor denfac.
            //     Compute potential main stem leaf growth, assuming that the main
            //  stem leaf is initiated at the same time as leaf (k,l,0).
            MainStemLeaf &main_stem_leaf = state.vegetative_branches[k].fruiting_branches[l].main_stem_leaf;
            if (main_stem_leaf.leaf_area <= 0)
            {
                main_stem_leaf.potential_growth_for_leaf_area = 0;
                main_stem_leaf.potential_growth_for_leaf_weight = 0;
                main_stem_leaf.potential_growth_for_petiole_weight = 0;
            }
            else
            {
                int lp1 = l + 1;
                smax = VarPar[5] + VarPar[6] * lp1 * (VarPar[7] - lp1);
                smax = smax * denfac;
                if (smax < VarPar[4])
                    smax = VarPar[4];
                c = vpotlf[10] + lp1 * vpotlf[11];
                if (state.vegetative_branches[k].fruiting_branches[l].nodes[0].leaf.age > 70)
                    rate = 0;
                else
                    rate = smax * c * p * exp(-c * pow(state.vegetative_branches[k].fruiting_branches[l].nodes[0].leaf.age, p)) * pow(state.vegetative_branches[k].fruiting_branches[l].nodes[0].leaf.age, (p - 1));
                //     Add leaf and petiole weight potential growth to SPDWL and SPDWP.
                if (rate >= 1e-12)
                {
                    main_stem_leaf.potential_growth_for_leaf_area = rate * wstrlf * TemperatureOnLeafGrowthRate(state.average_temperature);
                    main_stem_leaf.potential_growth_for_leaf_weight = main_stem_leaf.potential_growth_for_leaf_area * state.leaf_weight_area_ratio;
                    main_stem_leaf.potential_growth_for_petiole_weight = main_stem_leaf.potential_growth_for_leaf_area * state.leaf_weight_area_ratio * vpotlf[13];
                    PotGroAllLeaves += main_stem_leaf.potential_growth_for_leaf_weight;
                    PotGroAllPetioles += main_stem_leaf.potential_growth_for_petiole_weight;
                }
            }
            //     Assign smax value of this main stem leaf to smaxx, c to cc.
            //     Loop over the nodes of this fruiting branch.
            double smaxx = smax;                                                                                 // value of smax for the corresponding main stem leaf.
            double cc = c;                                                                                       // value of c for the corresponding main stem leaf.
            for (int m = 0; m < state.vegetative_branches[k].fruiting_branches[l].number_of_fruiting_nodes; m++) // loop of nodes on a fruiting branch
            {
                FruitingSite &site = state.vegetative_branches[k].fruiting_branches[l].nodes[m];
                if (site.leaf.area <= 0)
                {
                    site.leaf.potential_growth = 0;
                    site.petiole.potential_growth = 0;
                }
                //     Compute potential growth of leaf area and leaf weight for leaf
                //  on fruiting branch node (k,l,m).
                //    Add leaf and petiole weight potential growth to spdwl and spdwp.
                else
                {
                    int mp1 = m + 1;
                    //     smax and c are reduced as a function of node number on this fruiting branch.
                    smax = smaxx * (1 - VarPar[8] * mp1);
                    c = cc * (1 - VarPar[8] * mp1);
                    //     Compute potential growth for the leaves on fruiting branches.
                    if (site.leaf.age > 70)
                        rate = 0;
                    else
                        rate = smax * c * p * exp(-c * pow(site.leaf.age, p)) * pow(site.leaf.age, (p - 1));
                    if (rate >= 1e-12)
                    {
                        //     Growth rate is modified by water stress. Potential growth is computed as a function of average temperature.
                        site.leaf.potential_growth = rate * wstrlf * TemperatureOnLeafGrowthRate(state.average_temperature);
                        site.petiole.potential_growth = site.leaf.potential_growth * state.leaf_weight_area_ratio * vpotlf[13];
                        PotGroAllLeaves += site.leaf.potential_growth * state.leaf_weight_area_ratio;
                        PotGroAllPetioles += site.petiole.potential_growth;
                    } //if rate
                }     // if LeafAreaNodes
            }         //loop nnid
        }             //loop nbrch
    }                 //loop NumVegBranches
}
