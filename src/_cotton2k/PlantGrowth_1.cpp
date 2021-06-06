//  PlantGrowth_1.cpp
//
//   functions in this file:
// PhysiologicalAge()
// LeafWaterPotential()
// LeafResistance()
// PlantGrowth()
//
#include <cmath>
#include "global.h"
#include "exceptions.h"
#include "GeneralFunctions.h"
#include "RootGrowth.h"
#include "Simulation.hpp"

extern "C"
{
    double LeafResistance(double);
    double PotentialStemGrowth(double, int, Stage, double, double, double, double, double, double, double, double);
    double AddPlantHeight(double, double, uint32_t, Stage, double, double, double, double, double, double, double, double, double, double, double, double, double, double);
}

// PlantGrowth_2
void PotentialLeafGrowth(State &, double, double[61]);

void PotentialFruitGrowth(State &, double[61]);

// PlantGrowth_3
void DryMatterBalance(State &, double &, double &, double &, double &, double);

void ActualFruitGrowth(State &);

void ActualLeafGrowth(State &);

//////////////////////////
void LeafWaterPotential(State &state, double row_space)
//     This function simulates the leaf water potential of cotton plants. It has been
//  adapted from the model of Moshe Meron (The relation of cotton leaf water potential to
//  soil water content in the irrigated management range. PhD dissertation, UC Davis, 1984).
//     It is called from Stress(). It calls wcond() and LeafResistance().
//
//     The following global variables are referenced here:
//       AverageSoilPsi, vanGenuchtenBeta, dl, Kday, LeafAge, NumFruitBranches,
//       NumNodes, NumVegBranches,
//       pi, PlantHeight, PoreSpace, ReferenceETP, RootColNumLeft, RootColNumRight,
//       RootWtCapblUptake, SaturatedHydCond, SoilPsi, thad, thts, VolWaterContent, wk.
//     The following global variables are set here:
//       LwpMin, LwpMax.
//
{
    //     Constant parameters used:
    const double cmg = 3200;     // length in cm per g dry weight of roots, based on an average
                                 // root diameter of 0.06 cm, and a specific weight of 0.11 g dw per cubic cm.
    const double psild0 = -1.32; // maximum values of LwpMin
    const double psiln0 = -0.40; // maximum values of LwpMax.
    const double rtdiam = 0.06;  // average root diameter in cm.
    const double vpsil[] = {0.48, -5.0, 27000., 4000., 9200., 920., 0.000012,
                            -0.15, -1.70, -3.5, 0.1e-9, 0.025, 0.80};
    //     Leaf water potential is not computed during 10 days after
    //  emergence. Constant values are assumed for this period.
    if (state.kday <= 10)
    {
        LwpMax = psiln0;
        LwpMin = psild0;
        return;
    }
    //     Compute shoot resistance (rshoot) as a function of plant height.
    double rshoot; // shoot resistance, Mpa hours per cm.
    rshoot = vpsil[0] * state.plant_height / 100;
    //     Assign zero to summation variables
    double psinum = 0;  // sum of RootWtCapblUptake for all soil cells with roots.
    double rootvol = 0; // sum of volume of all soil cells with roots.
    double rrlsum = 0;  // weighted sum of reciprocals of rrl.
    double rroot = 0;   // root resistance, Mpa hours per cm.
    double sumlv = 0;   // weighted sum of root length, cm, for all soil cells with roots.
    double vh2sum = 0;  // weighted sum of soil water content, for all soil cells with roots.
                        //     Loop over all soil cells with roots. Check if RootWtCapblUptake is
                        // greater than vpsil[10].
                        //     All average values computed for the root zone, are weighted by RootWtCapblUptake
                        // (root weight capable of uptake), but the weight assigned will not be greater than vpsil[11].
    double rrl;         // root resistance per g of active roots.
    for (int l = 0; l < state.soil.number_of_layers_with_root; l++)
        for (int k = state.soil.layers[l].number_of_left_columns_with_root; k < state.soil.layers[l].number_of_right_columns_with_root; k++)
        {
            if (state.soil.cells[l][k].root.weight_capable_uptake >= vpsil[10])
            {
                psinum += std::min(state.soil.cells[l][k].root.weight_capable_uptake, vpsil[11]);
                sumlv += std::min(state.soil.cells[l][k].root.weight_capable_uptake, vpsil[11]) * cmg;
                rootvol += dl(l) * wk(k, row_space);
                if (SoilPsi[l][k] <= vpsil[1])
                    rrl = vpsil[2] / cmg;
                else
                    rrl = (vpsil[3] - SoilPsi[l][k] * (vpsil[4] + vpsil[5] * SoilPsi[l][k])) / cmg;
                rrlsum += std::min(state.soil.cells[l][k].root.weight_capable_uptake, vpsil[11]) / rrl;
                vh2sum += VolWaterContent[l][k] * std::min(state.soil.cells[l][k].root.weight_capable_uptake, vpsil[11]);
            }
        }
    //     Compute average root resistance (rroot) and average soil water content (vh2).
    double dumyrs; // intermediate variable for computing cond.
    double vh2;    // average of soil water content, for all soil soil cells with roots.
    if (psinum > 0 && sumlv > 0)
    {
        rroot = psinum / rrlsum;
        vh2 = vh2sum / psinum;
        dumyrs = sqrt(1 / (pi * sumlv / rootvol)) / rtdiam;
        if (dumyrs < 1.001)
            dumyrs = 1.001;
    }
    else
    {
        rroot = 0;
        vh2 = thad[0];
        dumyrs = 1.001;
    }
    //     Compute hydraulic conductivity (cond), and soil resistance near
    //  the root surface  (rsoil).
    double cond; // soil hydraulic conductivity near the root surface.
    cond = wcond(vh2, thad[0], thts[0], vanGenuchtenBeta[0], SaturatedHydCond[0], PoreSpace[0]) / 24;
    cond = cond * 2 * sumlv / rootvol / log(dumyrs);
    if (cond < vpsil[6])
        cond = vpsil[6];
    double rsoil = 0.0001 / (2 * pi * cond); // soil resistance, Mpa hours per cm.
                                             //     Compute leaf resistance (LeafResistance) as the average of the resistances
                                             // of all existing leaves. The resistance of an individual leaf is a
                                             // function of its age. Function LeafResistance is called to compute it. This
                                             // is executed for all the leaves of the plant.
    int numl = 0;                            // number of leaves.
    double sumrl = 0;                        // sum of leaf resistances for all the plant.
    for (int j = 0; j < state.number_of_pre_fruiting_nodes; j++) // loop prefruiting nodes
    {
        numl++;
        sumrl += LeafResistance(state.age_of_pre_fruiting_nodes[j]);
    }
    //
    for (int k = 0; k < state.number_of_vegetative_branches; k++) // loop for all other nodes
        for (int l = 0; l < state.vegetative_branches[k].number_of_fruiting_branches; l++)
            for (int m = 0; m < state.vegetative_branches[k].fruiting_branches[l].number_of_fruiting_nodes; m++)
            {
                numl++;
                sumrl += LeafResistance(state.vegetative_branches[k].fruiting_branches[l].nodes[m].leaf.age);
            }
    double rleaf = sumrl / numl; // leaf resistance, Mpa hours per cm.

    double rtotal = rsoil + rroot + rshoot + rleaf; // The total resistance to transpiration, MPa hours per cm, (rtotal) is computed.
    //     Compute maximum (early morning) leaf water potential, LwpMax,
    //  from soil water potential (AverageSoilPsi, converted from bars to MPa).
    //     Check for minimum and maximum values.
    LwpMax = vpsil[7] + 0.1 * AverageSoilPsi;
    if (LwpMax < vpsil[8])
        LwpMax = vpsil[8];
    if (LwpMax > psiln0)
        LwpMax = psiln0;
    //     Compute minimum (at time of maximum transpiration rate) leaf water potential, LwpMin, from
    //  maximum transpiration rate (etmax) and total resistance to transpiration (rtotal).
    double etmax = 0;                  // the maximum hourly rate of evapotranspiration for this day.
    for (int ihr = 0; ihr < 24; ihr++) //  hourly loop
    {
        if (state.hours[ihr].ref_et > etmax)
            etmax = state.hours[ihr].ref_et;
    }
    LwpMin = LwpMax - 0.1 * std::max(etmax, vpsil[12]) * rtotal;
    //     Check for minimum and maximum values.
    if (LwpMin < vpsil[9])
        LwpMin = vpsil[9];
    if (LwpMin > psild0)
        LwpMin = psild0;
}

////////////////////////////////////////////////////////////////////////////
void PlantGrowth(State &state, double density_factor, double per_plant_area, double row_space, double cultivar_parameters[61], int NumRootAgeGroups, unsigned int day_emerge, unsigned int day_topping, unsigned int first_square, unsigned int plant_row_column)
//     This function simulates the potential and actual growth of cotton plants.
//  It is called from SimulateThisDay(), and it calls the following functions:
//    ActualFruitGrowth(), ActualLeafGrowth(), ActualRootGrowth(), AddPlantHeight(),
//    DryMatterBalance(), PotentialFruitGrowth(), PotentialLeafGrowth(),
//    PotentialRootGrowth(), PotentialStemGrowth().
//
//     The following global variables are referenced here:
//        ActualStemGrowth, DayInc, FirstSquare, FruitingCode, Kday,
//        RowSpace, WaterStressStem.
//
//     The following global variables are set here:
//        LeafAreaIndex, PlantHeight, PotGroAllRoots, PotGroStem, StemWeight,
//        TotalPetioleWeight.
{
    //     Call PotentialLeafGrowth() to compute potential growth rate of leaves.
    PotentialLeafGrowth(state, density_factor, cultivar_parameters);
    //     If it is after first square, call PotentialFruitGrowth() to compute potential
    //  growth rate of squares and bolls.
    if (state.vegetative_branches[0].fruiting_branches[0].nodes[0].stage != Stage::NotYetFormed)
        PotentialFruitGrowth(state, cultivar_parameters);
    //     Active stem tissue (stemnew) is the difference between state.stem_weight
    //  and the value of StemWeight(kkday).
    int voldstm = 32;           // constant parameter (days for stem tissue to become "old")
    int kkday = state.kday - voldstm; // age of young stem tissue
    if (kkday < 1)
        kkday = 1;
    double stemnew = state.stem_weight - StemWeight[kkday]; // dry weight of active stem tissue.
                                                          //     Call PotentialStemGrowth() to compute PotGroStem, potential growth rate of stems.
                                                          //  The effect of temperature is introduced, by multiplying potential growth rate by DayInc.
                                                          //  Stem growth is also affected by water stress (WaterStressStem). PotGroStem is limited by (maxstmgr * per_plant_area) g per plant per day.
    PotGroStem = PotentialStemGrowth(stemnew, state.kday, state.vegetative_branches[0].fruiting_branches[2].nodes[0].stage, density_factor, cultivar_parameters[12], cultivar_parameters[13], cultivar_parameters[14], cultivar_parameters[15], cultivar_parameters[16], cultivar_parameters[17], cultivar_parameters[18]) * state.day_inc * state.water_stress_stem;
    double maxstmgr = 0.067; // maximum posible potential stem growth, g dm-2 day-1.
    if (PotGroStem > maxstmgr * per_plant_area)
        PotGroStem = maxstmgr * per_plant_area;
    //	   Call PotentialRootGrowth() to compute potential growth rate of roots.
    double sumpdr; // total potential growth rate of roots in g per slab. this is
    // computed in PotentialRootGrowth() and used in ActualRootGrowth().
    sumpdr = PotentialRootGrowth(state.soil.cells, NumRootAgeGroups, state.soil.number_of_layers_with_root, per_plant_area);
    //     Total potential growth rate of roots is converted from g per
    //  slab (sumpdr) to g per plant (PotGroAllRoots).
    PotGroAllRoots = sumpdr * 100 * per_plant_area / row_space;
    //     Limit PotGroAllRoots to (maxrtgr* per_plant_area) g per plant per day.
    double maxrtgr = 0.045; // maximum possible potential root growth, g dm-2 day-1.
    if (PotGroAllRoots > maxrtgr * per_plant_area)
        PotGroAllRoots = maxrtgr * per_plant_area;
    //     Call DryMatterBalance() to compute carbon balance, allocation of carbon to
    //  plant parts, and carbon stress. DryMatterBalance() also computes and returns the
    //  values of the following arguments:
    //     cdleaf is carbohydrate requirement for leaf growth, g per plant per day.
    //     cdpet is carbohydrate requirement for petiole growth, g per plant per day.
    //     cdroot is carbohydrate requirement for root growth, g per plant per day.
    //     cdstem is carbohydrate requirement for stem growth, g per plant per day.
    double cdstem, cdleaf, cdpet, cdroot;
    DryMatterBalance(state ,cdstem, cdleaf, cdpet, cdroot, per_plant_area);
    //     If it is after first square, call ActualFruitGrowth() to compute actual
    //  growth rate of squares and bolls.
    if (state.vegetative_branches[0].fruiting_branches[0].nodes[0].stage != Stage::NotYetFormed)
        ActualFruitGrowth(state);
    //     Initialize state.leaf_weight. It is assumed that cotyledons fall off
    //  at time of first square. Also initialize state.leaf_area and TotalPetioleWeight.
    if (first_square > 0)
    {
        state.leaf_weight = 0;
        state.leaf_area = 0;
    }
    else
    {
        double cotylwt = 0.20; // weight of cotyledons dry matter.
        state.leaf_weight = cotylwt;
        state.leaf_area = 0.6 * cotylwt;
    }
    TotalPetioleWeight = 0;
    //     Call ActualLeafGrowth to compute actual growth rate of leaves and compute leaf area index.
    ActualLeafGrowth(state);
    state.leaf_area_index = state.leaf_area / per_plant_area;
    //     Add ActualStemGrowth to state.stem_weight, and define StemWeight(Kday) for this day.
    state.stem_weight += ActualStemGrowth;
    StemWeight[state.kday] = state.stem_weight;
    //     Plant density affects growth in height of tall plants.
    double htdenf = 55; // minimum plant height for plant density affecting growth in height.
    double z1;          // intermediate variable to compute denf2.
    z1 = (state.plant_height - htdenf) / htdenf;
    if (z1 < 0)
        z1 = 0;
    if (z1 > 1)
        z1 = 1;
    double denf2; // effect of plant density on plant growth in height.
    denf2 = 1 + z1 * (density_factor - 1);
    //     Call AddPlantHeight to compute PlantHeight.
    int l, l1, l2; // node numbers of top three nodes.
    l = state.vegetative_branches[0].number_of_fruiting_branches - 1;
    l1 = l - 1;
    if (l < 1)
        l1 = 0;
    l2 = l - 2;
    if (l < 2)
        l2 = 0;
    double agetop; // average physiological age of top three nodes.
    agetop = (state.vegetative_branches[0].fruiting_branches[l].nodes[0].age + state.vegetative_branches[0].fruiting_branches[l1].nodes[0].age + state.vegetative_branches[0].fruiting_branches[l2].nodes[0].age) / 3;
    if (day_topping <= 0 || state.daynum < day_topping)
        state.plant_height += AddPlantHeight(denf2, state.day_inc, state.number_of_pre_fruiting_nodes, state.vegetative_branches[0].fruiting_branches[1].nodes[0].stage, state.age_of_pre_fruiting_nodes[state.number_of_pre_fruiting_nodes - 1], state.age_of_pre_fruiting_nodes[state.number_of_pre_fruiting_nodes - 2], agetop, state.water_stress_stem, state.carbon_stress, state.nitrogen_stress_vegetative, cultivar_parameters[19], cultivar_parameters[20], cultivar_parameters[21], cultivar_parameters[22], cultivar_parameters[23], cultivar_parameters[24], cultivar_parameters[25], cultivar_parameters[26]);
    //     Call ActualRootGrowth() to compute actual root growth.
    ComputeActualRootGrowth(state, sumpdr, row_space, per_plant_area, NumRootAgeGroups, day_emerge, plant_row_column);
}
