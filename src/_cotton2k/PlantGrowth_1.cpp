//  PlantGrowth_1.cpp
//
//   functions in this file:
// PhysiologicalAge()
// LeafWaterPotential()
// LeafResistance()
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
}

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
//       RootWtCapblUptake, SaturatedHydCond, SoilPsi, thad, thts.
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
                vh2sum += state.soil.cells[l][k].water_content * std::min(state.soil.cells[l][k].root.weight_capable_uptake, vpsil[11]);
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
