//  PlantGrowth_1.cpp
//
//   functions in this file:
// PhysiologicalAge()
// LeafWaterPotential()
// LeafResistance()
// PlantGrowth()
//
#include <math.h>

#include "CottonSimulation.h"
#include "GeneralFunctions.h"

//
#ifdef _DEBUG
#define new DEBUG_NEW
#endif
//////////////////////////////////////////////////
double PhysiologicalAge()  // computes physiological age
//     This function returns the daily 'physiological age' increment,
//  based on hourly temperatures. It is called each day by SimulateThisDay().
//     The following global variable is used here:
//        AirTemp[] = array of hourly temperatures.
{
    //     The following constant Parameters are used in this function:
    const double p1 = 12;  // threshold temperature, C
    const double p2 =
        14;  // temperature, C, above p1, for one physiological day.
    const double p3 = 1.5;  // maximum value of a physiological day.
    //     The threshold value is assumed to be 12 C (p1). One physiological day
    //     is
    //  equivalent to a day with an average temperature of 26 C, and therefore
    //  the heat units are divided by 14 (p2).
    //     A linear relationship is assumed between temperature and heat unit
    //  accumulation in the range of 12 C (p1) to 33 C (p2*p3+p1). the effect of
    //  temperatures higher than 33 C is assumed to be equivalent to that of 33
    //  C.
    double dayfd =
        0;  // the daily contribution to physiological age (return value).
    for (int ihr = 0; ihr < 24; ihr++) {
        double tfd = (AirTemp[ihr] - p1) /
                     p2;  // the hourly contribution to physiological age.
        if (tfd < 0) tfd = 0;
        if (tfd > p3) tfd = p3;
        dayfd += tfd;
    }
    return dayfd / 24;
}
//////////////////////////
void LeafWaterPotential()
//     This function simulates the leaf water potential of cotton plants. It has
//     been
//  adapted from the model of Moshe Meron (The relation of cotton leaf water
//  potential to soil water content in the irrigated management range. PhD
//  dissertation, UC Davis, 1984).
//     It is called from Stress(). It calls wcond() and LeafResistance().
//
//     The following global variables are referenced here:
//       AgeOfPreFruNode, AverageSoilPsi, beta, dl, Kday, LeafAge,
//       NumFruitBranches, NumLayersWithRoots, NumNodes, NumPreFruNodes,
//       NumVegBranches, pi, PlantHeight, PoreSpace, ReferenceETP,
//       RootColNumLeft, RootColNumRight, RootWtCapblUptake, SaturatedHydCond,
//       SoilPsi, thad, thts, VolWaterContent, wk.
//     The following global variables are set here:
//       LwpMin, LwpMax.
//
{
    //     Constant parameters used:
    const double cmg =
        3200;  // length in cm per g dry weight of roots, based on an average
    //               root diameter of 0.06 cm, and a specific weight of 0.11 g
    //               dw per cubic cm.
    const double psild0 = -1.32;  // maximum values of LwpMin
    const double psiln0 = -0.40;  // maximum values of LwpMax.
    const double rtdiam = 0.06;   // average root diameter in cm.
    const double vpsil[] = {0.48,  -5.0,  27000., 4000.,  9200., 920., 0.000012,
                            -0.15, -1.70, -3.5,   0.1e-9, 0.025, 0.80};
    //     Leaf water potential is not computed during 10 days after
    //  emergence. Constant values are assumed for this period.
    if (Kday <= 10) {
        LwpMax = psiln0;
        LwpMin = psild0;
        return;
    }
    //     Compute shoot resistance (rshoot) as a function of plant height.
    double rshoot;  // shoot resistance, Mpa hours per cm.
    rshoot = vpsil[0] * PlantHeight / 100;
    //     Assign zero to summation variables
    double psinum =
        0;  // sum of RootWtCapblUptake for all soil cells with roots.
    double rootvol = 0;  // sum of volume of all soil cells with roots.
    double rrlsum = 0;   // weighted sum of reciprocals of rrl.
    double rroot = 0;    // root resistance, Mpa hours per cm.
    double sumlv =
        0;  // weighted sum of root length, cm, for all soil cells with roots.
    double vh2sum = 0;  // weighted sum of soil water content, for all soil
                        // cells with roots.
    //     Loop over all soil cells with roots. Check if RootWtCapblUptake is
    //  greater than vpsil[10].
    //     All average values computed for the root zone, are weighted by
    //     RootWtCapblUptake
    //  (root weight capable of uptake), but the weight assigned will not be
    //  greater than vpsil[11].
    double rrl;  // root resistance per g of active roots.
    for (int l = 0; l < NumLayersWithRoots; l++)
        for (int k = RootColNumLeft[l]; k < RootColNumRight[l]; k++) {
            if (RootWtCapblUptake[l][k] >= vpsil[10]) {
                psinum += fmin(RootWtCapblUptake[l][k], vpsil[11]);
                sumlv += fmin(RootWtCapblUptake[l][k], vpsil[11]) * cmg;
                rootvol += dl[l] * wk[k];
                if (SoilPsi[l][k] <= vpsil[1])
                    rrl = vpsil[2] / cmg;
                else
                    rrl =
                        (vpsil[3] - SoilPsi[l][k] *
                                        (vpsil[4] + vpsil[5] * SoilPsi[l][k])) /
                        cmg;
                rrlsum += fmin(RootWtCapblUptake[l][k], vpsil[11]) / rrl;
                vh2sum += VolWaterContent[l][k] *
                          fmin(RootWtCapblUptake[l][k], vpsil[11]);
            }
        }
    //     Compute average root resistance (rroot) and average soil water
    //     content (vh2).
    double dumyrs;  // intermediate variable for computing cond.
    double vh2;  // average of soil water content, for all soil soil cells with
                 // roots.
    if (psinum > 0 && sumlv > 0) {
        rroot = psinum / rrlsum;
        vh2 = vh2sum / psinum;
        dumyrs = sqrt(1 / (pi * sumlv / rootvol)) / rtdiam;
        if (dumyrs < 1.001) dumyrs = 1.001;
    } else {
        rroot = 0;
        vh2 = thad[0];
        dumyrs = 1.001;
    }
    //     Compute hydraulic conductivity (cond), and soil resistance near
    //  the root surface  (rsoil).
    double cond;  // soil hydraulic conductivity near the root surface.
    cond = wcond(vh2, thad[0], thts[0], beta[0], SaturatedHydCond[0],
                 PoreSpace[0]) /
           24;
    cond = cond * 2 * sumlv / rootvol / log(dumyrs);
    if (cond < vpsil[6]) cond = vpsil[6];
    double rsoil =
        0.0001 / (2 * pi * cond);  // soil resistance, Mpa hours per cm.
    //     Compute leaf resistance (LeafResistance) as the average of the
    //     resistances
    //  of all existing leaves. The resistance of an individual leaf is a
    //  function of its age. Function LeafResistance is called to compute it.
    //  This is executed for all the leaves of the plant.
    int numl = 0;      // number of leaves.
    double sumrl = 0;  // sum of leaf resistances for all the plant.
    for (int j = 0; j < NumPreFruNodes; j++)  // loop prefruiting nodes
    {
        numl++;
        sumrl += LeafResistance(AgeOfPreFruNode[j]);
    }
    //
    int nbrch;  // number of fruiting branches on a vegetative branch.
    int nnid;   // number of nodes on a fruiting branch.
    for (int k = 0; k < NumVegBranches; k++)  // loop for all other nodes
    {
        nbrch = NumFruitBranches[k];
        for (int l = 0; l < nbrch; l++) {
            nnid = NumNodes[k][l];
            for (int m = 0; m < nnid; m++) {
                numl++;
                sumrl += LeafResistance(LeafAge[k][l][m]);
            }
        }
    }
    double rleaf = sumrl / numl;  // leaf resistance, Mpa hours per cm.
    //     The total resistance to transpiration, MPa hours per cm, (rtotal) is
    //     computed.
    double rtotal = rsoil + rroot + rshoot + rleaf;
    //     Compute maximum (early morning) leaf water potential, LwpMax,
    //  from soil water potential (AverageSoilPsi, converted from bars to MPa).
    //     Check for minimum and maximum values.
    LwpMax = vpsil[7] + 0.1 * AverageSoilPsi;
    if (LwpMax < vpsil[8]) LwpMax = vpsil[8];
    if (LwpMax > psiln0) LwpMax = psiln0;
    //     Compute minimum (at time of maximum transpiration rate) leaf water
    //     potential, LwpMin, from
    //  maximum transpiration rate (etmax) and total resistance to transpiration
    //  (rtotal).
    double etmax =
        0;  // the maximum hourly rate of evapotranspiration for this day.
    for (int ihr = 0; ihr < 24; ihr++)  //  hourly loop
    {
        if (ReferenceETP[ihr] > etmax) etmax = ReferenceETP[ihr];
    }
    LwpMin = LwpMax - 0.1 * fmax(etmax, vpsil[12]) * rtotal;
    //     Check for minimum and maximum values.
    if (LwpMin < vpsil[9]) LwpMin = vpsil[9];
    if (LwpMin > psild0) LwpMin = psild0;
}
//////////////////////////////////
double LeafResistance(double agel)
//     This function computes and returns the resistance of leaves of cotton
// plants to transpiration. It is assumed to be a function of leaf age.
// It is called from LeafWaterPotential().
//     The input argument (agel) is leaf age in physiological days.
{
    //     The following constant parameters are used:
    const double afac = 160;   // factor used for computing leaf resistance.
    const double agehi = 94;   // higher limit for leaf age.
    const double agelo = 48;   // lower limit for leaf age.
    const double rlmin = 0.5;  // minimum leaf resistance.
                               //
    double leafResistance;
    if (agel <= agelo)
        leafResistance = rlmin;
    else if (agel >= agehi)
        leafResistance = rlmin + (agehi - agelo) * (agehi - agelo) / afac;
    else {
        double ax = 2 * agehi - agelo;  // intermediate variable
        leafResistance = rlmin + (agel - agelo) * (ax - agel) / afac;
    }
    return leafResistance;
}
////////////////////////////////////////////////////////////////////////////
void PlantGrowth()
//     This function simulates the potential and actual growth of cotton plants.
//  It is called from SimulateThisDay(), and it calls the following functions:
//    ActualFruitGrowth(), ActualLeafGrowth(), ActualRootGrowth(),
//    AddPlantHeight(), DryMatterBalance(), PotentialFruitGrowth(),
//    PotentialLeafGrowth(), PotentialRootGrowth(), PotentialStemGrowth().
//
//     The following global variables are referenced here:
//        ActualStemGrowth, DayInc, FirstSquare, FruitingCode, Kday, pixdz,
//        PerPlantArea, RowSpace, WaterStressStem.
//
//     The following global variables are set here:
//        LeafAreaIndex, PlantHeight, PotGroAllRoots, PotGroStem, StemWeight,
//        TotalLeafArea, TotalLeafWeight, TotalPetioleWeight, TotalStemWeight.
{
    //     Call PotentialLeafGrowth() to compute potential growth rate of
    //     leaves.
    PotentialLeafGrowth();
    //     If it is after first square, call PotentialFruitGrowth() to compute
    //     potential
    //  growth rate of squares and bolls.
    if (FruitingCode[0][0][0] > 0) PotentialFruitGrowth();
    //     Active stem tissue (stemnew) is the difference between
    //     TotalStemWeight
    //  and the value of StemWeight(kkday).
    int voldstm =
        32;  // constant parameter (days for stem tissue to become "old")
    int kkday = Kday - voldstm;  // age of young stem tissue
    if (kkday < 1) kkday = 1;
    double stemnew = TotalStemWeight -
                     StemWeight[kkday];  // dry weight of active stem tissue.
    //     Call PotentialStemGrowth() to compute PotGroStem, potential growth
    //     rate of stems.
    //  The effect of temperature is introduced, by multiplying potential growth
    //  rate by DayInc. Stem growth is also affected by water stress
    //  (WaterStressStem) and possible PIX application (pixdz).   PotGroStem is
    //  limited by (maxstmgr * PerPlantArea) g per plant per day.
    PotGroStem =
        PotentialStemGrowth(stemnew) * DayInc * WaterStressStem * pixdz;
    double maxstmgr =
        0.067;  // maximum posible potential stem growth, g dm-2 day-1.
    if (PotGroStem > maxstmgr * PerPlantArea)
        PotGroStem = maxstmgr * PerPlantArea;
    //	   Call PotentialRootGrowth() to compute potential growth rate of roots.
    double sumpdr;  // total potential growth rate of roots in g per slab. this
                    // is computed in PotentialRootGrowth() and used in
                    // ActualRootGrowth().
    sumpdr = PotentialRootGrowth();
    //     Total potential growth rate of roots is converted from g per
    //  slab (sumpdr) to g per plant (PotGroAllRoots).
    PotGroAllRoots = sumpdr * 100 * PerPlantArea / RowSpace;
    //     Limit PotGroAllRoots to (maxrtgr*PerPlantArea) g per plant per day.
    double maxrtgr =
        0.045;  // maximum possible potential root growth, g dm-2 day-1.
    if (PotGroAllRoots > maxrtgr * PerPlantArea)
        PotGroAllRoots = maxrtgr * PerPlantArea;
    //     Call DryMatterBalance() to compute carbon balance, allocation of
    //     carbon to
    //  plant parts, and carbon stress. DryMatterBalance() also computes and
    //  returns the values of the following arguments:
    //     cdleaf is carbohydrate requirement for leaf growth, g per plant per
    //     day. cdpet is carbohydrate requirement for petiole growth, g per
    //     plant per day. cdroot is carbohydrate requirement for root growth, g
    //     per plant per day. cdstem is carbohydrate requirement for stem
    //     growth, g per plant per day.
    double cdstem, cdleaf, cdpet, cdroot;
    DryMatterBalance(cdstem, cdleaf, cdpet, cdroot);
    //     If it is after first square, call ActualFruitGrowth() to compute
    //     actual
    //  growth rate of squares and bolls.
    if (FruitingCode[0][0][0] > 0) ActualFruitGrowth();
    //     Initialize TotalLeafWeight. It is assumed that cotyledons fall off
    //  at time of first square. Also initialize TotalLeafArea and
    //  TotalPetioleWeight.
    TotalPetioleWeight = 0;
    //     Call ActualLeafGrowth to compute actual growth rate of leaves and
    //     compute leaf area index.
    ActualLeafGrowth();
    LeafAreaIndex = TotalLeafArea() / PerPlantArea;
    //     Add ActualStemGrowth to TotalStemWeight, and define StemWeight(Kday)
    //     for this day.
    TotalStemWeight += ActualStemGrowth;
    StemWeight[Kday] = TotalStemWeight;
    //     Plant density affects growth in height of tall plants.
    double htdenf = 55;  // minimum plant height for plant density affecting
                         // growth in height.
    double z1;           // intermediate variable to compute denf2.
    z1 = (PlantHeight - htdenf) / htdenf;
    if (z1 < 0) z1 = 0;
    if (z1 > 1) z1 = 1;
    double denf2;  // effect of plant density on plant growth in height.
    denf2 = 1 + z1 * (DensityFactor - 1);
    //     Call AddPlantHeight to compute PlantHeight.
    PlantHeight += AddPlantHeight(denf2);
    //     Call ActualRootGrowth() to compute actual root growth.
    ComputeActualRootGrowth(sumpdr);
}
