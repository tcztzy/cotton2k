//  File  PlantNitrogen.cpp
//
//   functions in this file:
// PlantNitrogen()
// NitrogenRequirement()
// NitrogenSupply()
// PetioleNitrateN()
// NitrogenAllocation()
// ExtraNitrogenAllocation()
// PlantNitrogenContent()
// GetNitrogenStress()
// NitrogenUptakeRequirement()
// PlantNitrogenBal()
//
#include "CottonSimulation.h"
#include "GeneralFunctions.h"
//
#ifdef _DEBUG
#define new DEBUG_NEW
#endif
//  File scope variables
double burres;  // reserve N in burrs, in g per plant.
double leafrs;  // reserve N in leaves, in g per plant.
double petrs;   // reserve N in petioles, in g per plant.
double stemrs;  // reserve N in stems, in g per plant.
double rootrs;  // reserve N in roots, in g per plant.
double rqnbur;  // nitrogen requirement for burr growth.
double rqnlef;  // nitrogen requirement for leaf growth.
double rqnpet;  // nitrogen requirement for petiole growth.
double rqnrut;  // nitrogen requirement for root growth.
double rqnsed;  // nitrogen requirement for seed growth.
double rqnsqr;  // nitrogen requirement for square growth.
double rqnstm;  // nitrogen requirement for stem growth.
double reqv;    // nitrogen requirement for vegetative shoot growth.
double reqf;    // nitrogen requirement for fruit growth.
double reqtot;  // total nitrogen requirement for plant growth.
double npool;   // total nitrogen available for growth.
double uptn;    // nitrogen uptake from the soil, g per plant.
double xtran;   // amount of nitrogen not used for growth of plant parts.
double addnf;   // daily added nitrogen to fruit, g per plant.
double addnr;   // daily added nitrogen to root, g per plant.
double addnv;   // daily added nitrogen to vegetative shoot, g per plant.
//////////////////////////////////////////////////
void PlantNitrogen()
//     This function simulates the nitrogen accumulation and distribution in
//     cotton plants,
//  and computes nitrogen stresses. It is called from SimulateThisDay().
//
//     The maximum and minimum N concentrations for the various organs are
//     modified from
//  those reported by: Jones et. al. (1974): Development of a nitrogen balance
//  for cotton growth models: a first approximation. Crop Sci. 14:541-546.
//
//     The following parameters are used in all plant N routines:
//            Growth requirement   Uptake requirement  Minimum content
//  leaves      lefcn0    = .064    vnreqlef  = .042    vlfnmin   = .018
//  petioles    petcn0    = .032    vnreqpet  = .036    vpetnmin  = .005
//  stems       stmcn0    = .036    vnreqstm  = .012    vstmnmin  = .006
//  roots       rootcn0   = .018    vnreqrt   = .010    vrtnmin   = .010
//  burrs       burcn0    = .026    vnreqbur  = .012    vburnmin  = .006
//  seeds       seedcn0   = .036    seedcn1   = .045
//  squares     sqrcn0    = .024    vnreqsqr  = .024
//
//     PlantNitrogen() calls the following functions:
//       NitrogenSupply(), NitrogenRequirement(), NitrogenAllocation(),
//       ExtraNitrogenAllocation(), PlantNitrogenContent(), GetNitrogenStress(),
//       NitrogenUptakeRequirement().
//     The following global variables are referenced here:
//       BurrNConc, BurrNitrogen, Kday, LeafNConc, LeafNitrogen,
//       PetioleNConc, PetioleNConc, PetioleNO3NConc, PetioleNitrogen,
//       RootNConc, RootNitrogen, SeedNConc, SeedNitrogen, SquareNitrogen,
//       StemNConc, StemNitrogen.
//     The following global and file scope variables are set here:
//       addnf, addnr, addnv, burres, leafrs, npool, petrs, reqf, reqtot, reqv,
//       rootrs, rqnbur, rqnlef, rqnpet, rqnrut, rqnsed, rqnsqr, rqnstm, stemrs,
//       uptn, xtran.
{
    //     Assign zero to some variables.
    leafrs = 0;  // reserve N in leaves, in g per plant.
    petrs = 0;   // reserve N in petioles, in g per plant.
    stemrs = 0;  // reserve N in stems, in g per plant.
    rootrs = 0;  // reserve N in roots, in g per plant.
    reqf = 0;    // nitrogen requirement for fruit growth.
    reqtot = 0;  // total nitrogen requirement for plant growth.
    reqv = 0;    // nitrogen requirement for vegetative shoot growth.
    rqnbur = 0;  // nitrogen requirement for burr growth.
    rqnlef = 0;  // nitrogen requirement for leaf growth.
    rqnpet = 0;  // nitrogen requirement for petiole growth.
    rqnrut = 0;  // nitrogen requirement for root growth.
    rqnsed = 0;  // nitrogen requirement for seed growth.
    rqnsqr = 0;  // nitrogen requirement for square growth.
    rqnstm = 0;  // nitrogen requirement for stem growth.
    npool = 0;   // total nitrogen available for growth.
    uptn = 0;    // nitrogen uptake from the soil, g per plant.
    xtran = 0;   // amount of nitrogen not used for growth of plant parts.
    addnf = 0;   // daily added nitrogen to fruit, g per plant.
    addnr = 0;   // daily added nitrogen to root, g per plant.
    addnv = 0;   // daily added nitrogen to vegetative shoot, g per plant.
                 //   The following subroutines are now called:
    NitrogenRequirement();  //  computes the N requirements for growth.
    NitrogenSupply();  //  computes the supply of N from uptake and reserves.
    NitrogenAllocation();  //  computes the allocation of N in the plant.
    if (xtran > 0)
        ExtraNitrogenAllocation();  // computes the further allocation of N in
                                    // the plant
    PlantNitrogenContent();  // computes the concentrations of N in plant dry
                             // matter.
    GetNitrogenStress();     //  computes nitrogen stress factors.
    NitrogenUptakeRequirement();  // computes N requirements for uptake
}
//////////////////////////
void NitrogenRequirement()
//     This function computes the N requirements for growth. It is called from
//     PlantNitrogen{}.
//
//     The following global variables are referenced here:
//       ActualBollGrowth, ActualBurrGrowth, ActualSquareGrowth,
//       ActualStemGrowth, CarbonAllocatedForRootGrowth, CottonWeightGreenBolls,
//       DayEmerge, Daynum, ExtraCarbon, Kday,
//       SeedNitrogen, TotalActualLeafGrowth, TotalActualPetioleGrowth.
//     The following global and file scope variables are set in this function:
//       PetioleNConc, PetioleNO3NConc, reqf, reqtot, reqv, rqnbur,
//       rqnlef, rqnpet, rqnrut, rqnsed, rqnsqr, rqnstm.
{
    //     The following constant parameters are used:
    const double burcn0 = .026;   //  maximum N content for burrs
    const double lefcn0 = .064;   //  maximum N content for leaves
    const double petcn0 = .032;   //  maximum N content for petioles
    const double rootcn0 = .018;  //  maximum N content for roots
    const double seedcn0 = .036;  //  maximum N content for seeds
    const double seedcn1 =
        .045;  //  additional requirement of N for existing seed tissue.
    const double seedratio =
        .64;  //  ratio of seeds to seedcotton in green bolls.
    const double sqrcn0 = .024;  //  maximum N content for squares
    const double stmcn0 = .036;  //  maximum N content for stems
    //     On emergence, assign initial values to petiole N concentrations.
    if (Daynum <= DayEmerge) {
        PetioleNConc = petcn0;
        PetioleNO3NConc = petcn0;
    }
    //     Compute the nitrogen requirements for growth, by multiplying
    //  the daily added dry weight by the maximum N content of each organ.
    //  Nitrogen requirements are based on actual growth rates.
    //     These N requirements will be used to compute the allocation
    //  of N to plant parts and the nitrogen stress factors.
    //     All nitrogen requirement variables are in g N per plant.
    rqnlef = lefcn0 * TotalActualLeafGrowth;     //      for leaf blade
    rqnpet = petcn0 * TotalActualPetioleGrowth;  //      for petiole
    rqnstm = stmcn0 * ActualStemGrowth;          //      for stem
    //     Add ExtraCarbon to CarbonAllocatedForRootGrowth to compute the total
    //     supply of
    //  carbohydrates for root growth.
    rqnrut =
        rootcn0 * (CarbonAllocatedForRootGrowth + ExtraCarbon);  // for root
    rqnsqr = ActualSquareGrowth * sqrcn0;  //      for squares
    double rqnsed1, rqnsed2;               // components of seed N requirements.
    rqnsed1 = ActualBollGrowth * seedratio * seedcn0;  //   for seed growth
    //     The N required for replenishing the N content of existing seed
    //  tissue (rqnsed2) is added to seed growth requirement.
    if (CottonWeightGreenBolls > ActualBollGrowth) {
        double rseedn;  // existing ratio of N to dry matter in the seeds.
        rseedn = SeedNitrogen /
                 ((CottonWeightGreenBolls - ActualBollGrowth) * seedratio);
        rqnsed2 = (CottonWeightGreenBolls - ActualBollGrowth) * seedratio *
                  (seedcn1 - rseedn);
        if (rqnsed2 < 0) rqnsed2 = 0;
    } else
        rqnsed2 = 0;
    //
    rqnsed = rqnsed1 + rqnsed2;          //   total requirement for seeds
    rqnbur = ActualBurrGrowth * burcn0;  //   for burrs
    reqf = rqnsqr + rqnsed + rqnbur;     //    total for fruit
    reqv = rqnlef + rqnpet + rqnstm;     //    total for shoot
    reqtot = rqnrut + reqv + reqf;       //     total N requirement
}
//////////////////////////
void NitrogenSupply()
//     This function computes the supply of N by uptake from the soil reserves,
//  It is called from PlantNitrogen(). It calls function PetioleNitrateN().
//
//     The following global variables are referenced here:
//       BurrWeightGreenBolls, Kday, reqtot, SupplyNH4N, SupplyNO3N,
//       TotalLeafWeight, TotalPetioleWeight, TotalRootWeight, TotalStemWeight.
//     The following global and file scope variables are set here:
//       burres, BurrNitrogen, LeafNitrogen, leafrs, npool, PetioleNitrogen,
//       petrs, RootNitrogen, rootrs, StemNitrogen, stemrs, uptn, xtran.
{
    //     The following constant parameters are used:
    const double MobilizNFractionBurrs =
        0.08;  //  fraction of N mobilizable for burrs
    const double MobilizNFractionLeaves =
        0.09;  //  fraction of N mobilizable for leaves and petioles
    const double MobilizNFractionStemRoot =
        0.40;  //  fraction of N mobilizable for stems and roots
    const double vburnmin = .006;  //  minimum N contents of burrs
    const double vlfnmin = .018;   //  minimum N contents of leaves
    const double vpetnmin =
        .005;  //  minimum N contents of petioles, non-nitrate fraction.
    const double vpno3min =
        .002;  //  minimum N contents of petioles, nitrate fraction.
    const double vrtnmin = .010;   //  minimum N contents of roots
    const double vstmnmin = .006;  //  minimum N contents of stems
    //     uptn is the total supply of nitrogen to the plant by uptake of
    //     nitrate and ammonium.
    uptn = SupplyNO3N + SupplyNH4N;
    double resn;  // total reserve N, in g per plant.
    //     If total N requirement is less than the supply, define npool as the
    //     supply and
    //  assign zero to the N reserves in all organs.
    if (reqtot <= uptn) {
        npool = uptn;
        leafrs = 0;
        petrs = 0;
        stemrs = 0;
        rootrs = 0;
        burres = 0;
        xtran = 0;
        resn = 0;
    } else {
        //     If total N requirement exceeds the supply, compute the nitrogen
        //  reserves in the plant. The reserve N in an organ is defined as a
        //  fraction of the nitrogen content exceeding a minimum N content in
        //  it.
        //     The N reserves in leaves, petioles, stems, roots and burrs of
        //  green bolls are computed, and their N content updated.
        leafrs =
            (LeafNitrogen - vlfnmin * TotalLeafWeight()) * MobilizNFractionLeaves;
        if (leafrs < 0) leafrs = 0;
        LeafNitrogen -= leafrs;
        //     The petiole N content is subdivided to nitrate and non-nitrate.
        //     The
        //  nitrate ratio in the petiole N is computed by calling function
        //  PetioleNitrateN(). Note that the nitrate fraction is more available
        //  for redistribution.
        double rpetno3;  // ratio of NO3 N to total N in petioles.
        rpetno3 = PetioleNitrateN();
        double petrs1, petrs2;  // components of reserve N in petioles, for
                                // non-NO3 and NO3 origin, respectively.
        petrs1 =
            (PetioleNitrogen * (1 - rpetno3) - vpetnmin * TotalPetioleWeight) *
            MobilizNFractionLeaves;
        if (petrs1 < 0) petrs1 = 0;
        petrs2 = (PetioleNitrogen * rpetno3 - vpno3min * TotalPetioleWeight) *
                 MobilizNFractionLeaves;
        if (petrs2 < 0) petrs2 = 0;
        petrs = petrs1 + petrs2;
        PetioleNitrogen -= petrs;
        //  Stem N reserves.
        stemrs = (StemNitrogen - vstmnmin * TotalStemWeight) *
                 MobilizNFractionStemRoot;
        if (stemrs < 0) stemrs = 0;
        StemNitrogen -= stemrs;
        //  Root N reserves
        rootrs = (RootNitrogen - vrtnmin * TotalRootWeight) *
                 MobilizNFractionStemRoot;
        if (rootrs < 0) rootrs = 0;
        RootNitrogen -= rootrs;
        //  Burr N reserves
        if (BurrWeightGreenBolls > 0) {
            burres = (BurrNitrogen - vburnmin * BurrWeightGreenBolls) *
                     MobilizNFractionBurrs;
            if (burres < 0) burres = 0;
            BurrNitrogen -= burres;
        } else
            burres = 0;
        //     The total reserves, resn, are added to the amount taken up from
        //     the soil, for computing
        //  npool. Note that N of seeds or squares is not available for
        //  redistribution in the plant.
        resn = leafrs + petrs + stemrs + rootrs + burres;
        npool = uptn + resn;
    }
}
//////////////////////////
double PetioleNitrateN()
//     This function computes the ratio of NO3 nitrogen to total N in the
//     petioles.
//  It is called from NitrogenSupply() and PlantNitrogenContent().
//     The following global variables are referenced here:
//       AgeOfPreFruNode, LeafAge, NumFruitBranches, NumNodes,
//       NumPreFruNodes, NumVegBranches.
{
    //     The following constant parameters are used:
    const double p1 =
        0.96;  //  the maximum ratio (of NO3 to total N in petioles).
    const double p2 = 0.015;  //  the rate of decline of this ratio with age.
    const double p3 = 0.02;   //  the minimum ratio
    //     The ratio of NO3 to total N in each individual petiole is computed
    //  as a linear function of leaf age. It is assumed that this ratio is
    //  maximum for young leaves and is declining with leaf age.
    int numl = 0;        // number of petioles computed.
    double spetno3 = 0;  // sum of petno3r.
    double petno3r;      // ratio of NO3 to total N in an individual petiole.
                         //     Loop of prefruiting node leaves.
    for (int j = 0; j < NumPreFruNodes; j++) {
        numl++;
        petno3r = p1 - AgeOfPreFruNode[j] * p2;
        if (petno3r < p3) petno3r = p3;
        spetno3 += petno3r;
    }
    //     Loop of all the other leaves, with the same computations.
    int nbrch;  // number of fruiting branches on a vegetative stem.
    int nnid;   // number of fruiting nodes on a fruiting branch.
    for (int k = 0; k < NumVegBranches; k++)  // loop of vegetative branches
    {
        nbrch = NumFruitBranches[k];
        for (int l = 0; l < nbrch; l++)  // loop of fruiting branches
        {
            nnid = NumNodes[k][l];
            for (int m = 0; m < nnid;
                 m++)  // loop of nodes on a fruiting branch
            {
                numl++;
                petno3r = p1 - LeafAge[k][l][m] * p2;
                if (petno3r < p3) petno3r = p3;
                spetno3 += petno3r;
            }  // m
        }      // l
    }          // k
    //     The return value of the function is the average ratio of NO3 to
    //  total N for all the petioles in the plant.
    double AvrgNO3Ratio = spetno3 / numl;
    return AvrgNO3Ratio;
}
//////////////////////////////////////////////////
void NitrogenAllocation()
//     This function computes the allocation of supplied nitrogen to
//  the plant parts. It is called from PlantNitrogen().
//
//     The following global and file scope variables are set here:
//       addnf, addnr, addnv, BurrNitrogen, LeafNitrogen, npool,
//       PetioleNitrogen, RootNitrogen, SeedNitrogen, SquareNitrogen,
//       StemNitrogen, xtran.
//     The following global and file scope variables are referenced in this
//     function:
//       reqtot, rqnbur, rqnlef, rqnpet, rqnrut, rqnsed, rqnsqr, rqnstm
{
    //     The following constant parameters are used:
    const double vseednmax =
        0.70;  // maximum proportion of N pool that can be added to seeds
    const double vsqrnmax =
        0.65;  // maximum proportion of N pool that can be added to squares
    const double vburnmax =
        0.65;  // maximum proportion of N pool that can be added to burrs
    const double vlfnmax =
        0.90;  // maximum proportion of N pool that can be added to leaves
    const double vstmnmax =
        0.70;  // maximum proportion of N pool that can be added to stems
    const double vpetnmax =
        0.75;  // maximum proportion of N pool that can be added to petioles
    //     If total N requirement is less than npool, add N required for growth
    //     to the N in
    //  each organ, compute added N to vegetative parts, fruiting parts and
    //  roots, and compute xtran as the difference between npool and the total N
    //  requirements.
    if (reqtot <= npool) {
        LeafNitrogen += rqnlef;
        PetioleNitrogen += rqnpet;
        StemNitrogen += rqnstm;
        RootNitrogen += rqnrut;
        SquareNitrogen += rqnsqr;
        SeedNitrogen += rqnsed;
        BurrNitrogen += rqnbur;
        addnv = rqnlef + rqnstm + rqnpet;
        addnf = rqnsqr + rqnsed + rqnbur;
        addnr = rqnrut;
        xtran = npool - reqtot;
        return;
    }
    //     If N requirement is greater than npool, execute the following:
    //     First priority is nitrogen supply to the growing seeds. It is assumed
    //     that up to
    //  vseednmax = 0.70 of the supplied N can be used by the seeds. Update seed
    //  N and addnf by the amount of nitrogen used for seed growth, and decrease
    //  npool by this amount. The same procedure is used for each organ,
    //  consecutively.
    double useofn;  // amount of nitrogen used in growth of a plant organ.
    if (rqnsed > 0) {
        useofn = min(vseednmax * npool, rqnsed);
        SeedNitrogen += useofn;
        addnf += useofn;
        npool -= useofn;
    }
    //     Next priority is for burrs, which can use N up to vburnmax = 0.65 of
    //     the
    //  remaining N pool, and for squares, which can use N up to vsqrnmax = 0.65
    if (rqnbur > 0) {
        useofn = min(vburnmax * npool, rqnbur);
        BurrNitrogen += useofn;
        addnf += useofn;
        npool -= useofn;
    }
    if (rqnsqr > 0) {
        useofn = min(vsqrnmax * npool, rqnsqr);
        SquareNitrogen += useofn;
        addnf += useofn;
        npool -= useofn;
    }
    //     Next priority is for leaves, which can use N up to vlfnmax = 0.90
    //  of the remaining N pool, for stems, up to vstmnmax = 0.70, and for
    //  petioles, up to vpetnmax = 0.75
    if (rqnlef > 0) {
        useofn = min(vlfnmax * npool, rqnlef);
        LeafNitrogen += useofn;
        addnv += useofn;
        npool -= useofn;
    }
    if (rqnstm > 0) {
        useofn = min(vstmnmax * npool, rqnstm);
        StemNitrogen += useofn;
        addnv += useofn;
        npool -= useofn;
    }
    if (rqnpet > 0) {
        useofn = min(vpetnmax * npool, rqnpet);
        PetioleNitrogen += useofn;
        addnv += useofn;
        npool -= useofn;
    }
    //     The remaining npool goes to root growth. If any npool remains
    //  it is defined as xtran.
    if (rqnrut > 0) {
        useofn = min(npool, rqnrut);
        RootNitrogen += useofn;
        addnr += useofn;
        npool -= useofn;
    }
    xtran = npool;
}
//////////////////////////
void ExtraNitrogenAllocation()
//     This function computes the allocation of extra nitrogen to the plant
//  parts. It is called from PlantNitrogen() if there is a non-zero xtran.
//
//     The following global and file scope variables are referenced here:
//       burres, BurrWeightGreenBolls, leafrs, petrs, rootrs, stemrs,
//       TotalLeafWeight, TotalPetioleWeight, TotalRootWeight, TotalStemWeight,
//       xtran.
//     The following global variables are set here:
//       BurrNitrogen, PetioleNitrogen, RootNitrogen, LeafNitrogen,
//       StemNitrogen.
{
    //     If there are any N reserves in the plant, allocate remaining xtran in
    //     proportion to
    //  the N reserves in each of these organs. Note: all reserves are in g per
    //  plant units.
    double addbur;   // reserve N to be added to the burrs.
    double addlfn;   // reserve N to be added to the leaves.
    double addpetn;  // reserve N to be added to the petioles.
    double addrt;    // reserve N to be added to the roots.
    double addstm;   // reserve N to be added to the stem.
    double rsum;     // sum of existing reserve N in plant parts.
    rsum = leafrs + petrs + stemrs + rootrs + burres;
    if (rsum > 0) {
        addlfn = xtran * leafrs / rsum;
        addpetn = xtran * petrs / rsum;
        addstm = xtran * stemrs / rsum;
        addrt = xtran * rootrs / rsum;
        addbur = xtran * burres / rsum;
    } else {
        //     If there are no reserves, allocate xtran in proportion to the dry
        //  weights in each of these organs.
        double vegwt;  // weight of vegetative plant parts, plus burrs.
        vegwt = TotalLeafWeight() + TotalPetioleWeight + TotalStemWeight +
                TotalRootWeight + BurrWeightGreenBolls;
        addlfn = xtran * TotalLeafWeight() / vegwt;
        addpetn = xtran * TotalPetioleWeight / vegwt;
        addstm = xtran * TotalStemWeight / vegwt;
        addrt = xtran * TotalRootWeight / vegwt;
        addbur = xtran * BurrWeightGreenBolls / vegwt;
    }
    //     Update N content in these plant parts. Note that at this stage of
    //  nitrogen allocation, only vegetative parts and burrs are updated (not
    //  seeds or squares).
    LeafNitrogen += addlfn;
    PetioleNitrogen += addpetn;
    StemNitrogen += addstm;
    RootNitrogen += addrt;
    BurrNitrogen += addbur;
}
//////////////////////////
void PlantNitrogenContent()
//     This function computes the concentrations of nitrogen in the dry
//  matter of the plant parts. It is called from PlantNitrogen(). It calls the
//  function PetioleNitrateN().
//
//     The following global variables are referenced here:
//       BurrNitrogen, BurrWeightGreenBolls, BurrWeightOpenBolls,
//       CottonWeightGreenBolls, CottonWeightOpenBolls, Gintot, LeafNitrogen,
//       PetioleNitrogen, RootNitrogen, SeedNitrogen, SquareNitrogen,
//       StemNitrogen, TotalLeafWeight, TotalPetioleWeight, TotalRootWeight,
//       TotalSquareWeight, TotalStemWeight.
//     The following global variables are set here:
//       BurrNConc, LeafNConc, PetioleNConc, PetioleNO3NConc, RootNConc,
//       SeedNConc, SquareNConc, StemNConc.
{
    //     The following constant parameter is used:
    const double seedratio = 0.64;
    //     Compute N concentration in plant organs as the ratio of N content to
    //     weight of dry matter.
    if (TotalLeafWeight() > 0.00001) LeafNConc = LeafNitrogen / TotalLeafWeight();
    if (TotalPetioleWeight > 0.00001) {
        PetioleNConc = PetioleNitrogen / TotalPetioleWeight;
        PetioleNO3NConc = PetioleNConc * PetioleNitrateN();
    }
    if (TotalStemWeight > 0) StemNConc = StemNitrogen / TotalStemWeight;
    if (TotalRootWeight > 0) RootNConc = RootNitrogen / TotalRootWeight;
    if (TotalSquareWeight > 0) SquareNConc = SquareNitrogen / TotalSquareWeight;
    double xxseed;  // weight of seeds in green and mature bolls.
    xxseed = CottonWeightOpenBolls * (1 - Gintot) +
             CottonWeightGreenBolls * seedratio;
    if (xxseed > 0) SeedNConc = SeedNitrogen / xxseed;
    double xxbur;  // weight of burrs in green and mature bolls.
    xxbur = BurrWeightOpenBolls + BurrWeightGreenBolls;
    if (xxbur > 0) BurrNConc = BurrNitrogen / xxbur;
}
//////////////////////////
void GetNitrogenStress()
//     This function computes the nitrogen stress factors. It is called from
//     PlantNitrogen().
//
//     The following global variables are set here:
//       NStressFruiting, NStressRoots, NStressVeg, NitrogenStress.
//     The following file scope variables are referenced in this function:
//       addnf, addnr, addnv, reqf, reqv, rqnrut
{
    //     Set the default values for the nitrogen stress coefficients to 1.
    NStressVeg = 1;
    NStressRoots = 1;
    NStressFruiting = 1;
    NitrogenStress = 1;
    //     Compute the nitrogen stress coefficients. NStressFruiting is the
    //     ratio of
    //  N added actually to the fruits, to their N requirements. NStressVeg is
    //  the same for vegetative shoot growth, and NStressRoots for roots. Also,
    //  an average stress coefficient for vegetative and reproductive organs is
    //  computed as NitrogenStress.
    //     Each stress coefficient has a value between 0 and 1.
    if (reqf > 0) {
        NStressFruiting = addnf / reqf;
        if (NStressFruiting > 1) NStressFruiting = 1;
        if (NStressFruiting < 0) NStressFruiting = 0;
    }
    if (reqv > 0) {
        NStressVeg = addnv / reqv;
        if (NStressVeg > 1) NStressVeg = 1;
        if (NStressVeg < 0) NStressVeg = 0;
    }
    if (rqnrut > 0) {
        NStressRoots = addnr / rqnrut;
        if (NStressRoots > 1) NStressRoots = 1;
        if (NStressRoots < 0) NStressRoots = 0;
    }
    if ((reqf + reqv) > 0) {
        NitrogenStress = (addnf + addnv) / (reqf + reqv);
        if (NitrogenStress > 1) NitrogenStress = 1;
        if (NitrogenStress < 0) NitrogenStress = 0;
    }
}
//////////////////////////
void NitrogenUptakeRequirement()
//     This function computes TotalRequiredN, the nitrogen requirements of the
//     plant -
//  to be used for simulating the N uptake from the soil (in function
//  NitrogenUptake() ) in the next day. It is called from PlantNitrogen().
//
//     The following global variables is set here:      TotalRequiredN
//     The following global variables are referenced here:
//       BurrNConc, BurrWeightGreenBolls, CottonWeightGreenBolls, Kday,
//       LeafNConc, PetioleNConc, reqtot, RootNConc, SeedNConc, SquareNConc,
//       StemNConc, StemWeight, TotalLeafWeight, TotalPetioleWeight,
//       TotalRootWeight, TotalSquareWeight, TotalStemWeight.
{
    //     The following constant parameters are used:
    const double seedcn1 =
        .045;  //   further requirement for existing seed tissue.
    const double seedratio =
        0.64;  // the ratio of seeds to seedcotton in green bolls.
    const double vnreqlef =
        .042;  //  coefficient for computing N uptake requirements of leaves
    const double vnreqpet =
        .036;  //  coefficient for computing N uptake requirements of petioles
    const double vnreqstm =
        .012;  //  coefficient for computing N uptake requirements of stems
    const double vnreqrt =
        .010;  //  coefficient for computing N uptake requirements of roots
    const double vnreqsqr =
        .024;  //  coefficient for computing N uptake requirements of squares
    const double vnreqbur =
        .012;  //  coefficient for computing N uptake requirements of burrs
    const int voldstm = 32;  // active stem tissue growth period, calendar days.
                             //
    TotalRequiredN = reqtot;
    //     After the requirements of today's growth are supplied, N is also
    //  required for supplying necessary functions in other active plant
    //  tissues.
    //    Add nitrogen uptake required for leaf and petiole tissue to
    //    TotalRequiredN.
    if (LeafNConc < vnreqlef)
        TotalRequiredN += TotalLeafWeight() * (vnreqlef - LeafNConc);
    if (PetioleNConc < vnreqpet)
        TotalRequiredN += TotalPetioleWeight * (vnreqpet - PetioleNConc);
    //    The active stem tissue is the stem formed during the last voldstm
    //  days (32 calendar days). add stem requirement to TotalRequiredN.
    double grstmwt;  // weight of actively growing stems.
    int kkday;  // day (from emergence) of oldest actively growing stem tissue.
    kkday = Kday - voldstm;
    if (kkday < 1)
        grstmwt = TotalStemWeight;
    else
        grstmwt = TotalStemWeight - StemWeight[kkday];
    if (StemNConc < vnreqstm)
        TotalRequiredN += grstmwt * (vnreqstm - StemNConc);
    //     Compute nitrogen uptake requirement for existing tissues of roots,
    //  squares, and seeds and burrs of green bolls. Add it to TotalRequiredN.
    if (RootNConc < vnreqrt)
        TotalRequiredN += TotalRootWeight * (vnreqrt - RootNConc);
    if (SquareNConc < vnreqsqr)
        TotalRequiredN += TotalSquareWeight * (vnreqsqr - SquareNConc);
    if (SeedNConc < seedcn1)
        TotalRequiredN +=
            CottonWeightGreenBolls * seedratio * (seedcn1 - SeedNConc);
    if (BurrNConc < vnreqbur)
        TotalRequiredN += BurrWeightGreenBolls * (vnreqbur - BurrNConc);
}
///////////////////////////////////////////////////////////////////////////////
void PlantNitrogenBal()
//     This function calculates the nitrogen balance in the cotton
//  plant, for diagnostic purposes. It is called from SimulateThisDay().
//     Units are g per plant.
//
//     The following global variables are referenced here:
//       BurrNitrogen, CumPlantNLoss, Kday, LeafNitrogen, PetioleNitrogen,
//       RootNitrogen, SeedNitrogen, SquareNitrogen, StemNitrogen, SupplyNH4N,
//       SupplyNO3N.
//
{
    //     addn is the cumulative supply of nitrogen from the soil,
    //  including N content of the seedling at emergence.
    static double addn;
    if (Kday == 1) addn = RootNitrogen + StemNitrogen + LeafNitrogen;
    addn += SupplyNO3N + SupplyNH4N;
    double plantn;  // total nitrogen in plant.
    plantn = RootNitrogen + StemNitrogen + LeafNitrogen + PetioleNitrogen +
             SeedNitrogen + BurrNitrogen + SquareNitrogen;
    // CumPlantNLoss is the amount lost in abscised fruit and leaves and dead
    // roots.
    double balpn;  // the plant nitrogen balance, which should be zero.
    balpn = addn - plantn - CumPlantNLoss;
}
