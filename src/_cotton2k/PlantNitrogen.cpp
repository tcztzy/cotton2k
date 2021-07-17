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
//
#include <string>
#include "global.h"
#include "Simulation.hpp"

using namespace std;

void NitrogenRequirement(State &, const int &, const int &, double);

void NitrogenSupply(State &);

double PetioleNitrateN(State &);

void NitrogenAllocation(State &);

void ExtraNitrogenAllocation(State &);

void PlantNitrogenContent(State &);

void GetNitrogenStress(State &);

void NitrogenUptakeRequirement(State &);

//  File scope variables
double burres; // reserve N in burrs, in g per plant.
double leafrs; // reserve N in leaves, in g per plant.
double petrs;  // reserve N in petioles, in g per plant.
double stemrs; // reserve N in stems, in g per plant.
double rootrs; // reserve N in roots, in g per plant.
double rqnbur; // nitrogen requirement for burr growth.
double rqnlef; // nitrogen requirement for leaf growth.
double rqnpet; // nitrogen requirement for petiole growth.
double rqnrut; // nitrogen requirement for root growth.
double rqnsed; // nitrogen requirement for seed growth.
double rqnsqr; // nitrogen requirement for square growth.
double rqnstm; // nitrogen requirement for stem growth.
double reqv;   // nitrogen requirement for vegetative shoot growth.
double reqf;   // nitrogen requirement for fruit growth.
double reqtot; // total nitrogen requirement for plant growth.
double npool;  // total nitrogen available for growth.
double uptn;   // nitrogen uptake from the soil, g per plant.
double xtran;  // amount of nitrogen not used for growth of plant parts.
double addnf;  // daily added nitrogen to fruit, g per plant.
double addnr;  // daily added nitrogen to root, g per plant.
double addnv;  // daily added nitrogen to vegetative shoot, g per plant.
//////////////////////////////////////////////////
void PlantNitrogen(Simulation &sim, uint32_t u)
//     This function simulates the nitrogen accumulation and distribution in cotton plants,
//  and computes nitrogen stresses. It is called from SimulateThisDay().
//
//     The maximum and minimum N concentrations for the various organs are modified from
//  those reported by: Jones et. al. (1974): Development of a nitrogen balance for cotton
//  growth models: a first approximation. Crop Sci. 14:541-546.
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
//       BurrNConc, Kday, LeafNConc,
//       PetioleNConc, PetioleNConc, RootNConc,
//       RootNitrogen, SeedNConc, StemNConc.
//     The following global and file scope variables are set here:
//       addnf, addnr, addnv, burres, leafrs, npool, petrs, reqf, reqtot, reqv, rootrs,
//       rqnbur, rqnlef, rqnpet, rqnrut, rqnsed, rqnsqr, rqnstm, stemrs, uptn, xtran.
{
    State &state = sim.states[u];
    //     Assign zero to some variables.
    leafrs = 0; // reserve N in leaves, in g per plant.
    petrs = 0;  // reserve N in petioles, in g per plant.
    stemrs = 0; // reserve N in stems, in g per plant.
    rootrs = 0; // reserve N in roots, in g per plant.
    reqf = 0;   // nitrogen requirement for fruit growth.
    reqtot = 0; // total nitrogen requirement for plant growth.
    reqv = 0;   // nitrogen requirement for vegetative shoot growth.
    rqnbur = 0; // nitrogen requirement for burr growth.
    rqnlef = 0; // nitrogen requirement for leaf growth.
    rqnpet = 0; // nitrogen requirement for petiole growth.
    rqnrut = 0; // nitrogen requirement for root growth.
    rqnsed = 0; // nitrogen requirement for seed growth.
    rqnsqr = 0; // nitrogen requirement for square growth.
    rqnstm = 0; // nitrogen requirement for stem growth.
    npool = 0;  // total nitrogen available for growth.
    uptn = 0;   // nitrogen uptake from the soil, g per plant.
    xtran = 0;  // amount of nitrogen not used for growth of plant parts.
    addnf = 0;  // daily added nitrogen to fruit, g per plant.
    addnr = 0;  // daily added nitrogen to root, g per plant.
    addnv = 0;  // daily added nitrogen to vegetative shoot, g per plant.

    //   The following subroutines are now called:
    NitrogenRequirement(state, state.daynum, sim.day_emerge, state.extra_carbon); //  computes the N requirements for growth.
    NitrogenSupply(state);                                                      //  computes the supply of N from uptake and reserves.
    NitrogenAllocation(state);                                                                         //  computes the allocation of N in the plant.
    if (xtran > 0)
        ExtraNitrogenAllocation(state); // computes the further allocation of N in the plant
    PlantNitrogenContent(state);   // computes the concentrations of N in plant dry matter.
    GetNitrogenStress(state);           //  computes nitrogen stress factors.
    NitrogenUptakeRequirement(state);   // computes N requirements for uptake
}

//////////////////////////
void NitrogenRequirement(State &state, const int &Daynum, const int &DayEmerge, double ExtraCarbon)
//     This function computes the N requirements for growth. It is called from PlantNitrogen{}.
//
//     The following global and file scope variables are set in this function:
//       PetioleNConc, reqf, reqtot, reqv, rqnbur,
//       rqnlef, rqnpet, rqnrut, rqnsed, rqnsqr, rqnstm.
{
    //     The following constant parameters are used:
    const double burcn0 = .026;   //  maximum N content for burrs
    const double lefcn0 = .064;   //  maximum N content for leaves
    const double petcn0 = .032;   //  maximum N content for petioles
    const double rootcn0 = .018;  //  maximum N content for roots
    const double seedcn0 = .036;  //  maximum N content for seeds
    const double seedcn1 = .045;  //  additional requirement of N for existing seed tissue.
    const double seedratio = .64; //  ratio of seeds to seedcotton in green bolls.
    const double sqrcn0 = .024;   //  maximum N content for squares
    const double stmcn0 = .036;   //  maximum N content for stems
                                  //     On emergence, assign initial values to petiole N concentrations.
    if (Daynum <= DayEmerge)
    {
        state.petiole_nitrogen_concentration = petcn0;
        state.petiole_nitrate_nitrogen_concentration = petcn0;
    }
    //     Compute the nitrogen requirements for growth, by multiplying
    //  the daily added dry weight by the maximum N content of each organ.
    //  Nitrogen requirements are based on actual growth rates.
    //     These N requirements will be used to compute the allocation
    //  of N to plant parts and the nitrogen stress factors.
    //     All nitrogen requirement variables are in g N per plant.
    rqnlef = lefcn0 * state.total_actual_leaf_growth;                //      for leaf blade
    rqnpet = petcn0 * state.total_actual_petiole_growth;             //      for petiole
    rqnstm = stmcn0 * state.actual_stem_growth;                              //      for stem
    // Add ExtraCarbon to carbon_allocated_for_root_growth to compute the total supply of carbohydrates for root growth.
    rqnrut = rootcn0 * (state.carbon_allocated_for_root_growth + ExtraCarbon); // for root
    rqnsqr = state.actual_square_growth * sqrcn0;                            //      for squares
    double rqnsed1, rqnsed2;                                         // components of seed N requirements.
    rqnsed1 = state.actual_boll_growth * seedratio * seedcn0;                //   for seed growth
                                                                     //     The N required for replenishing the N content of existing seed
                                                                     //  tissue (rqnsed2) is added to seed growth requirement.
    if (state.green_bolls_weight > state.actual_boll_growth)
    {
        double rseedn; // existing ratio of N to dry matter in the seeds.
        rseedn = state.seed_nitrogen / ((state.green_bolls_weight - state.actual_boll_growth) * seedratio);
        rqnsed2 = (state.green_bolls_weight - state.actual_boll_growth) * seedratio * (seedcn1 - rseedn);
        if (rqnsed2 < 0)
            rqnsed2 = 0;
    }
    else
        rqnsed2 = 0;
    //
    rqnsed = rqnsed1 + rqnsed2;         //   total requirement for seeds
    rqnbur = state.actual_burr_growth * burcn0; //   for burrs
    reqf = rqnsqr + rqnsed + rqnbur;    //    total for fruit
    reqv = rqnlef + rqnpet + rqnstm;    //    total for shoot
    reqtot = rqnrut + reqv + reqf;      //     total N requirement
}

//////////////////////////
void NitrogenSupply(State &state)
//     This function computes the supply of N by uptake from the soil reserves,
//  It is called from PlantNitrogen(). It calls function PetioleNitrateN().
//
//     The following global variables are referenced here:
//       reqtot.
//     The following global and file scope variables are set here:
//       burres, leafrs, npool, petrs,
//       RootNitrogen, rootrs, stemrs, uptn, xtran.
{
    //     The following constant parameters are used:
    const double MobilizNFractionBurrs = 0.08;    //  fraction of N mobilizable for burrs
    const double MobilizNFractionLeaves = 0.09;   //  fraction of N mobilizable for leaves and petioles
    const double MobilizNFractionStemRoot = 0.40; //  fraction of N mobilizable for stems and roots
    const double vburnmin = .006;                 //  minimum N contents of burrs
    const double vlfnmin = .018;                  //  minimum N contents of leaves
    const double vpetnmin = .005;                 //  minimum N contents of petioles, non-nitrate fraction.
    const double vpno3min = .002;                 //  minimum N contents of petioles, nitrate fraction.
    const double vrtnmin = .010;                  //  minimum N contents of roots
    const double vstmnmin = .006;                 //  minimum N contents of stems
                                                  //     uptn is the total supply of nitrogen to the plant by uptake of nitrate and ammonium.
    uptn = state.supplied_nitrate_nitrogen + state.supplied_ammonium_nitrogen;
    double resn; // total reserve N, in g per plant.
                 //     If total N requirement is less than the supply, define npool as the supply and
                 //  assign zero to the N reserves in all organs.
    if (reqtot <= uptn)
    {
        npool = uptn;
        leafrs = 0;
        petrs = 0;
        stemrs = 0;
        rootrs = 0;
        burres = 0;
        xtran = 0;
        resn = 0;
    }
    else
    {
        //     If total N requirement exceeds the supply, compute the nitrogen
        //  reserves in the plant. The reserve N in an organ is defined as a fraction
        //  of the nitrogen content exceeding a minimum N content in it.
        //     The N reserves in leaves, petioles, stems, roots and burrs of
        //  green bolls are computed, and their N content updated.
        leafrs = (state.leaf_nitrogen - vlfnmin * state.leaf_weight) * MobilizNFractionLeaves;
        if (leafrs < 0)
            leafrs = 0;
        state.leaf_nitrogen -= leafrs;
        //     The petiole N content is subdivided to nitrate and non-nitrate. The
        //  nitrate ratio in the petiole N is computed by calling function PetioleNitrateN().
        //  Note that the nitrate fraction is more available for redistribution.
        double rpetno3; // ratio of NO3 N to total N in petioles.
        rpetno3 = PetioleNitrateN(state);
        double petrs1, petrs2; // components of reserve N in petioles, for non-NO3 and NO3 origin, respectively.
        petrs1 = (state.petiole_nitrogen * (1 - rpetno3) - vpetnmin * state.petiole_weight) * MobilizNFractionLeaves;
        if (petrs1 < 0)
            petrs1 = 0;
        petrs2 = (state.petiole_nitrogen * rpetno3 - vpno3min * state.petiole_weight) * MobilizNFractionLeaves;
        if (petrs2 < 0)
            petrs2 = 0;
        petrs = petrs1 + petrs2;
        state.petiole_nitrogen -= petrs;
        //  Stem N reserves.
        stemrs = (state.stem_nitrogen - vstmnmin * state.stem_weight) * MobilizNFractionStemRoot;
        if (stemrs < 0)
            stemrs = 0;
        state.stem_nitrogen -= stemrs;
        //  Root N reserves
        rootrs = (RootNitrogen - vrtnmin * state.root_weight) * MobilizNFractionStemRoot;
        if (rootrs < 0)
            rootrs = 0;
        RootNitrogen -= rootrs;
        //  Burr N reserves
        if (state.green_bolls_burr_weight > 0)
        {
            burres = (state.burr_nitrogen - vburnmin * state.green_bolls_burr_weight) * MobilizNFractionBurrs;
            if (burres < 0)
                burres = 0;
            state.burr_nitrogen -= burres;
        }
        else
            burres = 0;
        //     The total reserves, resn, are added to the amount taken up from the soil, for computing
        //  npool. Note that N of seeds or squares is not available for redistribution in the plant.
        resn = leafrs + petrs + stemrs + rootrs + burres;
        npool = uptn + resn;
    }
}

//////////////////////////
double PetioleNitrateN(State &state)
//     This function computes the ratio of NO3 nitrogen to total N in the petioles.
//  It is called from NitrogenSupply() and PlantNitrogenContent().
//     The following global variables are referenced here:
//       LeafAge, NumFruitBranches, NumNodes, NumVegBranches.
{
    //     The following constant parameters are used:
    const double p1 = 0.96;  //  the maximum ratio (of NO3 to total N in petioles).
    const double p2 = 0.015; //  the rate of decline of this ratio with age.
    const double p3 = 0.02;  //  the minimum ratio
                             //     The ratio of NO3 to total N in each individual petiole is computed
                             //  as a linear function of leaf age. It is assumed that this ratio is
                             //  maximum for young leaves and is declining with leaf age.
    double spetno3 = 0;      // sum of petno3r.
    double petno3r;          // ratio of NO3 to total N in an individual petiole.
                             //     Loop of prefruiting node leaves.
    for (int j = 0; j < state.number_of_pre_fruiting_nodes; j++)
    {
        petno3r = p1 - state.age_of_pre_fruiting_nodes[j] * p2;
        if (petno3r < p3)
            petno3r = p3;
        spetno3 += petno3r;
    }
    //     Loop of all the other leaves, with the same computations.
    int numl = state.number_of_pre_fruiting_nodes; // number of petioles computed.
    int nbrch;                 // number of fruiting branches on a vegetative stem.
    int nnid;                  // number of fruiting nodes on a fruiting branch.
    for (int k = 0; k < state.number_of_vegetative_branches; k++)
        for (int l = 0; l < state.vegetative_branches[k].number_of_fruiting_branches; l++)
        {
            numl += state.vegetative_branches[k].fruiting_branches[l].number_of_fruiting_nodes;
            for (int m = 0; m < state.vegetative_branches[k].fruiting_branches[l].number_of_fruiting_nodes; m++)
            {
                petno3r = p1 - state.vegetative_branches[k].fruiting_branches[l].nodes[m].leaf.age * p2;
                if (petno3r < p3)
                    petno3r = p3;
                spetno3 += petno3r;
            }
        }
    //     The return value of the function is the average ratio of NO3 to
    //  total N for all the petioles in the plant.
    double AvrgNO3Ratio = spetno3 / numl;
    return AvrgNO3Ratio;
}

//////////////////////////////////////////////////
void NitrogenAllocation(State &state)
//     This function computes the allocation of supplied nitrogen to
//  the plant parts. It is called from PlantNitrogen().
//
//     The following global and file scope variables are set here:
//       addnf, addnr, addnv, npool, RootNitrogen, xtran.
//     The following global and file scope variables are referenced in this function:
//       reqtot, rqnbur, rqnlef, rqnpet, rqnrut, rqnsed, rqnsqr, rqnstm
{
    //     The following constant parameters are used:
    const double vseednmax = 0.70; // maximum proportion of N pool that can be added to seeds
    const double vsqrnmax = 0.65;  // maximum proportion of N pool that can be added to squares
    const double vburnmax = 0.65;  // maximum proportion of N pool that can be added to burrs
    const double vlfnmax = 0.90;   // maximum proportion of N pool that can be added to leaves
    const double vstmnmax = 0.70;  // maximum proportion of N pool that can be added to stems
    const double vpetnmax = 0.75;  // maximum proportion of N pool that can be added to petioles
                                   //     If total N requirement is less than npool, add N required for growth to the N in
                                   //  each organ, compute added N to vegetative parts, fruiting parts and roots, and compute
                                   //  xtran as the difference between npool and the total N requirements.
    if (reqtot <= npool)
    {
        state.leaf_nitrogen += rqnlef;
        state.petiole_nitrogen += rqnpet;
        state.stem_nitrogen += rqnstm;
        RootNitrogen += rqnrut;
        state.square_nitrogen += rqnsqr;
        state.seed_nitrogen += rqnsed;
        state.burr_nitrogen += rqnbur;
        addnv = rqnlef + rqnstm + rqnpet;
        addnf = rqnsqr + rqnsed + rqnbur;
        addnr = rqnrut;
        xtran = npool - reqtot;
        return;
    }
    //     If N requirement is greater than npool, execute the following:
    //     First priority is nitrogen supply to the growing seeds. It is assumed that up to
    //  vseednmax = 0.70 of the supplied N can be used by the seeds. Update seed N and
    //  addnf by the amount of nitrogen used for seed growth, and decrease npool by this amount.
    //  The same procedure is used for each organ, consecutively.
    double useofn; // amount of nitrogen used in growth of a plant organ.
    if (rqnsed > 0)
    {
        useofn = std::min(vseednmax * npool, rqnsed);
        state.seed_nitrogen += useofn;
        addnf += useofn;
        npool -= useofn;
    }
    //     Next priority is for burrs, which can use N up to vburnmax = 0.65 of the
    //  remaining N pool, and for squares, which can use N up to vsqrnmax = 0.65
    if (rqnbur > 0)
    {
        useofn = std::min(vburnmax * npool, rqnbur);
        state.burr_nitrogen += useofn;
        addnf += useofn;
        npool -= useofn;
    }
    if (rqnsqr > 0)
    {
        useofn = std::min(vsqrnmax * npool, rqnsqr);
        state.square_nitrogen += useofn;
        addnf += useofn;
        npool -= useofn;
    }
    //     Next priority is for leaves, which can use N up to vlfnmax = 0.90
    //  of the remaining N pool, for stems, up to vstmnmax = 0.70, and for
    //  petioles, up to vpetnmax = 0.75
    if (rqnlef > 0)
    {
        useofn = std::min(vlfnmax * npool, rqnlef);
        state.leaf_nitrogen += useofn;
        addnv += useofn;
        npool -= useofn;
    }
    if (rqnstm > 0)
    {
        useofn = std::min(vstmnmax * npool, rqnstm);
        state.stem_nitrogen += useofn;
        addnv += useofn;
        npool -= useofn;
    }
    if (rqnpet > 0)
    {
        useofn = std::min(vpetnmax * npool, rqnpet);
        state.petiole_nitrogen += useofn;
        addnv += useofn;
        npool -= useofn;
    }
    //     The remaining npool goes to root growth. If any npool remains
    //  it is defined as xtran.
    if (rqnrut > 0)
    {
        useofn = std::min(npool, rqnrut);
        RootNitrogen += useofn;
        addnr += useofn;
        npool -= useofn;
    }
    xtran = npool;
}

//////////////////////////
void ExtraNitrogenAllocation(State &state)
//     This function computes the allocation of extra nitrogen to the plant
//  parts. It is called from PlantNitrogen() if there is a non-zero xtran.
//
//     The following global and file scope variables are referenced here:
//       burres, leafrs, petrs, rootrs, stemrs, xtran.
//     The following global variables are set here:
//       RootNitrogen.
{
    //     If there are any N reserves in the plant, allocate remaining xtran in proportion to
    //  the N reserves in each of these organs. Note: all reserves are in g per plant units.
    double addbur;  // reserve N to be added to the burrs.
    double addlfn;  // reserve N to be added to the leaves.
    double addpetn; // reserve N to be added to the petioles.
    double addrt;   // reserve N to be added to the roots.
    double addstm;  // reserve N to be added to the stem.
    double rsum;    // sum of existing reserve N in plant parts.
    rsum = leafrs + petrs + stemrs + rootrs + burres;
    if (rsum > 0)
    {
        addlfn = xtran * leafrs / rsum;
        addpetn = xtran * petrs / rsum;
        addstm = xtran * stemrs / rsum;
        addrt = xtran * rootrs / rsum;
        addbur = xtran * burres / rsum;
    }
    else
    {
        //     If there are no reserves, allocate xtran in proportion to the dry
        //  weights in each of these organs.
        double vegwt; // weight of vegetative plant parts, plus burrs.
        vegwt = state.leaf_weight + state.petiole_weight + state.stem_weight + state.root_weight + state.green_bolls_burr_weight;
        addlfn = xtran * state.leaf_weight / vegwt;
        addpetn = xtran * state.petiole_weight / vegwt;
        addstm = xtran * state.stem_weight / vegwt;
        addrt = xtran * state.root_weight / vegwt;
        addbur = xtran * state.green_bolls_burr_weight / vegwt;
    }
    //     Update N content in these plant parts. Note that at this stage of
    //  nitrogen allocation, only vegetative parts and burrs are updated (not
    //  seeds or squares).
    state.leaf_nitrogen += addlfn;
    state.petiole_nitrogen += addpetn;
    state.stem_nitrogen += addstm;
    RootNitrogen += addrt;
    state.burr_nitrogen += addbur;
}

//////////////////////////
void PlantNitrogenContent(State &state)
//     This function computes the concentrations of nitrogen in the dry
//  matter of the plant parts. It is called from PlantNitrogen(). It calls the
//  function PetioleNitrateN().
//
//     The following global variables are referenced here:
//       RootNitrogen.
//     The following global variables are set here:
//       BurrNConc, LeafNConc, PetioleNConc, RootNConc,
//       SeedNConc, StemNConc.
{
    //     The following constant parameter is used:
    const double seedratio = 0.64;
    //     Compute N concentration in plant organs as the ratio of N content to weight of dry matter.
    if (state.leaf_weight > 0.00001)
        state.leaf_nitrogen_concentration = state.leaf_nitrogen / state.leaf_weight;
    if (state.petiole_weight > 0.00001)
    {
        state.petiole_nitrogen_concentration = state.petiole_nitrogen / state.petiole_weight;
        state.petiole_nitrate_nitrogen_concentration = state.petiole_nitrogen_concentration * PetioleNitrateN(state);
    }
    if (state.stem_weight > 0)
        StemNConc = state.stem_nitrogen / state.stem_weight;
    if (state.root_weight > 0)
        state.root_nitrogen_concentration = RootNitrogen / state.root_weight;
    if (state.square_weight > 0)
        state.square_nitrogen_concentration = state.square_nitrogen / state.square_weight;
    double xxseed; // weight of seeds in green and mature bolls.
    xxseed = state.open_bolls_weight * (1 - state.ginning_percent) + state.green_bolls_weight * seedratio;
    if (xxseed > 0)
        state.seed_nitrogen_concentration = state.seed_nitrogen / xxseed;
    double xxbur; // weight of burrs in green and mature bolls.
    xxbur = state.open_bolls_burr_weight + state.green_bolls_burr_weight;
    if (xxbur > 0)
        BurrNConc = state.burr_nitrogen / xxbur;
}

//////////////////////////
void GetNitrogenStress(State &state)
//     This function computes the nitrogen stress factors. It is called from PlantNitrogen().
//     The following file scope variables are referenced in this function:
//       addnf, addnr, addnv, reqf, reqv, rqnrut
{
    //     Set the default values for the nitrogen stress coefficients to 1.
    state.nitrogen_stress_vegetative = 1;
    state.nitrogen_stress_root = 1;
    state.nitrogen_stress_fruiting = 1;
    state.nitrogen_stress = 1;
    //     Compute the nitrogen stress coefficients. state.nitrogen_stress_fruiting is the ratio of
    //  N added actually to the fruits, to their N requirements. state.nitrogen_stress_vegetative is the
    //  same for vegetative shoot growth, and state.nitrogen_stress_root for roots. Also, an average
    //  stress coefficient for vegetative and reproductive organs is computed as NitrogenStress.
    //     Each stress coefficient has a value between 0 and 1.
    if (reqf > 0)
    {
        state.nitrogen_stress_fruiting = addnf / reqf;
        if (state.nitrogen_stress_fruiting > 1)
            state.nitrogen_stress_fruiting = 1;
        if (state.nitrogen_stress_fruiting < 0)
            state.nitrogen_stress_fruiting = 0;
    }
    if (reqv > 0)
    {
        state.nitrogen_stress_vegetative = addnv / reqv;
        if (state.nitrogen_stress_vegetative > 1)
            state.nitrogen_stress_vegetative = 1;
        if (state.nitrogen_stress_vegetative < 0)
            state.nitrogen_stress_vegetative = 0;
    }
    if (rqnrut > 0)
    {
        state.nitrogen_stress_root = addnr / rqnrut;
        if (state.nitrogen_stress_root > 1)
            state.nitrogen_stress_root = 1;
        if (state.nitrogen_stress_root < 0)
            state.nitrogen_stress_root = 0;
    }
    if ((reqf + reqv) > 0)
    {
        state.nitrogen_stress = (addnf + addnv) / (reqf + reqv);
        if (state.nitrogen_stress > 1)
            state.nitrogen_stress = 1;
        if (state.nitrogen_stress < 0)
            state.nitrogen_stress = 0;
    }
}

//////////////////////////
void NitrogenUptakeRequirement(State &state)
//     This function computes TotalRequiredN, the nitrogen requirements of the plant -
//  to be used for simulating the N uptake from the soil (in function NitrogenUptake() )
//  in the next day. It is called from PlantNitrogen().
//
//     The following global variables is set here:      TotalRequiredN
//     The following global variables are referenced here:
//       BurrNConc, LeafNConc, PetioleNConc, reqtot,
//       StemNConc, StemWeight.
{
    //     The following constant parameters are used:
    const double seedcn1 = .045;   //   further requirement for existing seed tissue.
    const double seedratio = 0.64; // the ratio of seeds to seedcotton in green bolls.
    const double vnreqlef = .042;  //  coefficient for computing N uptake requirements of leaves
    const double vnreqpet = .036;  //  coefficient for computing N uptake requirements of petioles
    const double vnreqstm = .012;  //  coefficient for computing N uptake requirements of stems
    const double vnreqrt = .010;   //  coefficient for computing N uptake requirements of roots
    const double vnreqsqr = .024;  //  coefficient for computing N uptake requirements of squares
    const double vnreqbur = .012;  //  coefficient for computing N uptake requirements of burrs
    const int voldstm = 32;        // active stem tissue growth period, calendar days.
                                   //
    state.total_required_nitrogen = reqtot;
    //     After the requirements of today's growth are supplied, N is also
    //  required for supplying necessary functions in other active plant tissues.
    //    Add nitrogen uptake required for leaf and petiole tissue to TotalRequiredN.
    if (state.leaf_nitrogen_concentration < vnreqlef)
        state.total_required_nitrogen += state.leaf_weight * (vnreqlef - state.leaf_nitrogen_concentration);
    if (state.petiole_nitrogen_concentration < vnreqpet)
        state.total_required_nitrogen += state.petiole_weight * (vnreqpet - state.petiole_nitrogen_concentration);
    //    The active stem tissue is the stem formed during the last voldstm
    //  days (32 calendar days). add stem requirement to TotalRequiredN.
    double grstmwt; // weight of actively growing stems.
    int kkday;      // day (from emergence) of oldest actively growing stem tissue.
    kkday = state.kday - voldstm;
    if (kkday < 1)
        grstmwt = state.stem_weight;
    else
        grstmwt = state.stem_weight - StemWeight[kkday];
    if (StemNConc < vnreqstm)
        state.total_required_nitrogen += grstmwt * (vnreqstm - StemNConc);
    //     Compute nitrogen uptake requirement for existing tissues of roots,
    //  squares, and seeds and burrs of green bolls. Add it to TotalRequiredN.
    if (state.root_nitrogen_concentration < vnreqrt)
        state.total_required_nitrogen += state.root_weight * (vnreqrt - state.root_nitrogen_concentration);
    if (state.square_nitrogen_concentration < vnreqsqr)
        state.total_required_nitrogen += state.square_weight * (vnreqsqr - state.square_nitrogen_concentration);
    if (state.seed_nitrogen_concentration < seedcn1)
        state.total_required_nitrogen += state.green_bolls_weight * seedratio * (seedcn1 - state.seed_nitrogen_concentration);
    if (BurrNConc < vnreqbur)
        state.total_required_nitrogen += state.green_bolls_burr_weight * (vnreqbur - BurrNConc);
}
