// File  PlantGrowth_2.cpp
//
//   functions in this file:
// PotentialStemGrowth()
// PotentialLeafGrowth()
// TemperatureOnLeafGrowthRate()
// PotentialFruitGrowth()
// TemperatureOnFruitGrowthRate()
//
#include "global.h"
#include "State.h"

extern "C" {
    double TemperatureOnLeafGrowthRate(double);
    double TemperatureOnFruitGrowthRate(double);
}

//////////////////////////
void PotentialLeafGrowth(State &state)
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
//        AgeOfPreFruNode, AvrgDailyTemp, DensityFactor, LeafAreaNodes, 
//        LeafAreaMainStem, LeafAreaPreFru, NumFruitBranches, NumNodes, NumPreFruNodes,
//        NumVegBranches, pixda, VarPar, WaterStress.
//     The following global variables are set here:
//        LeafWeightAreaRatio, PotGroAllLeaves, PotGroAllPetioles,  
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
//     Calculate wtfstrs, the effect of leaf water stress onLeafWeightAreaRatio (the ratio
//  of leaf dry weight to leaf area). This has also been empirically calibrated in COTTON2K.
    double wtfstrs = vpotlf[1] + vpotlf[2] * (1 - wstrlf);
//     Compute the ratio of leaf dry weight increment to leaf area increment (g per dm2), 
//  as a function of average daily temperature and water stress. Parameters for the effect 
//  of temperature are adapted from GOSSYM.
    double tdday = AvrgDailyTemp;  // limited value of today's average temperature.
    if (tdday < vpotlf[3])
        tdday = vpotlf[3];
    LeafWeightAreaRatio = wtfstrs / (vpotlf[4] + tdday * (vpotlf[5] - tdday * vpotlf[6]));
//     Assign zero to total potential growth of leaf and petiole.
    PotGroAllLeaves = 0;
    PotGroAllPetioles = 0;
    double c = 0; // parameter of the leaf growth rate equation.
    double smax = 0; // maximum possible leaf area, a parameter of the leaf growth rate equation.
    double rate; // growth rate of area of a leaf.
//     Compute the potential growth rate of prefruiting leaves. 
//  smax and c are functions of prefruiting node number.
    for (int j = 0; j < NumPreFruNodes; j++) // loop by prefruiting node.
    {
        if (LeafAreaPreFru[j] <= 0) {
            PotGroLeafAreaPreFru[j] = 0;
            PotGroLeafWeightPreFru[j] = 0;
            PotGroPetioleWeightPreFru[j] = 0;
        } else {
            int jp1 = j + 1;
            smax = jp1 * (VarPar[2] - VarPar[3] * jp1);
            if (smax < VarPar[4])
                smax = VarPar[4];
            c = vpotlf[7] + vpotlf[8] * jp1 * (jp1 - vpotlf[9]);
            rate = smax * c * p * exp(-c * pow(AgeOfPreFruNode[j], p))
                   * pow(AgeOfPreFruNode[j], (p - 1));
//     Growth rate is modified by water stress, pix and a function of average temperature.  
//     Compute potential growth of leaf area, leaf weight and petiole 
//  weight for leaf on node j. Add leaf weight potential growth to PotGroAllLeaves.
//  Add potential growth of petiole weight to PotGroAllPetioles.
            if (rate >= 1e-12) {
                PotGroLeafAreaPreFru[j] = rate * wstrlf * pixda
                                          * TemperatureOnLeafGrowthRate(AvrgDailyTemp);
                PotGroLeafWeightPreFru[j] = PotGroLeafAreaPreFru[j] * LeafWeightAreaRatio;
                PotGroPetioleWeightPreFru[j] = PotGroLeafAreaPreFru[j]
                                               * LeafWeightAreaRatio * vpotlf[13];
                PotGroAllLeaves += PotGroLeafWeightPreFru[j];
                PotGroAllPetioles += PotGroPetioleWeightPreFru[j];
            }  // rate
        }  // LeafAreaPreFru
    }    // NumPreFruNodes
//     denfac is the effect of plant density on leaf growth rate.
    double denfac = 1 - vpotlf[12] * (1 - DensityFactor);
    int nbrch;  // number of fruiting branches on a vegetative stem.
    int nnid;   // number of nodes on a fruiting branch.
//    Loop for all vegetative stems.
    for (int k = 0; k < NumVegBranches; k++) // loop of vegetative branches
    {
        nbrch = NumFruitBranches[k];
        for (int l = 0; l < nbrch; l++) // loop of fruiting branches
        {
//     smax and c are  functions of fruiting branch number.
//     smax is modified by plant density, using the density factor denfac.
//     Compute potential main stem leaf growth, assuming that the main
//  stem leaf is initiated at the same time as leaf (k,l,0).
            if (LeafAreaMainStem[k][l] <= 0) {
                PotGroLeafAreaMainStem[k][l] = 0;
                PotGroLeafWeightMainStem[k][l] = 0;
                PotGroPetioleWeightMainStem[k][l] = 0;
            } else {
                int lp1 = l + 1;
                smax = VarPar[5] + VarPar[6] * lp1 * (VarPar[7] - lp1);
                smax = smax * denfac;
                if (smax < VarPar[4])
                    smax = VarPar[4];
                c = vpotlf[10] + lp1 * vpotlf[11];
                if (state.site[k][l][0].leaf.age > 70)
                    rate = 0;
                else
                    rate = smax * c * p * exp(-c * pow(state.site[k][l][0].leaf.age, p)) * pow(state.site[k][l][0].leaf.age, (p - 1));
//     Add leaf and petiole weight potential growth to SPDWL and SPDWP.
                if (rate >= 1e-12) {
                    PotGroLeafAreaMainStem[k][l] = rate * wstrlf * pixda
                                                   * TemperatureOnLeafGrowthRate(AvrgDailyTemp);
                    PotGroLeafWeightMainStem[k][l] = PotGroLeafAreaMainStem[k][l]
                                                     * LeafWeightAreaRatio;
                    PotGroPetioleWeightMainStem[k][l] = PotGroLeafAreaMainStem[k][l]
                                                        * LeafWeightAreaRatio * vpotlf[13];
                    PotGroAllLeaves += PotGroLeafWeightMainStem[k][l];
                    PotGroAllPetioles += PotGroPetioleWeightMainStem[k][l];
                } // rate
            } // LeafAreaMainStem
//     Assign smax value of this main stem leaf to smaxx, c to cc.
//     Loop over the nodes of this fruiting branch.
            double smaxx = smax; // value of smax for the corresponding main stem leaf.
            double cc = c;       // value of c for the corresponding main stem leaf.
            nnid = NumNodes[k][l];
            for (int m = 0; m < nnid; m++)// loop of nodes on a fruiting branch
            {
                if (LeafAreaNodes[k][l][m] <= 0) {
                    PotGroLeafAreaNodes[k][l][m] = 0;
                    PotGroLeafWeightNodes[k][l][m] = 0;
                    PotGroPetioleWeightNodes[k][l][m] = 0;
                }
//     Compute potential growth of leaf area and leaf weight for leaf
//  on fruiting branch node (k,l,m).
//    Add leaf and petiole weight potential growth to spdwl and spdwp.
                else {
                    int mp1 = m + 1;
//     smax and c are reduced as a function of node number on this fruiting branch.
                    smax = smaxx * (1 - VarPar[8] * mp1);
                    c = cc * (1 - VarPar[8] * mp1);
//     Compute potential growth for the leaves on fruiting branches.
                    if (state.site[k][l][m].leaf.age > 70)
                        rate = 0;
                    else
                        rate = smax * c * p * exp(-c * pow(state.site[k][l][m].leaf.age, p)) * pow(state.site[k][l][m].leaf.age, (p - 1));
                    if (rate >= 1e-12) {
//     Growth rate is modified by water stress and pix. Potential growth
//  is computed as a function of average temperature.
                        PotGroLeafAreaNodes[k][l][m] = rate * wstrlf * pixda
                                                       * TemperatureOnLeafGrowthRate(AvrgDailyTemp);
                        PotGroLeafWeightNodes[k][l][m] = PotGroLeafAreaNodes[k][l][m]
                                                         * LeafWeightAreaRatio;
                        PotGroPetioleWeightNodes[k][l][m] = PotGroLeafAreaNodes[k][l][m]
                                                            * LeafWeightAreaRatio * vpotlf[13];
                        PotGroAllLeaves += PotGroLeafWeightNodes[k][l][m];
                        PotGroAllPetioles += PotGroPetioleWeightNodes[k][l][m];
                    } //if rate
                } // if LeafAreaNodes
            } //loop nnid
        } //loop nbrch
    }  //loop NumVegBranches
}

////////////////////////////////////////////////////////////////////////
void PotentialFruitGrowth(State &state, const double &DayLength)
//     This function simulates the potential growth of fruiting sites of cotton plants.
//  It is called from PlantGrowth(). It calls TemperatureOnFruitGrowthRate()
//
//     The following gobal variables are referenced here:
//        AgeOfBoll, AgeOfSite, DayLength, FruitingCode, FruitFraction, NumFruitBranches, 
//        NumNodes,  NumVegBranches, DayTimeTemp, NightTimeTemp, VarPar, WaterStress.
//     The following global variables are set here:
//        PotGroAllBolls, PotGroAllBurrs, PotGroAllSquares, PotGroBolls, PotGroBurrs, PotGroSquares.
{
//     The constant parameters used:
    const double vpotfrt[5] = {0.72, 0.30, 3.875, 0.125, 0.17};
//    Compute tfrt for the effect of temperature on boll and burr growth rates. Function 
// TemperatureOnFruitGrowthRate() is used (with parameters derived from GOSSYM), for day time
// and night time temperatures, weighted by day and night lengths.
    double tfrt;  // the effect of temperature on rate of boll, burr or square growth.
    tfrt = (DayLength * TemperatureOnFruitGrowthRate(DayTimeTemp)
            + (24 - DayLength) * TemperatureOnFruitGrowthRate(NightTimeTemp)) / 24;
//     Assign zero to sums of potential growth of squares, bolls and burrs.
    PotGroAllSquares = 0;
    PotGroAllBolls = 0;
    PotGroAllBurrs = 0;
//     Assign values for the boll growth equation parameters. These are cultivar - specific.
    double agemax = VarPar[9]; // maximum boll growth period (physiological days).
    double rbmax = VarPar[10]; // maximum rate of boll (seed and lint) growth, g per boll per physiological day.
    double wbmax = VarPar[11]; // maximum possible boll (seed and lint) weight, g per boll.
//     Loop for all vegetative stems.
    for (int k = 0; k < NumVegBranches; k++) // loop of vegetative stems
    {
        int nbrch = NumFruitBranches[k]; // number of fruiting branches on a vegetative stem.
//     Loop for all fruiting branches on this vegetative stem.
        for (int l = 0; l < nbrch; l++) // loop of fruiting branches
        {
            int nnid = NumNodes[k][l]; // number of nodes on a fruiting branch.
            for (int m = 0; m < nnid; m++) // loop for nodes on a fruiting branch
            {
//     Calculate potential square growth for node (k,l,m).
//     Sum potential growth rates of squares as PotGroAllSquares.
                if (FruitingCode[k][l][m] == 1) {
//     ratesqr is the rate of square growth, g per square per day.
//  The routine for this is derived from GOSSYM, and so are the parameters used.
                    double ratesqr = tfrt * vpotfrt[3]
                                     * exp(-vpotfrt[2] + vpotfrt[3] * state.site[k][l][m].age);
                    PotGroSquares[k][l][m] = ratesqr * FruitFraction[k][l][m];
                    PotGroAllSquares += PotGroSquares[k][l][m];
                }
//     Growth of seedcotton is simulated separately from the growth of burrs.
//  The logistic function is used to simulate growth of seedcotton. The constants
//  of this function for cultivar 'Acala-SJ2', are based on the data of Marani (1979); 
//  they are derived from calibration for other cultivars
//     agemax is the age of the boll (in physiological days after
//  bloom) at the time when the boll growth rate is maximal.
//     rbmax is the potential maximum rate of boll growth (g seeds
//  plus lint dry weight per physiological day) at this age.
//     wbmax is the maximum potential weight of seed plus lint (g dry
//  weight per boll).
//     The auxiliary variable pex is computed as
//            pex = exp(-4 * rbmax * (t - agemax) / wbmax)
//  where t is the physiological age of the boll after bloom (= agebol).
//     Boll weight (seed plus lint) at age T, according to the
//  logistic function is:
//            wbol = wbmax / (1 + pex)
//  and the potential boll growth rate at this age will be the
//  derivative of this function:
//            ratebol = 4 * rbmax * pex / (1. + pex)**2
                else if (FruitingCode[k][l][m] == 2 || FruitingCode[k][l][m] == 7) {
//     pex is an intermediate variable to compute boll growth.
                    double pex = exp(-4 * rbmax * (AgeOfBoll[k][l][m] - agemax) / wbmax);
//  ratebol is the rate of boll (seed and lint) growth, g per boll per day.
                    double ratebol = 4 * rbmax * pex / pow((1 + pex), 2);
                    ratebol = ratebol * tfrt;
//     Potential growth rate of the burrs is assumed to be constant (vpotfrt[4] g dry weight 
//  per day) until the boll reaches its final volume. This occurs at the age of 22
//  physiological days in 'Acala-SJ2'. Both ratebol and ratebur are modified by 
//  temperature (tfrt) and ratebur is also affected by water stress (wfdb).
//     Compute wfdb for the effect of water stress on burr growth rate.
//  wfdb is the effect of water stress on rate of burr growth.
                    double wfdb = vpotfrt[0] + vpotfrt[1] * state.water_stress;
                    if (wfdb < 0)
                        wfdb = 0;
                    if (wfdb > 1)
                        wfdb = 1;
                    double ratebur; // rate of burr growth, g per boll per day.
                    if (AgeOfBoll[k][l][m] >= 22)
                        ratebur = 0;
                    else
                        ratebur = vpotfrt[4] * tfrt * wfdb;
//     Potential boll (seeds and lint) growth rate (ratebol) and
//  potential burr growth rate (ratebur) are multiplied by FruitFraction to
//  compute PotGroBolls and PotGroBurrs for node (k,l,m).
                    PotGroBolls[k][l][m] = ratebol * FruitFraction[k][l][m];
                    PotGroBurrs[k][l][m] = ratebur * FruitFraction[k][l][m];
//     Sum potential growth rates of bolls and burrs as PotGroAllBolls and
//  PotGroAllBurrs, respectively.
                    PotGroAllBolls += PotGroBolls[k][l][m];
                    PotGroAllBurrs += PotGroBurrs[k][l][m];
                }
//     If these are not green bolls, their potential growth is 0. End loop.
                else {
                    PotGroBolls[k][l][m] = 0;
                    PotGroBurrs[k][l][m] = 0;
                } // if FruitingCode
            } // for m
        } // for l
    } // for k
}
/*     References:
    Marani, A. 1979. Growth rate of cotton bolls and their
 components. Field Crops Res. 2:169-175.
    Marani, A., Phene, C.J. and Cardon, G.E. 1992. CALGOS, a
 version of GOSSYM adapted for irrigated cotton.  III. leaf and boll
 growth routines. Beltwide Cotton Grow, Res. Conf. 1992:1361-1363.  
*/
