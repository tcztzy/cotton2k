//  RootGrowth_1.cpp
//
//   functions in this file:
// PotentialRootGrowth()
// RootImpedance()
// SoilTemOnRootGrowth()
// SoilMechanicResistance()
// SoilAirOnRootGrowth()
// SoilNitrateOnRootGrowth()
// SoilWaterOnRootGrowth()
// ComputeActualRootGrowth()
//
#include "global.h"

using namespace std;

void RootImpedance(const int &);

extern "C" {
    double SoilTemOnRootGrowth(double);
    double SoilAirOnRootGrowth(double, double, double);
    double SoilNitrateOnRootGrowth(double);
    double SoilWaterOnRootGrowth(double);
    double SoilMechanicResistance(double);
}

tuple<int> RedistRootNewGrowth(int, int, double, int, double[40][20]);

tuple<int> TapRootGrowth(const int &, int, double[40][20][3], double[40][20]);

void InitiateLateralRoots();

void LateralRootGrowthLeft(int, const int &, double[40][20][3], double[40][20]);

void LateralRootGrowthRight(int, const int &, double[40][20][3], double[40][20]);

void RootAging(int, int, double[40][20][3], double[40][20]);

double RootDeath(int, int, double, double[40][20][3], const double[40][20]);

double RootCultivation(int, const int &, double, double[40][20][3]);

void RootSummation(const string &, const int &, const int &, const int &, double[40][20][3]);

//////////////////////////////////////////////////
//                   THE COTTON ROOT SUB-MODEL.
//     The following is a documentation of the root sub-model used in COTTON2K. It is 
//  derived from the principles of RHIZOS, as implemented in GOSSYM and in GLYCIM, 
//  and from some principles of ROOTSIMU (Hoogenboom and Huck, 1986). It is devised
//  to be generally applicable, and may be used with root systems of different crops by
//  redefining the parameters, which are set here as constants, and some of them are
//  set in function InitializeRootData(). These parameters are of course specific for the
//  crop species, and perhaps also for cultivars or cultivar groups.
//
//     This is a two-dimensional model and it may be used with soil cells of different 
//  sizes. The grid can be defined by the modeler. The maximum numbers of layers 
//  and columns are given by the parameters maxl and maxk, respectively. These are set 
//  to 40 and 20, in this version of COTTON2K. The grid is set in function InitializeGrid().
//
//     The whole slab is being simulated. Thus, non-symmetrical processes (such as 
//  side-dressing of fertilizers or drip-irrigation) can be handled. The plant is assumed to 
//  be situated at the center of the soil slab, or off-center for skip-rows. Adjoining soil
//  slabs are considered as mirror-images of each other. Alternate-row drip systems (or any 
//  other agricultural input similarly situated) are located at one edge of the slab.
//
//     The root mass in each cell is made up of NumRootAgeGroups classes, whose number is to be 
//  defined by the modeler. The maximum number of classes is 3 in this version of COTTON2K.
//
//     The following functions account for root growth morphology: TapRootGrowth() describes  
//  growth of the taproot, and LateralRootGrowth() describes growth of the lateral roots.
//
//     The calling sequence of the root submodel modules is as follows:
//     InitializeRootData() and ReadSoilImpedance() are called from ReadInput() 
//  at the start of the simulation (see their code in file gettinginput_2.cpp) .
//     PotentialRootGrowth() and ActualRootGrowth() are called each day from PlantGrowth().
//     PotentialRootGrowth() calls RootImpedance(), SoilMechanicResistance(), SoilAirOnRootGrowth(), 
//  SoilNitrateOnRootGrowth(), SoilTemOnRootGrowth(), SoilWaterOnRootGrowth().
//     ActualRootGrowth() calls RedistRootNewGrowth(), TapRootGrowth(), LateralRootGrowth(),
//  RootAging(), RootDeath(), RootCultivation(), RootSummation().
///////////////////////////////////////////////////////////////////////////////////
double PotentialRootGrowth(const int &NumRootAgeGroups, const int &NumLayersWithRoots, const int &ncurve,
                           const double RootWeight[40][20][3], const double RootAge[40][20])
//     This function calculates the potential root growth rate.  The return value 
//  is the sum of potential root growth rates for the whole slab (sumpdr).It is called from 
//  PlantGrowth(). It calls: RootImpedance(), SoilNitrateOnRootGrowth(), SoilAirOnRootGrowth(), 
//  SoilMechanicResistance(), SoilTemOnRootGrowth() and SoilWaterOnRootGrowth().
//
//     The following global variables are referenced here:
//       cgind, NumRootAgeGroups, nk, OutIndex, 
//       PerPlantArea, PoreSpace, RootAge, RootWeight. SoilPsi, SoilTempDailyAvrg, 
//       VolNo3NContent, VolWaterContent.
//     The following global variables are set here:    PotGroRoots, RootGroFactor
{
//     The following constant parameter is used:
    const double rgfac = 0.36; // potential relative growth rate of the roots (g/g/day).
//     Initialize to zero the PotGroRoots array.
    for (int l = 0; l < NumLayersWithRoots; l++)
        for (int k = 0; k < nk; k++)
            PotGroRoots[l][k] = 0;
    RootImpedance(ncurve);
    double sumpdr = 0;  // sum of potential root growth rate for the whole slab
    for (int l = 0; l < NumLayersWithRoots; l++)
        for (int k = 0; k < nk; k++) {
//     Check if this soil cell contains roots (if RootAge is greater
//  than 0), and execute the following if this is true.
//     In each soil cell with roots, the root weight capable of growth rtwtcg is computed as
//  the sum of RootWeight[l][k][i] * cgind[i] for all root classes.
            if (RootAge[l][k] > 0) {
                double rtwtcg = 0; // root weight capable of growth in a soil soil cell.
                for (int i = 0; i < NumRootAgeGroups; i++)
                    rtwtcg += RootWeight[l][k][i] * cgind[i];
//     Compute the temperature factor for root growth by calling function
//  SoilTemOnRootGrowth() for this layer.
                double stday; // soil temperature, C, this day's average for this cell.
                stday = SoilTempDailyAvrg[l][k] - 273.161;
                double temprg; // effect of soil temperature on root growth.
                temprg = SoilTemOnRootGrowth(stday);
//     Compute soil mechanical resistance for each soil cell by
//  calling SoilMechanicResistance{}, the effect of soil aeration on root growth by 
//  calling SoilAirOnRootGrowth(), and the effect of soil nitrate on root growth 
//  by calling SoilNitrateOnRootGrowth().
                double rtpct; // effect of soil mechanical resistance on root growth (returned from SoilMechanicResistance).
                //
                int lp1;      // layer below l.
                if (l == nl - 1) // last layer
                    lp1 = l;
                else
                    lp1 = l + 1;
            //
                int km1, kp1; // columns to the left and to the right of k.
                if (k == nk - 1) // last column
                    kp1 = k;
                else
                    kp1 = k + 1;
            //
                if (k == 0) // first column
                    km1 = 0;
                else
                    km1 = k - 1;
            //
                double rtimpd0 = RootImpede[l][k];
                double rtimpdkm1 = RootImpede[l][km1];
                double rtimpdkp1 = RootImpede[l][kp1];
                double rtimpdlp1 = RootImpede[lp1][k];
                double rtimpdmin = rtimpd0;      // minimum value of rtimpd
                if (rtimpdkm1 < rtimpdmin)
                    rtimpdmin = rtimpdkm1;
                if (rtimpdkp1 < rtimpdmin)
                    rtimpdmin = rtimpdkp1;
                if (rtimpdlp1 < rtimpdmin)
                    rtimpdmin = rtimpdlp1;
                rtpct = SoilMechanicResistance(rtimpdmin);
                double rtrdo; // effect of oxygen deficiency on root growth (returned from SoilAirOnRootGrowth).
                rtrdo = SoilAirOnRootGrowth(SoilPsi[l][k], PoreSpace[l], VolWaterContent[l][k]);
                double rtrdn; // effect of nitrate deficiency on root growth (returned from SoilNitrateOnRootGrowth).
                rtrdn = SoilNitrateOnRootGrowth(VolNo3NContent[l][k]);
//     The root growth resistance factor RootGroFactor(l,k), which can take a
//  value between 0 and 1, is computed as the minimum of these
//  resistance factors. It is further modified by multiplying it by 
//  the soil moisture function SoilWaterOnRootGrowth().
//     Potential root growth PotGroRoots(l,k) in each cell is computed as a
//  product of rtwtcg, rgfac, the temperature function temprg, and
//  RootGroFactor(l,k). It is also multiplied by PerPlantArea / 19.6, for the effect
//  of plant population density on root growth: it is made comparable
//  to a population of 5 plants per m in 38" rows.
//     The sum of the potential growth for the whole slab is computed
//  as sumpdr.
                double minres = min(rtpct, rtrdo);
                if (rtrdn < minres)
                    minres = rtrdn;
                double rtpsi = SoilWaterOnRootGrowth(SoilPsi[l][k]);
                RootGroFactor[l][k] = rtpsi * minres;
                PotGroRoots[l][k] = rtwtcg * rgfac * temprg * RootGroFactor[l][k] * PerPlantArea / 19.6;
                sumpdr += PotGroRoots[l][k];
            } // if RootAge
        } // end loop l & k
    return sumpdr;
}

//////////////////////////
void RootImpedance(const int &ncurve)
//     This function calculates soil mechanical impedance to root growth, rtimpd(l,k), 
//  for all soil cells. It is called from PotentialRootGrowth(). The impedance is a function  
//  of bulk density and water content in each soil soil cell. No changes have been made 
//  in the original GOSSYM code.
//
//     The following global variables are referenced here:
//  BulkDensity, gh2oc, impede, inrim, ncurve, nk, nl, 
//  SoilHorizonNum, tstbd, VolWaterContent.
//     The following global variables are set here:    RootImpede.
{
    for (int l = 0; l < nl; l++) {
        int j = SoilHorizonNum[l];
        double Bd = BulkDensity[j]; // bulk density for this layer
//
        int jj;
        for (jj = 0; jj < inrim; jj++) {
            if (Bd <= tstbd[jj][0])
                break;
        }
        int j1 = jj;
        if (j1 > inrim - 1)
            j1 = inrim - 1;
        int j0 = jj - 1;
        if (j0 < 0)
            j0 = 0;
//
        for (int k = 0; k < nk; k++) {
            double Vh2o = VolWaterContent[l][k] / Bd;
            int ik;
            for (ik = 0; ik < ncurve; ik++) {
                if (Vh2o <= gh2oc[ik])
                    break;
            }
            int i1 = ik;
            if (i1 > ncurve - 1)
                i1 = ncurve - 1;
            int i0 = ik - 1;
            if (i0 < 0)
                i0 = 0;
//
            if (j1 == 0) {
                if (i1 == 0 || Vh2o <= gh2oc[i1])
                    RootImpede[l][k] = impede[j1][i1];
                else
                    RootImpede[l][k] = impede[j1][i0] - (impede[j1][i0] - impede[j1][i1])
                                                        * (Vh2o - gh2oc[i0]) / (gh2oc[i1] - gh2oc[i0]);
            } else {
                if (i1 == 0 || Vh2o <= gh2oc[i1])
                    RootImpede[l][k] = impede[j0][i1] - (impede[j0][i1] - impede[j1][i1])
                                                        * (tstbd[j0][i1] - Bd) / (tstbd[j0][i1] - tstbd[j1][i1]);
                else {
                    double temp1 = impede[j0][i1] - (impede[j0][i1] - impede[j1][i1])
                                                    * (tstbd[j0][i1] - Bd) / (tstbd[j0][i1] - tstbd[j1][i1]);
                    double temp2 = impede[j0][i0] - (impede[j0][i0] - impede[j1][i1])
                                                    * (tstbd[j0][i0] - Bd) / (tstbd[j0][i0] - tstbd[j1][i0]);
                    RootImpede[l][k] = temp2 + (temp1 - temp2)
                                               * (Vh2o - gh2oc[i0]) / (gh2oc[i1] - gh2oc[i0]);
                }
            }
        }
    }
//
}

//////////////////////////
tuple<int>
ComputeActualRootGrowth(double sumpdr, const string &ProfileName, const int &Daynum, const int &DayOfSimulation,
                        const int &DayEmerge, int NumLayersWithRoots, const int &NumRootAgeGroups,
                        double RootWeight[40][20][3], double RootAge[40][20])
//     This function calculates the actual root growth rate. It is called from function 
//  PlantGrowth(). It calls the following functions:  InitiateLateralRoots(), 
//  LateralRootGrowthLeft(), LateralRootGrowthRight(), RedistRootNewGrowth(), RootAging(), 
//	RootCultivation(), RootDeath(), RootSummation(), TapRootGrowth().
//
//     The following global variables are referenced here:
//       CarbonAllocatedForRootGrowth, cgind, CultivationDate, Daynum, DepthLastRootLayer, 
//       dl, ExtraCarbon, LateralRootFlag, LastTaprootLayer, nk, nl, NumLayersWithRoots,
//       NumRootAgeGroups, PerPlantArea, pixcon, PlantRowColumn, PotGroRoots, RootAge, 
//       RootNConc, RowSpace, TapRootLength, wk.
//     The following global variables are set here:
//       ActualRootGrowth, CumPlantNLoss, DailyRootLoss, PixInPlants, 
//       RootNitrogen, RootWeight, RootWeightLoss.
//     The following argument is used:
//       sumpdr - Sum of potential root growth rate for the whole slab.
{
//     The following constant parameters are used:
//     The index for the relative partitioning of root mass produced by new growth to class i.
    const double RootGrowthIndex[3] = {1.0, 0.0, 0.0};
    const double rtminc = 0.0000001; // the threshold ratio of root mass capable of growth
    // to soil cell volume (g/cm3); when this threshold is reached, a
    // part of root growth in this cell may be extended to adjoining cells.
//
    static double pavail; // residual available carbon for root growth from previous day.
//     Assign zero to pavail if this is the day of emergence.
    if (Daynum <= DayEmerge)
        pavail = 0;
    double adwr1[maxl][maxk]; // actual growth rate from roots existing in this soil cell.
//     Assign zero to the arrays of actual root growth rate.
    for (int l = 0; l < nl; l++)
        for (int k = 0; k < nk; k++) {
            ActualRootGrowth[l][k] = 0;
            adwr1[l][k] = 0;
        }
//     The amount of carbon allocated for root growth is calculated from  
//  CarbonAllocatedForRootGrowth, converted to g dry matter per slab, and added to previously
//  allocated carbon that has not been used for growth (pavail). if there is no potential root
//  growth, this will be stored in pavail. Otherwise, zero is assigned to pavail.
    if (sumpdr <= 0) {
        pavail += CarbonAllocatedForRootGrowth * 0.01 * RowSpace / PerPlantArea;
        return make_tuple(NumLayersWithRoots);
    }
    double actgf; // actual growth factor (ratio of available C to potential growth).
//     The ratio of available C to potential root growth (actgf) is calculated. 
//  pavail (if not zero) is used here, and zeroed after being used.
    actgf = (pavail + CarbonAllocatedForRootGrowth * 0.01 * RowSpace / PerPlantArea)
            / sumpdr;
    pavail = 0;
//
    for (int l = 0; l < NumLayersWithRoots; l++)
        for (int k = 0; k < nk; k++) {
//     adwr1(l,k), is proportional to the potential growth rate of roots in this cell.
            if (RootAge[l][k] > 0)
                adwr1[l][k] = PotGroRoots[l][k] * actgf;
        }
//     If extra carbon is available, it is assumed to be added to the taproot. 
    if (ExtraCarbon > 0) {
        double availt; // available carbon for taproot growth, in g dry matter per slab.
//  ExtraCarbon is converted to availt (g dry matter per slab).
        availt = ExtraCarbon * 0.01 * RowSpace / PerPlantArea;
        double sdl = TapRootLength - DepthLastRootLayer;; // distance from the tip of the taproot, cm.
        double tpwt[maxl][2]; // proportionality factors for allocating added dry matter among taproot soil cells.
        double sumwt = 0; // sum of the tpwt array.
//     Extra Carbon (availt) is added to soil cells with roots in the columns immediately to the
//  left and to the right of the location of the plant row.
        for (int l = LastTaprootLayer; l >= 0; l--) {
//     The weighting factors for allocating the carbon (tpwt) are proportional to the volume 
//  of each soil cell and its distance (sdl) from the tip of the taproot. 
            sdl += dl[l];
            tpwt[l][0] = sdl * dl[l] * wk[PlantRowColumn];
            tpwt[l][1] = sdl * dl[l] * wk[PlantRowColumn + 1];
            sumwt += tpwt[l][0] + tpwt[l][1];
        }
//     The proportional amount of mass is added to the mass of the last (inactive) 
//  root class in each soil cell.
        for (int l = 0; l <= LastTaprootLayer; l++) {
            RootWeight[l][PlantRowColumn][NumRootAgeGroups - 1] += availt * tpwt[l][0] / sumwt;
            RootWeight[l][PlantRowColumn + 1][NumRootAgeGroups - 1] += availt * tpwt[l][1] / sumwt;
        }
    }
//     Check each cell if the ratio of root weight capable of growth to cell volume (rtconc) 
//  exceeds the threshold rtminc, and call RedistRootNewGrowth() for this cell. Otherwise, all 
//  new growth is contained in the same cell, and the actual growth in this cell, 
//  ActualRootGrowth(l,k) will be equal to adwr1(l,k).
    for (int l = 0; l < nl; l++)
        for (int k = 0; k < nk; k++)
            if (RootAge[l][k] > 0) {
                double rtconc = 0; // ratio of root weight capable of growth to cell volume.
                for (int i = 0; i < NumRootAgeGroups; i++)
                    rtconc += RootWeight[l][k][i] * cgind[i];
                rtconc = rtconc / (dl[l] * wk[k]);
                if (rtconc > rtminc)
                    tie(NumLayersWithRoots) = RedistRootNewGrowth(l, k, adwr1[l][k], NumLayersWithRoots, RootAge);
                else
                    ActualRootGrowth[l][k] += adwr1[l][k];
            }
//     The new actual growth ActualRootGrowth(l,k) in each cell is partitioned among the root 
//  classes in it in proportion to the parameters RootGrowthIndex(i), and the previous values of 
//  RootWeight(k,l,i), and added to RootWeight(k,l,i).
    double sumind = 0; // sum of the growth index for all classes in a cell.
    for (int i = 0; i < NumRootAgeGroups; i++)
        sumind += RootGrowthIndex[i];
//
    for (int l = 0; l < NumLayersWithRoots; l++)
        for (int k = 0; k < nk; k++) {
            if (RootAge[l][k] > 0) {
                double sumgr = 0; // sum of growth index multiplied by root weight, for all classes in a cell.
                for (int i = 0; i < NumRootAgeGroups; i++)
                    sumgr += RootGrowthIndex[i] * RootWeight[l][k][i];
                for (int i = 0; i < NumRootAgeGroups; i++) {
                    if (sumgr > 0)
                        RootWeight[l][k][i] +=
                                ActualRootGrowth[l][k] * RootGrowthIndex[i] * RootWeight[l][k][i] / sumgr;
                    else
                        RootWeight[l][k][i] += ActualRootGrowth[l][k] * RootGrowthIndex[i] / sumind;
                }
            }
        }
//     Call function TapRootGrowth() for taproot elongation, if the taproot
//  has not already reached the bottom of the slab.
    if (LastTaprootLayer < nl - 1 || TapRootLength < DepthLastRootLayer)
        tie(NumLayersWithRoots) = TapRootGrowth(NumRootAgeGroups, NumLayersWithRoots, RootWeight, RootAge);
//     Call functions for growth of lateral roots
    InitiateLateralRoots();
    for (int l = 0; l < LastTaprootLayer; l++) {
        if (LateralRootFlag[l] == 2) {
            LateralRootGrowthLeft(l, NumRootAgeGroups, RootWeight, RootAge);
            LateralRootGrowthRight(l, NumRootAgeGroups, RootWeight, RootAge);
        }
    }
//     Initialize DailyRootLoss (weight of sloughed roots) for this day.
    double DailyRootLoss = 0; // total weight of sloughed roots, g per plant per day.
    for (int l = 0; l < NumLayersWithRoots; l++)
        for (int k = 0; k < nk; k++) {
//     Check RootAge to determine if this soil cell contains roots, and then compute root 
//  aging and root death by calling RootAging() and RootDeath() for each soil cell with roots.
            if (RootAge[l][k] > 0) {
                RootAging(l, k, RootWeight, RootAge);
                DailyRootLoss = RootDeath(l, k, DailyRootLoss, RootWeight, RootAge);
            }
        }
//     Check if cultivation is executed in this day and call RootCultivation().
    for (int j = 0; j < 5; j++)
        if (CultivationDate[j] == Daynum)
            DailyRootLoss = RootCultivation(j, NumRootAgeGroups, DailyRootLoss, RootWeight);
//     Convert DailyRootLoss to g per plant units and add it to RootWeightLoss.
    DailyRootLoss = DailyRootLoss * 100. * PerPlantArea / RowSpace;
    RootWeightLoss += DailyRootLoss;
//     Adjust RootNitrogen (root N content) and PixInPlants (plant Pix content)
//  for loss by death of roots.
    RootNitrogen -= DailyRootLoss * RootNConc;
    CumPlantNLoss += DailyRootLoss * RootNConc;
    PixInPlants -= DailyRootLoss * pixcon;
//     Call function RootSummation().
    RootSummation(ProfileName, DayOfSimulation, NumRootAgeGroups, NumLayersWithRoots, RootWeight);
    return make_tuple(NumLayersWithRoots);
}
