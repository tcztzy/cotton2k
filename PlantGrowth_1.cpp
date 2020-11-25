//  PlantGrowth_1.cpp
//
//   functions in this file:
// PhysiologicalAge()
// Stress()
// LeafWaterPotential()
// LeafResistance()
// GetNetPhotosynthesis()
// PlantGrowth()
//
#include "CottonSimulation.h"
#include "GeneralFunctions.h"
#include "resource.h"
#include <math.h>
//
#ifdef _DEBUG
#define new DEBUG_NEW
#endif
//
double ptsred; // The effect of moisture stress on the photosynthetic rate
//////////////////////////////////////////////////
double PhysiologicalAge()  // computes physiological age
//     This function returns the daily 'physiological age' increment, 
//  based on hourly temperatures. It is called each day by SimulateThisDay().
//     The following global variable is used here:    
//        AirTemp[] = array of hourly temperatures.
{
//     The following constant Parameters are used in this function:
      const double p1 = 12; // threshold temperature, C
	  const double p2 = 14; // temperature, C, above p1, for one physiological day.
	  const double p3 = 1.5; // maximum value of a physiological day.
//     The threshold value is assumed to be 12 C (p1). One physiological day is 
//  equivalent to a day with an average temperature of 26 C, and therefore the 
//  heat units are divided by 14 (p2).
//     A linear relationship is assumed between temperature and heat unit
//  accumulation in the range of 12 C (p1) to 33 C (p2*p3+p1). the effect of
//  temperatures higher than 33 C is assumed to be equivalent to that of 33 C.
      double dayfd = 0; // the daily contribution to physiological age (return value).
      for (int ihr = 0; ihr < 24; ihr++)
	  {
         double tfd = (AirTemp[ihr] - p1) / p2; // the hourly contribution to physiological age.
         if (tfd < 0)
			 tfd = 0;
         if (tfd > p3)
			 tfd = p3;
         dayfd += tfd;
      }
      return dayfd / 24;
}
/////////////////////////////////////////////////////////////////////////
void Stress()                      
//     This function computes the water stress variables affecting
// the cotton plants. It is called by SimulateThisDay() and calls LeafWaterPotential().
//
//     The following global variables are referenced here:
//        Kday, LwpMin, LwpMax.
//     The following global variables are set here:
//        AverageLwp, AverageLwpMin, LwpMinX, LwpX, ptsred, WaterStress, WaterStressStem.
//
{
//     The following constant parameters are used:
	  const double vstrs[9] = { -3.0, 3.229, 1.907, 0.321, -0.10, 1.230, 0.340, 0.30, 0.05 };
//     Call LeafWaterPotential() to compute leaf water potentials.
      LeafWaterPotential();
//     The running averages, for the last three days, are computed:
//  AverageLwpMin is the average of LwpMin, and AverageLwp of LwpMin + LwpMax.
      AverageLwpMin += (LwpMin - LwpMinX[2]) / 3;
      AverageLwp += (LwpMin + LwpMax - LwpX[2]) / 3;
      for (int i = 2; i > 0; i--)
	  {
         LwpMinX[i] = LwpMinX[i-1];
         LwpX[i] = LwpX[i-1];
      }
      LwpMinX[0] = LwpMin;
      LwpX[0] = LwpMin + LwpMax;
//     No stress effects before 5th day after emergence.
      if ( Kday < 5 ) 
	  {
         ptsred = 1;
         WaterStress = 1;
         WaterStressStem = 1;
         return;
      }
//     The computation of ptsred, the effect of moisture stress on
//  the photosynthetic rate, is based on the following work: Ephrath,
//  J.E., Marani, A., Bravdo, B.A., 1990. Effects of moisture stress on
//  stomatal resistance and photosynthetic rate in cotton (Gossypium
//  hirsutum) 1. Controlled levels of stress. Field Crops Res.23:117-131.
//  It is a function of AverageLwpMin (average LwpMin for the last three days).
      if ( AverageLwpMin < vstrs[0] ) 
		   AverageLwpMin = vstrs[0];
      ptsred = vstrs[1] + AverageLwpMin * (vstrs[2] + vstrs[3] * AverageLwpMin);
      if ( ptsred > 1 )
		   ptsred = 1;
//     The general moisture stress factor (WaterStress) is computed as an
// empirical function of AverageLwp. psilim, the value of AverageLwp at the
// maximum value of the function, is used for truncating it.
//     The minimum value of WaterStress is 0.05, and the maximum is 1.
      double psilim; // limiting value of AverageLwp.
      psilim = -0.5 * vstrs[5] / vstrs[6];
      if ( AverageLwp > psilim ) 
         WaterStress = 1;
      else
	  {
         WaterStress = vstrs[4] - AverageLwp * (vstrs[5] + vstrs[6] * AverageLwp);
         if ( WaterStress > 1 )
			  WaterStress = 1;
         if ( WaterStress < 0.05 )
			  WaterStress = 0.05;
      }
//     Water stress affecting plant height and stem growth (WaterStressStem) is assumed
//  to be more severe than WaterStress, especially at low WaterStress values.
      WaterStressStem =  WaterStress * (1 + vstrs[7] * (2 - WaterStress)) - vstrs[7];
      if ( WaterStressStem < vstrs[8] )
		   WaterStressStem = vstrs[8];
}
//////////////////////////
void LeafWaterPotential()
//     This function simulates the leaf water potential of cotton plants. It has been 
//  adapted from the model of Moshe Meron (The relation of cotton leaf water potential to 
//  soil water content in the irrigated management range. PhD dissertation, UC Davis, 1984).
//     It is called from Stress(). It calls wcond() and LeafResistance().
//
//     The following global variables are referenced here:
//       AgeOfPreFruNode, AverageSoilPsi, beta, dl, Kday, LeafAge, NumFruitBranches, 
//       NumLayersWithRoots, NumNodes, NumPreFruNodes, NumVegBranches, OutIndex, 
//       pi, PlantHeight, PoreSpace, ReferenceETP, RootColNumLeft, RootColNumRight, 
//       RootWtCapblUptake, SaturatedHydCond, SoilPsi, thad, thts, VolWaterContent, wk.
//     The following global variables are set here:
//       LwpMin, LwpMax.
//
{
//     Constant parameters used:
      const double cmg = 3200; // length in cm per g dry weight of roots, based on an average
//               root diameter of 0.06 cm, and a specific weight of 0.11 g dw per cubic cm.
      const double psild0 = -1.32; // maximum values of LwpMin 
	  const double psiln0 = -0.40; // maximum values of LwpMax. 
      const double rtdiam = 0.06; // average root diameter in cm.
      const double vpsil[] = { 0.48, -5.0, 27000., 4000., 9200., 920., 0.000012,
                              -0.15, -1.70, -3.5, 0.1e-9, 0.025, 0.80 };
//     Leaf water potential is not computed during 10 days after
//  emergence. Constant values are assumed for this period.
      if (Kday <= 10) 
	  {
         LwpMax = psiln0;
         LwpMin = psild0;
         return;
      }
//     Compute shoot resistance (rshoot) as a function of plant height.
      double rshoot; // shoot resistance, Mpa hours per cm.
      rshoot = vpsil[0] * PlantHeight / 100;
//     Assign zero to summation variables
      double psinum = 0;  // sum of RootWtCapblUptake for all soil cells with roots.
      double rootvol = 0; // sum of volume of all soil cells with roots.
      double rrlsum = 0;  // weighted sum of reciprocals of rrl.
      double rroot = 0;   // root resistance, Mpa hours per cm.
      double sumlv = 0;   // weighted sum of root length, cm, for all soil cells with roots.
      double vh2sum = 0;  // weighted sum of soil water content, for all soil cells with roots.
//     Loop over all soil cells with roots. Check if RootWtCapblUptake is
//  greater than vpsil[10].
//     All average values computed for the root zone, are weighted by RootWtCapblUptake 
//  (root weight capable of uptake), but the weight assigned will not be greater than vpsil[11].
      double rrl; // root resistance per g of active roots.
      for (int l = 0; l < NumLayersWithRoots; l++)
         for (int k = RootColNumLeft[l]; k < RootColNumRight[l]; k++)
		 {
            if ( RootWtCapblUptake[l][k] >= vpsil[10] ) 
			{
               psinum += min(RootWtCapblUptake[l][k],vpsil[11]);
               sumlv += min(RootWtCapblUptake[l][k],vpsil[11]) * cmg;
               rootvol += dl[l] * wk[k];
               if ( SoilPsi[l][k] <= vpsil[1] ) 
                  rrl = vpsil[2] / cmg;
               else
                  rrl = (vpsil[3] - SoilPsi[l][k] * (vpsil[4] + vpsil[5] * SoilPsi[l][k])) / cmg;
               rrlsum += min(RootWtCapblUptake[l][k],vpsil[11]) / rrl;
               vh2sum += VolWaterContent[l][k] * min(RootWtCapblUptake[l][k],vpsil[11]);
            }
         }
//     Compute average root resistance (rroot) and average soil water content (vh2).
      double dumyrs; // intermediate variable for computing cond.
      double vh2; // average of soil water content, for all soil soil cells with roots.
      if ( psinum > 0 && sumlv > 0 ) 
	  {
         rroot = psinum / rrlsum;
         vh2 = vh2sum / psinum;
         dumyrs = sqrt(1 / (pi * sumlv / rootvol)) / rtdiam;
         if ( dumyrs < 1.001 )
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
      cond = wcond(vh2,thad[0],thts[0],beta[0],SaturatedHydCond[0],PoreSpace[0]) / 24;
      cond = cond * 2 * sumlv / rootvol / log(dumyrs);
      if ( cond < vpsil[6] )
		   cond = vpsil[6];
      double rsoil =  0.0001 / ( 2 * pi * cond ); // soil resistance, Mpa hours per cm.
//     Compute leaf resistance (LeafResistance) as the average of the resistances
//  of all existing leaves. The resistance of an individual leaf is a
//  function of its age. Function LeafResistance is called to compute it. This
//  is executed for all the leaves of the plant.
      int numl = 0; // number of leaves.
      double sumrl = 0; // sum of leaf resistances for all the plant.
      for (int j = 0; j < NumPreFruNodes; j++) // loop prefruiting nodes
	  {
         numl++;
         sumrl += LeafResistance( AgeOfPreFruNode[j] );
      }
//
      int nbrch; // number of fruiting branches on a vegetative branch.
      int nnid; // number of nodes on a fruiting branch.
      for (int k = 0; k < NumVegBranches; k++)  // loop for all other nodes
	  {
         nbrch = NumFruitBranches[k];
         for (int l = 0; l < nbrch; l++)
		 {
            nnid = NumNodes[k][l];
            for (int m = 0; m < nnid; m++)
			{
               numl++;
               sumrl += LeafResistance( LeafAge[k][l][m] );
			}
		 }
	  }
      double rleaf = sumrl / numl; // leaf resistance, Mpa hours per cm.
//     The total resistance to transpiration, MPa hours per cm, (rtotal) is computed.
      double rtotal = rsoil + rroot + rshoot + rleaf; 
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
      double etmax = 0; // the maximum hourly rate of evapotranspiration for this day.
      for (int ihr = 0; ihr < 24; ihr++)  //  hourly loop
	  {
         if (ReferenceETP[ihr] > etmax) 
             etmax = ReferenceETP[ihr];
      }
      LwpMin = LwpMax - 0.1 * max( etmax, vpsil[12] ) * rtotal;
//     Check for minimum and maximum values.
      if ( LwpMin < vpsil[9] )
		   LwpMin = vpsil[9];
      if ( LwpMin > psild0 )
		   LwpMin = psild0;
//     If the output flag is non-zero, write daily values to file 44 (*.LWP)
      if (OutIndex[19] > 0) 
	  {
         ofstream File44("Output\\" + ProfileName + ".LWP", ios::app);
		 File44.unsetf(ios::left);
         File44.width(5);
		 File44 << Kday;
		 File44 << "   ";
		 File44.setf(ios::fixed);
	     File44.precision(2);
         File44.width(9);
		 File44 << AverageSoilPsi;
		 File44 << " ";
         File44.width(9);
	     File44.precision(3);
		 File44.setf(ios::scientific);
		 File44 << cond;
		 File44.setf(ios::fixed);
		 File44 << " ";
	     File44.precision(2);
         File44.width(9);
		 File44 << rsoil;
		 File44 << " ";
         File44.width(9);
		 File44 << rroot;
		 File44 << " ";
	     File44.precision(3);
         File44.width(9);
		 File44 << rshoot;
		 File44 << " ";
         File44.width(9);
		 File44 << rleaf;
		 File44 << " ";
         File44.width(9);
		 File44 << LwpMax;
		 File44 << " ";
         File44.width(9);
		 File44 << LwpMin;
		 File44 << endl;
	  }
}
//////////////////////////////////
double LeafResistance( double agel )
//     This function computes and returns the resistance of leaves of cotton
// plants to transpiration. It is assumed to be a function of leaf age. 
// It is called from LeafWaterPotential().
//     The input argument (agel) is leaf age in physiological days.
{
//     The following constant parameters are used:
      const double afac = 160;  // factor used for computing leaf resistance. 
      const double agehi = 94;  // higher limit for leaf age.
      const double agelo = 48;  // lower limit for leaf age.
      const double rlmin = 0.5; // minimum leaf resistance.
//
      double leafResistance;
      if ( agel <= agelo ) 
         leafResistance = rlmin;
      else if ( agel >= agehi ) 
         leafResistance = rlmin + (agehi - agelo)*(agehi - agelo) / afac;
      else
	  {
         double ax = 2 * agehi - agelo; // intermediate variable
         leafResistance = rlmin + (agel - agelo) * (ax - agel) / afac;
      }
      return leafResistance;
}
////////////////////////////////
void GetNetPhotosynthesis()            // computes net photosynthesis.
//     This function simulates the net photosynthesis of cotton  plants. It is called 
// daily by SimulateThisDay(). This is essentially the routine of GOSSYM with minor changes.
//     The following global and file scope variables are referenced here:
//       BurrWeightOpenBolls, CO2EnrichmentFactor, CottonWeightOpenBolls, 
//       DayLength, Daynum, DayEmerge, DayEndCO2, DayStartCO2, 
//       DayTimeTemp, iyear, Kday, LeafNConc, LightIntercept, PerPlantArea, 
//       PlantWeight, ptsred, StemWeight, TotalLeafWeight. 
//     The following global variables are set here:      
//       bEnd, CumNetPhotosynth, NetPhotosynthesis.
{
//  constants:
      const double gsubr = 0.375; // the growth resiration factor.
      const double rsubo = 0.0032; // maintenance respiration factor.
	  const double vpnet[4] = { 1.30, 0.034, 0.010, 0.32 }; 
      const double co2parm[45] =  // parameters used to correct photosynthesis for ambient CO2 concentration.
           {1.0235, 1.0264, 1.0285, 1.0321, 1.0335, 1.0353, 1.0385, 1.0403, 1.0431, 1.0485,
            1.0538, 1.0595, 1.0627, 1.0663, 1.0716, 1.0752, 1.0784, 1.0823, 1.0880, 1.0923,
            1.0968, 1.1019, 1.1087, 1.1172, 1.1208, 1.1243, 1.1311, 1.1379, 1.1435, 1.1490, 
            1.1545, 1.1601, 1.1656, 1.1712, 1.1767, 1.1823, 1.1878, 1.1934, 1.1990, 1.2045, 
            1.2101, 1.2156, 1.2212, 1.2267, 1.2323};
//     Note: co2parm is for icrease in ambient CO2 concentration changes from 1959 (308 ppm). 
//  The first 28 values (up to 1987) are from GOSSYM. The other values (up to 2004) 
//  are derived from data of the Carbon Dioxide Information Analysis Center (CDIAC).
//
//     Exit the function and end simulation if there are no leaves.
      if ( TotalLeafWeight <= 0 ) 
      {
          bEnd = true;
		  return;
      }
//     If this is the first time the function is executed, get the ambient CO2 correction.
      static double AmbientCO2Factor; // correction factor for ambient CO2 in air
      if ( Daynum <= DayEmerge )
	  {
         int co2indx = iyear - 1960; // count of years from 1960.
         if ( co2indx < 0 )
			  AmbientCO2Factor = 1;
         else if ( co2indx < 45 )  // for years 1960 to 2004
             AmbientCO2Factor = co2parm[co2indx];
         else  // extrapolate for years after 2004
             AmbientCO2Factor = 1.2323 + 0.004864 * (co2indx - 45);
      }
//     Get the CO2 correction factor (pnetcor) for photosynthesis,
//  using AmbientCO2Factor and a factor that may be variety specific (vpnet[0]).
//     CO2EnrichmentFactor is used for CO2 enrichment simulations, between DOY
//  dates DayStartCO2 and DayEndCO2.
      double pnetcor; // correction factor for gross photosynthesis.
      if ( Daynum >= DayStartCO2 && Daynum <= DayEndCO2 && CO2EnrichmentFactor > 1) 
         pnetcor = AmbientCO2Factor * vpnet[0] * CO2EnrichmentFactor;
      else
         pnetcor = AmbientCO2Factor * vpnet[0];
//     Compute ptnfac, the effect of leaf N concentration on
//  photosynthesis, using an empirical relationship.
      double ptnfac; // correction factor for low nitrogen content in leaves.
      ptnfac = vpnet[3] + (LeafNConc - vpnet[2]) *
              (1 - vpnet[3]) / (vpnet[1] - vpnet[2]);
      if ( ptnfac > 1 )
		   ptnfac = 1;
      if ( ptnfac < vpnet[3] )
		   ptnfac = vpnet[3];
//     Convert the average daily short wave radiation from langley per 
//  day, to Watts per square meter (wattsm).
      double wattsm; // average daily global radiation, W m-2.
      wattsm = GetFromClim("rad", Daynum) * 697.45 / (DayLength * 60);
//     Compute pstand as an empirical function of wattsm (based on Baker et al., 1972).
      double pstand; // gross photosynthesis for a non-stressed full canopy.
      pstand = 2.3908 + wattsm * (1.37379 - wattsm * 0.00054136);
//     Convert it to gross photosynthesis per plant (pplant), using
//  PerPlantArea and corrections for light interception by canopy, ambient CO2
//  concentration, water stress and low N in the leaves.
      double pplant; // actual gross photosynthetic rate, g per plant per day.
      pplant = 0.001 * pstand * LightIntercept * PerPlantArea * ptsred * pnetcor * ptnfac;
//     Compute the photorespiration factor (rsubl) as a linear
//  function af average day time temperature.
      double rsubl = 0.0032125 + 0.0066875 * DayTimeTemp; // photorespiration factor.
//     Photorespiration (lytres) is computed as a proportion of gross
//  photosynthetic rate.
      double lytres; // rate of photorespiration, g per plant per day.
      lytres = rsubl * pplant;
//     Old stems are those more than voldstm = 32 calendar days old.
//     Maintenance respiration is computed on the basis of plant dry
//  weight, minus the old stems and the dry tissue of opened bolls.
      double oldstmwt; // weight of old stems.
      int voldstm = 32; 
      int kkday = Kday - voldstm; // day of least recent actively growing stems.
      if ( kkday < 1 ) 
         oldstmwt = 0;
      else
         oldstmwt = StemWeight[kkday];
      double bmain; // maintenance respiration, g per plant per day.
      bmain = (PlantWeight - CottonWeightOpenBolls - BurrWeightOpenBolls - oldstmwt) * rsubo;
//     Net photosynthesis is computed by substracting photo-respiration and maintenance
//  respiration from the gross rate of photosynthesis. To avoid computational problems, 
//  make sure that pts is positive and non-zero.
      double pts; // intermediate computation of NetPhotosynthesis.
      pts = pplant - lytres - bmain;
      if (pts < 0.00001)
		  pts = 0.00001;
//     The growth respiration (gsubr) supplies energy for converting the supplied carbohydrates 
//  to plant tissue dry matter. 0.68182 converts CO2 to CH2O. NetPhotosynthesis is the computed  
//  net photosynthesis, in g per plant per day.
      NetPhotosynthesis = pts / (1 + gsubr) * 0.68182;
//     CumNetPhotosynth is the cumulative value of NetPhotosynthesis, from day of emergence.
      CumNetPhotosynth += NetPhotosynthesis;
}
/*     References:
     Baker et. al. (1972). Simulation of Growth and Yield in
  Cotton: I. Gross photosynthesis, respiration and growth. Crop Sci.
  12:431-435.
     Harper et. al. (1973) Carbon dioxide and the photosynthesis of
  field crops.  A metered carbon dioxide release in cotton under
  field conditions.  Agron. J. 65:7-11.
     Baker (1965)  Effects of certain environmental factors on net
  assimilation in cotton.  Crop Sci. 5:53-56 (Fig 5).  
*/
////////////////////////////////////////////////////////////////////////////
void PlantGrowth()                 
//     This function simulates the potential and actual growth of cotton plants. 
//  It is called from SimulateThisDay(), and it calls the following functions:
//    ActualFruitGrowth(), ActualLeafGrowth(), ActualRootGrowth(), AddPlantHeight(),
//    DryMatterBalance(), PotentialFruitGrowth(), PotentialLeafGrowth(), 
//    PotentialRootGrowth(), PotentialStemGrowth().
//
//     The following global variables are referenced here:
//        ActualStemGrowth, DayInc, FirstSquare, FruitingCode, Kday, pixdz, 
//        PerPlantArea, RowSpace, WaterStressStem.
//
//     The following global variables are set here:
//        LeafAreaIndex, PlantHeight, PotGroAllRoots, PotGroStem, StemWeight, 
//        TotalLeafArea, TotalLeafWeight, TotalPetioleWeight, TotalStemWeight. 
{
//     Call PotentialLeafGrowth() to compute potential growth rate of leaves.
      PotentialLeafGrowth();
//     If it is after first square, call PotentialFruitGrowth() to compute potential
//  growth rate of squares and bolls.
      if ( FruitingCode[0][0][0] > 0 ) 
		  PotentialFruitGrowth();
//     Active stem tissue (stemnew) is the difference between TotalStemWeight 
//  and the value of StemWeight(kkday). 
      int voldstm = 32; // constant parameter (days for stem tissue to become "old")
      int kkday = Kday - voldstm; // age of young stem tissue
      if ( kkday < 1 )
		   kkday = 1;
      double stemnew = TotalStemWeight - StemWeight[kkday]; // dry weight of active stem tissue.
//     Call PotentialStemGrowth() to compute PotGroStem, potential growth rate of stems.
//  The effect of temperature is introduced, by multiplying potential growth rate by DayInc. 
//  Stem growth is also affected by water stress (WaterStressStem) and possible PIX application 
//  (pixdz).   PotGroStem is limited by (maxstmgr * PerPlantArea) g per plant per day.
      PotGroStem = PotentialStemGrowth (stemnew) * DayInc * WaterStressStem * pixdz;
      double maxstmgr = 0.067; // maximum posible potential stem growth, g dm-2 day-1.
      if ( PotGroStem > maxstmgr * PerPlantArea )  
		   PotGroStem = maxstmgr * PerPlantArea;
//	   Call PotentialRootGrowth() to compute potential growth rate of roots.
      double sumpdr; // total potential growth rate of roots in g per slab. this is 
	                 // computed in PotentialRootGrowth() and used in ActualRootGrowth().
      sumpdr = PotentialRootGrowth();
//     Total potential growth rate of roots is converted from g per
//  slab (sumpdr) to g per plant (PotGroAllRoots).
      PotGroAllRoots = sumpdr * 100 * PerPlantArea / RowSpace;
//     Limit PotGroAllRoots to (maxrtgr*PerPlantArea) g per plant per day.
      double maxrtgr  = 0.045; // maximum possible potential root growth, g dm-2 day-1. 
      if ( PotGroAllRoots > maxrtgr * PerPlantArea )
		   PotGroAllRoots = maxrtgr * PerPlantArea;
//     Call DryMatterBalance() to compute carbon balance, allocation of carbon to
//  plant parts, and carbon stress. DryMatterBalance() also computes and returns the
//  values of the following arguments:
//     cdleaf is carbohydrate requirement for leaf growth, g per plant per day.
//     cdpet is carbohydrate requirement for petiole growth, g per plant per day.
//     cdroot is carbohydrate requirement for root growth, g per plant per day.
//     cdstem is carbohydrate requirement for stem growth, g per plant per day.
      double cdstem, cdleaf, cdpet, cdroot;
      DryMatterBalance(cdstem, cdleaf, cdpet, cdroot);
//     If it is after first square, call ActualFruitGrowth() to compute actual
//  growth rate of squares and bolls.
      if ( FruitingCode[0][0][0] > 0 ) 
		   ActualFruitGrowth();
//     Initialize TotalLeafWeight. It is assumed that cotyledons fall off
//  at time of first square. Also initialize TotalLeafArea and TotalPetioleWeight.
      if (FirstSquare > 0) 
	  {
         TotalLeafWeight = 0;
         TotalLeafArea = 0;
	  }
      else
	  {
	     double cotylwt = 0.20; // weight of cotyledons dry matter.
         TotalLeafWeight = cotylwt;
         TotalLeafArea = 0.6 * cotylwt;
      }
      TotalPetioleWeight = 0;
//     Call ActualLeafGrowth to compute actual growth rate of leaves and compute leaf area index.
      ActualLeafGrowth();
      LeafAreaIndex = TotalLeafArea / PerPlantArea;
//     Add ActualStemGrowth to TotalStemWeight, and define StemWeight(Kday) for this day.
      TotalStemWeight += ActualStemGrowth;
      StemWeight[Kday] = TotalStemWeight;
//     Plant density affects growth in height of tall plants.
      double htdenf = 55; // minimum plant height for plant density affecting growth in height.
      double z1;          // intermediate variable to compute denf2.
      z1 = (PlantHeight - htdenf) / htdenf;
      if ( z1 < 0 )
		   z1 = 0;
      if ( z1 > 1 ) 
		   z1 = 1;
      double denf2; // effect of plant density on plant growth in height.
      denf2 = 1 + z1 * (DensityFactor - 1);
//     Call AddPlantHeight to compute PlantHeight.
      PlantHeight += AddPlantHeight(denf2);
//     Call ActualRootGrowth() to compute actual root growth.
      ComputeActualRootGrowth( sumpdr );
//     Output data to file *.CHB
      if (OutIndex[18] > 0) 
	  {
         ofstream File36("Output\\" + ProfileName + ".CHB", ios::app);
		 File36.unsetf(ios::left);
         File36.width(11);
         File36 << Date;
		 File36.setf(ios::fixed);
	     File36.precision(5);
         File36.width(10);
         File36 << cdstem;
         File36.width(10);
         File36 << cdleaf;
         File36.width(10);
         File36 << cdpet;
         File36.width(10);
         File36 << cdroot;
         File36.width(10);
         File36 << CarbonStress;
		 File36.width(10);
         File36 << TotalStemWeight;
         File36.width(10);
         File36 << TotalLeafWeight;
         File36.width(10);
         File36 << TotalPetioleWeight;
         File36.width(10);
         File36 << TotalRootWeight;
         File36.width(10);
         File36 << sumpdr;
         File36.width(10);
         File36 << NStressRoots;
         File36 << endl;
	  }
}
