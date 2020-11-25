// File SoilNitrogen.cpp
//
//   functions in this file:
// SoilNitrogen()
// UreaHydrolysis()
// SoilWaterEffect();
// MineralizeNitrogen()
// SoilTemperatureEffect()
// Nitrification()
// Denitrification()
// SoilNitrogenBal()
// SoilNitrogenAverage()
//
#include "CottonSimulation.h"
#include "GeneralFunctions.h"
#include "resource.h"
#include <math.h>  

//////////////////////////
/*                References for soil nitrogen routines:
                ======================================
           Godwin, D.C. and Jones, C.A. 1991. Nitrogen dynamics
  in soil - plant systems. In: J. Hanks and J.T. Ritchie (ed.)
  Modeling Plant and Soil Systems, American Society of Agronomy,
  Madison, WI, USA, pp 287-321.
           Quemada, M., and Cabrera, M.L. 1995. CERES-N model predictions
  of nitrogen mineralized from cover crop residues. Soil  Sci. Soc.
  Am. J. 59:1059-1065.
           Rolston, D.E., Sharpley, A.N., Toy, D.W., Hoffman, D.L., and
  Broadbent, F.E. 1980. Denitrification as affected by irrigation
  frequency of a field soil. EPA-600/2-80-06. U.S. Environmental
  Protection Agency, Ada, OK.
           Vigil, M.F., and Kissel, D.E. 1995. Rate of nitrogen mineralized
  from incorporated crop residues as influenced by temperarure. Soil
  Sci. Soc. Am. J. 59:1636-1644.
           Vigil, M.F., Kissel, D.E., and Smith, S.J. 1991. Field crop
  recovery and modeling of nitrogen mineralized from labeled sorghum
  residues. Soil Sci. Soc. Am. J. 55:1031-1037.
************************************************************************/
void SoilNitrogen()  
//     This function computes the transformations of the nitrogen
// compounds in the soil. It is called each day from SimulateThisDay().
//     It calls UreaHydrolysis(), MineralizeNitrogen{}, Nitrification(),
// Denitrification().
//
//     The following global variables are referenced here:
//       Daynum, DayStart, dl, FieldCapacity, nk, nl, SoilTempDailyAvrg, 
//       VolWaterContent, VolNh4NContent, VolNo3NContent, VolUreaNContent.
{
      static double depth[maxl]; // depth to the end of each layer.
//     At start compute depth[l] as the depth to the bottom of each layer, cm.
      if ( Daynum <= DayStart ) 
	  {
         double sumdl = 0; // sum of layer thicknesses.
         for (int l = 0; l < nl; l++)
		 {
           sumdl += dl[l];
           depth[l] = sumdl;
		 }
	  }
//     For each soil cell: call functions UreaHydrolysis(), MineralizeNitrogen(), 
//  Nitrification() and Denitrification().
      for (int l = 0; l < nl; l++)
         for (int k = 0; k < nk; k++)
		 {
           if (VolUreaNContent[l][k] > 0)    
                          UreaHydrolysis( l, k );
           MineralizeNitrogen( l, k );
           if ( VolNh4NContent[l][k] > 0.00001 )  
                          Nitrification( l, k, depth[l] );
//     Denitrification() is called if there are enough water and nitrates in the
//  soil cell. cparmin is the minimum temperature C for denitrification.
           const double cparmin = 5; 
           if ( VolNo3NContent[l][k] > 0.001 && VolWaterContent[l][k] > FieldCapacity[l]
                && SoilTempDailyAvrg[l][k] >= (cparmin + 273.161) )
                          Denitrification( l, k );
		 }
}
//////////////////////////
void UreaHydrolysis(int l, int k)
//     This function computes the hydrolysis of urea to ammonium in the soil.
//  It is called by function SoilNitrogen(). It calls the function SoilWaterEffect().
//     The following procedure is based on the CERES routine, as
//  documented by Godwin and Jones (1991).
//
//     The following global variables are referenced here:
//       BulkDensity, FieldCapacity, FreshOrganicMatter, HumusOrganicMatter, 
//       SoilHorizonNum, SoilTempDailyAvrg, thetar, thts, VolWaterContent.
//     The following global variables are set here: 
//       VolNh4NContent, VolUreaNContent.
//     The arguments (k, l) are soil column and layer numbers.
{
//     The following constant parameters are used:
      const double ak0 = 0.25; // minimal value of ak.
      const double cak1 = 0.3416;
	  const double cak2 = 0.0776; // constant parameters for computing ak from organic carbon.
      const double stf1 = 40.0; 
	  const double stf2 = 0.20; // constant parameters for computing stf.
      const double swf1 = 0.20; // constant parameter for computing swf.
//     Compute the organic carbon in the soil (converted from mg / cm3 to % by weight) for the
//  sum of stable and fresh organic matter, assuming 0.4 carbon content in soil organic matter.
      int j = SoilHorizonNum[l]; // profile horizon number for this soil layer.
      double oc; // organic carbon in the soil (% by weight).
      oc = 0.4 * (FreshOrganicMatter[l][k] + HumusOrganicMatter[l][k]) * 0.1 / BulkDensity[j]; 
//     Compute the potential rate of hydrolysis of urea. It is assumed
//  that the potential rate will not be lower than ak0 = 0.25 .
      double ak;  // potential rate of urea hydrolysis (day-1).
      ak = cak1 + cak2 * oc; 
      if (ak < ak0)
		  ak = ak0;
//     Compute the effect of soil moisture using function SoilWaterEffect on the rate of urea 
//  hydrolysis. The constant swf1 is added to the soil moisture function for mineralization, 
      double swf; // soil moisture effect on rate of urea hydrolysis.
      swf = SoilWaterEffect(l, k,  0.5) + swf1;
      if (swf < 0)
		  swf = 0;
      if (swf > 1)
		  swf = 1;
//     Compute the effect of soil temperature. The following parameters are used for the
//  temperature function: stf1, stf2.
      double stf; // soil temperature effect on rate of urea hydrolysis.
      stf = (SoilTempDailyAvrg[l][k] - 273.161) / stf1 + stf2;
      if (stf > 1)
		  stf = 1;
      if (stf < 0)
		  stf = 0;
//    Compute the actual amount of urea hydrolized, and update VolUreaNContent
// and VolNh4NContent.
      double hydrur; // amount of urea hydrolized, mg N cm-3 day-1.
      hydrur = ak * swf * stf * VolUreaNContent[l][k];
      if (hydrur > VolUreaNContent[l][k])
		  hydrur = VolUreaNContent[l][k];
      VolUreaNContent[l][k] -= hydrur;
      VolNh4NContent[l][k] += hydrur;
/*********************************************************************
     Note: Since COTTON2K does not require soil pH in the input, the
  CERES rate equation was modified as follows: 
     ak is a function of soil organic matter, with two site-dependent 
  parameters cak1 and cak2. Their values are functions of the prevalent pH:
               cak1 = -1.12 + 0.203 * pH
               cak2 = 0.524 - 0.062 * pH
    Some examples of these values:
            pH      cak1     cak2
           6.8      .2604    .1024
           7.2      .3416    .0776
           7.6      .4228    .0528
   The values for pH 7.2 are used.
*/
}
/////////////////////////
double SoilWaterEffect (int l, int k, double xx)
//     This function computes the effect of soil moisture on the rate of mineralization of 
//  organic mineralizable nitrogen, and on the rates of urea hydrolysis and nitrification.
//     It is based on Godwin and Jones (1991).
//     The following global variables are referenced:
//       FieldCapacity, thetar, thts, VolWaterContent.
//     The argument xx is 0.5 when used for mineralization and urea hydrolysis,
//  or 1.0 when used for nitrification.
//     l, k are layer and column of this cell.
//
{
      double wf; // the effect of soil moisture on process rate.
      if ( VolWaterContent[l][k] <= FieldCapacity[l]) 
//     Soil water content less than field capacity:
         wf = (VolWaterContent[l][k] - thetar[l]) / (FieldCapacity[l] - thetar[l]);
      else
//     Soil water content more than field capacity:
         wf = 1 - xx * (VolWaterContent[l][k] - FieldCapacity[l])
                     / (thts[l] - FieldCapacity[l]);
//
	  if (wf < 0)
		  wf = 0;
      return wf;
}
///////////////////////////////////////////
void MineralizeNitrogen(int l, int k)
//     This function computes the mineralization of organic nitrogen in the soil, and the 
//  immobilization of mineral nitrogen by soil microorganisms. It is called by function 
//  SoilNitrogen(). 
//     It calls the following functions: SoilTemperatureEffect(), SoilWaterEffect().
//     The procedure is based on the CERES routines, as documented by Godwin and Jones (1991).
//     Note: CERES routines assume freshly incorporated organic matter consists of 20% 
//  carbohydrates, 70% cellulose, and 10% lignin, with maximum decay rates of 0.2, 0.05 
//  and 0.0095, respectively. Quemada and Cabrera (1995) suggested decay rates of 0.14, 
//  0.0034, and 0.00095 per day for carbohydrates, cellulose and lignin, respectively.
//     Assuming cotton stalks consist of 20% carbohydrates, 50% cellulose and 30% lignin - 
//  the average maximum decay rate of FreshOrganicMatter decay rate = 0.03 will be used here.
//
//     The following global variables are referenced here:
//  Daynum, DayStart, dl, SoilTempDailyAvrg, wk.
//     The following global variables are set here:
//  FreshOrganicMatter, FreshOrganicNitrogen, HumusNitrogen, HumusOrganicMatter, 
//  MineralizedOrganicN, VolNh4NContent, VolNo3NContent
//     The arguments (k, l) are soil column and layer numbers.
{
//     The following constant parameters are used:
      const double cnfresh = 25; // C/N ratio in fresh organic matter.
      const double cnhum = 10;   // C/N ratio in stabilized organic matter (humus).
      const double cnmax = 13;   // C/N ratio higher than this reduces rate of mineralization.
      const double cparcnrf = 0.693; // constant parameter for computing cnRatioEffect.
      const double cparHumusN = 0.20;  // ratio of N released from fresh OM incorporated in the humus.
	  const double cparMinNH4 = 0.00025; // mimimum NH4 N remaining after mineralization.
      const double decayRateFresh = 0.03; // decay rate constant for fresh organic matter.
	  const double decayRateHumus = 0.000083;// decay rate constant for humic organic matter.
//     On the first day of simulation set initial values for N in fresh organic 
//  matter and in humus, assuming C/N ratios of cnfresh = 25 and cnhum = 10, 
//  respectively. Carbon in soil organic matter is 0.4 of its dry weight.
      if ( Daynum <= DayStart ) 
	  {
         FreshOrganicNitrogen[l][k] = FreshOrganicMatter[l][k] * 0.4 / cnfresh;
         HumusNitrogen[l][k] = HumusOrganicMatter[l][k] * 0.4 / cnhum;
      }
//     This function will not be executed for soil cells with no organic matter in them.
      if ( FreshOrganicMatter[l][k] <= 0 && HumusOrganicMatter[l][k] <= 0 )
		  return;
//
// **  C/N ratio in soil **
//     The C/N ratio (cnRatio) is computed for the fresh organic matter and the nitrate and 
//  ammonium nitrogen in the soil. It is assumed that C/N ratios higher than cnmax reduce the 
// rate of mineralization. Following the findings of Vigil et al. (1991) the value of cnmax 
// is set to 13. 
      double cnRatio = 1000;    // C/N ratio in fresh organic matter and mineral N in soil.
      double cnRatioEffect = 1; // the effect of C/N ratio on rate of mineralization.
      double totalSoilN;        // total N in the soil cell, excluding the stable humus fraction, mg/cm3
      totalSoilN = FreshOrganicNitrogen[l][k] + VolNo3NContent[l][k] + VolNh4NContent[l][k];
      if ( totalSoilN > 0) 
	  {
         cnRatio = FreshOrganicMatter[l][k] * 0.4  / totalSoilN;
         if (cnRatio >= 1000) 
             cnRatioEffect = 0;
         else if (cnRatio > cnmax) 
             cnRatioEffect = exp(-cparcnrf * (cnRatio - cnmax) / cnmax);
		 else
             cnRatioEffect = 1;
      }
//
// **  Mineralization of fresh organic matter **
//     The effects of soil moisture (wf) and of soil temperature (tfac) are computed.
      double wf = SoilWaterEffect(l, k, 0.5);
      double tfac = SoilTemperatureEffect(SoilTempDailyAvrg[l][k] - 273.161 );
//     The gross release of dry weight and of N from decomposition of fresh organic matter is computed.
      double grossReleaseN; // gross release of N from decomposition, mg/cm3
      double immobilizationRateN; // immobilization rate of N associated with decay of residues, mg/cm3 .
      if ( FreshOrganicMatter[l][k] > 0.00001 ) 
	  {
//     The decayRateFresh constant (= 0.03) is modified by soil temperature, soil moisture, 
//  and the C/N ratio effect.
         double g1; // the actual decay rate of fresh organic matter, day-1.
         g1 = tfac * wf * cnRatioEffect * decayRateFresh;
         double grossReleaseDW; // the gross release of dry weight from decomposition, mg/cm3
         grossReleaseDW = g1 * FreshOrganicMatter[l][k];
         grossReleaseN = g1 * FreshOrganicNitrogen[l][k];
//     The amount of N required for microbial decay of a unit of fresh organic matter suggested 
//  in CERES is 0.02 (derived from:  C fraction in FreshOrganicMatter (=0.4) * biological 
//  efficiency of C turnover by microbes (=0.4) * N/C ratio in microbes (=0.125) ). However,
//  Vigil et al. (1991) suggested that this value is 0.0165.
	     const double cparnreq = 0.0165; // The amount of N required for decay of fresh organic matter
//     Substract from this the N ratio in the decaying FreshOrganicMatter, and multiply
//  by grossReleaseDW to get the amount needed (immobilizationRateN) in mg cm-3. Negative value
//  indicates that there is enough N for microbial decay.
         immobilizationRateN =  grossReleaseDW 
                * (cparnreq - FreshOrganicNitrogen[l][k]/FreshOrganicMatter[l][k]);
//     All computations assume that the amounts of VolNh4NContent and VNO3C will
//  each not become lower than cparMinNH4 (= 0.00025) .
         double rnac1; // the maximum possible value of immobilizationRateN, mg/cm3 .
         rnac1 = VolNh4NContent[l][k] + VolNo3NContent[l][k] - 2 * cparMinNH4;
         if (immobilizationRateN > rnac1)
			 immobilizationRateN = rnac1;
         if (immobilizationRateN < 0)
			 immobilizationRateN = 0;
//     FreshOrganicMatter and FreshOrganicNitrogen (the N in it) are now updated.
         FreshOrganicMatter[l][k] -= grossReleaseDW;
         FreshOrganicNitrogen[l][k] += immobilizationRateN - grossReleaseN;
	  }
      else
	  {
		  grossReleaseN = 0;
          immobilizationRateN = 0;
      }
//
// **  Mineralization of humic organic matter **
//     The mineralization of the humic fraction (rhmin) is now computed. decayRateHumus = 0.000083
//  is the humic fraction decay rate (day-1). It is modified by soil temperature and soil moisture.
      double rhmin; // N mineralized from the stable humic fraction, mg/cm3 .
      rhmin = HumusNitrogen[l][k] * decayRateHumus * tfac * wf;
//     rhmin is substacted from HumusNitrogen, and a corresponding amount of dry matter is
//  substracted from HumusOrganicMatter (assuming C/N = cnhum = 10).
//     It is assumed that 20% (=cparHumusN) of the N released from the fresh
//  organic matter is incorporated in the humus, and a parallel amount of
//  dry matter is also incorporated in it (assuming C/N = cnfresh = 25).
      HumusNitrogen[l][k] -= rhmin + cparHumusN * grossReleaseN;
      HumusOrganicMatter[l][k] -= cnhum * rhmin / 0.4 
                                + cparHumusN * cnfresh * grossReleaseN / 0.4;
//     80% (1 - cparHumusN) of the N released from the fresh organic matter , the N released
//  from the decay of the humus, and the immobilized N are used to compute netNReleased.
//     Negative value of netNReleased indicates net N immobilization.
      double netNReleased; // the net N released from all organic sources (mg/cm3).
      netNReleased = (1 - cparHumusN) * grossReleaseN + rhmin - immobilizationRateN;
//     If the net N released is positive, it is added to the NH4 fraction. MineralizedOrganicN, 
//  the accumulated nitrogen released by mineralization in the slab, is updated.
      if ( netNReleased > 0 ) 
	  {
         VolNh4NContent[l][k] += netNReleased;
         MineralizedOrganicN += netNReleased * dl[l] * wk[k];
	  }
//     If net N released is negative (net immobilization), the NH4 fraction
//  is reduced, but at least 0.25 ppm (=cparMinNH4 in mg cm-3) of NH4 N
//  should remain. A matching amount of N is added to the organic N fraction.
//     MineralizedOrganicN, the accumulated nitrogen released by mineralization in the
//  slab is updated.
      else
	  {
         double addvnc = 0; // immobilised N added to the organic fraction.
         double nnom1 = 0; // temporary storage of netNReleased (if N is also
                           // immobilized from NO3).
         if ( VolNh4NContent[l][k] > cparMinNH4 ) 
		 {
            if ( fabs(netNReleased) < (VolNh4NContent[l][k] - cparMinNH4) ) 
               addvnc = -netNReleased;
            else
               addvnc = VolNh4NContent[l][k] - cparMinNH4;
            VolNh4NContent[l][k] -= addvnc;
            MineralizedOrganicN -= addvnc * dl[l] * wk[k];
            FreshOrganicNitrogen[l][k] += addvnc;
            nnom1 = netNReleased + addvnc;
		 } 
//     If immobilization is larger than the use of NH4 nitrogen, the NO3
//  fraction is reduced in a similar procedure.
         if ( nnom1 < 0 && VolNo3NContent[l][k] > cparMinNH4 ) 
		 {
            if ( fabs(nnom1) < (VolNo3NContent[l][k] - cparMinNH4) ) 
               addvnc = -nnom1;
            else
               addvnc = VolNo3NContent[l][k] - cparMinNH4;
            VolNo3NContent[l][k] -= addvnc;
            FreshOrganicNitrogen[l][k] += addvnc;
            MineralizedOrganicN -= addvnc * dl[l] * wk[k];
         }
      }  
}
/////////////////////////
double SoilTemperatureEffect (double tt )
//     This function computes the effect of temperature on the rate
//  of mineralization of organic mineralizable nitrogen. It is based on
//  GODWIN and JONES (1991).
//     The following argument is used:  tt - soil temperature (C).
{
//     The following constant parameters are used:
      const double tfpar1 =  0.010645;
      const double tfpar2 =  0.12979;
//     The temperature function of CERES is replaced by the function
//  suggested by Vigil and Kissel (1995):
//               tfm = 0.010645 * exp(0.12979 * tt)
//     Note: tfm = 0.5 for 29.66 C, tfm = 1 for 35 C, tfm = 2 for 40.34 C.
	  double tfm;
      tfm = tfpar1 * exp( tfpar2 * tt );
      if (tfm < 0)
		  tfm = 0;
      if (tfm > 2)
		  tfm = 2;
      return tfm;
}
//////////////////////////////////
void Nitrification(int l, int k, double DepthOfLayer)
//     This function computes the transformation of soil ammonia nitrogen to nitrate. 
//  It is called by SoilNitrogen(). It calls the function SoilWaterEffect()
//
//     The following global variable is referenced here:    SoilTempDailyAvrg
//     The following global variables are set here:   VolNh4NContent, VolNo3NContent
//     The following arguments are used:
//  DepthOfLayer - depth to the bottom of this layer, cm.
//  k, l - soil column and layer numbers.
//
{
//     The following constant parameters are used:
      const double cpardepth = 0.45;
	  const double cparnit1  = 24.635;
	  const double cparnit2  = 8227;
	  const double cparsanc  = 204; // this constant parameter is modified from kg/ha units in CERES
                              // to mg/cm3 units of VolNh4NContent (assuming 15 cm layers)
      double sanc; // effect of NH4 N in the soil on nitrification rate (0 to 1).
      if (VolNh4NContent[l][k] <  0.1) 
	     sanc = 1 - exp(-cparsanc * VolNh4NContent[l][k]);
	  else
  	     sanc = 1;
//     The rate of nitrification, con1, is a function of soil temperature. It is slightly 
//  modified from GOSSYM. it is transformed from immediate rate to a daily time step ratenit.
//     The rate is modified by soil depth, assuming that for an increment
//  of 30 cm depth, the rate is decreased by 55% (multiply by a power of
//  cpardepth). It is also multiplied by the environmental limiting
//  factors (sanc, SoilWaterEffect) to get the actual rate of nitrification.
//     The maximum rate is assumed not higher than 10%.
	  double con1;    // rate of nitrification as a function of temperature.
      con1 = exp( cparnit1 - cparnit2 / SoilTempDailyAvrg[l][k] );
      double ratenit; // actual rate of nitrification (day-1).
      ratenit = 1 - exp(-con1);
      double tff; // effect of soil depth on nitrification rate.
      tff = (DepthOfLayer - 30) / 30;
      if (tff < 0) 
		  tff = 0;
//     Add the effects of NH4 in soil, soil water content, and depth of soil layer.
      ratenit = ratenit * sanc * SoilWaterEffect(l, k, 1) * pow(cpardepth, tff);
      if (ratenit < 0) 
		  ratenit = 0;
      if (ratenit > 0.10)
		  ratenit = 0.10;
// Compute the actual amount of N nitrified, and update VolNh4NContent and VolNo3NContent.
      double dnit; // actual nitrification (mg n cm-3 day-1).
      dnit = ratenit * VolNh4NContent[l][k];
      VolNh4NContent[l][k] -= dnit;
      VolNo3NContent[l][k] += dnit;
}
/////////////////////////
void Denitrification(int l, int k)
//     This function computes the denitrification of nitrate N in the soil. 
//     It is called by function SoilNitrogen().
//     The procedure is based on the CERES routine, as documented by Godwin and Jones (1991).
//
//     The following global variables are referenced here:
//       dl, FieldCapacity, HumusOrganicMatter, SoilTempDailyAvrg, thts, VolWaterContent, wk
//    The following global variables are set here:     SoilNitrogenLoss, VolNo3NContent
{
//    The following constant parameters are used:
      const double cpar01 = 24.5;
	  const double cpar02 = 3.1;
	  const double cpardenit = 0.00006;
	  const double cparft = 0.046;
	  const double cparhum = 0.58;
      const double vno3min = 0.00025; 
//
      double soilc; // soil carbon content, mg/cm3.
//     soilc is calculated as 0.58 (cparhum) of the stable humic fraction
//  (following CERES), and cw is estimated following Rolston et al. (1980).
      soilc = cparhum * HumusOrganicMatter[l][k];
      double cw;    // water soluble carbon content of soil, ppm.
      cw = cpar01 + cpar02 * soilc;
//     The effects of soil moisture (fw) and soil temperature (ft) are computed as 0 to 1 factors.
      double fw; // effect of soil moisture on denitrification rate.
      fw = (VolWaterContent[l][k] - FieldCapacity[l]) / (thts[l] - FieldCapacity[l]);
      if (fw < 0)
		  fw = 0;
      double ft; // effect of soil temperature on denitrification rate.
      ft = 0.1 * exp(cparft * (SoilTempDailyAvrg[l][k] - 273.161));
      if (ft > 1)
		  ft = 1;
//     The actual rate of denitrification is calculated. The equation is modified from CERES to 
//  units of mg/cm3/day. 
      double dnrate; // actual rate of denitrification, mg N per cm3 of soil per day.
      dnrate = cpardenit * cw * VolNo3NContent[l][k] * fw * ft;
//     Make sure that a minimal amount of nitrate will remain after denitrification.
      if ( dnrate > (VolNo3NContent[l][k] - vno3min) )
           dnrate = VolNo3NContent[l][k] - vno3min;
      if ( dnrate < 0 )
		   dnrate = 0;
//     Update VolNo3NContent, and add the amount of nitrogen lost to SoilNitrogenLoss.
      VolNo3NContent[l][k] -= dnrate;
      SoilNitrogenLoss += dnrate * dl[l] * wk[k];
}
/////////////////////////
void SoilNitrogenBal(const string& ProfileName)  
//     This function computes the nitrogen balances in the soil,
//  for diagnostic purposes. It is called by SimulateThisDay().
//     The following global variables are referenced here:
//       CumFertilizerN, CumNitrogenUptake, Kday, MineralizedOrganicN, SoilNitrogenAtStart, 
//       SoilNitrogenLoss, TotalSoilNitrogen.
{
      double balsn; // soil nitrogen balance. It should be zero.
//     Calculate soil nitrogen balance, in units of mg per slab.
//     The "plus side" is the amount of mineral N in the soil at the beginning of simulation, 
//  the amount added by fertilizer, and the amount generated by mineralization of organic matter. 
//     The "minus side" is the amount taken up by the plants, present amount of mineral
//  N in the soil, and the amount lost by drainage.
      balsn = SoilNitrogenAtStart + CumFertilizerN + MineralizedOrganicN
            - CumNitrogenUptake - TotalSoilNitrogen - SoilNitrogenLoss;
//     Output to file NB1
      ofstream File47("Output\\" + ProfileName + ".NB1", ios::app);
	  File47.width(4);
	  File47 << Kday;
	  File47.setf(ios::fixed);
	  File47.precision(2);
	  File47.width(9);
	  File47 << balsn;
	  File47.width(9);
	  File47 << CumFertilizerN;
	  File47.width(9);
	  File47 << MineralizedOrganicN;
	  File47.width(9);
	  File47 << CumNitrogenUptake;
	  File47.width(9);
	  File47 << TotalSoilNitrogen;
	  File47.width(9);
	  File47 << SoilNitrogenLoss << endl;
}
//////////////////////////
void SoilNitrogenAverage(const string& ProfileName)
//     This function computes average values of soil N. It is called by SimulateThisDay().
//     The following global variables are referenced here:
//  dl, Kday, nk, VolNh4NContent, VolNo3NContent.
//
{
      double avno30 = 0; // Aaverage NO3-N content in four successive 30 cm layers of soil.
      double avno60 = 0;
      double avno90 = 0;
      double avno120 = 0;
      double avnh30 = 0; // Aaverage NH4-N content in four successive 30 cm layers of soil.
      double avnh60 = 0;
      double avnh90 = 0;
      double avnh120 = 0;
//
//     Compute average soil N content by layers 0-30, 30-60, 60-90,
//  and 90-120 cm, in ppm per volume, and write it to file NB0.
      for (int k = 0; k < nk; k++)
	  {
         for (int l = 0; l < 8; l++)
		 {
		    avno30 += VolNo3NContent[l][k] * dl[l];
            avnh30 += VolNh4NContent[l][k] * dl[l];
         }
         for (int l = 8; l < 14; l++)
		 {
            avno60 += VolNo3NContent[l][k] * dl[l];
            avnh60 += VolNh4NContent[l][k] * dl[l];
         }
         for (int l = 14; l < 20; l++)
		 {
            avno90 += VolNo3NContent[l][k] * dl[l];
            avnh90 += VolNh4NContent[l][k] * dl[l];
         }
         for (int l = 20; l < 26; l++)
		 {
            avno120 += VolNo3NContent[l][k] * dl[l];
            avnh120 += VolNh4NContent[l][k] * dl[l];
		 }
      }
//
      avno30 = 1000 * avno30 / (30 * nk);
      avnh30 = 1000 * avnh30 / (30 * nk);
      avno60 = 1000 * avno60 / (30 * nk);
      avnh60 = 1000 * avnh60 / (30 * nk);
      avno90 = 1000 * avno90 / (30 * nk);
      avnh90 = 1000 * avnh90 / (30 * nk);
      avno120 = 1000 *  avno120 / (30 * nk);
      avnh120 = 1000 *  avnh120 / (30 * nk);
//
      ofstream File35("Output\\" + ProfileName + ".NB0", ios::app);
	  File35.width(4);
	  File35 << Kday;
	  File35.setf(ios::fixed);
	  File35.precision(2);
	  File35.width(9);
	  File35 << avno30;
	  File35.width(9);
	  File35 << avno60;
	  File35.width(9);
	  File35 << avno90;
	  File35.width(9);
	  File35 << avno120;
	  File35.width(9);
	  File35 << avnh30;
	  File35.width(9);
	  File35 << avnh60;
	  File35.width(9);
	  File35 << avnh90;
	  File35.width(9);
	  File35 << avnh120 << endl;
}
