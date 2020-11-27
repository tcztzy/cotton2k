//  FruitAbscission.cpp
//
//   functions in this file:
//       FruitingSitesAbscission()
//       SiteAbscissionRatio()
//       SquareAbscission()
//       BollAbscission()
//       AdjustAbscission()
//       AdjustSquareAbscission()
//       AdjustYoungBollAbscission()
//       AdjustSetBollAbscission()
//       AdjustBollAbscission()
//       ComputeSiteNumbers()
//
#include "CottonSimulation.h"
#include "GeneralFunctions.h"
#include <math.h>

//////////////////////////////////////////////////
void FruitingSitesAbscission()
//     This function simulates the abscission of squares and bolls.
//  It is called from function CottonPhenology().  It calls SiteAbscissionRatio(), 
//	SquareAbscission(), BollAbscission(), AdjustAbscission() and ComputeSiteNumbers()
//
//     The following global variables are referenced here:
//  CarbonStress, DayInc, Daynum, FruitingCode, ginp, Gintot, Kday, KdayAdjust, NitrogenStress, 
//  NumAdjustDays, NumFruitBranches, NumNodes, NumVegBranches, VarPar, WaterStress.
//
//     The following global variable are set here:
//  AbscissionLag, NumSheddingTags, ShedByCarbonStress, ShedByNitrogenStress, ShedByWaterStress.
//
{
//     The following constant parameters are used:
	  const double vabsfr[9] = { 21.0, 0.42, 30.0, 0.05, 6.0, 2.25, 0.60, 5.0, 0.20 };
//
//     Update tags for shedding: Increment NumSheddingTags by 1, and move the array members 
// of ShedByCarbonStress, ShedByNitrogenStress, ShedByWaterStress, and AbscissionLag.
      NumSheddingTags++;
      if ( NumSheddingTags > 1 ) 
         for (int lt = NumSheddingTags-1; lt > 0; lt--)
		 {
			 int ltm1 = lt - 1;
             ShedByCarbonStress[lt] = ShedByCarbonStress[ltm1];
             ShedByNitrogenStress[lt] = ShedByNitrogenStress[ltm1];
             ShedByWaterStress[lt] = ShedByWaterStress[ltm1];
             AbscissionLag[lt] = AbscissionLag[ltm1];
		 }
//     Calculate shedding intensity: The shedding intensity due to stresses of this day is assigned
//  to the first members of the arrays ShedByCarbonStress, ShedByNitrogenStress, and ShedByWaterStress.
      if ( CarbonStress < VarPar[43] )
           ShedByCarbonStress[0] = (VarPar[43] - CarbonStress) / VarPar[43];
	  else
           ShedByCarbonStress[0] = 0;
      if ( NitrogenStress < vabsfr[1] )
           ShedByNitrogenStress[0] = (vabsfr[1] - NitrogenStress) / vabsfr[1];
	  else
           ShedByNitrogenStress[0] = 0;
      if ( WaterStress < VarPar[44] )
           ShedByWaterStress[0] = (VarPar[44] - WaterStress) / VarPar[44];
	  else
           ShedByWaterStress[0] = 0;
//     Assign 0.01 to the first member of AbscissionLag.
      AbscissionLag[0] = 0.01;
//     Updating age of tags for shedding: Each member of array AbscissionLag is 
//  incremented by physiological age of today. It is further increased (i.e., shedding 
//  will occur sooner) when maximum temperatures are high.
      double tmax = GetFromClim("tmax", Daynum);
      for (int lt = 0; lt < NumSheddingTags; lt++)
	  {
         AbscissionLag[lt] += max(DayInc, 0.40);
         if ( tmax > vabsfr[2] )
            AbscissionLag[lt] += (tmax - vabsfr[2]) * vabsfr[3];
	  }
//     Assign zero to idecr, and start do loop over all days since
//  first relevant tagging. If AbscissionLag reaches a value of vabsfr[4],
//  calculate actual shedding of each site:
      int idecr = 0; // decrease in NumSheddingTags after shedding has been executed.
      for (int lt = 0; lt < NumSheddingTags; lt++)
	  {
         if ( AbscissionLag[lt] >= vabsfr[4] || lt >= 20 ) 
		 {
            double gin1;
            if ( Gintot > 0 ) 
               gin1 = Gintot;
            else
               gin1 = ginp;
//      Start loop over all possible fruiting sites. The abscission functions
//  will be called for sites that are squares or green bolls.
            for (int k = 0; k < NumVegBranches; k++)
			{
               int nbrch = NumFruitBranches[k]; // fruiting branch number.
               for (int l = 0; l < nbrch; l++)
			   {
                  int nnid = NumNodes[k][l]; // node number on fruiting branch.
                  for (int m = 0; m < nnid; m++)
				  {
                     if ( FruitingCode[k][l][m] == 1 || FruitingCode[k][l][m] == 2 || FruitingCode[k][l][m] == 7 )   
					 {
                        double abscissionRatio; // ratio of abscission for a fruiting site.
                        abscissionRatio = SiteAbscissionRatio(k, l, m, lt);
						if (abscissionRatio > 0)
						{
                           if (FruitingCode[k][l][m] == 1) 
                              SquareAbscission(k, l, m, abscissionRatio);
                           else
                              BollAbscission(k,l,m, abscissionRatio, gin1);
						}
                     }
                  }  // for m
			   }  // for l
            }  // for k
//      Assign zero to the array members for this day.
            ShedByCarbonStress[lt] = 0;
            ShedByNitrogenStress[lt] = 0;
            ShedByWaterStress[lt] = 0;
            AbscissionLag[lt] = 0;
            idecr++;
         } // if AbscissionLag
	  } // for lt
//  Decrease NumSheddingTags. If plantmap adjustments are necessary for square 
// number, or green boll number, or open boll number - call AdjustAbscission().
      NumSheddingTags = NumSheddingTags - idecr;
//
      if (Kday > KdayAdjust && Kday <= (KdayAdjust + NumAdjustDays)) 
          AdjustAbscission();
//
      ComputeSiteNumbers();
}
/////////////////////////
double SiteAbscissionRatio(int k, int l, int m, int lt)
//     This function computes and returns the probability of abscission of a single 
//  site (k, l, m). It is called from function FruitingSitesAbscission(). 
//
//     The following global variables are referenced here:   
//        AgeOfSite, AgeOfBoll, FruitingCode, ShedByCarbonStress, 
//        ShedByNitrogenStress, ShedByWaterStress, VarPar
//
//     The following arguments are used here:
//        k, l, m - indices defining position of this site.
//        lt - lag index for this node.
//
{
//     The following constant parameters are used:
	  const double vabsc[5] = { 21.0, 2.25, 0.60, 5.0, 0.20 };
//
//     For each site, compute the probability of its abscission (pabs) 
//  as afunction of site age, and the total shedding ratio (shedt) as a
//  function of plant stresses that occurred when abscission was
//  triggered.
      double pabs = 0;  // probability of abscission of a fruiting site.
      double shedt = 0; // total shedding ratio, caused by various stresses.
//     (1) Squares (FruitingCode = 1).  
      if ( FruitingCode[k][l][m] == 1 ) 
	  {
         if ( AgeOfSite[k][l][m] < vabsc[3] ) 
            pabs = 0; // No abscission of very young squares (AgeOfSite less than vabsc(3))
         else
		 {
            double xsqage; // square age after becoming susceptible to shedding.
            xsqage = AgeOfSite[k][l][m] - vabsc[3];
            if ( xsqage >= vabsc[0] ) 
               pabs = VarPar[46]; // Old squares have a constant probability of shedding.
            else
//     Between these limits, pabs is a function of xsqage.
               pabs = VarPar[46] + (VarPar[45] - VarPar[46]) 
			          * pow( ((vabsc[0] - xsqage) / vabsc[0]), vabsc[1] );
         }
//     Total shedding ratio (shedt) is a product of the effects of
//  carbohydrate stress and nitrogen stress.
         shedt = 1 - (1 - ShedByCarbonStress[lt]) * (1 - ShedByNitrogenStress[lt]);
      }
//     (2) Very young bolls (FruitingCode = 7, and AgeOfBoll less than VarPar[47]). 
      else if ( FruitingCode[k][l][m] == 7 && AgeOfBoll[k][l][m] <= VarPar[47] ) 
	  {
//     There is a constant probability of shedding (VarPar[48]), and shedt is a product 
//  of the effects carbohydrate, and nitrogen stresses. Note that nitrogen stress has only a
//  partial effect in this case, as modified by vabsc[2].
         pabs = VarPar[48];
         shedt = 1 - (1 - ShedByCarbonStress[lt]) * (1 - vabsc[2] * ShedByNitrogenStress[lt]);
      }
//     (3) Medium age bolls (AgeOfBoll between VarPar[47] and VarPar[47] + VarPar[49]). 
      else if ( AgeOfBoll[k][l][m] > VarPar[47] && AgeOfBoll[k][l][m] <= (VarPar[47] + VarPar[49]) ) 
	  {
//     pabs is linearly decreasing with age, and shedt is a product of the effects 
//  carbohydrate, nitrogen and water stresses.  Note that nitrogen stress has only 
//  a partial effect in this case, as modified by vabsc[4].
         pabs = VarPar[48] - (VarPar[48] - VarPar[50]) * (AgeOfBoll[k][l][m] - VarPar[47]) / VarPar[49];
         shedt = 1 - (1 - ShedByCarbonStress[lt]) * (1 - vabsc[4] * ShedByNitrogenStress[lt]) * (1 - ShedByWaterStress[lt]);
      }
//     (4) Older bolls (AgeOfBoll between VarPar[47] + VarPar[49] and VarPar[47] + 2*VarPar[49]). 
      else if (AgeOfBoll[k][l][m] > (VarPar[47] + VarPar[49]) && AgeOfBoll[k][l][m] <= (VarPar[47] + 2*VarPar[49])) 
	  {
//     pabs is linearly decreasing with age, and shedt is affected only by water stress.
         pabs = VarPar[50] / VarPar[49] * (VarPar[47] + 2 * VarPar[49] - AgeOfBoll[k][l][m]);
         shedt = ShedByWaterStress[lt];
      }
//     (5) bolls older than VarPar[47] + 2*VarPar[49]
      else if ( AgeOfBoll[k][l][m] > (VarPar[47] + 2*VarPar[49]))  
		  pabs = 0; // no abscission
//      Actual abscission of tagged sites (abscissionRatio) is a product of pabs,
//  shedt and DayInc for this day. It can not be greater than 1.
      double abscissionRatio = pabs * shedt * DayInc;
      if ( abscissionRatio > 1 )
		   abscissionRatio = 1;
      return abscissionRatio;
}
//////////////////////////////////////////////////////////////////
void SquareAbscission(int k, int l, int m, double abscissionRatio)
//     This function simulates the abscission of a single square
//  at site (k, l, m). It is called from function FruitingSitesAbscission() 
//  if this site is a square. 
//
//     The following global variable is referenced here:   SquareNConc
//
//     The following global variable are set here:
//   BloomWeightLoss, CumPlantNLoss, FruitingCode, FruitFraction, SquareNitrogen, 
//   SquareWeight, TotalSquareWeight.
//
//     The following arguments are used in this function:
//        abscissionRatio - ratio of abscission of a fruiting site.
//        k, l, m - indices defining position of this site.
//
{
//     Compute the square weight lost by shedding (wtlos) as a proportion of SquareWeight
//  of this site. Update SquareNitrogen, CumPlantNLoss, SquareWeight[k][l][m], BloomWeightLoss, 
//  TotalSquareWeight, and FruitFraction[k][l][m]. 
      double wtlos = SquareWeight[k][l][m] * abscissionRatio; // weight lost by shedding at this site.
      SquareNitrogen -= wtlos * SquareNConc;
      CumPlantNLoss += wtlos * SquareNConc;
      SquareWeight[k][l][m] -= wtlos;
      BloomWeightLoss += wtlos;
      TotalSquareWeight -= wtlos;
      FruitFraction[k][l][m] = FruitFraction[k][l][m] * (1 - abscissionRatio);
//     If FruitFraction[k][l][m] is less than 0.001 make it zero, and update
//  SquareNitrogen, CumPlantNLoss, BloomWeightLoss, TotalSquareWeight, SquareWeight[k][l][m], 
//  and assign 5 to FruitingCode.
      if (FruitFraction[k][l][m] <= 0.001) 
	  {
          FruitFraction[k][l][m] = 0;
          SquareNitrogen -= SquareWeight[k][l][m] * SquareNConc;
          CumPlantNLoss += SquareWeight[k][l][m] * SquareNConc;
          BloomWeightLoss += SquareWeight[k][l][m];
          TotalSquareWeight -= SquareWeight[k][l][m];
          SquareWeight[k][l][m] = 0;
          FruitingCode[k][l][m] = 5;
      }
}
////////////////////////
void BollAbscission(int k, int l, int m, double abscissionRatio, double gin1)
//     This function simulates the abscission of a single green boll
//  at site (k, l, m). It is called from function FruitingSitesAbscission() if this site
//  is a green boll. 
//
//     The following global variables are referenced here:
//   BurrNConc, pixcon, SeedNConc
//
//     The following global variable are set here:
//   BollWeight, BurrNitrogen, BurrWeight, BurrWeightGreenBolls, CottonWeightGreenBolls,
//   CumPlantNLoss, FruitFraction, FruitingCode, FruitFraction, GreenBollsLost, 
//   PixInPlants, SeedNitrogen, 
//
//     The following arguments are used in this function:
//        abscissionRatio - ratio of abscission of a fruiting site.
//        gin1 - percent of seeds in seedcotton, used to compute lost nitrogen.
//        k, l, m - location of this site on the plant
//
{
//     Update SeedNitrogen, BurrNitrogen, CumPlantNLoss, PixInPlants, GreenBollsLost, CottonWeightGreenBolls, BurrWeightGreenBolls, 
//  BollWeight[k][l][m], BurrWeight[k][l][m], and FruitFraction[k][l][m].
      SeedNitrogen -= BollWeight[k][l][m] * abscissionRatio * (1 - gin1) * SeedNConc;
      BurrNitrogen -= BurrWeight[k][l][m] * abscissionRatio * BurrNConc;
      CumPlantNLoss += BollWeight[k][l][m] * abscissionRatio * (1.- gin1) * SeedNConc;
      CumPlantNLoss += BurrWeight[k][l][m] * abscissionRatio * BurrNConc;
      PixInPlants -= (BollWeight[k][l][m] + BurrWeight[k][l][m]) * abscissionRatio * pixcon;
      GreenBollsLost += (BollWeight[k][l][m] + BurrWeight[k][l][m]) * abscissionRatio;
      CottonWeightGreenBolls -= BollWeight[k][l][m] * abscissionRatio;
      BurrWeightGreenBolls -= BurrWeight[k][l][m] * abscissionRatio;
      BollWeight[k][l][m] -= BollWeight[k][l][m] * abscissionRatio;
      BurrWeight[k][l][m] -= BurrWeight[k][l][m] * abscissionRatio;
      FruitFraction[k][l][m] -= FruitFraction[k][l][m] * abscissionRatio;
//
//     If FruitFraction[k][l][m] is less than 0.001 make it zero, update SeedNitrogen,
//  BurrNitrogen, CumPlantNLoss, PixInPlants, CottonWeightGreenBolls, BurrWeightGreenBolls, GreenBollsLost,
//  BollWeight[k][l][m], BurrWeight[k][l][m], and assign 4 to FruitingCode.
//
      if ( FruitFraction[k][l][m] <= 0.001 ) 
	  {
          FruitingCode[k][l][m] = 4;
          SeedNitrogen -= BollWeight[k][l][m] * (1 - gin1)  * SeedNConc;
          BurrNitrogen -= BurrWeight[k][l][m] * BurrNConc;
          CumPlantNLoss += BollWeight[k][l][m] * (1 - gin1)  * SeedNConc;
          CumPlantNLoss += BurrWeight[k][l][m] * BurrNConc;
          PixInPlants -= (BollWeight[k][l][m] + BurrWeight[k][l][m]) * pixcon;
          FruitFraction[k][l][m] = 0;
          CottonWeightGreenBolls -= BollWeight[k][l][m];
          BurrWeightGreenBolls -= BurrWeight[k][l][m];
          GreenBollsLost += BollWeight[k][l][m] + BurrWeight[k][l][m];
          BollWeight[k][l][m] = 0;
          BurrWeight[k][l][m] = 0;
      }
}
////////////////////////
void AdjustAbscission()
//     This function adjusts the abscission of squares and bolls, to make their 
//  numbers approximate the recorded plant map data. It is called from function
//  FruitingSitesAbscission() when plant adjustment is needed.
//
//     The following global variables are referenced here:
//        FruitingCode, FruitFraction, AdjGreenBollAbsc, ginp, Gintot, NumFruitBranches,
//        NumNodes, NumOpenBolls, NumVegBranches, AdjSquareAbsc.
//
{
      int jx[2]; // if jx[0] = 1 squares are adjusted, if jx[1] = 1 green bolls, 
	  for (int i = 0; i < 2; i++)
	       jx[i] = 0;
//
      double abscsq; // adjusted rate of abscission of a square.
//     Compute rate of square abscission for adjusting square number
      if (nadj[3] && AdjSquareAbsc > 0) 
	  {
         jx[0] = 1;
         abscsq = AdjSquareAbsc;
      }
//     Compute all green bolls susceptible to shedding
      int nbrch; // fruiting branch number.
      int nnid; // node number on fruiting branch.
      double abscgb; // adjusted rate of abscission of a green boll.
      if (nadj[4] && AdjGreenBollAbsc > 0) 
	  {
         jx[1] = 1;
         double gbolnum = 0; // number of bolls susceptible to shedding per plant.
         for (int k = 0; k < NumVegBranches; k++)
		 {
            nbrch = NumFruitBranches[k];
            for (int l = 0; l < nbrch; l++)
			{
               nnid = NumNodes[k][l];
               for (int m = 0; m < nnid; m++)
			   {
                  if ( FruitingCode[k][l][m] == 7 )
                      gbolnum += FruitFraction[k][l][m];
			   }
			}
		 }
//     Compute rate of young boll abscission for adjusting green boll number
	     if (gbolnum > 0) 
            abscgb = AdjGreenBollAbsc * NumGreenBolls / gbolnum;
	     else
            abscgb = 0;
	  }
//     Compute ginning percentage for computing seed weights.
      double gin1; // ginning percent, used to compute lost nitrogen.
      if (jx[1] == 1) 
	  {
         if ( Gintot > 0 ) 
            gin1 = Gintot;
         else
            gin1 = ginp;
      }
//     Start a loop for all fruiting sites.
      for (int k = 0; k < NumVegBranches; k++)
	  {
         nbrch = NumFruitBranches[k];
         for (int l = 0; l < nbrch; l++)
		 {
            nnid = NumNodes[k][l];
            for (int m = 0; m < nnid; m++)
			{
               if ( jx[0] == 1 && FruitingCode[k][l][m] == 1 ) 
                       AdjustSquareAbscission(k, l, m, abscsq);       
               if ( jx[1] == 1 && FruitingCode[k][l][m] == 7 ) 
                       AdjustYoungBollAbscission(k, l, m, abscgb, gin1);
			}
		}
	 }
}
////////////////////////
void AdjustSquareAbscission(int k, int l, int m, double abscsq)
//     This function adjusts the abscission of squares, to
//  make their numbers approximate the recorded plant map data.
//  It is called from function AdjustAbscission().
//
//     The following global variables are referenced here: SquareNConc.
//
//     The following global variable are set here:
//  BloomWeightLoss, CumPlantNLoss, FruitingCode, FruitFraction, SquareNitrogen, 
//  SquareWeight, TotalSquareWeight.
//
//     The following arguments are used:
//        abscsq - adjusted rate of abscission of a square.
//        k, l, m - site location
//
{
//     If this is a square induce abscission ratio abscsq. Compute changes
//  in dry weight and nitrogen content of squares, and associated losses.
//  Also, compute FruitFraction for this site.
      double wtlos; // weight lost by shedding at this site.
      wtlos = SquareWeight[k][l][m] * abscsq;
      SquareNitrogen -= wtlos * SquareNConc;
      CumPlantNLoss += wtlos * SquareNConc;
      SquareWeight[k][l][m] -= wtlos;
      BloomWeightLoss += wtlos;
      TotalSquareWeight -= wtlos;
      FruitFraction[k][l][m] = FruitFraction[k][l][m] * (1 - abscsq);
//     If FruitFraction is very small, make it 0 and change FruitingCode to 5.
      if (FruitFraction[k][l][m] <= 0.001) 
	  {
          FruitFraction[k][l][m] = 0;
          SquareNitrogen -= SquareWeight[k][l][m] * SquareNConc;
          CumPlantNLoss += SquareWeight[k][l][m] * SquareNConc;
          BloomWeightLoss += SquareWeight[k][l][m];
          TotalSquareWeight -= SquareWeight[k][l][m];
          SquareWeight[k][l][m] = 0;
          FruitingCode[k][l][m] = 5;
      }
//     If FruitFraction is greater than 1, make it 1 and compute associated changes
//  in dry weight and nitrogen and their losses.
      if (FruitFraction[k][l][m] > 1) 
	  {
          wtlos = SquareWeight[k][l][m] * (1 - 1 / FruitFraction[k][l][m]);
          FruitFraction[k][l][m] = 1;
          SquareNitrogen -= wtlos * SquareNConc;
          CumPlantNLoss += wtlos * SquareNConc;
          SquareWeight[k][l][m] -= wtlos;
          BloomWeightLoss += wtlos;
          TotalSquareWeight -= wtlos;
      }
}
/////////////////////
void AdjustYoungBollAbscission(int k, int l, int m, double abscgb, double gin1)
//     This function adjusts the abscission of young green bolls, to
//  make their numbers near the recorded plant map data.
//     It is called from function AdjustAbscission().
//
//     The following global variables are referenced here:
//        BurrNConc, pixcon, SeedNConc.
//
//     The following global variable are set here:
//  BollWeight, BurrNitrogen, BurrWeight, BurrWeightGreenBolls, CottonWeightGreenBolls, 
//  CumPlantNLoss, FruitFraction, GreenBollsLost, 
//        PixInPlants, SeedNitrogen.
//
//     The following arguments are used:
//        abscgb - adjusted rate of abscission of a green boll.
//        gin1 - ginning percent, used to compute lost nitrogen.
//        k, l, m - loop index.
//
{
//     If this is a green boll induce abscission ratio abscgb. Compute
//  changes in dry weight and nitrogen content of seedcotton and burrs,
//  and associated losses. Also compute FruitFraction for this site.
//
      SeedNitrogen -= BollWeight[k][l][m] * abscgb * (1 - gin1) * SeedNConc;
      BurrNitrogen -= BurrWeight[k][l][m] * abscgb * BurrNConc;
      CumPlantNLoss += BollWeight[k][l][m] * abscgb * (1 - gin1) * SeedNConc;
      CumPlantNLoss += BurrWeight[k][l][m] * abscgb * BurrNConc;
      double pixlos; // amount of pix lost by abscission, g per plant.
      pixlos = (BollWeight[k][l][m] + BurrWeight[k][l][m]) * abscgb * pixcon;
      PixInPlants = PixInPlants - pixlos;
      GreenBollsLost += (BollWeight[k][l][m] + BurrWeight[k][l][m]) * abscgb;
      CottonWeightGreenBolls -= BollWeight[k][l][m] * abscgb;
      BurrWeightGreenBolls -= BurrWeight[k][l][m] * abscgb;
      BollWeight[k][l][m] = BollWeight[k][l][m] * (1 - abscgb);
      BurrWeight[k][l][m] = BurrWeight[k][l][m] * (1 - abscgb);
      FruitFraction[k][l][m] = FruitFraction[k][l][m] * (1 - abscgb);
//
      int jx = 2; // if jx = 2 green bolls, and if jx = 3 open bolls.
      AdjustBollAbscission(k, l, m, jx, gin1);
}
/////////////////////
void AdjustBollAbscission(int k, int l, int m, int jx, double gin1)
//     This function adjusts the abscission of squares and bolls, to make their numbers 
//  near the recorded plant map data. It is called from functions 
//  AdjustYoungBollAbscission() and AdjustSetBollAbscission().
//
//     The following global variables are referenced here:
//  BurrNConc, pixcon, SeedNConc.
//
//     The following global variable are set here:
//  BollWeight, BurrNitrogen, BurrWeight, BurrWeightGreenBolls, BurrWeightOpenBolls, 
//  CottonWeightGreenBolls, CottonWeightOpenBolls, CumPlantNLoss, 
//  FruitingCode, FruitFraction, GreenBollsLost, PixInPlants, SeedNitrogen.
//
//     The following arguments are used:
//        gin1 - ginning percent, used to compute lost nitrogen.
//        jx - if jx = 2 these are green bolls, and if jx = 3 open bolls.
//        k, l, m - site index.
//
{
//     The following section is activated for both green and open bolls.
//  If FruitFraction is very small, make it 0 and change FruitingCode to 4.
      if ( FruitFraction[k][l][m] <= 0.001 ) 
	  {
         SeedNitrogen -= BollWeight[k][l][m] * (1 - gin1) * SeedNConc;
         BurrNitrogen -= BurrWeight[k][l][m] * BurrNConc;
         CumPlantNLoss += BollWeight[k][l][m] * (1 - gin1) * SeedNConc;
         CumPlantNLoss += BurrWeight[k][l][m] * BurrNConc;
         PixInPlants -= (BollWeight[k][l][m] + BurrWeight[k][l][m]) * pixcon;
         GreenBollsLost += BollWeight[k][l][m] + BurrWeight[k][l][m];
//
	     if (jx == 2) 
		 {
            CottonWeightGreenBolls -= BollWeight[k][l][m];
            BurrWeightGreenBolls -= BurrWeight[k][l][m];
		 }
	     else if (jx == 3) 
		 {
            CottonWeightOpenBolls -= BollWeight[k][l][m];
            BurrWeightOpenBolls -= BurrWeight[k][l][m];
		 }
//
         FruitFraction[k][l][m] = 0;
         BollWeight[k][l][m] = 0;
         BurrWeight[k][l][m] = 0;
         FruitingCode[k][l][m] = 4;
      }
//     If FruitFraction is greater than 1, make it 1 and compute associated changes
//  in dry weight and nitrogen and their losses.
      if (FruitFraction[k][l][m] > 1) 
	  {
         double bolwlos; // weight of seed cotton lost from a boll.
         double burwlos; // weight of burrs lost from a boll.
         bolwlos = BollWeight[k][l][m] * (1 - 1 / FruitFraction[k][l][m]);
         burwlos = BurrWeight[k][l][m] * (1 - 1 / FruitFraction[k][l][m]);
         FruitFraction[k][l][m] = 1;
         SeedNitrogen -= bolwlos * (1 - gin1) * SeedNConc;
         BurrNitrogen -= burwlos * BurrNConc;
         CumPlantNLoss += bolwlos * (1 - gin1) * SeedNConc;
         CumPlantNLoss += burwlos * BurrNConc;
         PixInPlants -= (bolwlos + burwlos) * pixcon;
         GreenBollsLost += bolwlos + burwlos;
         BollWeight[k][l][m] -= bolwlos;
         BurrWeight[k][l][m] -= burwlos;
//
	     if (jx == 2) 
		 {
            CottonWeightGreenBolls -= bolwlos;
            BurrWeightGreenBolls -= burwlos;
		 }
	     else if (jx == 3) 
		 {
            CottonWeightOpenBolls -= bolwlos;
            BurrWeightOpenBolls -= burwlos;
		 }
      }
}
////////////////////
void ComputeSiteNumbers()
//     This function calculates square, green boll, open boll, and abscised site numbers 
//  (NumSquares, NumGreenBolls, NumOpenBolls, and AbscisedFruitSites, respectively), as 
//  the sums of FruitFraction in all sites with appropriate FruitingCode.
//  It is called from function FruitingSitesAbscission(). 
//
//     The following global variables are referenced here:
//   FruitFraction, FruitingCode, NumFruitBranches, NumFruitSites, NumNodes, NumVegBranches.
//
//     The following global variable are set here:
//   AbscisedFruitSites,  GreenBollsLost, NumGreenBolls, NumOpenBolls, NumSquares, .
//
{
      NumSquares = 0;
      NumGreenBolls = 0;
      NumOpenBolls = 0;
      for (int k = 0; k < NumVegBranches; k++)
	  {
         int nbrch = NumFruitBranches[k];
         for (int l = 0; l < nbrch; l++)
		 {
            int nnid = NumNodes[k][l];
            for (int m = 0; m < nnid; m++)
			{
               if ( FruitingCode[k][l][m] == 1 ) 
			        NumSquares += FruitFraction[k][l][m];
               else if ( FruitingCode[k][l][m] == 7 || FruitingCode[k][l][m] == 2 ) 
                    NumGreenBolls += FruitFraction[k][l][m];
               else if ( FruitingCode[k][l][m] == 3 ) 
                    NumOpenBolls += FruitFraction[k][l][m];
			}
		 }
	  }
//
      AbscisedFruitSites = NumFruitSites - NumSquares - NumGreenBolls - NumOpenBolls;
}
