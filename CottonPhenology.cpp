//  CottonPhenology.cpp
//
//   Functions in this file:
// CottonPhenology()
// PreFruitingNode()
// DaysToFirstSquare()
// CreateFirstSquare()
// AddVegetativeBranch()
// AddFruitingBranch()
// AddFruitingNode()
// FruitingSite{}
// NewBollFormation()
// BollOpening()
//
#include "global.h"
#include "LeafAbscission.h"
#include "FruitAbscission.h"
#include "GeneralFunctions.h"

void PreFruitingNode(double, const double&);
double DaysToFirstSquare(const int&, const int&, const double&);
tuple<double> CreateFirstSquare(double, double);
void AddVegetativeBranch(double, double, double);
void AddFruitingBranch(int, double, double, const double&);
void AddFruitingNode(int, int, double, double, const double& WaterStress);
void FruitingSite(int, int, int, int&, const int&, const double&, const double&);
void NewBollFormation(int, int, int);
void BollOpening(int, int, int, double, const int&);

//   Declaration of file-scope variables:  
   double FibLength;           // fiber length
   double FibStrength;         // fiber strength
   double PhenDelayByNStress;  // phenological delay caused by vegetative nitrogen stress.
//////////////////////////////////////////////////
//      The following documentation describes the implementation of
//  the simulation of the cotton plant phenology and the abscission of
//  leaves and of fruiting sites in Cotton2K.
//
//      The calling sequence is as follows:
//      LeafAbscission() calls PreFruitLeafAbscission(), MainStemLeafAbscission(),
//          FruitNodeLeafAbscission(), DefoliationLeafAbscission()
//          === see file LeafAbscission.cpp
//      FruitingSitesAbscission() calls SiteAbscissionRatio(), SquareAbscission(), BollAbscission(),
//	        AdjustAbscission() and ComputeSiteNumbers()
//          === see file FruitAbscission.cpp
//////////////////////////////////////////////////
tuple<int, double, double> CottonPhenology(const int& Daynum, const int& DayEmerge, int FirstSquare, const double& DayInc, const double& WaterStress, double AbscisedLeafWeight)
//     This is is the main function for simulating events of phenology and abscission
//  in the cotton plant. It is called each day from DailySimulation().
//     CottonPhenology() calls PreFruitingNode(), DaysToFirstSquare(), CreateFirstSquare(),
//  AddVegetativeBranch(), AddFruitingBranch(), AddFruitingNode(), FruitingSite(),
//  LeafAbscission(), FruitingSitesAbscission().
//     The following global variables are referenced here:
//        CarbonStress, DensityFactor, Kday, NStressVeg, NumFruitBranches, 
//        NumNodes, NumVegBranches, PerPlantArea, StemNitrogen, TotalStemWeight, VarPar. 
//     The following global variable are set here:
//        FirstSquare, NumFruitSites.
//
{
//     The following constant parameters are used:
	  double vpheno[8] = { 0.65, -0.83, -1.67, -0.25, -0.75, 10.0, 15.0, 7.10}; 
//     
      static int nwfl = 0; // the node of the most recent white flower. Note: this variable
//  is not used. It is kept for compatibility with previous versions, and may be use in future versions.
      NumFruitSites  = 0;
      double stemNRatio; // the ratio of N to dry matter in the stems.
      stemNRatio = StemNitrogen / TotalStemWeight;
//     Compute the phenological delays:
//     PhenDelayByNStress, the delay caused by nitrogen stress, is assumed to be
//  a function of the vegetative nitrogen stress..
      PhenDelayByNStress = vpheno[0] * (1 - NStressVeg);  // a file-scope variable
      if ( PhenDelayByNStress > 1 )
		   PhenDelayByNStress = 1;
      else if ( PhenDelayByNStress < 0 )
		   PhenDelayByNStress = 0;
//  
      double delayVegByCStress; // delay in formation of new fruiting branches caused by carbon stress.
      delayVegByCStress = VarPar[27] + CarbonStress * (vpheno[3] + vpheno[4] * CarbonStress);
      if ( delayVegByCStress > 1 )
		   delayVegByCStress = 1;
      else if ( delayVegByCStress < 0 ) 
		   delayVegByCStress = 0;
//
      double delayFrtByCStress; // delay in formation of new fruiting sites caused by carbon stress.
      delayFrtByCStress = VarPar[28] + CarbonStress * (vpheno[1] + vpheno[2] * CarbonStress);
      if ( delayFrtByCStress > VarPar[29] )
		   delayFrtByCStress = VarPar[29];
      if ( delayFrtByCStress < 0 ) 
		   delayFrtByCStress = 0;
//
      static double DaysTo1stSqare; // number of days from emergence to 1st square
//      The following section is executed if the first square has not yet been
//  formed. Function DaysToFirstSquare() is called to compute the  number of days
//  to 1st square, and function PreFruitingNode() is called to simulate the 
//  formation of prefruiting nodes.
      if (FirstSquare <= 0)
      {
         DaysTo1stSqare = DaysToFirstSquare(Daynum, DayEmerge, WaterStress);
         PreFruitingNode(stemNRatio, DayInc);
//      When first square is formed, FirstSquare is assigned the day of year.
//  Function CreateFirstSquare() is called for formation of first square.
         if ( Kday >= (int) DaysTo1stSqare ) 
         {
            FirstSquare = Daynum;
            tie(AbscisedLeafWeight) = CreateFirstSquare(stemNRatio, AbscisedLeafWeight);
         }
//      if a first square has not been formed, call LeafAbscission() and exit.
         else
         {
            tie(AbscisedLeafWeight) = LeafAbscission(Daynum, FirstSquare, DayInc, AbscisedLeafWeight);
	        return make_tuple(FirstSquare, 0, AbscisedLeafWeight);
         }
      }
//     The following is executed after the appearance of the first square.
//     If there are only one or two vegetative branches, and if plant
//  population allows it, call AddVegetativeBranch() to decide if a new vegetative
//  branch is to be added. Note that dense plant populations (large
//  PerPlantArea) prevent new vegetative branch formation.
      if ( NumVegBranches == 1 && PerPlantArea >= vpheno[5] )
           AddVegetativeBranch(delayVegByCStress,stemNRatio, DaysTo1stSqare);
      if ( NumVegBranches == 2 && PerPlantArea >= vpheno[6] )
           AddVegetativeBranch(delayVegByCStress,stemNRatio, DaysTo1stSqare);
//     The maximum number of nodes per fruiting branch (nidmax) is
//  affected by plant density. It is computed as a function of DensityFactor.
      int nidmax; // maximum number of nodes per fruiting branch.
      nidmax = (int) (vpheno[7] * DensityFactor + 0.5);
      if ( nidmax > 5 ) 
		   nidmax = 5;
//     Start loop over all existing vegetative branches.
//     Call AddFruitingBranch() to decide if a new node (and a new fruiting
//  branch) is to be added on this stem.
      for (int k = 0; k < NumVegBranches; k++)
	  {
         if ( NumFruitBranches[k] < 30 ) 
            AddFruitingBranch( k, delayVegByCStress, stemNRatio, WaterStress );
//     Loop over all existing fruiting branches, and call AddFruitingNode() to 
//  decide if a new node on this fruiting branch is to be added.
         for (int l = 0; l < NumFruitBranches[k]; l++)
		 {
            if ( NumNodes[k][l] < nidmax ) 
                 AddFruitingNode( k, l, delayFrtByCStress, stemNRatio, WaterStress );
//     Loop over all existing fruiting nodes, and call FruitingSite() to
//  simulate the condition of each fruiting node.
            for (int m = 0; m < NumNodes[k][l]; m++)
               FruitingSite( k, l, m, nwfl, Daynum, DayInc, WaterStress);
		 }
	  }
//     Call FruitingSitesAbscission() to simulate the abscission of fruiting parts.
      double AbscisedFruitSites;
      tie(AbscisedFruitSites) = FruitingSitesAbscission(Daynum, DayInc, WaterStress);
//     Call LeafAbscission() to simulate the abscission of leaves.
      tie(AbscisedLeafWeight) = LeafAbscission(Daynum, FirstSquare, DayInc, AbscisedLeafWeight);
      return make_tuple(FirstSquare, AbscisedFruitSites, AbscisedLeafWeight);
}
//////////////////////////
void PreFruitingNode(double stemNRatio, const double& DayInc)
//     This function checks if a new prefruiting node is to be added, and then sets it. 
//  It is called from function CottonPhenology().
//     The following global variables are referenced here:
//        DayInc, LeafWeightAreaRatio, VarPar.
//     The following global variable are set here:
//        AgeOfPreFruNode, LeafAreaPreFru, LeafNitrogen, LeafWeightPreFru, NumPreFruNodes, 
//        StemNitrogen, TotalLeafWeight, TotalStemWeight.
//     The following argument is used:
//        stemNRatio - the ratio of N to dry matter in the stems.
//
{
//     The following constant parameter is used:
      const double MaxAgePreFrNode = 66; // maximum age of a prefruiting node (constant)
//     When the age of the last prefruiting node exceeds MaxAgePreFrNode,
//  this function is not activated.
      if ( AgeOfPreFruNode[NumPreFruNodes - 1] > MaxAgePreFrNode ) 
		  return;
//      Loop over all existing prefruiting nodes.
//      Increment the age of each prefruiting node in physiological days.
      for (int j = 0; j < NumPreFruNodes; j++)
         AgeOfPreFruNode[j] += DayInc;
//      For the last prefruiting node (if there are less than 9
//  prefruiting nodes): The period (timeToNextPreFruNode) until the formation of the next
//  node is VarPar(31), but it is modified for the first three nodes. If
//  the physiological age of the last prefruiting node is more than timeToNextPreFruNode,
//  form a new prefruiting node - increase NumPreFruNodes, assign the initial
//  average temperature for the new node, and initiate a new leaf on
//  this node.
      if (NumPreFruNodes >= 9 ) 
          return;
      double timeToNextPreFruNode; // time, in physiological days, for the next prefruiting node to be formed.
      timeToNextPreFruNode = VarPar[31]; 
      if ( NumPreFruNodes <= 2 )
		  timeToNextPreFruNode *= VarPar[32];
      else if ( NumPreFruNodes == 3 )
		  timeToNextPreFruNode *= VarPar[33];
//
      if ( AgeOfPreFruNode[NumPreFruNodes-1] >= timeToNextPreFruNode )
	  {
          NumPreFruNodes++;
          LeafAreaPreFru[NumPreFruNodes-1] = VarPar[34];
          LeafWeightPreFru[NumPreFruNodes-1] = LeafAreaPreFru[NumPreFruNodes-1] * LeafWeightAreaRatio;
          TotalLeafWeight += LeafWeightPreFru[NumPreFruNodes-1];
          TotalStemWeight -= LeafWeightPreFru[NumPreFruNodes-1];
          LeafNitrogen += LeafWeightPreFru[NumPreFruNodes-1] * stemNRatio;
          StemNitrogen  -= LeafWeightPreFru[NumPreFruNodes-1] * stemNRatio;
	  }
}
////////////////////////////////////////////////////////////////////////////////
double DaysToFirstSquare(const int& Daynum, const int& DayEmerge, const double& WaterStress)
//     This function computes and returns tsq1, the number of days from emergence to first
//  square. It is called from CottonPhenology().
//     The following global variables are referenced here:
//        AvrgDailyTemp, Kday, NStressVeg, VarPar, WaterStress.
//
{
//     The following constant parameters are used:
      double p1 = 34;
      double p2 = 132.2;
      double p3 = -7;
      double p4 = 0.125;
      double p5 = 0.08;
      double p6 = 0.30;
//     On day of emergence assign initial values to some local variables.
	  static double avtemp; // average temperature from day of emergence.
      static double sumstrs; // cumulative effect of water and N stresses on date of first square.
      if ( Daynum <= DayEmerge ) 
	  {
         avtemp = AvrgDailyTemp;
         sumstrs = 0;
      }
//      Compute the average temperature from the day of emergence to
//  now. Compute the number of days from emergence to first
//  square, as a function of this average temperature. This
//  relationships is derived from data of K. R. Reddy et al.
//  (unpublished), CSRU, for Delta cultivars.
      avtemp = ((Kday-1) * avtemp + AvrgDailyTemp ) / Kday;
      if ( avtemp > p1 ) 
		   avtemp = p1;
      double tsq1 = p2 + avtemp * ( p3 + avtemp * p4 );
      tsq1 = tsq1 * VarPar[30];
//      The cumulative effect of water stress and vegetative N stress
//  is computed. This effect is substacted from TSQ, assuming that the
//  first square appears earlier under water or N stressed conditions.
      sumstrs += p5 * (1 - WaterStress) + p6 * (1 - NStressVeg);
      tsq1 -= sumstrs;
//
      return tsq1;
}
//////////////////////////////////////////////////////////////////////////////
tuple<double> CreateFirstSquare(double stemNRatio, double AbscisedLeafWeight)
//     This function initiates the first square. It is called from function CottonPhenology().
//     The following global variables are referenced here:
//        AvrgDailyTemp, Kday, LeafWeightAreaRatio, pixcon, VarPar.
//     The following global variable are set here:
//        AbscisedLeafWeight, AvrgNodeTemper, CumPlantNLoss, FirstSquare, FruitFraction, 
//        FruitGrowthRatio, FruitingCode, LeafAreaNodes, LeafNitrogen, LeafWeightNodes,
//        NumFruitBranches, NumNodes, PixInPlants, StemNitrogen, TotalLeafWeight, TotalStemWeight.
//     Argument used:
//        stemNRatio - the ratio of N to dry matter in the stems.
//
{
//     FruitFraction and FruitingCode are assigned 1 for the first fruiting site.
      FruitingCode[0][0][0] = 1;
      FruitFraction[0][0][0] = 1;
//     Initialize a new leaf at this position. define its initial weight and area. 
//  VarPar[34] is the initial area of a new leaf. The mass and nitrogen of the new leaf 
//  are substacted from the stem.
      LeafAreaNodes[0][0][0] = VarPar[34];
      LeafWeightNodes[0][0][0] = VarPar[34] * LeafWeightAreaRatio;
      TotalStemWeight -= LeafWeightNodes[0][0][0];
      TotalLeafWeight += LeafWeightNodes[0][0][0];
      LeafNitrogen += LeafWeightNodes[0][0][0] * stemNRatio;
      StemNitrogen  -= LeafWeightNodes[0][0][0] * stemNRatio;
//      Define the initial values of NumFruitBranches, NumNodes, FruitGrowthRatio, and AvrgNodeTemper.
      NumFruitBranches[0] = 1;
      NumNodes[0][0] = 1;
      FruitGrowthRatio = 1;
      AvrgNodeTemper[0][0][0] = AvrgDailyTemp;
//     It is assumed that the cotyledons are dropped at time of first square. compute changes  
//  in AbscisedLeafWeight, TotalLeafWeight, LeafNitrogen, CumPlantNLoss, and PixInPlants caused
//  by the abscission of the cotyledons.
	  double cotylwt = 0.20; // cotylwt is the leaf weight of the cotyledons. 
      AbscisedLeafWeight += cotylwt;
      TotalLeafWeight -= cotylwt;
      CumPlantNLoss  += cotylwt * LeafNitrogen / TotalLeafWeight;
      LeafNitrogen -= cotylwt * LeafNitrogen / TotalLeafWeight;
      PixInPlants -= cotylwt * pixcon;
      return make_tuple(AbscisedLeafWeight);
}
///////////////////////////////////////////////////////////////////////////////////
void AddVegetativeBranch(double delayVegByCStress, double stemNRatio, double DaysTo1stSqare)
//     This function decides whether a new vegetative branch is to be added, and
//  then forms it. It is called from CottonPhenology().
//     The following global variables are referenced here:
//        AgeOfSite, LeafWeightAreaRatio, AvrgDailyTemp, VarPar.
//     The following global variable are set here:
//        AvrgNodeTemper, FruitFraction, FruitingCode, LeafAreaMainStem, LeafAreaNodes, 
//        LeafNitrogen, LeafWeightMainStem, LeafWeightNodes, NumFruitBranches, NumNodes, 
//        NumVegBranches, StemNitrogen, TotalLeafWeight, TotalStemWeight.
//     Arguments used in this function:
//        delayVegByCStress - delay in formation of new fruiting branches caused by carbohydrate stress.
//        stemNRatio - the ratio of N to dry matter in the stems.
//        DaysTo1stSqare - days to 1st square
//
{
//     The following constant parameters are used:
      const double vpvegb[3] = { 13.39, -0.696, 0.012 }; 
//      TimeToNextVegBranch is computed as a function of this average temperature.
      double TimeToNextVegBranch; // time, in physiological days, for the next vegetative branch to be formed.
      TimeToNextVegBranch = vpvegb[0] + AvrgNodeTemper[NumVegBranches-1][0][0] * (vpvegb[1] 
                          + AvrgNodeTemper[NumVegBranches-1][0][0] * vpvegb[2]);
//      Compare the age of the first fruiting site of the last formed
//  vegetative branch with TimeToNextVegBranch plus DaysTo1stSqare and the delays caused by
//  stresses, in order to decide if a new vegetative branch is to be formed.
      if ( AgeOfSite[NumVegBranches-1][0][0] 
           < (TimeToNextVegBranch + delayVegByCStress + PhenDelayByNStress + DaysTo1stSqare) ) 
          return;
//      When a new vegetative branch is formed, increase NumVegBranches by 1.
      NumVegBranches++;
	  if (NumVegBranches > 3) 
      {
	        NumVegBranches = 3;
	        return;
      }
//      Assign 1 to FruitFraction and FruitingCode of the first site of this branch.
      FruitFraction[NumVegBranches-1][0][0] = 1;
      FruitingCode[NumVegBranches-1][0][0] = 1;
//      Add a new leaf to the first site of this branch.
      LeafAreaNodes[NumVegBranches-1][0][0] = VarPar[34];
      LeafWeightNodes[NumVegBranches-1][0][0] = VarPar[34] * LeafWeightAreaRatio;
//      Add a new mainstem leaf to the first node of this branch.
      LeafAreaMainStem[NumVegBranches-1][0]  = VarPar[34];
      LeafWeightMainStem[NumVegBranches-1][0]  = LeafAreaMainStem[NumVegBranches-1][0] * LeafWeightAreaRatio;
//      The initial mass and nitrogen in the new leaves are
//  substracted from the stem.
      TotalStemWeight -= ( LeafWeightNodes[NumVegBranches-1][0][0] +  LeafWeightMainStem[NumVegBranches-1][0] );
      TotalLeafWeight += LeafWeightNodes[NumVegBranches-1][0][0] +  LeafWeightMainStem[NumVegBranches-1][0];
      double addlfn; // nitrogen moved to new leaves from stem.
      addlfn = (LeafWeightNodes[NumVegBranches-1][0][0] + LeafWeightMainStem[NumVegBranches-1][0]) * stemNRatio;
      LeafNitrogen += addlfn;
      StemNitrogen  -= addlfn;
//      Assign the initial value of the average temperature of the first site.
//      Define initial NumFruitBranches and NumNodes for the new vegetative branch.
      AvrgNodeTemper[NumVegBranches-1][0][0] = AvrgDailyTemp;
      NumFruitBranches[NumVegBranches-1] = 1;
      NumNodes[NumVegBranches-1][0] = 1;
}
/////////////////////////////////////////////////////////////////////////////////////////////
void AddFruitingBranch (int k, double delayVegByCStress, double stemNRatio, const double& WaterStress)
//     This function decides if a new fruiting branch is to be added to a vegetative
//  branch, and forms it. It is called from function CottonPhenology().
//     The following global variables are referenced here:
//  AdjAddMSNodesRate, AgeOfSite, DensityFactor, LeafWeightAreaRatio, AvrgDailyTemp, 
//  VarPar, WaterStress.
//     The following global variable are set here:
//        AvrgNodeTemper, DelayNewFruBranch, FruitFraction, FruitingCode, LeafAreaMainStem,
//        LeafAreaNodes, LeafNitrogen, LeafWeightMainStem, LeafWeightNodes, NumFruitBranches, 
//        NumNodes, StemNitrogen, TotalLeafWeight, TotalStemWeight,
//     The following arguments are used in this function:
//        delayVegByCStress - delay in formation of new fruiting branches caused by
//                 carbohydrate stress.
//        k - index of this vegetative branch.
//        stemNRatio - the ratio of N to dm in the stems.
//
{
//     The following constant parameters are used:
	  const double vfrtbr[8] = { 0.8, 0.95, 33.0, 4.461, -0.1912, 0.00265, 1.8, -1.32 };
//      Compute the cumulative delay for the appearance of the next caused by carbohydrate,
//  nitrogen, and water stresses.
      DelayNewFruBranch[k] += delayVegByCStress + vfrtbr[0] * PhenDelayByNStress;
      DelayNewFruBranch[k] += vfrtbr[1] * (1 - WaterStress);
//     Define nbrch and compute TimeToNextFruBranch, the time (in physiological days) needed 
//  for the formation of each successive fruiting branch, as a function of 
//  the average temperature. This function is derived from data of K. R.
//  Reddy, CSRU, adjusted for age expressed in physiological days.  It
//  is different for the main stem (k = 0) than for the other vegetative 
//  branches. TimeToNextFruNode is modified for plant density. Add DelayNewFruBranch to TimeToNextFruNode.
      int nbrch = NumFruitBranches[k] - 1; // index of new fruiting branch on this vegetative branch.
      double tav = AvrgNodeTemper[k][nbrch][0]; // modified average daily temperature.
      if ( tav > vfrtbr[2] ) 
		   tav = vfrtbr[2];
//     TimeToNextFruBranch is the time, in physiological days, for the next fruiting branch to be formed.
      double TimeToNextFruBranch = VarPar[35] + tav * (vfrtbr[3] + tav * (vfrtbr[4] 
                                 + tav * vfrtbr[5]));
      if ( k > 0 )
		   TimeToNextFruBranch =  TimeToNextFruBranch * vfrtbr[6];
      TimeToNextFruBranch = TimeToNextFruBranch * (1 + vfrtbr[7] * (1 - DensityFactor))
                          + DelayNewFruBranch[k];
//     Apply adjustment to ageinc if plant map data indicate it.
      if (Kday > KdayAdjust && Kday <= (KdayAdjust + NumAdjustDays)) 
	  {
         if (nadj[0]) 
         {
            if (AdjAddMSNodesRate == 0)
                TimeToNextFruBranch = 100;
            else
                TimeToNextFruBranch = TimeToNextFruBranch / AdjAddMSNodesRate; 
         }
      }
//     Check if the the age of the last fruiting branch exceeds TimeToNextFruBranch. If so, form the new fruiting branch:
      if ( AgeOfSite[k][nbrch][0] < TimeToNextFruBranch ) 
          return;
//     Increment NumFruitBranches, define newbr, and assign 1 to NumNodes, FruitFraction and FruitingCode.
      NumFruitBranches[k]++;
	  if (NumFruitBranches[k] > 30) 
      {
	        NumFruitBranches[k] = 30;
	        return;
      }
      int newbr; // the index number of the new fruiting branch on this vegetative branch, after a new branch has been added.
      newbr = NumFruitBranches[k] - 1;
      NumNodes[k][newbr] = 1;
      FruitFraction[k][newbr][0] = 1;
      FruitingCode[k][newbr][0] = 1;
//     Initiate new leaves at the first node of the new fruiting branch, and at the
//  corresponding main stem node. The mass and nitrogen in the new leaves is substacted 
//  from the stem.
      LeafAreaNodes[k][newbr][0] = VarPar[34];
      LeafWeightNodes[k][newbr][0] = VarPar[34] * LeafWeightAreaRatio;
      LeafAreaMainStem[k][newbr] = VarPar[34];
      LeafWeightMainStem[k][newbr]  = LeafAreaMainStem[k][newbr] * LeafWeightAreaRatio;
      TotalStemWeight -= (LeafWeightMainStem[k][newbr] + LeafWeightNodes[k][newbr][0]);
      TotalLeafWeight +=  LeafWeightMainStem[k][newbr] +  LeafWeightNodes[k][newbr][0];
// addlfn is the nitrogen added to new leaves from stem.
      double addlfn = ( LeafWeightMainStem[k][newbr] + LeafWeightNodes[k][newbr][0] ) * stemNRatio;
      LeafNitrogen += addlfn;
      StemNitrogen  -= addlfn;
//      Begin computing AvrgNodeTemper of the new node and assign zero to DelayNewFruBranch.
      AvrgNodeTemper[k][newbr][0]  = AvrgDailyTemp;
      DelayNewFruBranch[k] = 0;
}
/////////////////////////////////////////////////////////////////////////////////////////
void AddFruitingNode(int k, int l, double delayFrtByCStress, double stemNRatio, const double& WaterStress)
//     Function AddFruitingNode() decides if a new node is to be added to a fruiting branch, 
//  and forms it. It is called from function CottonPhenology().
//     The following global variables are referenced here:
//  AdjAddSitesRate, AgeOfSite, AvrgDailyTemp, DensityFactor, LeafWeightAreaRatio, 
//  pixdn, VarPar, WaterStress,
//     The following global variable are set here:
//        AvrgNodeTemper, DelayNewNode, FruitFraction, FruitingCode, LeafAreaNodes, LeafWeightNodes, 
//        NumNodes, LeafNitrogen, StemNitrogen, TotalLeafWeight, TotalStemWeight.
//     The following arguments are used in this function:
//        delayFrtByCStress - delay caused by carbohydrate and nitrogen stresses.
//        k, l - indices of this vegetative branch and fruiting branch.
//        stemNRatio - the ratio of n to dm in the stems.
//
{
//     The following constant parameters are used:
	  const double vfrtnod[6] = { 1.32, 0.90, 33.0, 7.6725, -0.3297, 0.004657 };
//      Compute the cumulative delay for the appearance of the next
//  node on the fruiting branch, caused by carbohydrate, nitrogen, and
//  water stresses.
      DelayNewNode[k][l] += (delayFrtByCStress + vfrtnod[0] * PhenDelayByNStress) / pixdn;
      DelayNewNode[k][l] += vfrtnod[1] * (1 - WaterStress);
//     Define nnid, and compute the average temperature of the last 
//  node of this fruiting branch, from the time it was formed.
      int nnid = NumNodes[k][l] - 1; // the number of the last node on this fruiting branche.
      double tav = AvrgNodeTemper[k][l][nnid]; // modified daily average temperature.
      if ( tav > vfrtnod[2] )
		   tav = vfrtnod[2];
//     Compute TimeToNextFruNode, the time (in physiological days) needed for the
//  formation of each successive node on the fruiting branch. This is
//  a function of temperature, derived from data of K. R. Reddy, CSRU,
//  adjusted for age in physiological days. It is modified for plant density.
      double TimeToNextFruNode; // time, in physiological days, for the next node on the fruiting branch to be formed 
      TimeToNextFruNode = VarPar[36] 
                        + tav * (vfrtnod[3] + tav * (vfrtnod[4] + tav * vfrtnod[5]));
      TimeToNextFruNode = TimeToNextFruNode 
                        * (1 + VarPar[37] * (1 - DensityFactor)) + DelayNewNode[k][l];
//     Check if plant adjustment by map data is needed
      if (Kday > KdayAdjust && Kday <= (KdayAdjust + NumAdjustDays)) 
	  {
         if (nadj[2]) 
         {
             if ( AdjAddSitesRate == 0)
                  TimeToNextFruNode = 100;
             else
                  TimeToNextFruNode = TimeToNextFruNode / AdjAddSitesRate; 
         }
      }
//     Check if the the age of the last node on the fruiting branch exceeds TimeToNextFruNode.
//  If so, form the new node:
      if ( AgeOfSite[k][l][nnid] < TimeToNextFruNode )
           return;
//     Increment NumNodes, define newnod, and assign 1 to FruitFraction and FruitingCode.
      NumNodes[k][l]++;
      if ( NumNodes[k][l] > 5 ) 
      {
	     NumNodes[k][l] = 5;
	     return;
      }
      int newnod = nnid + 1; // the number of the new node on this fruiting branche.
      FruitFraction[k][l][newnod] = 1;
      FruitingCode[k][l][newnod] = 1;
//     Initiate a new leaf at the new node. The mass and nitrogen in
//  the new leaf is substacted from the stem.
      LeafAreaNodes[k][l][newnod] = VarPar[34];
      LeafWeightNodes[k][l][newnod] = VarPar[34] * LeafWeightAreaRatio;
      TotalStemWeight -= LeafWeightNodes[k][l][newnod];
      TotalLeafWeight += LeafWeightNodes[k][l][newnod];
      LeafNitrogen += LeafWeightNodes[k][l][newnod] * stemNRatio;
      StemNitrogen  -= LeafWeightNodes[k][l][newnod] * stemNRatio;
//     Begin computing AvrgNodeTemper of the new node, and assign zero to DelayNewNode.
      AvrgNodeTemper[k][l][newnod] = AvrgDailyTemp;
      DelayNewNode[k][l] = 0;
}
//////////////////////////////////////////////////
void FruitingSite(int k, int l, int m, int & NodeRecentWhiteFlower, const int& Daynum, const double& DayInc, const double& WaterStress)
//     Function FruitingSite() simulates the development of each fruiting site. 
//  It is called from function CottonPhenology().
//     The following global variables are referenced here:
//        AvrgDailyTemp, CottonWeightGreenBolls, 
//        DayFirstDef, DayInc, Kday, KdayAdjust, LeafAreaIndex, nadj, 
//        NStressVeg, NumAdjustDays, NumFruitBranches, NStressFruiting, WaterStress, VarPar.
//     The following global variable are set here:
//        AgeOfSite, AgeOfBoll, AvrgNodeTemper, BollWeight, BurrWeight, 
//        FirstBloom, FruitingCode, LeafAge, NumFruitSites, 
//     The following arguments are used in this function:
//        k, l, m - indices of vegetative branch, fruiting branch, and
//                  node on fruiting branch for this site.
//        NodeRecentWhiteFlower - the node of the most recent white flower is computed in this
//                  function, although it is not used in this version.
//
{
//     The following constant parameters are used:
	  const double vfrsite[15] = { 0.60, 0.40, 12.25, 0.40, 33.0, 0.20, 0.04, 0.45, 
                    26.10, 9.0, 0.10, 3.0, 1.129, 0.043, 0.26 };
//      FruitingCode = 0 indicates that this node has not yet been formed.
//  In this case, assign zero to boltmp and return.
      static double boltmp[3][30][5];    // cumulative boll temperature.
      if ( FruitingCode[k][l][m] <= 0 ) 
	  {
         boltmp[k][l][m] = 0;
         return;
      }
//      Assign zero to FibLength and FibStrength before any sites have been formed.
      if ( NumFruitSites <= 0 ) 
	  {
         FibLength = 0;
         FibStrength = 0;
      }
      NumFruitSites++; //      Increment site number.
//     LeafAge(k,l,m) is the age of the leaf at this site. it is updated
//  by adding the physiological age of this day, the effect of water
//  and nitrogen stresses (agefac).
      double agefac; // effect of water and nitrogen stresses on leaf aging.
      agefac = (1 - WaterStress) * vfrsite[0] + (1 - NStressVeg) * vfrsite[1];
      LeafAge[k][l][m] += DayInc + agefac;
//  After the application of defoliation, add the effect of defoliation on leaf age.
      if ( DayFirstDef > 0 && Daynum > DayFirstDef )
                 LeafAge[k][l][m] += VarPar[38];
//     FruitingCode = 3, 4, 5 or 6 indicates that this node has an open boll,
//  or has been completely abscised. Return in this case.
      if ( FruitingCode[k][l][m] >= 3 && FruitingCode[k][l][m] <= 6 )
		  return;
//     Age of node is modified for low minimum temperatures and for high
//  maximum temperatures.
      double tmin = GetFromClim("tmin", Daynum);
      double tmax = GetFromClim("tmax", Daynum);
      double ageinc = DayInc; // daily addition to site age.
//     Adjust leaf aging for low minimum twmperatures.
      if ( tmin < vfrsite[2] ) 
		   ageinc += vfrsite[3] * (vfrsite[2] - tmin);
//     Adjust leaf aging for high maximum twmperatures.
      if ( tmax > vfrsite[4] )
	  {
		  double adjust = vfrsite[6] * (tmax - vfrsite[4]);  // vfrsite [4] = 33.0  [5] 0.20  [6] 0.04
		  if (adjust > vfrsite[5])
			  adjust = vfrsite[5];
		  ageinc -= adjust;
	  }
//
      if ( ageinc < vfrsite[7] ) 
		   ageinc = vfrsite[7];
//     Compute average temperature of this site since formation.
      AvrgNodeTemper[k][l][m] = ( AvrgNodeTemper[k][l][m] * AgeOfSite[k][l][m] + AvrgDailyTemp * ageinc )
                              / ( AgeOfSite[k][l][m] + ageinc );
//     Update the age of this node, AgeOfSite(k,l,m), by adding ageinc.
      AgeOfSite[k][l][m] += ageinc;
//     The following is executed if this node is a square (FruitingCode =  1):
//     If square is old enough, make it a green boll: initialize the
//  computations of average boll temperature (boltmp) and boll age
//  (AgeOfBoll). FruitingCode will now be 7.
      if ( FruitingCode[k][l][m] == 1 ) 
	  {
         if ( AgeOfSite[k][l][m] >= vfrsite[8] ) 
		 {
            boltmp[k][l][m] = AvrgDailyTemp;
            AgeOfBoll[k][l][m] = DayInc;
            FruitingCode[k][l][m] = 7;
            NewBollFormation(k,l,m);
//     If this is the first flower, define FirstBloom.
            if ( CottonWeightGreenBolls > 0 && FirstBloom <= 1 )
				 FirstBloom = Daynum;
//     Determine node of most recent white flower.
            if ( k == 0 && m == 0 )
				 NodeRecentWhiteFlower = max(NodeRecentWhiteFlower,l);
         }
	     return;
      }
//     If there is a boll at this site:
//     Calculate average boll temperature (boltmp), and boll age
//  (AgeOfBoll) which is its physiological age, modified by water stress.
//  If leaf area index is low, dum is calculated as an intermediate
//  variable. It is used to increase boll temperature and to accelerate
//  boll aging when leaf cover is decreased. Boll age is also modified
//  by nitrogen stress (NStressFruiting).
      if ( BollWeight[k][l][m] > 0 ) 
	  {
         double dum; // effect of leaf area index on boll temperature and age.
         if ( LeafAreaIndex <= vfrsite[11] && Kday > 100 )
              dum = vfrsite[12] - vfrsite[13] * LeafAreaIndex;
		 else
			  dum = 1;
         double dagebol; // added physiological age of boll on this day.
         dagebol = DayInc * dum + vfrsite[14] * (1 - WaterStress) + vfrsite[10] * (1 - NStressFruiting);
         boltmp[k][l][m] = (boltmp[k][l][m] * AgeOfBoll[k][l][m] + AvrgDailyTemp * dagebol)
                         / (AgeOfBoll[k][l][m] + dagebol);
         AgeOfBoll[k][l][m] += dagebol;
      }
//     If this node is a young green boll (FruitingCode = 7):
//     Check boll age and after a fixed age convert it to an "old"
//  green boll (FruitingCode = 2).
      if ( FruitingCode[k][l][m] == 7 ) 
	  {
         if ( AgeOfBoll[k][l][m] >= vfrsite[9] )
			  FruitingCode[k][l][m] = 2;
         return;
      }
//     If this node is an older green boll (FruitingCode = 2):
      if ( FruitingCode[k][l][m] == 2 ) 
	      BollOpening(k, l, m, boltmp[k][l][m], Daynum);
}
/////////////////////////
void NewBollFormation(int k, int l, int m)
//     Function NewBollFormation() simulates the formation of a new boll at a 
//   fruiting site. It is called from function FruitingSite().
//
//     The following global variables are referenced here:
//        bPollinSwitch, SquareNConc
//     The following global variable are set here:
//        BloomWeightLoss, BollWeight, BurrNitrogen, BurrWeight, CottonWeightGreenBolls, 
//        BurrWeightGreenBolls, CumPlantNLoss, FruitFraction, FruitingCode, 
//        SeedNitrogen, SquareNitrogen, SquareWeight, TotalSquareWeight.
//     The following arguments are used:
//        k, l, m - indices of vegetative branch, fruiting branch, and
//                  node on fruiting branch for this site.
{
//     The following constant parameters are used:
      const double seedratio = 0.64; // ratio of seeds in seedcotton weight.
      const double vnewboll[2] = { 0.31, 0.02 };
//     If bPollinSwitch is false accumulate number of blooms to be dropped,
//  and define FruitingCode as 6.
      if (! bPollinSwitch ) 
	  {
         FruitingCode[k][l][m] = 6;
         FruitFraction[k][l][m] = 0;
         BloomWeightLoss += SquareWeight[k][l][m];
         SquareWeight[k][l][m] = 0;
         return;
	  }
//     The initial weight of the new boll (BollWeight) and new burr (BurrWeight)
//  will be a fraction of the square weight, and the rest will be added
//  to BloomWeightLoss. 80% of the initial weight will be in the burr.
//     The nitrogen in the square is partitioned in the same proportions. The nitrogen
//  that was in the square is transferred to the burrs. Update CottonWeightGreenBolls, 
//  BurrWeightGreenBolls and TotalSquareWeight. assign zero to SquareWeight at this site.
      double bolinit; // initial weight of boll after flowering.
      bolinit = vnewboll[0] * SquareWeight[k][l][m];
      BollWeight[k][l][m] = 0.2 * bolinit;
      BurrWeight[k][l][m] = bolinit - BollWeight[k][l][m];
      BloomWeightLoss += SquareWeight[k][l][m] - bolinit;
//
      double sqr1n; // the nitrogen content of one square before flowering.
      sqr1n = SquareNConc * SquareWeight[k][l][m];
      SquareNitrogen -= sqr1n;
      CumPlantNLoss += sqr1n * (1 - vnewboll[0]);
      sqr1n = sqr1n * vnewboll[0];
//
      double seed1n; // the nitrogen content of seeds in a new boll on flowering.
      seed1n = BollWeight[k][l][m] * seedratio * vnewboll[1];
      if ( seed1n > sqr1n ) 
		   seed1n = sqr1n;
      SeedNitrogen += seed1n;
      BurrNitrogen += sqr1n - seed1n;
//
      CottonWeightGreenBolls += BollWeight[k][l][m];
      BurrWeightGreenBolls += BurrWeight[k][l][m];
      TotalSquareWeight -= SquareWeight[k][l][m];
      SquareWeight[k][l][m] = 0;
}
/////////////////////////
void BollOpening(int k, int l, int m, double tmpboll, const int& Daynum)
//     Function BollOpening() simulates the transition of each fruiting site 
//  from green to dehissed (open) boll. It is called from FruitingSite().
//     The following global variables are referenced here:
//        AgeOfBoll, BollWeight, BurrWeight, DayFirstDef, FruitFraction, 
//        LeafAreaIndex, PlantPopulation, VarPar 
//     The following global variable are set here:
//        BurrWeightGreenBolls, BurrWeightOpenBolls, CottonWeightGreenBolls, CottonWeightOpenBolls,
//        FibLength, FruitingCode, FibStrength, ginp, Gintot, NumOpenBolls, LintYield.
//     The following arguments are used in this function:
//        k, l, m - indices of vegetative branch, fruiting branch, and
//                  node on fruiting branch for this site.
//        tmpboll - average temperature of this boll.
//
{
//     The following constant parameters are used:
      double ddpar1 = 1;
	  double ddpar2 = 0.8; // constant parameters for computing fdhslai.
	  double vboldhs[11] = { 30.0, 41.189, -1.6057, 0.020743, 70.0,
                      0.994, 56.603, -2.921, 0.059, 1.219, 0.0065 };
//     Assign atn as the average boll temperature (tmpboll), and check that it
//  is not higher than a maximum value.
	  double atn; // modified average temperature of this boll.
      atn = tmpboll;
      if ( atn > vboldhs[0] )
		   atn = vboldhs[0];
//     Compute dehiss as a function of boll temperature. 
      double dehiss; // days from flowering to boll opening.
      dehiss = VarPar[39] + atn * (vboldhs[1] + atn * (vboldhs[2] + atn * vboldhs[3]));
      dehiss = dehiss * VarPar[40];
      if ( dehiss > vboldhs[4] )
		   dehiss = vboldhs[4];
//     Dehiss is decreased after a defoliation.
      if ( DayFirstDef > 0 && Daynum > DayFirstDef )
           dehiss = dehiss * pow(vboldhs[5], (Daynum - DayFirstDef));
//     If leaf area index is less than dpar1, decrease dehiss.
      if ( LeafAreaIndex < ddpar1 ) 
	  {
         double fdhslai; // effect of small lai on dehiss
         fdhslai = ddpar2 + LeafAreaIndex * (1 - ddpar2) / ddpar1;
         if ( fdhslai < 0 ) 
			  fdhslai = 0;
         if ( fdhslai > 1 )
			  fdhslai = 1;
         dehiss = dehiss * fdhslai;
      }
      if ( AgeOfBoll[k][l][m] < dehiss ) 
          return;
//     If green boll is old enough (AgeOfBoll greater than dehiss), make
//  it an open boll, set FruitingCode to 3, and update CottonWeightOpenBolls, BurrWeightOpenBolls,
//  CottonWeightGreenBolls, BurrWeightGreenBolls.
      FruitingCode[k][l][m] = 3;
      CottonWeightOpenBolls  += BollWeight[k][l][m];
      BurrWeightOpenBolls  += BurrWeight[k][l][m];
      CottonWeightGreenBolls -= BollWeight[k][l][m];
      BurrWeightGreenBolls -= BurrWeight[k][l][m];
//     Compute the ginning percentage as a function of boll temperature.
//     Compute the average ginning percentage of all the bolls opened
//  until now (Gintot).
      ginp = (VarPar[41] - VarPar[42] * atn) / 100;
      Gintot = (Gintot * NumOpenBolls + ginp * FruitFraction[k][l][m]) /
               (NumOpenBolls + FruitFraction[k][l][m]);
//     Cumulative lint yield (LintYield) is computed in kg per ha.
      LintYield += ginp * BollWeight[k][l][m] * PlantPopulation * .001;
//     Note: computation of fiber properties is as in GOSSYM, it is
//  not used in COTTON2K, and it has not been tested. It is included here
//  for compatibility, and it may be developed in future versions.
//     fsx (fiber strength in g / tex at 1/8 inch) is computed, and
//  averaged (as FibStrength) for all open bolls.
//     flx (fiber length in inches, 2.5% span) is computed, and
//  averaged (as FibLength) for all open bolls.
      double flx; // fiber length (inches, 2.5% span) of this boll.
      double fsx; // fiber strength (g / tex at 1/8 inch) of this boll.
      fsx = vboldhs[6] + atn * (vboldhs[7] + vboldhs[8] * atn);
      flx = vboldhs[9] - vboldhs[10] * atn;
      FibStrength = (FibStrength * NumOpenBolls + fsx * FruitFraction[k][l][m]) / (NumOpenBolls + FruitFraction[k][l][m]);
      FibLength = (FibLength * NumOpenBolls + flx * FruitFraction[k][l][m]) / (NumOpenBolls + FruitFraction[k][l][m]);
//     Update the number of open bolls per plant (nopen).
      NumOpenBolls += FruitFraction[k][l][m];
}
