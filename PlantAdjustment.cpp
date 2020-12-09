// File   PlantAdjustment.cpp
//
//   List of functions in this file:
//       WriteStateVariables()
//       PlantAdjustments()
//       GoBack()
//
#include "global.h"
#include "GeneralFunctions.h"

tuple<string, int, int, double> GoBack(const string&, const int&, int, double);

///////////////////////////////////////////////////////////////////////////////////////////////////
void WriteStateVariables(bool bAdjusting, const string& Date, const int& Daynum, const int& NumLayersWithRoots, const double& PlantHeight)
//     This function stores all state or rate variables, needed for output, or for rerunning
//  plant adjustments, in the structure Scratch21. It is called from DailySimulation(),
//  DoAdjustments(), and DailyOutput().
//     Each record is an array cell, starting from day of start of simulation. 
//     If the argument bAdjusting is true, the variables needed for plant adjustment rerun
//  are stored. The variables needed for output are always stored.
//
{
//     Variables needed for output:
         Scratch21[DayOfSimulation-1].kday = Kday;
//
         Scratch21[DayOfSimulation-1].abscisedFruitSites = AbscisedFruitSites;
         Scratch21[DayOfSimulation-1].averageSoilPsi = AverageSoilPsi;
         Scratch21[DayOfSimulation-1].avrgDailyTemp = AvrgDailyTemp;
         Scratch21[DayOfSimulation-1].burrNConc = BurrNConc;
         Scratch21[DayOfSimulation-1].burrWeightOpenBolls = BurrWeightOpenBolls;
         Scratch21[DayOfSimulation-1].carbonStress = CarbonStress;
         Scratch21[DayOfSimulation-1].cottonWeightGreenBolls = CottonWeightGreenBolls;
         Scratch21[DayOfSimulation-1].cottonWeightOpenBolls = CottonWeightOpenBolls;
         Scratch21[DayOfSimulation-1].cumEvaporation = CumEvaporation;
         Scratch21[DayOfSimulation-1].cumFertilizerN = CumFertilizerN;
         Scratch21[DayOfSimulation-1].cumNetPhotosynth = CumNetPhotosynth;
         Scratch21[DayOfSimulation-1].cumNitrogenUptake = CumNitrogenUptake;
         Scratch21[DayOfSimulation-1].cumTranspiration = CumTranspiration;
         Scratch21[DayOfSimulation-1].cumWaterAdded = CumWaterAdded;
         Scratch21[DayOfSimulation-1].cumWaterDrained = CumWaterDrained;
         Scratch21[DayOfSimulation-1].deadwt = AbscisedLeafWeight + BloomWeightLoss + GreenBollsLost + RootWeightLoss;
         Scratch21[DayOfSimulation-1].date = Date;
         Scratch21[DayOfSimulation-1].daynum = Daynum;
         Scratch21[DayOfSimulation-1].dayTimeTemp = DayTimeTemp;
         Scratch21[DayOfSimulation-1].gbw = CottonWeightGreenBolls + BurrWeightGreenBolls;
//     h2obal is computed as the water balance in mm. It should always be zero.
//  The "positive" amount is the initial water in the soil slab, plus water
//  added by rain and irrigation, and also water added from the water-table.
//  The "negative" is the present total soil water in the soil slab, and cumulative
//  amounts lost by transpiration, evaporation and drainage.
         Scratch21[DayOfSimulation-1].h2obal = InitialTotalSoilWater + CumWaterAdded 
                                             + addwtbl - TotalSoilWater - CumTranspiration 
                                             - CumEvaporation - CumWaterDrained;
         Scratch21[DayOfSimulation-1].leafAreaIndex = LeafAreaIndex;
         Scratch21[DayOfSimulation-1].leafNConc = LeafNConc;
         Scratch21[DayOfSimulation-1].lightIntercept = LightIntercept;
         Scratch21[DayOfSimulation-1].lintYield = LintYield;
         Scratch21[DayOfSimulation-1].lwpMin = LwpMin;
         Scratch21[DayOfSimulation-1].mainStemNodes = MainStemNodes;
         Scratch21[DayOfSimulation-1].mineralizedOrganicN = MineralizedOrganicN;
         Scratch21[DayOfSimulation-1].netPhotosynthesis = NetPhotosynthesis;
         Scratch21[DayOfSimulation-1].nightTimeTemp = NightTimeTemp;
         Scratch21[DayOfSimulation-1].nitrogenStress = NitrogenStress;
         Scratch21[DayOfSimulation-1].nStressFruiting = NStressFruiting;
         Scratch21[DayOfSimulation-1].nStressVeg = NStressVeg;
         Scratch21[DayOfSimulation-1].numGreenBolls = NumGreenBolls;
         Scratch21[DayOfSimulation-1].numOpenBolls = NumOpenBolls;
         Scratch21[DayOfSimulation-1].numSquares = NumSquares;
         Scratch21[DayOfSimulation-1].petioleNConc = PetioleNConc;
         Scratch21[DayOfSimulation-1].petioleNO3NConc = PetioleNO3NConc;
         Scratch21[DayOfSimulation-1].plantHeight = PlantHeight;
         Scratch21[DayOfSimulation-1].plantWeight = PlantWeight;
         Scratch21[DayOfSimulation-1].rad = GetFromClim("rad", Daynum);
         Scratch21[DayOfSimulation-1].rain = GetFromClim("rain", Daynum);
         Scratch21[DayOfSimulation-1].reserveC = ReserveC;
		 Scratch21[DayOfSimulation-1].rn = Rn;
         Scratch21[DayOfSimulation-1].rootNConc = RootNConc;
         Scratch21[DayOfSimulation-1].seedNConc = SeedNConc;
         Scratch21[DayOfSimulation-1].soilNitrogenLoss = SoilNitrogenLoss;
         Scratch21[DayOfSimulation-1].stemNConc = StemNConc;
         Scratch21[DayOfSimulation-1].sumNO3N90 = SumNO3N90;
         Scratch21[DayOfSimulation-1].tmax = GetFromClim("tmax", Daynum);
         Scratch21[DayOfSimulation-1].tmin = GetFromClim("tmin", Daynum);
         Scratch21[DayOfSimulation-1].totalLeafWeight = TotalLeafWeight;
         Scratch21[DayOfSimulation-1].totalPetioleWeight = TotalPetioleWeight;
         Scratch21[DayOfSimulation-1].totalRootWeight = TotalRootWeight;
         Scratch21[DayOfSimulation-1].totalSquareWeight = TotalSquareWeight;
         Scratch21[DayOfSimulation-1].totalSoilWater = TotalSoilWater;
         Scratch21[DayOfSimulation-1].totalStemWeight = TotalStemWeight;
         Scratch21[DayOfSimulation-1].waterStress = WaterStress;
         Scratch21[DayOfSimulation-1].waterStressStem = WaterStressStem;
         Scratch21[DayOfSimulation-1].wind = GetFromClim("wind", Daynum);
//
         for (int l = 0; l < maxl; l++)
            for (int k = 0; k < maxk; k++)
            {
                Scratch21[DayOfSimulation-1].rootWtCapblUptake[l][k] = RootWtCapblUptake[l][k];
                Scratch21[DayOfSimulation-1].soilPsi[l][k] = SoilPsi[l][k];
                Scratch21[DayOfSimulation-1].soilTempDailyAvrg[l][k] = SoilTempDailyAvrg[l][k];
                Scratch21[DayOfSimulation-1].volNh4NContent[l][k] = VolNh4NContent[l][k];
                Scratch21[DayOfSimulation-1].volNo3NContent[l][k] = VolNo3NContent[l][k];
                Scratch21[DayOfSimulation-1].volWaterContent[l][k] = VolWaterContent[l][k];
            }
//
     if (bAdjusting)
     {
//     Variables needed for plant adjustment only:
         Scratch21[DayOfSimulation-1].abscisedLeafWeight = AbscisedLeafWeight;
         Scratch21[DayOfSimulation-1].averageLwp = AverageLwp;
         Scratch21[DayOfSimulation-1].averageLwpMin = AverageLwpMin;
         Scratch21[DayOfSimulation-1].bloomWeightLoss = BloomWeightLoss;
         Scratch21[DayOfSimulation-1].burrNitrogen = BurrNitrogen;
         Scratch21[DayOfSimulation-1].burrWeightGreenBolls = BurrWeightGreenBolls;

         Scratch21[DayOfSimulation-1].cumPlantNLoss = CumPlantNLoss;
         Scratch21[DayOfSimulation-1].deepSoilTemperature = DeepSoilTemperature;
         Scratch21[DayOfSimulation-1].extraCarbon = ExtraCarbon;

         Scratch21[DayOfSimulation-1].firstBloom = FirstBloom;
         Scratch21[DayOfSimulation-1].firstSquare = FirstSquare;
         Scratch21[DayOfSimulation-1].fruitGrowthRatio = FruitGrowthRatio;
         Scratch21[DayOfSimulation-1].gintot = Gintot;
         Scratch21[DayOfSimulation-1].greenBollsLost = GreenBollsLost;
         Scratch21[DayOfSimulation-1].lastTaprootLayer = LastTaprootLayer;
         Scratch21[DayOfSimulation-1].leafNitrogen = LeafNitrogen;

         Scratch21[DayOfSimulation-1].nStressRoots = NStressRoots;
         Scratch21[DayOfSimulation-1].numAbscisedLeaves = NumAbscisedLeaves;
         Scratch21[DayOfSimulation-1].numFruitSites = NumFruitSites;
         Scratch21[DayOfSimulation-1].numLayersWithRoots = NumLayersWithRoots;
         Scratch21[DayOfSimulation-1].numPreFruNodes = NumPreFruNodes;
         Scratch21[DayOfSimulation-1].numSheddingTags = NumSheddingTags;
         Scratch21[DayOfSimulation-1].numVegBranches = NumVegBranches;

         Scratch21[DayOfSimulation-1].petioleNitrogen = PetioleNitrogen;
         Scratch21[DayOfSimulation-1].PixInPlants = PixInPlants;

         Scratch21[DayOfSimulation-1].rootNitrogen = RootNitrogen;
         Scratch21[DayOfSimulation-1].rootWeightLoss = RootWeightLoss;

         Scratch21[DayOfSimulation-1].seedNitrogen = SeedNitrogen;
         Scratch21[DayOfSimulation-1].squareNConc = SquareNConc;
         Scratch21[DayOfSimulation-1].squareNitrogen = SquareNitrogen;
         Scratch21[DayOfSimulation-1].stemNitrogen = StemNitrogen;
         Scratch21[DayOfSimulation-1].supplyNH4N = SupplyNH4N;
         Scratch21[DayOfSimulation-1].supplyNO3N = SupplyNO3N;

         Scratch21[DayOfSimulation-1].totalRequiredN = TotalRequiredN;
         Scratch21[DayOfSimulation-1].totalSoilNh4N = TotalSoilNh4N;
         Scratch21[DayOfSimulation-1].totalSoilNo3N = TotalSoilNo3N;
         Scratch21[DayOfSimulation-1].totalSoilUreaN = TotalSoilUreaN;

         for (int i = 0; i < 9; i++)
         {
            Scratch21[DayOfSimulation-1].ageOfPreFruNode[i] = AgeOfPreFruNode[i];
            Scratch21[DayOfSimulation-1].leafAreaPreFru[i] = LeafAreaPreFru[i];
            Scratch21[DayOfSimulation-1].leafWeightPreFru[i] = LeafWeightPreFru[i];
            Scratch21[DayOfSimulation-1].petioleWeightPreFru[i] = PetioleWeightPreFru[i];
         }
         for (int i = 0; i < 20; i++)
         {
            Scratch21[DayOfSimulation-1].abscissionLag[i] = AbscissionLag[i];
            Scratch21[DayOfSimulation-1].shedByCarbonStress[i] = ShedByCarbonStress[i];
            Scratch21[DayOfSimulation-1].shedByNitrogenStress[i] = ShedByNitrogenStress[i];
            Scratch21[DayOfSimulation-1].shedByWaterStress[i] = ShedByWaterStress[i];
         }
         for (int k = 0; k < 3; k++)
         {
            Scratch21[DayOfSimulation-1].delayNewFruBranch[k] = DelayNewFruBranch[k];
            Scratch21[DayOfSimulation-1].lwpMinX[k] = LwpMinX[k];
            Scratch21[DayOfSimulation-1].lwpX[k] = LwpX[k];
            Scratch21[DayOfSimulation-1].numFruitBranches[k] = NumFruitBranches[k];
            for (int l = 0; l < 30; l ++)
            {
               Scratch21[DayOfSimulation-1].delayNewNode[k][l] = DelayNewNode[k][l];
               Scratch21[DayOfSimulation-1].leafAreaMainStem[k][l] = LeafAreaMainStem[k][l];
               Scratch21[DayOfSimulation-1].leafWeightMainStem[k][l] = LeafWeightMainStem[k][l];
               Scratch21[DayOfSimulation-1].numNodes[k][l] = NumNodes[k][l];
               Scratch21[DayOfSimulation-1].petioleWeightMainStem[k][l] = PetioleWeightMainStem[k][l];
               Scratch21[DayOfSimulation-1].potGroLeafAreaMainStem[k][l] = PotGroLeafAreaMainStem[k][l];
               Scratch21[DayOfSimulation-1].potGroLeafWeightMainStem[k][l] = PotGroLeafWeightMainStem[k][l];
               Scratch21[DayOfSimulation-1].potGroPetioleWeightMainStem[k][l] = PotGroPetioleWeightMainStem[k][l];
               for (int m = 0; m < 5; m ++)
               {
                  Scratch21[DayOfSimulation-1].ageOfSite[k][l][m] = AgeOfSite[k][l][m];
                  Scratch21[DayOfSimulation-1].ageOfBoll[k][l][m] = AgeOfBoll[k][l][m];
                  Scratch21[DayOfSimulation-1].avrgNodeTemper[k][l][m] = AvrgNodeTemper[k][l][m];
                  Scratch21[DayOfSimulation-1].bollWeight[k][l][m] = BollWeight[k][l][m];
                  Scratch21[DayOfSimulation-1].burrWeight[k][l][m] = BurrWeight[k][l][m];
                  Scratch21[DayOfSimulation-1].fruitFraction[k][l][m] = FruitFraction[k][l][m];
                  Scratch21[DayOfSimulation-1].fruitingCode[k][l][m] = FruitingCode[k][l][m];
                  Scratch21[DayOfSimulation-1].leafAge[k][l][m] = LeafAge[k][l][m];
                  Scratch21[DayOfSimulation-1].leafAreaNodes[k][l][m] = LeafAreaNodes[k][l][m];
                  Scratch21[DayOfSimulation-1].leafWeightNodes[k][l][m] = LeafWeightNodes[k][l][m];
                  Scratch21[DayOfSimulation-1].petioleWeightNodes[k][l][m] = PetioleWeightNodes[k][l][m];
                  Scratch21[DayOfSimulation-1].squareWeight[k][l][m] = SquareWeight[k][l][m];
               }
            }
         }
         for (int l = 0; l < maxl; l++)
         {
            Scratch21[DayOfSimulation-1].lateralRootFlag[l] = LateralRootFlag[l];
            Scratch21[DayOfSimulation-1].rootColNumLeft[l] = RootColNumLeft[l];
            Scratch21[DayOfSimulation-1].rootColNumRight[l] = RootColNumRight[l];
            for (int k = 0; k < maxk; k++)
            {
                if (l == 0)
                {
                   Scratch21[DayOfSimulation-1].foliageTemp[k] = FoliageTemp[k];
                   Scratch21[DayOfSimulation-1].mulchTemp[k] = MulchTemp[k];
                }
                Scratch21[DayOfSimulation-1].rootAge[l][k] = RootAge[l][k];
                Scratch21[DayOfSimulation-1].freshOrganicMatter[l][k] = FreshOrganicMatter[l][k];
                Scratch21[DayOfSimulation-1].humusOrganicMatter[l][k] = HumusOrganicMatter[l][k];
                Scratch21[DayOfSimulation-1].nhum[l][k] = HumusNitrogen[l][k];
                Scratch21[DayOfSimulation-1].soilTemp[l][k] = SoilTemp[l][k];
                Scratch21[DayOfSimulation-1].freshOrganicNitrogen[l][k] = FreshOrganicNitrogen[l][k];
                Scratch21[DayOfSimulation-1].volUreaNContent[l][k] = VolUreaNContent[l][k];
                for (int i = 0; i < 3; i++)
                    Scratch21[DayOfSimulation-1].rootWeight[l][k][i] = RootWeight[l][k][i];
            }
         }
     }
}
//////////////////////////
tuple<string, int, int, double> PlantAdjustments(int i, int jj, const string& ProfileName, const string& Date, const int& daynum, int NumLayersWithRoots, double PlantHeight)
//     This function adjusts plant height and plant fruiting map, when data for such 
//  adjustments are available.
//     This function is called from DoAdjustments(). it calls GoBack().
//     The following global variables are referenced here:
//  MapDataAllSiteNum, MapDataGreenBollNum, MapDataMainStemNodes, 
//  MapDataPlantHeight, MapDataSquareNum, NumAdjustDays, NumFruitBranches, NumFruitSites, 
//  NumGreenBolls, NumOpenBolls, NumPreFruNodes, NumSquares, PlantHeight.
//     The following global variables are set:
//  AdjAddHeightRate, AdjAddMSNodesRate, AdjAddSitesRate, AdjGreenBollAbsc, 
//  AdjSquareAbsc, FirstSquare, nadj.
//     The following arguments are used:
//  i - the number of this adjustment, as read from the *.MAP file.
//  jj - the type of this adjustment.
{
//     Define nadj(jj) as true or false, where jj is:
//  0 - for main stem nodes.
//  1 - for height.
//  2 - for total site number.
//  3 - for square number.
//  4 - for green boll number.
//  5 - for open boll number.
     string date = Date;
     int Daynum = daynum;
     switch (jj)
     {
     case 0:   // adjust the number of main stem nodes
         if (MapDataMainStemNodes[i] <= 0)
                return make_tuple(date, Daynum, NumLayersWithRoots, PlantHeight);
//     Compute targeted number of adjusted number of fruiting branches (mntarget).
//  The difference between the measured and the simulated number of nodes should be more than
//  one node. Otherwise no adjustments are made.
         int mnsim;     // simulated number of main stem node count.
         mnsim = Scratch21[DayOfSimulation - 1].mainStemNodes;
         double mntarget; //  targeted number of main stem nodes.
         if (fabs(MapDataMainStemNodes[i] - mnsim) > 1) 
         {
               mntarget = MapDataMainStemNodes[i];
               nadj[0] = true;
         }
         else
         {
               nadj[0] = false;
	           return make_tuple(date, Daynum, NumLayersWithRoots, PlantHeight);
         }
//     Compute the ratio used to adjust the rate of formation of mainstem nodes.
         int mn0;   // simulated number of main stem nodes, at the start of plant adjustment period.
         mn0 = Scratch21[DayOfSimulation - NumAdjustDays - 1].mainStemNodes;
         AdjAddMSNodesRate = 1;
         if (mnsim != mn0) 
         {
               AdjAddMSNodesRate = (mntarget - mn0) / (double) (mnsim - mn0);
               if (AdjAddMSNodesRate > 0.98 && AdjAddMSNodesRate < 1.02) 
                   AdjAddMSNodesRate = 1;
	           if (AdjAddMSNodesRate < 0)
                   AdjAddMSNodesRate = 0;
         }
//
         if (AdjAddMSNodesRate == 1)
             nadj[0] = false;
         else
         {
             nadj[0] = true;
//     AdjAddMSNodesRate will be used in function AddFruitingBranch()
             ofstream File46(fs::path("output") / (ProfileName + ".F01"), ios::app);
			 File46 << " Apply plant adjustment for main stem nodes to date " << date << endl;
             tie(date, Daynum, NumLayersWithRoots, PlantHeight) = GoBack(date, Daynum, NumLayersWithRoots, PlantHeight);
         }
         return make_tuple(date, Daynum, NumLayersWithRoots, PlantHeight);
//
      case 1: // Plant stem height
         if (MapDataPlantHeight[i] <= 0)  
                return make_tuple(date, Daynum, NumLayersWithRoots, PlantHeight);
//     Compute targeted adjusted plant height, if difference is larger than
//  5% of the average.  Note that adjustment is by 90% of the difference.
         double zsim; //   simulated plant height.
         zsim = PlantHeight;
         double ddif; // the difference between simulated and plant adjustment values.
         ddif = MapDataPlantHeight[i] - PlantHeight;
         double ztarget; //  targeted plant height.
         double pHeight; //  plant height at start of adjustment
//     The difference should be at least 5% of the height, otherwise no adjustments made.
         if (fabs(ddif) > (MapDataPlantHeight[i] + PlantHeight) / 40) 
         {
             ztarget = PlantHeight + 0.9 * ddif;
             nadj[1] = true;
         }
	     else
         {
             nadj[1] = false;
             return make_tuple(date, Daynum, NumLayersWithRoots, PlantHeight);
         }
//     Compute the ratio to adjust rate of growth of plant height.
         AdjAddHeightRate = 1;
         pHeight = Scratch21[DayOfSimulation - NumAdjustDays - 1].plantHeight;
         if (fabs(zsim - pHeight) > 0)
             AdjAddHeightRate = (ztarget - pHeight) / (zsim - pHeight);
         if (AdjAddHeightRate < 0)  //  no negative growth rates for plant height.
             AdjAddHeightRate = 0;
//     If AdjAddHeightRate value is near 1, no adjustments are made.
         if (AdjAddHeightRate > 0.98 && AdjAddHeightRate < 1.02)
             nadj[1] = false;
         else
         {
             nadj[1] = true;
//     AdjAddHeightRate will be used in function AddPlantHeight()
             ofstream File46(fs::path("output") / (ProfileName + ".F01"), ios::app);
		     File46 << " Apply plant adjustment for stem height to date " << date << endl;
             tie(date, Daynum, NumLayersWithRoots, PlantHeight) = GoBack(date, Daynum, NumLayersWithRoots, PlantHeight);
         }
         return make_tuple(date, Daynum, NumLayersWithRoots, PlantHeight);
//
      case 2: // total number of fruiting sites
         if ( MapDataAllSiteNum[i] <= 0)
                return make_tuple(date, Daynum, NumLayersWithRoots, PlantHeight);
//     If first square has not yet been simulated, but map data shows squares,
//  adjust the day of first square.
		 if (FirstSquare <= 0) 
         {
	         FirstSquare = (int) (Daynum - 3 * MapDataAllSiteNum[i]);
             nadj[2] = true;
             return make_tuple(date, Daynum, NumLayersWithRoots, PlantHeight);
         }
         double sitesim;      // simulated number of total fruiting sites.
         sitesim = NumFruitSites;
//     Compute targeted adjusted total site number. 
         double sitarget; // targeted number of total fruiting sites.
         if (fabs(MapDataAllSiteNum[i] - NumFruitSites) > 2) 
         {
               sitarget = MapDataAllSiteNum[i];
               nadj[2] = true;
         }
         else
         {
             nadj[2] = false;
             return make_tuple(date, Daynum, NumLayersWithRoots, PlantHeight);
         }
//     Compute the ratio to adjust the rate of formation of sites on fruiting branches. 
         if (sitarget <= 0) 
	         AdjAddSitesRate = 1;
         int pNumFruitSites; // number of sites at start of adjustment
         pNumFruitSites = Scratch21[DayOfSimulation - NumAdjustDays - 1].numFruitSites;
         if (sitesim != pNumFruitSites) 
         {
             AdjAddSitesRate = (sitarget - pNumFruitSites) / (sitesim - pNumFruitSites);
             if (AdjAddSitesRate > 0.98 && AdjAddSitesRate < 1.02)
                 AdjAddSitesRate = 1;
             if (AdjAddSitesRate < 0) 
                 AdjAddSitesRate = 0;
         }
	     if (AdjAddSitesRate == 1) 
             nadj[2] = false;
         else
         {
             nadj[2] = true;
//     AdjAddSitesRate will be used in function AddFruitingNode()
             ofstream File46(fs::path("output") / (ProfileName + ".F01"), ios::app);
			 File46 << " Apply plant adjustment for total number of sites to date " << date << endl;
             tie(date, Daynum, NumLayersWithRoots, PlantHeight) = GoBack(date, Daynum, NumLayersWithRoots, PlantHeight);
         }
         return make_tuple(date, Daynum, NumLayersWithRoots, PlantHeight);
//
      case 3: // number of squares
         if (MapDataSquareNum[i] <= 0)
            return make_tuple(date, Daynum, NumLayersWithRoots, PlantHeight);
//     Adjust actual square numbers. Amount to add each day is the
//  difference between target and actual simulated value on adjustment date
//  divided by the length of the  adjustment period. 
         if ((NumSquares - MapDataSquareNum[i]) > 1 && NumAdjustDays > 0) 
         {
             nadj[3] = true;
	         AdjSquareAbsc = 1 - pow( (MapDataSquareNum[i] / NumSquares), (1 / (double) NumAdjustDays) );
//     AdjSquareAbsc will be used in function AdjustAbscission()
             ofstream File46(fs::path("output") / (ProfileName + ".F01"), ios::app);
			 File46 << " Apply plant adjustment for number of squares to date " << date << endl;
             tie(date, Daynum, NumLayersWithRoots, PlantHeight) = GoBack(date, Daynum, NumLayersWithRoots, PlantHeight);
         }
         else
             nadj[3] = false;
         return make_tuple(date, Daynum, NumLayersWithRoots, PlantHeight);
//
      case 4: // number of green bolls
         if (MapDataGreenBollNum[i] <= 0)
                return make_tuple(date, Daynum, NumLayersWithRoots, PlantHeight);
//     Adjust actual green boll numbers using the same method as for squares. 
         if ((NumGreenBolls - MapDataGreenBollNum[i]) > 1 && NumAdjustDays > 0) 
         {
             AdjGreenBollAbsc = 1 - pow( (MapDataGreenBollNum[i] / NumGreenBolls), (1 / (double) NumAdjustDays) );
             nadj[4] = true;
//     AdjGreenBollAbsc will be used in function AdjustAbscission()
             ofstream File46(fs::path("output") / (ProfileName + ".F01"), ios::app);
			 File46 << " Apply plant adjustment for number of green bolls to date " << date << endl;
             tie(date, Daynum, NumLayersWithRoots, PlantHeight) = GoBack(date, Daynum, NumLayersWithRoots, PlantHeight);
         }
         else
             nadj[4] = false;
         return make_tuple(date, Daynum, NumLayersWithRoots, PlantHeight);
     }  // end switch
     return make_tuple(date, Daynum, NumLayersWithRoots, PlantHeight);
} 
///////////////////////////////////////////////////////////////////////
tuple<string, int, int, double> GoBack(const string& Date, const int& daynum, int NumLayersWithRoots, double PlantHeight)
//     This function reads state variables retroactively NumAdjustDays days earlier. 
//  Thus, all required global state variables will assume their values
//  at beginning of adjustment period.
//     Note that retroactive adjustment begins at (Kday - NumAdjustDays) days earlier.
//
{
         string date = Date;
         int Daynum = daynum;
	     if (Kday <= NumAdjustDays)
	            return make_tuple(date, Daynum, NumLayersWithRoots, PlantHeight);  // before emergence
         int irec; // the record number in the structure containing the state variables.
         irec = DayOfSimulation - NumAdjustDays - 1;
//     Get state variables from Scratch21 structure for the day defined by irec.
         Kday = Scratch21[irec].kday;
//
         AbscisedFruitSites = Scratch21[irec].abscisedFruitSites;
         AbscisedLeafWeight = Scratch21[irec].abscisedLeafWeight;
         AverageLwp = Scratch21[irec].averageLwp;
         AverageLwpMin = Scratch21[irec].averageLwpMin;
         AverageSoilPsi = Scratch21[irec].averageSoilPsi;
         AvrgDailyTemp = Scratch21[irec].avrgDailyTemp;
         BloomWeightLoss = Scratch21[irec].bloomWeightLoss;
         BurrNConc = Scratch21[irec].burrNConc;
         BurrNitrogen = Scratch21[irec].burrNitrogen;
         BurrWeightGreenBolls = Scratch21[irec].burrWeightGreenBolls;
         BurrWeightOpenBolls = Scratch21[irec].burrWeightOpenBolls;
         CarbonStress = Scratch21[irec].carbonStress;
         CottonWeightGreenBolls = Scratch21[irec].cottonWeightGreenBolls;
         CottonWeightOpenBolls = Scratch21[irec].cottonWeightOpenBolls;
         CumEvaporation = Scratch21[irec].cumEvaporation;
         CumFertilizerN = Scratch21[irec].cumFertilizerN;
         CumNetPhotosynth = Scratch21[irec].cumNetPhotosynth;
         CumNitrogenUptake = Scratch21[irec].cumNitrogenUptake;
         CumPlantNLoss = Scratch21[irec].cumPlantNLoss;
         CumTranspiration = Scratch21[irec].cumTranspiration;
         CumWaterAdded = Scratch21[irec].cumWaterAdded;
         CumWaterDrained = Scratch21[irec].cumWaterDrained;
         date = Scratch21[irec].date;
         Daynum = Scratch21[irec].daynum;
         DayTimeTemp = Scratch21[irec].dayTimeTemp;
         DeepSoilTemperature = Scratch21[irec].deepSoilTemperature;
         ExtraCarbon = Scratch21[irec].extraCarbon;
         FirstBloom = Scratch21[irec].firstBloom;
         FirstSquare = Scratch21[irec].firstSquare;
         FruitGrowthRatio = Scratch21[irec].fruitGrowthRatio;
         Gintot = Scratch21[irec].gintot;
         GreenBollsLost = Scratch21[irec].greenBollsLost;
         LastTaprootLayer = Scratch21[irec].lastTaprootLayer;
         LeafAreaIndex = Scratch21[irec].leafAreaIndex;
         LeafNConc = Scratch21[irec].leafNConc;
         LeafNitrogen = Scratch21[irec].leafNitrogen;
         LightIntercept = Scratch21[irec].lightIntercept;
         LintYield = Scratch21[irec].lintYield;
         LwpMin = Scratch21[irec].lwpMin;
         MainStemNodes = Scratch21[irec].mainStemNodes;
         MineralizedOrganicN = Scratch21[irec].mineralizedOrganicN;
         NetPhotosynthesis = Scratch21[irec].netPhotosynthesis;
         NightTimeTemp = Scratch21[irec].nightTimeTemp;
         NitrogenStress = Scratch21[irec].nitrogenStress;
         NStressFruiting = Scratch21[irec].nStressFruiting;
         NStressRoots = Scratch21[irec].nStressRoots;
         NStressVeg = Scratch21[irec].nStressVeg;
         NumAbscisedLeaves = Scratch21[irec].numAbscisedLeaves;
         NumFruitSites = Scratch21[irec].numFruitSites;
         NumGreenBolls = Scratch21[irec].numGreenBolls;
         NumLayersWithRoots = Scratch21[irec].numLayersWithRoots;
         NumOpenBolls = Scratch21[irec].numOpenBolls;
         NumPreFruNodes = Scratch21[irec].numPreFruNodes;
         NumSheddingTags = Scratch21[irec].numSheddingTags;
         NumSquares = Scratch21[irec].numSquares;
         NumVegBranches = Scratch21[irec].numVegBranches;
         PetioleNConc = Scratch21[irec].petioleNConc;
         PetioleNitrogen = Scratch21[irec].petioleNitrogen;
         PetioleNO3NConc = Scratch21[irec].petioleNO3NConc;
         PixInPlants = Scratch21[irec].PixInPlants;
         PlantHeight = Scratch21[irec].plantHeight;
         PlantWeight = Scratch21[irec].plantWeight;
         ReserveC = Scratch21[irec].reserveC;
		 Rn = Scratch21[irec].rn;
         RootNConc = Scratch21[irec].rootNConc;
         RootNitrogen = Scratch21[irec].rootNitrogen;
         RootWeightLoss = Scratch21[irec].rootWeightLoss;
         SeedNConc = Scratch21[irec].seedNConc;
         SeedNitrogen = Scratch21[irec].seedNitrogen;
         SoilNitrogenLoss = Scratch21[irec].soilNitrogenLoss;
         SquareNConc = Scratch21[irec].squareNConc;
         SquareNitrogen = Scratch21[irec].squareNitrogen;
         StemNConc = Scratch21[irec].stemNConc;
         StemNitrogen = Scratch21[irec].stemNitrogen;
         SumNO3N90 = Scratch21[irec].sumNO3N90;
         SupplyNO3N = Scratch21[irec].supplyNO3N;
         SupplyNH4N = Scratch21[irec].supplyNH4N;
         TotalLeafWeight = Scratch21[irec].totalLeafWeight;
         TotalPetioleWeight = Scratch21[irec].totalPetioleWeight;
         TotalRequiredN = Scratch21[irec].totalRequiredN;
         TotalRootWeight = Scratch21[irec].totalRootWeight;
         TotalSquareWeight = Scratch21[irec].totalSquareWeight;
         TotalSoilNh4N = Scratch21[irec].totalSoilNh4N;
         TotalSoilNo3N = Scratch21[irec].totalSoilNo3N;
         TotalSoilUreaN = Scratch21[irec].totalSoilUreaN;
         TotalSoilWater = Scratch21[irec].totalSoilWater;
         TotalStemWeight = Scratch21[irec].totalStemWeight;
         WaterStress = Scratch21[irec].waterStress;
         WaterStressStem = Scratch21[irec].waterStressStem;
//
         for (int i = 0; i < 9; i++)
         {
            AgeOfPreFruNode[i] = Scratch21[irec].ageOfPreFruNode[i];
            LeafAreaPreFru[i] = Scratch21[irec].leafAreaPreFru[i];
            LeafWeightPreFru[i] = Scratch21[irec].leafWeightPreFru[i];
            PetioleWeightPreFru[i] = Scratch21[irec].petioleWeightPreFru[i];
         }
         for (int i = 0; i < 20; i++)
         {
            AbscissionLag[i] = Scratch21[irec].abscissionLag[i];
            ShedByCarbonStress[i] = Scratch21[irec].shedByCarbonStress[i];
            ShedByNitrogenStress[i] = Scratch21[irec].shedByNitrogenStress[i];
            ShedByWaterStress[i] = Scratch21[irec].shedByWaterStress[i];
         }
         for (int k = 0; k < 3; k++)
         {
            DelayNewFruBranch[k] = Scratch21[irec].delayNewFruBranch[k];
            LwpMinX[k] = Scratch21[irec].lwpMinX[k];
            LwpX[k] = Scratch21[irec].lwpX[k];
            NumFruitBranches[k] = Scratch21[irec].numFruitBranches[k];
            for (int l = 0; l < 30; l ++)
            {
               DelayNewNode[k][l] = Scratch21[irec].delayNewNode[k][l];
               LeafAreaMainStem[k][l] = Scratch21[irec].leafAreaMainStem[k][l];
               LeafWeightMainStem[k][l] = Scratch21[irec].leafWeightMainStem[k][l];
               NumNodes[k][l] = Scratch21[irec].numNodes[k][l];
               PetioleWeightMainStem[k][l] = Scratch21[irec].petioleWeightMainStem[k][l];
               PotGroLeafAreaMainStem[k][l] = Scratch21[irec].potGroLeafAreaMainStem[k][l];
               PotGroLeafWeightMainStem[k][l] = Scratch21[irec].potGroLeafWeightMainStem[k][l];
               PotGroPetioleWeightMainStem[k][l] = Scratch21[irec].potGroPetioleWeightMainStem[k][l];
               for (int m = 0; m < 5; m ++)
               {
                  AgeOfBoll[k][l][m] = Scratch21[irec].ageOfBoll[k][l][m];
                  AgeOfSite[k][l][m] = Scratch21[irec].ageOfSite[k][l][m];
                  AvrgNodeTemper[k][l][m] = Scratch21[irec].avrgNodeTemper[k][l][m];
                  BollWeight[k][l][m] = Scratch21[irec].bollWeight[k][l][m];
                  BurrWeight[k][l][m] = Scratch21[irec].burrWeight[k][l][m];
                  FruitingCode[k][l][m] = Scratch21[irec].fruitingCode[k][l][m];
                  FruitFraction[k][l][m] = Scratch21[irec].fruitFraction[k][l][m];
                  LeafAge[k][l][m] = Scratch21[irec].leafAge[k][l][m];
                  LeafAreaNodes[k][l][m] = Scratch21[irec].leafAreaNodes[k][l][m];
                  LeafWeightNodes[k][l][m] = Scratch21[irec].leafWeightNodes[k][l][m];
                  PetioleWeightNodes[k][l][m] = Scratch21[irec].petioleWeightNodes[k][l][m];
                  SquareWeight[k][l][m] = Scratch21[irec].squareWeight[k][l][m];
               }
            }
         }
         for (int l = 0; l < maxl; l++)
         {
            LateralRootFlag[l] = Scratch21[irec].lateralRootFlag[l];
            RootColNumLeft[l] = Scratch21[irec].rootColNumLeft[l];
            RootColNumRight[l] = Scratch21[irec].rootColNumRight[l];
            for (int k = 0; k < maxk; k++)
            {
                if (l == 0)
                {
                   FoliageTemp[k] = Scratch21[irec].foliageTemp[k];
                   MulchTemp[k] = Scratch21[irec].mulchTemp[k];
                }
                FreshOrganicMatter[l][k] = Scratch21[irec].freshOrganicMatter[l][k];
                HumusOrganicMatter[l][k] = Scratch21[irec].humusOrganicMatter[l][k];
                FreshOrganicNitrogen[l][k] = Scratch21[irec].freshOrganicNitrogen[l][k];
                HumusNitrogen[l][k] = Scratch21[irec].nhum[l][k];
                RootAge[l][k] = Scratch21[irec].rootAge[l][k];
                RootWtCapblUptake[l][k] = Scratch21[irec].rootWtCapblUptake[l][k];
                SoilTemp[l][k] = Scratch21[irec].soilTemp[l][k];
                VolNh4NContent[l][k] = Scratch21[irec].volNh4NContent[l][k];
                VolNo3NContent[l][k] = Scratch21[irec].volNo3NContent[l][k];
                VolWaterContent[l][k] = Scratch21[irec].volWaterContent[l][k];
                VolUreaNContent[l][k] = Scratch21[irec].volUreaNContent[l][k];
                for (int i = 0; i < 3; i++)
                    RootWeight[l][k][i] = Scratch21[irec].rootWeight[l][k][i];
            }
         }
         return make_tuple(date, Daynum, NumLayersWithRoots, PlantHeight);
}
