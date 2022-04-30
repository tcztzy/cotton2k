// CottonSimulation.h : Defines the functions of the application.
//
#pragma once

#include "global.h"
//
//  definition of functions
//  =======================
// A
void ActualFruitGrowth();
void ActualLeafGrowth();
void AddFruitingBranch(int k, double delayVegByCStress, double stemNRatio);
void AddFruitingNode(int k, int l, double delayFrtByCStress, double stemNRatio);
double AddPlantHeight(double denf2);
void AddVegetativeBranch(double delayVegByCStress, double stemNRatio,
                         double DaysTo1stSqare);
void AdjustAbscission();
void AdjustBollAbscission(int k, int l, int m, int jx, double gin1);
//    void   AdjustSetBollAbscission(int k, int l, int m, double abscob, double
//    gin1);
void AdjustSquareAbscission(int k, int l, int m, double abscsq);
void AdjustYoungBollAbscission(int k, int l, int m, double abscgb, double gin1);
double AveragePsi();
// B
void BollAbscission(int k, int l, int m, double abscissionRatio, double gin1);
void BollOpening(int k, int l, int m, double tmpboll);
void bollsize();
// C
void CanopyBalance(int ihr, int k, double etp1, double rlzero, double rsv,
                   double c2, double sf, double so, double thet, double tm,
                   double &tv);
void CapillaryFlow();
void CheckDryMatterBal();
void ComputeActualRootGrowth(double sumpdr);
void ComputeIrrigation();
void ComputeSiteNumbers();
void cotplt(int nmap);
void CottonPhenology();
void CreateFirstSquare(double stemNRatio);
// D
double DaysToFirstSquare();
void Defoliate();
void DefoliationLeafAbscission();
void Denitrification(int l, int k);
double Drain();
void DryMatterBalance(double &cdstem, double &cdleaf, double &cdpet,
                      double &cdroot);
// E
void ExtraNitrogenAllocation();
// F
double form(double c0, double d0, double g0);
void FruitingSite(int k, int l, int m, int &NodeRecentWhiteFlower);
void FruitingSitesAbscission();
void FruitNodeLeafAbscission(int k, int l, int m, double droplf);
// G, H
void GetNitrogenStress();
double GetTargetStress();
void GoBack();
double TemperatureOnLeafGrowthRate(double t);
void HeatBalance(int nn);
// I
void InitiateLateralRoots();
// L
void LateralRootGrowthLeft(int l);
void LateralRootGrowthRight(int l);
void LeafAbscission();
// M
void MainStemLeafAbscission(int k, int l, double droplf);
void MineralizeNitrogen(int l, int k);
void MulchSurfaceBalance(int ihr, int k, double rlsp, double rls5, double rsm,
                         double sf, double hsgp, double hsgm, double so,
                         double thet, double &tm, double tv);
// N
void NewBollFormation(int k, int l, int m);
void Nitrification(int l, int k, double DepthOfLayer);
void NitrogenAllocation();
void NitrogenFlow(int nn, double q01[], double q1[], double dd[], double nit[],
                  double nur[]);
void NitrogenRequirement();
void NitrogenSupply();
void NitrogenUptake(int l, int k, double reqnc);
void NitrogenUptakeRequirement();
// P
double PetioleNitrateN();
void PlantNitrogen();
void PlantNitrogenBal();
void PlantNitrogenContent();
void PotentialFruitGrowth();
void PotentialLeafGrowth();
double PotentialStemGrowth(double stemnew);
void PredictEmergence(int hour);
void PredictDripIrrigation(double TargetStress);
void PredictSurfaceIrrigation(double TargetStress);
void PreFruitingNode(double stemNRatio);
void PreFruitLeafAbscission(double droplf);
double PsiOnTranspiration(double PsiAverage);
// R
int ReadSoilHydraulicData();
void RedistRootNewGrowth(int l, int k, double adwr1);
void RootAging(int l, int k);
void RootCultivation(int j);
void RootDeath(int l, int k);
void RootImpedance();
void RootsCapableOfUptake();
void RootSummation();
// S
double SensibleHeatTransfer(double tsf, double tenviron, double PlantHeight,
                            double wndcanp);
double SiteAbscissionRatio(int k, int l, int m, int lt);
void sitecode();
double SoilAirOnRootGrowth(double psislk, double poreSpace, double vh2olk);
void SoilHeatFlux(double dlt, int iv, int nn, int layer, int n0);
double SoilMechanicResistance(int l, int k);
bool SoilMulchBalance(int ihr, int k, double rlzero, double rsm, double rss,
                      double sf, double &so, double &so2, double &so3,
                      double thet, double &tm, double tv, double wndcanp);
double SoilNitrateOnRootGrowth(double vno3lk);
void SoilNitrogen();
void SoilNitrogenAverage();
void SoilNitrogenBal();
void SoilSurfaceBalance(int ihr, int k, double ess, double rlzero, double rss,
                        double sf, double hsg, double &so, double &so2,
                        double &so3, double thet, double tm, double tv);
void SoilSum();
double SoilTemOnRootGrowth(double t);
void SoilWaterDataInp();
double SoilWaterOnRootGrowth(double psislk);
void SortArray(int size, double InData[], int indexk[], int indexl[],
               int indexm[]);
void SquareAbscission(int k, int l, int m, double abscissionRatio);
// T
void TapRootGrowth();
double TemperatureOnFruitGrowthRate(double t);
double SoilTemperatureEffect(double tt);
double ThermalCondSoil(double q0, double t0, int l0);
// U, V
void UreaHydrolysis(int l, int k);
// W
void WaterBalance(double q1[], double qx[], double dd[], int nn);
void WaterFlux(double q1[], double psi1[], double dd[], double qr1[],
               double qs1[], double pp1[], int nn, int iv, int ll,
               long numiter);
void WaterUptake();  // UPTAKE
double SoilWaterEffect(int l, int k, double xx);
void WriteStateVariables(bool bAdjusting);