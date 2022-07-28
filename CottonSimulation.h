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
// F
double form(double c0, double d0, double g0);
void FruitingSite(int k, int l, int m, int &NodeRecentWhiteFlower);
void FruitingSitesAbscission();
void FruitNodeLeafAbscission(int k, int l, int m, double droplf);
// G, H
double GetTargetStress();
void GoBack();
double TemperatureOnLeafGrowthRate(double t);
// I
void InitiateLateralRoots();
// L
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
void NitrogenFlow(int nn, double q01[], double q1[], double dd[], double nit[],
                  double nur[]);
void NitrogenUptake(int l, int k, double reqnc);
// P
void PotentialLeafGrowth();
void PredictEmergence(int hour);
void PredictDripIrrigation(double TargetStress);
void PredictSurfaceIrrigation(double TargetStress);
void PreFruitingNode(double stemNRatio);
void PreFruitLeafAbscission(double droplf);
double PsiOnTranspiration(double PsiAverage);
// R
int ReadSoilHydraulicData();
void RootCultivation(int j);
void RootSummation();
// S
double SensibleHeatTransfer(double tsf, double tenviron, double PlantHeight,
                            double wndcanp);
double SiteAbscissionRatio(int k, int l, int m, int lt);
void sitecode();
double SoilAirOnRootGrowth(double psislk, double poreSpace, double vh2olk);
double SoilMechanicResistance(int l, int k);
double SoilNitrateOnRootGrowth(double vno3lk);
void SoilNitrogen();
void SoilNitrogenAverage();
void SoilNitrogenBal();
void SoilSurfaceBalance(int ihr, int k, double ess, double rlzero, double rss,
                        double sf, double hsg, double &so, double &so2,
                        double &so3, double thet, double tm, double tv);
void SoilSum();
void SoilWaterDataInp();
double SoilWaterOnRootGrowth(double psislk);
void SortArray(int size, double InData[], int indexk[], int indexl[],
               int indexm[]);
void SquareAbscission(int k, int l, int m, double abscissionRatio);
// T
double SoilTemperatureEffect(double tt);
double ThermalCondSoil(double q0, double t0, int l0);
// U, V
void UreaHydrolysis(int l, int k);
// W
void WaterBalance(double q1[], double qx[], double dd[], int nn);
void WaterFlux(double q1[], double psi1[], double dd[], double qr1[],
               double qs1[], double pp1[], int nn, int iv, int ll,
               long numiter);
double SoilWaterEffect(int l, int k, double xx);
void WriteStateVariables(bool bAdjusting);
