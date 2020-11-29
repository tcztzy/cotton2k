// CottonSimulation.h : Defines the functions of the application.
//
#pragma once
#include "global.h"
#include <fstream>  // Necessary for file I/O
#include <algorithm>
#include <tuple>
using namespace std;
#ifdef _DEBUG
#define new DEBUG_NEW
#endif
//
//  definition of functions
//  =======================
// A 
//    void   AdjustSetBollAbscission(int k, int l, int m, double abscob, double gin1); 
    void   ApplyFertilizer();   
    double AveragePsi();      
// B
	void   bollsize();
// C
    void   CanopyBalance (int ihr, int k, double etp1, double rlzero, double rsv, 
		   double c2, double sf, double so, double thet, double tm, double &tv);
    void   CapillaryFlow(const int& DayStart);
    double CellDistance(int l, int k, int l0, int k0);
    void   ColumnShading();      
    void   ComputeIrrigation(const string& ProfileName); 
    void   cotplt (int nmap, const string& ProfileName, const string& Date);
// D
    void   DailyOutput(const string& ProfileName, const string& Date);
    tuple<string>   DataOutput(const string& ProfileName, const string& Date, const int& DayStart);
    void   Denitrification(int l, int k);
    double Drain();
    void   DripFlow(double WaterApplied);       
// E
    void   EnergyBalance (int ihr, int k, bool bMulchon, double ess, double etp1);
    void   ExtraNitrogenAllocation();
// F
// G, H
    void   GetNitrogenStress();
    double GetTargetStress();
    void   GravityFlow(double WaterApplied);       
    void   HeatBalance( int nn );
// I
    void   InitiateLateralRoots();
// L
    void   LateralRootGrowthLeft(int l);
    void   LateralRootGrowthRight(int l);
// M
    void   MineralizeNitrogen(int l, int k, const int& DayStart);
    void   MulchSurfaceBalance(int ihr, int k, double rlsp, double rls5, double rsm,
           double sf, double hsgp, double hsgm, double so, double thet, double &tm, double tv);
// N
    void   Nitrification(int l, int k, double DepthOfLayer);
    void   NitrogenAllocation();
    void   NitrogenFlow(int nn, double q01[], double q1[], double dd[], double nit[], double nur[]);
    void   NitrogenRequirement (const string& ProfileName);
    void   NitrogenSupply (const string& ProfileName);
    void   NitrogenUptake(int l, int k, double reqnc);
    void   NitrogenUptakeRequirement();
// O
    void   output1(const string& ProfileName, const string& Date);
    void   output2(const string& ProfileName);
    void   output3(const string& ProfileName);
    void   output4(const string& ProfileName);
    void   output5(const string& ProfileName);
    void   output6(const string& ProfileName);
    void   output7(const string& ProfileName, const int& DayStart);
    void   OutputForSoilMaps(int irec, int igo, int nday, const string& ProfileName);
    void   outputplt(const string& ProfileName);
    void   OutputPredictedIrrigation(double AppliedWater, double TargetStress, const string& ProfileName); 
// P
    double PetioleNitrateN();
    void   PlantNitrogen(const string& ProfileName);     
	void   PlantNitrogenBal(const string& ProfileName);  
    void   PlantNitrogenContent(); 
    void   PredictEmergence(int hour, const string& ProfileName);
    void   PredictDripIrrigation(double TargetStress); 
    void   PredictSurfaceIrrigation(double TargetStress);
    double PsiOnTranspiration (double PsiAverage);
// R
    void   RedistRootNewGrowth(int l, int k, double adwr1);
    void   RootAging(int l, int k);
    void   RootCultivation(int j);
    void   RootDeath(int l, int k);
    void   RootImpedance();
    void   RootsCapableOfUptake();
    void   RootSummation(const string& ProfileName);
// S
    double SensibleHeatTransfer(double tsf, double tenviron, double PlantHeight, double wndcanp);
    void   sitecode();
    double SoilAirOnRootGrowth(double psislk, double poreSpace, double vh2olk);
    void   SoilHeatFlux( double dlt, int iv, int nn, int layer, int n0 );
    void   SoilInit();
    double SoilMechanicResistance(int l, int k);
    bool   SoilMulchBalance (int ihr, int k, double rlzero, double rsm, double rss, double sf,
	       double &so, double &so2, double &so3, double thet, double &tm, double tv, double wndcanp);
    double SoilNitrateOnRootGrowth(double vno3lk);
    void   SoilNitrogen(const int& DayStart);
	void   SoilNitrogenAverage(const string& ProfileName);
	void   SoilNitrogenBal(const string& ProfileName);   
    void   SoilProcedures(const string& ProfileName, const int& DayStart);
    void   SoilSurfaceBalance (int ihr, int k, double ess, double rlzero, double rss, double sf, 
		   double hsg, double &so, double &so2, double &so3, double thet, double tm, double tv);
    void   SoilSum();           
    double SoilTemOnRootGrowth(double t);
    void   SoilTemperature(const string& ProfileName, const int& DayStart);
    void   SoilTemperatureInit(int &jt1, int &jt2, const string& ProfileName, const int& DayStart);
    void   SoilWaterDataInp();
    double SoilWaterOnRootGrowth(double psislk);
// T
    void   TapRootGrowth();
    double SoilTemperatureEffect (double tt );
    double ThermalCondSoil (double q0, double t0, int l0 );
// U, V
	void   UreaHydrolysis(int l, int k);
// W
    void   WaterBalance ( double q1[], double qx[], double dd[], int nn );
    void   WaterFlux( double q1[], double psi1[], double dd[], double qr1[],
		   double qs1[], double pp1[], int nn, int iv, int ll, long numiter );
    void   WaterTable();        // WATERTBL
    void   WaterUptake();       // UPTAKE
    double SoilWaterEffect (int l, int k, double xx);
	void   WriteLine22(ofstream &File22, double i00, double i01, double i02);
