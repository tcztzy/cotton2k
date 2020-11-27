//  File SoilTemperature_2.cpp
//
//   List of functions in this file:
//       EnergyBalance()
//       SensibleHeatTransfer()
//       SoilSurfaceBalance()
//       SoilMulchBalance()
//
#include "CottonSimulation.h"
#include "GeneralFunctions.h"
#include <math.h>

//////////////////////////
void EnergyBalance (int ihr, int k, bool bMulchon, double ess, double etp1)
//     This function solves the energy balance equations at the soil surface, and at the
//  foliage / atmosphere interface. It computes the resulting temperatures of the soil
//  surface, plastic mulch (if exists) and the plant canopy.
//     Units for all energy fluxes are: cal cm-2 sec-1.
//     It is called from SoilTemperature(), on each hourly time step and for each soil column.
//     It calls functions clearskyemiss(), VaporPressure(), SensibleHeatTransfer(),
//  SoilMulchBalance(), SoilSurfaceBalance() and CanopyBalance.()
//
//     The following arguments are used in this function:
//       ihr - the time of day in hours.
//       k - soil column number.
//       bMulchon - is true if this column is covered with plastic mulch, false if not.
//       ess - evaporation from surface of a soil column (mm / sec).
//       etp1 - actual transpiration rate (mm / sec).
//     The following global variables are referenced here:
//       AirTemp, albedo, CloudCoverRatio, CloudTypeCorr, FieldCapacity, MulchTemp, MulchTranSW,
//       PlantHeight, Radiation, RelativeHumidity, rracol, SitePar, thad, VolWaterContent, WindSpeed.
//    The following global variables are set here:
//       bEnd, SoilTemp, FoliageTemp, MulchTemp.
{
//     Constants used:
	  const double stefa1 = 1.38e-12;  // Stefan-Boltsman constant.
	  const double wndfac = 0.60;      // Ratio of wind speed under partial canopy cover.
	  const double cswint = 0.75;      // proportion of short wave radiation (on fully 
                                       // shaded soil surface)intercepted by the canopy.
//     Set initial values
      double sf = 1- rracol[k];             // fraction of shaded soil area
      double thet = AirTemp[ihr] + 273.161; // air temperature, K
      double so =  SoilTemp[0][k];          // soil surface temperature, K
      double so2 = SoilTemp[1][k];          // 2nd soil layer temperature, K
      double so3 = SoilTemp[2][k];          // 3rd soil layer temperature, K
//     Compute soil surface albedo (based on Horton and Chung, 1991):
	  double ag;                    // albedo of the soil surface
      if (VolWaterContent[0][k] <= thad[0])
         ag = SitePar[15];
      else if (VolWaterContent[0][k] >= FieldCapacity[0])
         ag = SitePar[16];
      else
         ag = SitePar[16] +  ( SitePar[15] - SitePar[16] ) * (FieldCapacity[0] - VolWaterContent[0][k]) / (FieldCapacity[0] - thad[0]);
//  ****   SHORT WAVE RADIATION ENERGY BALANCE   ****
//     Division by 41880 (= 698 * 60) converts from Joules per sq m to
// langley (= calories per sq cm) Or: from Watt per sq m to langley per sec.
//     Modify incoming short wave radiation to mulched soil surface.
      double rzero = Radiation[ihr] / 41880;      //short wave (global) radiation (ly / sec).
      double rss0 = rzero * ( 1 - sf * cswint );  //global radiation after passing through canopy
      double tm;       // temperature of mulch (K)
      double rsup;     // global radiation reflected up to the vegetation
      double rsm = 0;  // global radiation absorbed by mulch
      double rss;      // global radiation absorbed by soil surface
      if (bMulchon)
	  {
	     tm = MulchTemp[k];
//     Assume all non transfered radiation is absorbed by mulch
         rss = rss0 * MulchTranSW * (1 - ag);               // absorbed by soil surface
	     rsm = rss0 * (1 - MulchTranSW)                     // absorbed by mulch
             + rss0 * MulchTranSW * ag * (1 - MulchTranSW); // reflected up from soil surface and absorbed by mulch
         rsup = rss0 * MulchTranSW * ag * MulchTranSW;      // reflected up from soil surface through mulch
	  }
      else
	  {
	     tm = 0;
         rss = rss0 * (1 - ag);  // absorbed by soil surface
         rsup = rss0 * ag;       // reflected up from soil surface
	  }
//   ****   LONG WAVE RADIATION EMITTED FROM SKY    ****
      double vp = 0.01 * RelativeHumidity[ihr] * VaporPressure(AirTemp[ihr]); //air vapor pressure, KPa.
      double ea0 = clearskyemiss(vp, thet);          //sky emissivity from clear portions of the sky.
      double rlzero;  // incoming long wave radiation (ly / sec).
      rlzero = (ea0 * (1 - CloudCoverRatio[ihr]) + CloudCoverRatio[ihr]) * stefa1 * pow(thet, 4)
		     - CloudTypeCorr[ihr] / 41880; // CloudTypeCorr converted from W m-2 to ly sec-1.
//
//     Set initial values of canopy temperature and air temperature in canopy. 
      double tv; // temperature of plant foliage (K)
      double tafk; // temperature (K) of air inside the canopy.
      if (sf < 0.05) // no vegetation
	  {
		  tv = thet;          
		  tafk = thet;        
	  }
//     Wind velocity in canopy is converted to cm / s.
      double wndhr;     // wind speed in cm /sec
      wndhr = WindSpeed[ihr] * 100;     
      double rocp; // air density * specific heat at constant pressure = 0.24 * 2 * 1013 / 5740
                   // divided by tafk.
      double c2;   // multiplier for sensible heat transfer (at plant surface).
      double  rsv; // global radiation absorbed by the vegetation
      if (sf >= 0.05)  // a shaded soil column
	  {
         tv = FoliageTemp[k];  // vegetation temperature
//     Short wave radiation intercepted by the canopy:
         rsv = rzero * ( 1 - albedo[ihr]) * sf * cswint   //  from above
              + rsup * ( 1 - albedo[ihr]) * sf * cswint;  //  reflected from soil surface
//     Air temperature inside canopy is the average of soil, air, and plant temperatures, 
//  weighted by 0.1, 0.3, and 0.6, respectively. In case of mulch, mulch temperature replaces 
//  soil temperature.
         if (bMulchon)
            tafk = (1 - sf) * thet + sf * (0.1 * tm + 0.3 * thet + 0.6 * tv);
	     else
            tafk = (1 - sf) * thet + sf * (0.1 * so + 0.3 * thet + 0.6 * tv);
//
//     Call SensibleHeatTransfer() to compute sensible heat transfer coefficient. Factor 2.2 
//  for sensible heat transfer: 2 sides of leaf plus stems and petioles.
         double varcc;   // sensible heat transfer coefficient for soil
         varcc = SensibleHeatTransfer(tv, tafk, PlantHeight, wndhr); // canopy to air
	     if(bEnd)
			 return;
         rocp = 0.08471 / tafk;
	     c2 = 2.2 * sf * rocp * varcc;
      }
      int menit = 0;      // counter of iteration number
      double soold = so;  // previous value of soil surface temperature
      double tvold = tv;  // previous value of vegetation temperature
      double tmold = tm;  // previous value of temperature of mulch (K)
//     Starting iterations for soil, mulch and canopy energy balance
      do
      {
          soold = so;
          double wndcanp = (1 - sf * (1 - wndfac)) * wndhr;  // estimated wind speed under canopy
          if (bMulchon)
	      {
//     This section executed for mulched columns: call SoilMulchBalance() for soil / mulch interface.
             tmold = tm;
             bEnd = SoilMulchBalance (ihr, k, rlzero, rsm, rss, sf, so, so2, so3, thet, tm, tv, wndcanp);
	         if(bEnd)
			    return;
	      }
	      else
	      {
//     This section executed for non-mulched columns
//     Call SensibleHeatTransfer() to compute sensible heat transfer for soil surface to air
             tafk = (1 - sf) * thet + sf * (0.1 * so + 0.3 * thet + 0.6 * tv);
             double varc; // sensible heat transfer coefficientS for soil
             varc = SensibleHeatTransfer(so, tafk, 0, wndcanp);
	         if(bEnd)
			     return;
             rocp = 0.08471 / tafk;
             double hsg = rocp * varc; // multiplier for computing sensible heat transfer soil to air. 
//     Call SoilSurfaceBalance() for energy balance in soil surface / air interface.
		     SoilSurfaceBalance (ihr, k, ess, rlzero, rss, sf, hsg, so, so2, so3, thet, 0, tv);
	         if(bEnd)
			     return;
	      }
//
          if (sf >= 0.05)
	      {
//     This section executed for shaded columns only.
             tvold = tv; 
//     Compute canopy energy balance for shaded columns
             CanopyBalance (ihr, k, etp1, rlzero, rsv, c2, sf, so, thet, tm, tv);
//     Increment the number of iterations.
             menit++;
             if (menit > 10)
		     {
//     The following is used to reduce fluctuations.
	           so = 0.5 * (so + soold); 
	           if (bMulchon)
			       tm = 0.5 * (tm + tmold);
	           tv = 0.5 * (tv + tvold);
		     }
             if (menit > 30)
		     {
//     If more than 30 iterations are needed - stop simulation.
			     bEnd = true;
	             return;
		     }
          }  // end if sf
      } while (fabs(tv - tvold) > 0.05 || fabs(so - soold) > 0.05 || fabs(tm - tmold) > 0.05);
//     After convergence - set global variables for the following temperatures:
      if (sf >= 0.05)
          FoliageTemp[k] = tv;
      SoilTemp[0][k] = so;
      SoilTemp[1][k] = so2;
      SoilTemp[2][k] = so3;
      if (bMulchon)
         MulchTemp[k] = tm;
      else
         MulchTemp[k] = 0;
}
/////////////////////////////////////////////////////////////////////////////
double SensibleHeatTransfer(double tsf, double tenviron, double height, double wndcanp)
//     This function computes the sensible heat transfer coefficient, using the friction 
//  potential (shear) temperature (thstar), and the surface friction (shear) velocity (ustar) 
//  at the atmospheric boundary. It is called three times from EnergyBalance(): for canopy,
//  soil surface, or mulch boundaries with their environment.
//
//     The following arguments are referenced in this function:
//       tenviron - temperature (K) of the environment - air at 200 cm height
//                  for columns with no canopy, or tafk when canopy is present .
//       tsf - surface temperature, K, of soil, mulch or canopy.
//       wndcanp - wind speed in the canopy (if present), cm s-1.
//       height - canopy height, cm, or zero for soil surface and mulch boundaries.
//     The return value computed in this function:
//       sensibleHeatTransfer - raw sensible heat transfer coefficient
//     Global variable set:   bEnd
{
//	    Constant values used:
      const double grav = 980;       // acceleration due to gravity (980 cm sec-2).
      const double s40 = 0.13;       // calibration constant.
      const double s42 = 0.63;       // calibration constant.
      const double stmin = 5;        // minimal value of ustar.
      const double vonkar = 0.40;    // Von-Karman constant (0.40).
      const double zalit1  = 0.0962; // parameter .....
//     Wind velocity not allowed to be less than 100 cm s-1.
      double u = wndcanp; // wind speed at 200 cm height, cm / s.
      if (u < 100)
          u = 100;
//     Assign initial values to z0 and gtop, and set dt.
      double z0 = s40 * height; // surface roughness parameter, cm.
      if ( z0 < 1 )
		   z0 = 1;
      double gtop = log( (200 - s42 * height) / z0); // logarithm of ratio of height of measurement to surface roughness parameter.
      double dt = tsf - tenviron; // temperature difference.
//     Set approximate initial values for ustar and thstar (to reduce iterations).
      double thstar; // friction potential (shear) temperature.
      double ustar;  // Surface friction (shear) velocity (cm sec-1).
	  if (dt >= 0)
	  {
         ustar = 1.873 + 0.570172 * dt + .07438568 * u;
         thstar = - 0.05573 * dt + 1.987 / u - 6.657 * dt / u;
	  }
      else
	  {
         ustar = -4.4017 + 1.067 * dt + 0.25957 * u - 0.001683 * dt * u;
	     if (ustar < 5)
             ustar = 5;
         thstar = -0.0096 - 0.1149 * dt + 0.0000377 * u + 0.0002367 * dt * u;
         if (thstar < 0.03)
             thstar = 0.03;
      }
      double tbot1 = tsf; // surface temperature corrected for friction (shear) potential temperature.
      long mtest = 0;     // count of iterations.
      double g1; // temporary derived variable
      double ug1chk; // previous value of ug1.
      double ug1; // ratio of ustar to g1.
      double ug1res; // residual value of ug1.
//     Start iterations.
      do
      {
         double thekz; // previous value of thstar.
         double uchek; // previous value of ustar.
         ug1chk = 0; // previous value of UG1.
         if (mtest > 0)
	     {
//     Assign values to tbot1, uchek, thekz, and ug1chk.
            tbot1 = tsf + zalit1 * thstar * pow( (ustar*z0/15), 0.45) / vonkar;
            uchek = ustar;
            thekz = thstar;
	        if (g1 != 0)
			    ug1chk = ustar / g1;
	     }
//     Compute air temperature at 1 cm,, and compute zl and lstar.
         double zl; // nondimensional height.
         if ( fabs(thstar) < 1e-30 )
             zl = 0;
         else
	     {
             double thetmn = (tenviron + tbot1) * 0.5;// mean temperature (K) of air and surface.
             double lstar = (thetmn * ustar * ustar) / (vonkar * grav * thstar);
             zl = (200 - s42 * height)  / lstar;
//     Added to decrease fluctuations in zl:
             if (zl <-5)
			     zl = -5;
             if (zl > 0.5)
			     zl = 0.5;
	     }
//     Compute g1u, and g2.
	     double g1u, g2; // temporary derived variables.
         if ( zl > 0 )
	     {
             g1u = -4.7 * zl;
             g2 = -6.35135 * zl;
             if (g2 < -1)
			     g2 = -1;
	     }
         else
	     {
             double tmp1 = pow( (1 - 15 * zl), 0.25); // intermediate variable.
             g1u = 2 * log( (1 + tmp1) / 2 ) + log( (1 + tmp1 * tmp1) / 2 )
                 - 2 * atan(tmp1 + 1.5708);
             g2 = 2 * log( (1 + sqrt(1 - 9 * zl)) / 2 );
	     }
         if (g2 > gtop)
		     g2 = gtop;
//     Compute ustar and check for minimum value.
         ustar = vonkar * u / ( gtop - g1u );
         if (ustar < stmin)
		     ustar = stmin;
//     Compute g1 and thstar.
         g1 = 0.74 * ( gtop - g2) + zalit1 * pow((ustar * z0 / 0.15), 0.45);
         thstar = -dt * vonkar / g1;
//     If more than 30 iterations, reduce fluctuations
         if (mtest > 30)
	     {
	         thstar = (thstar + thekz) * 0.5;
             ustar = (ustar + uchek) * 0.5;
	     }
         mtest++;
//     Stop simulation if no convergence after 100 iterations.
         if ( mtest > 100 )
	     {
             string msg = " Infinite loop in SensibleHeatTransfer(). Abnormal stop!! \n";
             char C1[12];
             sprintf(C1, "%10.3g", tenviron);
             msg += " tenviron = " + (string) C1 + "\n";
             sprintf(C1, "%10.3g", tsf);
             msg += " tsf      = " + (string) C1 + "\n";
             sprintf(C1, "%10.3g", height);
             msg += " PlantHeight = " + (string) C1 + "\n";
             sprintf(C1, "%10.3g", u);
             msg += " u = " + (string) C1 + "\n";
		 throw Cotton2KException(msg);
         }
//     Compute ug1 and  ug1res to check convergence
         ug1 = ustar / g1; 
         if ( fabs(ug1chk) <= 1.e-30 )
             ug1res = fabs(ug1);
         else
             ug1res = fabs( (ug1chk - ug1) / ug1chk );
//     If ug1 did not converge, go to next iteration.
      }  while (fabs(ug1 - ug1chk) > 0.05 && ug1res > 0.01);
//     No more iterations. Compute sensibleHeatTransfer.
      double sensibleHeatTransfer = ustar * vonkar / g1;
      return sensibleHeatTransfer;
}
/////////////////////////////////////////
void SoilSurfaceBalance (int ihr, int k, double ess, double rlzero, double rss, double sf,
		double hsg, double &so, double &so2, double &so3, double thet, double tm, double tv)
//     This function is called from EnergyBalance(). It calls function ThermalCondSoil().
//     It solves the energy balance equations at the soil surface, and
//  computes the resulting temperature of the soil surface.
//     Units for all energy fluxes are: cal cm-2 sec-1.
//
//     The following global variables are referenced here:
//       Daynum, dl, MulchTranLW, VolWaterContent.
//     The following global variable is set here:      bEnd.
//     The following arguments are set in this function:
//       so - temperature of soil surface.
//       so2 - temperature of soil 2nd layer
//       so3 - temperature of soil 3rd layer
//     The following arguments are referenced in this function:
//       ess - evaporation from soil surface (mm / sec).
//       hsg - multiplier for computing sensible heat transfer from soil to air
//             or from soil to mulch if this is a mulched column.
//       ihr - the time in hours.
//       k - soil column number.
//       rlzero - incoming long wave radiation (ly / sec).
//       rss - global radiation absorbed by soil surface
//       sf - fraction of shaded soil area
//       thet - air temperature (K).
//       tm - temperature of plastic mulch (K). When tm = 0 there is no plastic mulch.
//       tv - temperature of plant canopy (K).
//
{
// Constants:
      const double ef = 0.95; // emissivity of the foliage surface
      const double eg = 0.95; // emissivity of the soil surface
      const double stefa1 = 1.38e-12; // Stefan-Boltsman constant.
//     Long wave radiation reaching the soil:
      double rls1; // long wave energy reaching soil surface
      if (sf >= 0.05)                       // shaded column
	     rls1 = (1 - sf) * eg * rlzero               // from sky on unshaded soil surface
              + sf * eg * ef * stefa1 * pow(tv, 4);  // from foliage on shaded soil surface
      else
         rls1 = eg * rlzero;                         // from sky in unshaded column
	  if (tm > 0)     // modify by mulch transmissivity and add lw emitted from mulch
         rls1 = rls1 * MulchTranLW + eg * (1 - MulchTranLW) * stefa1 * pow(tm, 4);
//     rls4 is the multiplier of so**4 for emitted long wave radiation from soil
      double rls4 = eg * stefa1;
      double bbex;      // previous value of bbadjust.
      int mon = 0;      // count of iterations for soil surface energy balance
      double soex = so; //  previous value of so.
//     Start itrations for soil surface enegy balance.
      while (mon < 50)
	  {
//     Compute latent heat flux from soil evaporation: convert from mm sec-1 to
//  cal cm-2 sec-1. Compute derivative of hlat
//     hlat is the energy used for evaporation from soil surface (latent heat)
         double hlat = ( 75.5255 - 0.05752 * so ) * ess;
         double dhlat = -0.05752 * ess; // derivative of hlat
//     Compute the thermal conductivity of layers 1 to 3 by function ThermalCondSoil().
         double rosoil1; // heat conductivity of 1st soil layer in cal / (cm sec deg).
         double rosoil2; // heat conductivity of 2nd soil layer in cal / (cm sec deg).
         double rosoil3; // heat conductivity of 3rd soil layer in cal / (cm sec deg).
         rosoil1 = ThermalCondSoil(VolWaterContent[0][k], so, 1);
         rosoil2 = ThermalCondSoil(VolWaterContent[1][k], so2, 2);
         rosoil3 = ThermalCondSoil(VolWaterContent[2][k], so3, 3);
//     Compute average rosoil between layers 1 to 3,and heat transfer from
//  soil surface to 3rd soil layer.
         double rosoil; // multiplier for heat flux between 1st and 3rd soil layers.
         rosoil = (rosoil1 * dl[0] + rosoil2 * dl[1] + rosoil3 * dl[2])
                 / (dl[0] + dl[1] + dl[2])
                 / (.5 * dl[0] + dl[1] + .5 * dl[2]);
//     bbsoil is the heat energy transfer by conductance from soil surface to soil
	     double bbsoil = rosoil * (so - so3);
//     emtlw is emitted long wave radiation from soil surface
         double emtlw =  rls4 * pow(so, 4);
//     Sensible heat transfer and its derivative
         double dsenheat; // derivative of senheat
         double senheat;  // sensible heat transfer from soil surface
         if (tm > 0)
		 {
            senheat = hsg * (so - tm);
            dsenheat = hsg;
		 }
	     else
		 {
            double tafk;     // average air temperature above soil surface (K)
            tafk = (1 - sf) * thet + sf * (0.1 * so + 0.3 * thet + 0.6 * tv);
            senheat = hsg * (so - tafk);
            dsenheat = hsg * (1 - sf * 0.1);
		 }
//     Compute the energy balance bb. (positive direction is upward)
         double bb = emtlw         // (a) long wave radiation emitted from soil surface
                   - rls1          // long wave radiation reaching the soil surface
                   + bbsoil        // (b) heat transfer from soil surface to next soil layer
                   + hlat          // (c) latent heat transfer
                   - rss           // global radiation reaching the soil surface
                   + senheat;      // (d) heat transfer from soil surface to air
//
         if (fabs(bb) < 10e-6)
		    return; // end computation for so
//     If bb is not small enough, compute its derivative by so.
         double demtlw; // The derivative of emitted long wave radiation (emtlw)
         demtlw = 4 * rls4 * pow(so, 3);
//     Compute derivative of bbsoil
         double sop001 = so + 0.001; // soil surface temperature plus 0.001
         double rosoil1p;   // heat conductivity of 1st soil layer for so+0.001
         rosoil1p = ThermalCondSoil(VolWaterContent[0][k], sop001, 1);
         double rosoilp;    // rosoil for so+0.001
         rosoilp = (rosoil1p * dl[0] + rosoil2 * dl[1] + rosoil3 * dl[2])
                 / (dl[0] + dl[1] + dl[2])
                 / (.5 * dl[0] + dl[1] + .5 * dl[2]);
	     double drosoil = (rosoilp - rosoil) / 0.001; // derivative of rosoil
   	     double dbbsoil = rosoil + drosoil * (so - so3); // derivative of bbsoil
//     The derivative of the energy balance function
         double bbp = demtlw             // (a)
                    + dbbsoil            // (b)
                    + dhlat              // (c)
                    + dsenheat;          // (d)
//     Correct the upper soil temperature by the ratio of bb to bbp.
         double bbadjust; // the adjustment of soil surface temperature before next iteration
         bbadjust = bb / bbp;
//     If adjustment is small enough, no more iterations are needed.
	     if (fabs(bbadjust) < 0.002)
	         return;
//     If bbadjust is not the same sign as bbex, reduce fluctuations
         if (mon <= 1)
             bbex = 0;
         else if (mon >= 2)
		 {
            if(fabs(bbadjust + bbex) < fabs(bbadjust - bbex))
			{
              bbadjust = (bbadjust + bbex) / 2;
              so = (so + soex) / 2;
			}
		 }
//
	     if (bbadjust > 10)
			 bbadjust = 10;
	     if (bbadjust < -10)
			 bbadjust = -10;
//
         so = so - bbadjust;
	     so2 = so2 + (so - soex) / 2;
	     so3 = so3 + (so - soex) / 3;
         soex = so;
	     bbex = bbadjust;
         mon++;
	  } // end while loop
//     If (mon >= 50) send message on error and end simulation.
      string msg = " Infinite loop in SoilSurfaceBalance(). Abnormal stop!! \n";
      char C1[12];
      sprintf(C1, "%3d %3d %3d", Daynum, ihr, k);
      msg += " Daynum, ihr, k = " + (string) C1 + "\n";
      sprintf(C1, "%10.3g", so);
      msg += " so      = " + (string) C1 + "\n";
      sprintf(C1, "%10.3g", so2);
      msg += " so2 = " + (string) C1 + "\n";
      sprintf(C1, "%10.3g", so3);
      msg += " so3 = " + (string) C1 + "\n";
	throw Cotton2KException(msg);
}
/////////////////////////////////////////
bool SoilMulchBalance (int ihr, int k, double rlzero, double rsm, double rss, double sf,
	double &so, double &so2, double &so3, double thet, double &tm, double tv, double wndcanp)
//     This function solves the energy balance equations at the interface of the soil 
//  surface and the plastic mulch cover and computes the resulting temperatures of the 
//  soil surface and of the plastic mulch.  
//     it is called from EnergyBalance(), on each time step and for each soil column, if
//  this column is covered with a plastic mulch.  It calls functions SensibleHeatTransfer(),
//  SoilSurfaceBalance() and MulchSurfaceBalance().
//     Units for all energy fluxes are: cal cm-2 sec-1.
//     If the return value is true, it means there was an error and simulation will end.
//
//     The following global variables are referenced here:
//       bEnd, Daynum, MulchTranLW, .
//     The following arguments are set in this function:
//       so - temperature of soil surface.
//       so2 - temperature of soil 2nd layer
//       so3 - temperature of soil 3rd layer
//       tm - temperature of plastic mulch (K). When tm = 0 there is no plastic mulch.
//     The following arguments are referenced in this function:
//       ihr - the time in hours.
//       k - soil column number.
//       rlzero - incoming long wave radiation (ly / sec).
//       rsm - global radiation absorbed by mulch
//       rss - global radiation absorbed by soil surface
//       sf - fraction of shaded soil area
//       thet - air temperature (K).
//       tv - temperature of plant canopy (K).
//       wndcanp - estimated wind speed under canopy
//
{
//     Constant variables:
      const double ef = 0.95; // emissivity of the foliage surface.
      const double eg = 0.95; // emissivity of the soil surface.
      const double stefa1 = 1.38e-12; // stefan-boltsman constant.
//     Compute long wave radiation reaching the surface mulch from above, and
//  the air temperature above it.
      double rlsp0; // long wave radiation reaching the mulch from above.
      double tafk; //  temperature (K) of air inside the canopy.
      if (sf > 0.05)       // shaded column
	  {
	     rlsp0 = (1 - sf) * (1 - MulchTranLW) * rlzero              // from sky in unshaded segment
               + sf * (1 - MulchTranLW) * ef * stefa1 * pow(tv, 4); // from foliage in shaded segment
         tafk = (1 - sf) * thet
              + sf * (0.1 * tm + 0.3 * thet + 0.6 * tv);
	  }
      else                 // unshaded column
	  {
         rlsp0 = (1 - MulchTranLW) * rlzero;
   	     tafk = thet;
	  }
//     rls5 is the multiplier of tm**4 for emitted long wave radiation from mulch,
//  sum of upward and downward emittance.
      double rls5; // multiplier for emitted long wave radiation from mulch,
      rls5 = 2 * (1 - MulchTranLW) * stefa1;
//     Call SensibleHeatTransfer() to compute sensible heat transfer between plastic mulch and air
      double varcm; // sensible heat transfer coefficients for mulch to air (before multiplying by ROCP).
      varcm = SensibleHeatTransfer(tm, tafk, 0, wndcanp);
	  if (bEnd)
		  return true;
//
      double hsgm; // multiplier for computing sensible heat transfer from soil to mulch.
      double hsgp; // multiplier for computing sensible heat transfer from mulch to air.
      double rocp; // air density * specific heat at constant pressure [= 0.24 * 2 * 1013 / (5740 * tk) ]
      rocp = 0.08471 / tafk;
      hsgp = rocp * varcm;
//     Compute sensible heat transfer between plastic mulch and soil surface
      if (tm > 0)
         rocp = 0.08471 / tm;
      int mtnit = 0; // counter for numbet of iterations
      double soold1; // previous value of temperature of soil surface (k)
      double tmold1; // previous value of temperature of mulch (k)
//     Start iterations:
      do
      {
          soold1 = so; // previous value of temperature of soil surface (k)
          tmold1 = tm; // previous value of temperature of mulch (k)
//     Energy balance for soil surface (mulch interface)
          hsgm = 2 * rocp * fabs(so - tm);
	      SoilSurfaceBalance (ihr, k, 0, rlzero, rss, sf, hsgm, so, so2, so3, thet, tm, tv);
	      if (bEnd)
		      return true;
//     Add Long wave radiation reaching the mulch from the soil:
          double rlsp; // total long wave radiation reaching the mulch.
          rlsp = rlsp0 + (1 - MulchTranLW) * eg * stefa1 * pow(so, 4);
//     Energy balance for mulch (soil and air interface)
          hsgm = 2 * rocp * fabs(so - tm);
          MulchSurfaceBalance (ihr, k, rlsp, rls5, rsm, sf, hsgp, hsgm, so, thet, tm, tv);
	      if (bEnd)
		      return true;
//     Check number of iterations - do not exceed 30 iterations.
          mtnit++;
          if (mtnit > 8)
	      {
	          so = 0.5 * (so + soold1);
              tm = 0.5 * (tm + tmold1);
	      }
          if (mtnit > 30)
	      {
              string msg = " Infinite loop in SoilMulchBalance(). Abnormal stop!! \n";
              char C1[12];
              sprintf(C1, "%3d %3d %3d", Daynum, ihr, k);
              msg += " Daynum, ihr, k = " + (string) C1 + "\n";
              sprintf(C1, "%10.3g", so);
              msg += " so      = " + (string) C1 + "\n";
              sprintf(C1, "%10.3g", tm);
              msg += " tm = " + (string) C1 + "\n";
	        throw Cotton2KException(msg);
	      }
      } while (fabs(tm - tmold1) > 0.05 || fabs(so - soold1) > 0.05);
//
      return false;
}
