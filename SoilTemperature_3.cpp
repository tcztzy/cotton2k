//  File SoilTemperature_3.cpp
//
//   List of functions in this file:
//       CanopyBalance()
//       MulchSurfaceBalance()
//       ThermalCondSoil()
//       PredictEmergence()
//
#include <math.h>
#include <stdio.h>

#include "CottonSimulation.h"
#include "GeneralFunctions.h"

//
#ifdef _DEBUG
#define new DEBUG_NEW
#endif
/////////////////////////////////////////
void CanopyBalance(int ihr, int k, double etp1, double rlzero, double rsv,
                   double c2, double sf, double so, double thet, double tm,
                   double &tv)
//     This function solves the energy balance equations at the foliage / air
//     interface, and
//  computes the resulting temperature of the foliage. It is Called from
//  EnergyBalance().
//     Units for all energy fluxes are: cal cm-2 sec-1.
//
//     The following global variables are referenced here:    Daynum,
//     MulchTranLW. The following argument is set in this function:
//       tv - temperature of plant canopy.
//     The following arguments are referenced in this function:
//       c2 - multiplier for sensible heat transfer at plant surface.
//       etp1 - transpiration (mm / sec).
//       ihr - the time in hours.
//       k - soil column number.
//       rlzero - incoming long wave radiation (ly / sec).
//       rsv - global radiation absorbed by the vegetation
//       sf - fraction of shaded soil area
//       so - temperature of soil surface (k).
//       thet - air temperature (k).
//       tm - plastic mulch temperature (K), tm = 0 if column not mulched.
{
    //     Constants:
    const double ef = 0.95;          // emissivity of the foliage surface
    const double eg = 0.95;          // emissivity of the soil surface
    const double stefa1 = 1.38e-12;  // stefan-boltsman constant.
                                     //
    double rlv1;                     // long wave radiation reaching the canopy
    if (tm > 0)                      // column is mulched
        rlv1 = sf * ef * rlzero      // from sky
               + sf * ef * eg * stefa1 * MulchTranLW * pow(so, 4)  // from soil
               +
               sf * ef * (1 - MulchTranLW) * stefa1 * pow(tm, 4);  // from mulch
    else                                              // column is not mulched
        rlv1 = sf * ef * rlzero                       // from sky
               + sf * ef * eg * stefa1 * pow(so, 4);  // from soil
    //     rlv4 is the multiplier of tv**4 for emitted long wave radiation from
    //     vegetation, corrected
    //  for the amount reflected back from soil surface and absorbed by foliage.
    //  This is two-sided ( note that when eg = ef = 1, the coefficient corr
    //  will be 2)
    double corr = 1 + eg / (ef + eg - ef * eg);  // coefficient
    double rlv4 = stefa1 * sf * ef * corr;
    double dsenfheat = c2 * (1 - 0.6 * sf);  // derivative of senfheat
                                             //     Start iterations for tv:
    int mot = 0;  // count of iterations for plant surface
    while (mot < 50) {
        //     Latent heat flux (hvlat) is computed from the transpiration rate.
        double hvlat = (75.5255 - 0.05752 * tv) * etp1;
        double dhvlat = -0.05752 * etp1;  // derivative of hvlat
        //     Emitted long wave radiation from vegetation (cclwe)
        double cclwe = rlv4 * pow(tv, 4);
        //     Sensible heat transfer from vegetation
        double tvex = tv;  // previous value of tv
        double
            tafk;  // average air temperature above soil surface (k) in canopy
        if (tm > 0.)
            tafk = (1 - sf) * thet + sf * (0.1 * tm + 0.3 * thet + 0.6 * tv);
        else
            tafk = (1 - sf) * thet + sf * (0.1 * so + 0.3 * thet + 0.6 * tv);
        double senfheat;  // sensible heat transfer from foliage
        senfheat = c2 * (tv - tafk);
        //     Compute the energy balance at the plant surface (cc), and
        //  if it is small enough end the computation.
        double cc =
            cclwe        // (a) long wave emission from vegetation
            - rlv1       // long wave radiation reaching the vegetation
            + hvlat      // (b) latent heat transfer
            - rsv        // global radiation on vegetation
            + senfheat;  // (c) sensible heat transfer from vegetation to air
        if (fabs(cc) < 10e-6) {
            tv = 0.5 * (tvex + tv);
            return;  // end iterations for tv
        }
        //     If cc is not small enough, compute its derivative by tv (ccp).
        double dcclwe = 4 * rlv4 * pow(tv, 3);  // derivative of cclwe
        //     ccp is the derivative of energy balance at the plant surface (by
        //     tv)
        double ccp = dcclwe        // (a)
                     + dhvlat      // (b)
                     + dsenfheat;  // (c)
        //     Correct the canopy temperature by  the ratio of cc to ccp.
        double ccadjust = cc / ccp;  // adjustment of tv before next iteration
        //     If adjustment is small enough, no more iterations are needed.
        if (fabs(ccadjust) < 0.002) return;
        //     If ccadjust is not the same sign as ccadx, reduce fluctuations
        double ccadx;  // previous value of ccadjust
        if (mot >= 2) {
            if (fabs(ccadjust - ccadx) > fabs(ccadjust + ccadx)) {
                ccadjust = (ccadjust + ccadx) * 0.5;
                tv = (tv + tvex) * 0.5;
            }
        }
        if (ccadjust > 10) ccadjust = 10;
        if (ccadjust < -10) ccadjust = -10;
        tv = tv - ccadjust;
        tvex = tv;
        ccadx = ccadjust;
        mot++;
    }
    //     If reached 50 iterations there must be an error somewhere!
    fprintf(stderr, " Infinite loop in CanopyBalance(). Abnormal stop!! \n");
    fprintf(stderr, " Daynum, ihr, k = %3d %3d %3d\n", Daynum, ihr, k);
    fprintf(stderr, " so      = %10.3g\n", so);
    fprintf(stderr, " tv = %10.3g\n", tv);
    bEnd = true;
}
/////////////////////////////////////////
void MulchSurfaceBalance(int ihr, int k, double rlsp, double rls5, double rsm,
                         double sf, double hsgp, double hsgm, double so,
                         double thet, double &tm, double tv)
//     This function solves the energy balance equations for the transparent
//     (polyethylene) mulch,
//  and computes its resulting temperature (tm). It is called from
//  EnrgyBalance(), and uses function ThermalCondSoil(). Units for all energy
//  fluxes are: cal cm-2 sec-1.
//
//     The following global variables are referenced here:        Daynum,
//     MulchTranLW. The following argument is set in this function:
//       tm - temperature of surface mulch(K)
//     The following arguments are used in this function:
//       hsgp - multiplier for computing sensible heat transfer from mulch to
//       air hsgm - multiplier for computing sensible heat transfer from mulch
//       to soil ihr - the time in hours. k - soil column number. rls5 =
//       multiplier for long wave energy emitted from soil surface rlsp -
//       incoming long wave radiation (ly / sec) reaching plastic layer. rsm -
//       global radiation absorbed by mulch surface sf - fraction of shaded soil
//       area so - temperature of soil surface. thet - air temperature (K). tv -
//       temperature of plant canopy
{
    double dsenheat;  // derivative of sensible heat transfer
    if (sf <= 0.05)   // if this is a shaded column
        dsenheat = hsgp * (1 - sf * 0.1) + hsgm;
    else
        dsenheat = hsgp + hsgm;
    //     Start iterations for soil mulch temperature (tm)
    int mop = 0;  // count of iterations for soil mulch energy balance
    do {
        //     Emitted long wave radiation from soil mulch
        double emtlw;  // emitted long wave radiation from soil surface
        emtlw = rls5 * pow(tm, 4);
        //     Sensible heat transfer
        double tafk;     // average air temperature above soil surface and in
                         // canopy (K)
        double senheat;  // sensible heat transfer from soil surface
        if (sf >= 0.05)  // if this is a shaded column
            tafk = (1 - sf) * thet + sf * (0.1 * tm + 0.3 * thet + 0.6 * tv);
        else
            tafk = thet;
        senheat = hsgp * (tm - tafk) + hsgm * (tm - so);
        //     Compute the energy balance tmbal. (positive direction is upward)
        double tmbal;       // energy balance of the soil mulch
        tmbal = emtlw       // (a) long wave radiation emitted from soil mulch
                - rlsp      // long wave radiation reaching the soil mulch
                - rsm       // global radiation absorbed
                + senheat;  // (d) heat transfer from soil surface to air
        if (fabs(tmbal) < 10.e-6) return;  // end iterations for tm
        //     If tmbal is not small enough, compute its derivative by tm
        //     (tmbalp). The derivative of emitted long wave radiation
        double demtlw;  // derivative of emtlw
        demtlw = 4 * rls5 * pow(tm, 3);
        //	   The derivative of the energy balance function
        double
            tmbalp;  // derivative of energy balance at the soil mulch (by tm)
        tmbalp = demtlw       // (a)
                 + dsenheat;  // (d)
        //     Correct the upper soil temperature by  the ratio of tmbal to
        //     tmbalp.
        double tmbaladjust;  // adjustment of soil mulch temperature before next
                             // iteration
        tmbaladjust = tmbal / tmbalp;
        //     If adjustment is small enough, no more iterations are needed.
        if (fabs(tmbaladjust) < 0.002) return;
        //     If tmbaladjust is not the same sign as tmbalex, reduce
        //     fluctuations
        double tmex;     // previous value of mulch temperature (k)
        double tmbalex;  // previous value of tmbaladjust.
        if (mop >= 2)
            if (fabs(tmbaladjust + tmbalex) < fabs(tmbaladjust - tmbalex)) {
                tmbaladjust = (tmbaladjust + tmbalex) / 2;
                tm = (tm + tmex) / 2;
            }
        if (tmbaladjust > 10) tmbaladjust = 10;
        if (tmbaladjust < -10) tmbaladjust = -10;
        //
        tm = tm - tmbaladjust;
        tmex = tm;
        tmbalex = tmbaladjust;
        mop++;
        //     Go to next iteration
    } while (mop < 50);
    //     If reached 50 iterations there must be an error somewhere!
    fprintf(stderr,
            " Infinite loop in MulchSurfaceBalance(). Abnormal stop!! \n");
    fprintf(stderr, " Daynum, ihr, k = %3d %3d %3d\n", Daynum, ihr, k);
    fprintf(stderr, " so      = %10.3g\n", so);
    fprintf(stderr, " tv = %10.3g\n", tv);
    bEnd = true;
}
/////////////////////////////////////////
double ThermalCondSoil(double q0, double t0, int l0)
//     This function computes and returns the thermal conductivity of the soil
//  (cal cm-1 s-1 oC-1). It is based on the work of De Vries(1963).
//
//     The following global variables are referenced here:
//       ClayVolumeFraction, FieldCapacity, HeatCondDrySoil,
//       MarginalWaterContent, PoreSpace, SandVolumeFraction.
//     The following arguments are used in this function:
//       l0 - soil layer.
//       q0 - volumetric soil moisture content.
//       t0 - soil temperature (K).
{
    //     Constant parameters:
    const double bclay =
        7.0;  // heat conductivity of clay (= 7 mcal cm-1 s-1 oc-1).
    const double bsand =
        20.0;  // heat conductivity of sand (= 20 mcal cm-1 s-1 oc-1).
    const double cka =
        0.0615;  // heat conductivity of air (= 0.0615 mcal cm-1 s-1 oc-1).
    const double ckw =
        1.45;  // heat conductivity of water (= 1.45 mcal cm-1 s-1 oc-1).
               //     Convert soil temperature to degrees C.
    double tcel = t0 - 273.161;  // soil temperature, in C.
    //     Compute cpn, the apparent heat conductivity of air in soil pore
    //     spaces, when saturated with
    //  water vapor, using a function of soil temperature, which changes
    //  linearly between 36 and 40 C.
    double bb;  // effect of temperature on heat conductivity of air saturated
                // with water vapor.
    if (tcel <= 36)
        bb = 0.06188;
    else if (tcel > 36 && tcel <= 40)
        bb = 0.06188 + (tcel - 36) * (0.05790 - 0.06188) / (40 - 36);
    else
        bb = 0.05790;
    double cpn;  // apparent heat conductivity of air in soil pore spaces, when
                 // it is saturated with water vapor.
    cpn = cka + 0.05 * exp(bb * tcel);
    //     Compute xair, air content of soil per volume, from soil porosity and
    //     moisture content. Compute thermal conductivity (a) for wet soil (soil
    //     moisture higher than field capacity),
    //                                  (b) for less wet soil.
    //     In each case compute first ga, and then dair.
    double xair;  // air content of soil, per volume.
    xair = PoreSpace[l0] - q0;
    if (xair < 0) xair = 0;
    double dair;   // aggregation factor for air in soil pore spaces.
    double ga;     // shape factor for air in pore spaces.
    double hcond;  // computed heat conductivity of soil, mcal cm-1 s-1 oc-1.
    if (q0 >= FieldCapacity[l0]) {
        //     (a) Heat conductivity of soil wetter than field capacity.
        ga = 0.333 - 0.061 * xair / PoreSpace[l0];
        dair = form(cpn, ckw, ga);
        hcond = (q0 * ckw + dsand * bsand * SandVolumeFraction[l0] +
                 dclay * bclay * ClayVolumeFraction[l0] + dair * cpn * xair) /
                (q0 + dsand * SandVolumeFraction[l0] +
                 dclay * ClayVolumeFraction[l0] + dair * xair);
    } else {
        //     (b) For soil less wet than field capacity, compute also ckn (heat
        //     conductivity
        // of air in the soil pores).
        double qq;   // soil water content for computing ckn and ga.
        double ckn;  // heat conductivity of air in pores in soil.
        qq = q0;
        if (qq < MarginalWaterContent[l0]) qq = MarginalWaterContent[l0];
        ckn = cka + (cpn - cka) * qq / FieldCapacity[l0];
        ga = 0.041 + 0.244 * (qq - MarginalWaterContent[l0]) /
                         (FieldCapacity[l0] - MarginalWaterContent[l0]);
        dair = form(ckn, ckw, ga);
        hcond = (qq * ckw + dsand * bsand * SandVolumeFraction[l0] +
                 dclay * bclay * ClayVolumeFraction[l0] + dair * ckn * xair) /
                (qq + dsand * SandVolumeFraction[l0] +
                 dclay * ClayVolumeFraction[l0] + dair * xair);
        //     When soil moisture content is less than the limiting value
        //     MarginalWaterContent,
        //  modify the value of hcond.
        if (qq <= MarginalWaterContent[l0])
            hcond =
                (hcond - HeatCondDrySoil[l0]) * q0 / MarginalWaterContent[l0] +
                HeatCondDrySoil[l0];
    }  // q0
       //     The result is hcond converted from mcal to cal.
    double result = hcond / 1000;
    return result;
}
////////////////////////////////////////////////////////////////////////////////////
void PredictEmergence(int hour)
//     This function predicts date of emergence. It is called from
//     soil_thermology(). There is one referenced argument (hour).
//
//     The following global variables are referenced here:
//       Daynum, dl, iyear, DayPlant, PlantRowColumn, nl, SoilPsi,
//       SoilTemp.
//     The following global variables are set here:
//       DayEmerge, isw, Kday.
//
{
    const double dpl = 5;            // depth of planting, cm (assumed 5).
    static double DelayOfEmergence;  // effect of negative values of xt on
                                     // germination rate.
    static double HypocotylLength;   // length of hypocotyl, cm.
    static double
        SeedMoisture;       // moisture content of germinating seeds, percent.
    static int nSeedLayer;  // layer number where the seeds are located.
    //     Define some initial values on day of planting.
    if (Daynum == DayPlant && hour == 0) {
        DelayOfEmergence = 0;
        HypocotylLength = 0.3;
        SeedMoisture = 8;
        //     Compute soil layer number for seed depth.
        double sumdl = 0;  // depth to the bottom of a soil layer.
        for (int l = 0; l < nl; l++) {
            sumdl += dl[l];
            if (sumdl >= dpl) {
                nSeedLayer = l;
                break;
            }
        }
    }
    //     Compute matric soil moisture potential at seed location.
    //     Define te as soil temperature at seed location, C.
    double psi;  // matric soil moisture potential at seed location.
    double te;   // soil temperature at seed depth, C.
    psi = SoilPsi[nSeedLayer][PlantRowColumn];
    te = SoilTemp[nSeedLayer][PlantRowColumn] - 273.161;
    if (te < 10) te = 10;
    //
    //     Phase 1 of of germination - imbibition. This phase is executed when
    //     the moisture
    //  content of germinating seeds is not more than 80%.
    double dw;  // rate of moisture addition to germinating seeds, percent per
                // hour.
    double
        xkl;  // a function of temperature and moisture, used to calculate dw.
    if (SeedMoisture <= 80) {
        xkl = .0338 + .0000855 * te * te - 0.003479 * psi;
        if (xkl < 0) xkl = 0;
        //     Compute the rate of moisture addition to germinating seeds,
        //     percent per hour.
        dw = xkl * (80 - SeedMoisture);
        //     Compute delw, the marginal value of dw, as a function of soil
        //  temperature and soil water potential.
        double delw;  // marginal value of dw.
        if (te < 21.2)
            delw = -0.1133 + .000705 * te * te - .001348 * psi +
                   .001177 * psi * psi;
        else if (te < 26.66)
            delw =
                -.3584 + .001383 * te * te - .03509 * psi + .003507 * psi * psi;
        else if (te < 32.3)
            delw = -.6955 + .001962 * te * te - .08335 * psi +
                   .007627 * psi * psi - .006411 * psi * te;
        else
            delw = 3.3929 - .00197 * te * te - .36935 * psi +
                   .00865 * psi * psi + .007306 * psi * te;
        if (delw < 0.01) delw = 0.01;
        //    Add dw to tw, or if dw is less than delw assign 100% to tw.
        if (dw > delw)
            SeedMoisture += dw;
        else
            SeedMoisture = 100;
        return;
    }
    //
    //     Phase 2 of of germination - hypocotyl elongation.
    double xt;  // a function of temperature, used to calculate de.
    if (te > 39.9)
        xt = 0;
    else
        xt = 0.0853 - 0.0057 * (te - 34.44) * (te - 34.44) / (41.9 - te);
    //     At low soil temperatures, when negative values of xt occur,
    //  compute the delay in germination rate.
    if (xt < 0 && te < 14) {
        DelayOfEmergence += xt / 2;
        return;
    } else {
        if (DelayOfEmergence < 0) {
            if (DelayOfEmergence + xt < 0) {
                DelayOfEmergence += xt;
                return;
            } else {
                xt += DelayOfEmergence;
                DelayOfEmergence = 0;
            }
        }
    }
    //     Compute elongation rate of hypocotyl, de, as a sigmoid
    //  function of HypocotylLength. Add de to HypocotylLength.
    double de;  // rate of hypocotyl growth, cm per hour.
    de = 0.0567 * xt * HypocotylLength * (10 - HypocotylLength);
    HypocotylLength += de;
    //     Check for completion of emergence (when HypocotylLength exceeds
    //     planting
    //  depth) and report germination to output. isw is 2 after emergence.
    if (HypocotylLength > dpl) {
        isw = 2;
        DayEmerge = Daynum;
        Kday = 1;
    }
}
