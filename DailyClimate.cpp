// DailyClymate.cpp : daily climate simulation.
//
//    This file contains the C++ code, commented and referenced, for computing
// hourly values of weather data and components of potential evapotranspiration.
//    Note that some variables are defined as "global", meaning that they are
// computed and used in several functions of the program. They are defined in
// files "global.h" and "global.cpp" (which also contains a "Dictionary" of these
// variables).
//     Functions in this file:
//  DayClim()
//  ComputeDayLength()
//  dayrad()
//  daytmp()
//  tdewhour()
//  dayrh()
//  daywnd()
//  AverageAirTemperatures()
//  VaporPressure()
//  EvapoTranspiration()
//  clearskyemiss()
//  cloudcov()
//  clcor()
//  del()
//  gam()
//  refalbed()
//  sunangle()
//  SimulateRunoff()
//
#include "global.h"
#include "GeneralFunctions.h"
#include "DailyClimate.h"

extern "C"
{
    double dayrad(double, double, double, double);
    double dayrh(double, double);
    double refalbed(double, double, double, double);
    double del(double, double);
    double gam(double, double);
    double cloudcov(double, double, double);
    double daywnd(double, double, double, double, double, double);
    double daytmp(Simulation &, uint32_t, double, double, double, double, uint32_t, double, double);
    double clcor(uint8_t, double, double, double, double, double, double);
    void AverageAirTemperatures(State &, const double *, double &, double &, double &);
    void ComputeDayLength(uint32_t, int32_t, double, double, double &, double &, double &, double &, double &, double &);
    void sunangle(double, double, double, double, double &, double &);
    double tdewhour(Simulation &, uint32_t, uint32_t, double, double, double, double, double, double, double, double);
    double SimulateRunoff(Simulation &, uint32_t, double, double, uint32_t);
}

void EvapoTranspiration(Simulation &, uint32_t, int);

//     Definition of file scope variables:
double declination, // daily declination angle, in radians.
    SolarNoon,      // time of solar noon, hours.
    sunr,           // time of sunrise, hours.
    suns,           // time of sunset, hours.
    tmpisr;         // extraterrestrial radiation, W / m2.
//////////////////////////////////////////////////////////////////////////////
void DayClim(Simulation &sim, uint32_t u)
//     The function DayClim() is called daily from SimulateThisDay(). It calls the
//  the following functions:
//     ComputeDayLength(), SimulateRunoff(), AverageAirTemperatures(), dayrad(),
//     daytemp(), EvapoTranspiration(), tdewhour(), dayrh(), daywnd()
//         Global variables referenced:
//    DayStart, declination, LastDayWeatherData, Latitude,
//    OutIndex, pi, SitePar
//         Global variables set:
//    AirTemp, bPollinSwitch, DewPointTemp, Radiation, RelativeHumidity, WindSpeed
//
{
    State &state = sim.states[u];
    //     Compute day length and related variables:
    ComputeDayLength(sim.day_start + u, sim.year, sim.latitude, sim.longitude, declination, tmpisr, SolarNoon, sim.states[u].day_length, sunr, suns);
    //
    double xlat = sim.latitude * pi / 180;    // latitude converted to radians.
    double cd = cos(xlat) * cos(declination); // amplitude of the sine of the solar height.
    double sd = sin(xlat) * sin(declination); // seasonal offset of the sine of the solar height.
                                              //     The computation of the daily integral of global radiation (from
                                              //  sunrise to sunset) is based on Spitters et al. (1986).
    const double c11 = 0.4;                   // constant parameter.
    double radsum;                            // daily radiation integral.
    if (fabs(sd / cd) >= 1)
        radsum = 0; //  arctic circle
    else
    {
        //     dsbe is the integral of sinb * (1 + c11 * sinb) from sunrise to sunset,
        double dsbe = acos(-sd / cd) * 24 / pi * (sd + c11 * sd * sd + 0.5 * c11 * cd * cd) + 12 * (cd * (2 + 3 * c11 * sd)) * sqrt(1 - (sd / cd) * (sd / cd)) / pi;
        //     The daily radiation integral is computed for later use in function Radiation.
        //  Daily radiation intedral is converted from langleys to Watt m-2, and divided by dsbe.
        //      11.630287 = 1000000 / 3600 / 23.884
        radsum = sim.climate[u].Rad * 11.630287 / dsbe;
    }
    //     Set 'pollination switch' for rainy days (as in GOSSYM).
    double rainToday; // The amount of rain today, mm
    rainToday = sim.climate[u].Rain;
    if (rainToday >= 2.5)
        bPollinSwitch = false;
    else
        bPollinSwitch = true;
    //     Call SimulateRunoff() only if the daily rainfall is more than 2 mm.
    //     Note: this is modified from the original GOSSYM - RRUNOFF routine. It is called here
    //  for rainfall only, but it is not activated when irrigation is applied.
    double runoffToday = 0; // amount of runoff today, mm
    if (rainToday >= 2.0)
    {
        runoffToday = SimulateRunoff(sim, u, SandVolumeFraction[0], ClayVolumeFraction[0], NumIrrigations);
        if (runoffToday < rainToday)
            rainToday -= runoffToday;
        else
            rainToday = 0;
        sim.climate[u].Rain = rainToday;
    }
    state.runoff = runoffToday;
    //     Set period for detailed output of weather variables, if requested.
    static int j1 = 0;
    static int j2 = 0;
    if (OutIndex[15] > 0 && u <= 0)
    {
        // NOTE: Used to dialog input DayStart DayFinish
        j1 = sim.day_start;
        j2 = sim.day_finish;
    }
    int jtout; // jtout > 0 if output is required
    if (OutIndex[15] > 0 && sim.day_start + u >= j1 && sim.day_start + u <= j2)
        jtout = OutIndex[15];
    else
        jtout = 0;
    //     Parameters for the daily wind function are now computed:
    //     Note:  SitePar[] are site specific parameters.
    double t1 = sunr + SitePar[1];      // the hour at which wind begins to blow (SitePar(1) hours after sunrise).
    double t2 = SolarNoon + SitePar[2]; // the hour at which wind speed is maximum (SitePar(2) hours after solar noon).
    double t3 = suns + SitePar[3];      // the hour at which wind stops to blow (SitePar(3) hours after sunset).
    double wnytf = SitePar[4];          // used for estimating night time wind (from time t3 to time t1 next day).
                                        //
    for (int ihr = 0; ihr < 24; ihr++)  //  Start hourly loop.
    {
        double ti = ihr + 0.5;                                   // time in the middle of each hourly interval.
        double sinb = sd + cd * cos(pi * (ti - SolarNoon) / 12); // sine of the solar elevation.
                                                                 //     Compute hourly global radiation, using function dayrad.
        state.hours[ihr].radiation = dayrad(ti, radsum, sinb, c11);
        //     Compute hourly temperature, using function daytmp.
        AirTemp[ihr] = daytmp(sim, u, ti, sim.states[u].day_length, SolarNoon, SitePar[8], LastDayWeatherData, sunr, suns);
        //     Compute hourly dew point temperature, using function tdewhour.
        DewPointTemp[ihr] = tdewhour(sim, u, LastDayWeatherData, ti, AirTemp[ihr], sunr, SolarNoon, SitePar[8], SitePar[12], SitePar[13], SitePar[14]);
        //     Compute hourly relative humidity, using function dayrh.
        RelativeHumidity[ihr] = dayrh(AirTemp[ihr], DewPointTemp[ihr]);
        //     Compute hourly wind speed, using function daywnd, and daily sum of wind.
        WindSpeed[ihr] = daywnd(ti, sim.climate[u].Wind, t1, t2, t3, wnytf);
    }
    //     Write output file if requested.
    if (jtout > 1)
    {
        ofstream File18(fs::path("output") / (string(sim.profile_name) + ".TM2"), ios::app);
        File18 << endl;
        File18 << " hour     radiation    temperature             RelativeHumidity           wind" << endl;
        for (int ihr = 0; ihr < 24; ihr++) //  Start hourly loop.
        {
            File18.width(5);
            File18 << ihr;
            File18.precision(5);
            File18.width(15);
            File18 << state.hours[ihr].radiation;
            File18.width(15);
            File18 << AirTemp[ihr];
            File18.width(15);
            File18 << RelativeHumidity[ihr];
            File18.width(15);
            File18 << WindSpeed[ihr] << endl;
        }
    }
    //     Compute average daily temperature, using function AverageAirTemperatures.
    AverageAirTemperatures(state, AirTemp, AvrgDailyTemp, DayTimeTemp, NightTimeTemp);
    //     Write output file if requested.
    if (jtout > 1)
    {
        ofstream File18(fs::path("output") / (string(sim.profile_name) + ".TM2"), ios::app);
        File18 << endl
               << " DayTimeTemp, NightTimeTemp, AvrgDailyTemp =   ";
        File18.precision(2);
        File18.width(10);
        File18 << DayTimeTemp;
        File18.width(10);
        File18 << NightTimeTemp;
        File18.width(10);
        File18 << AvrgDailyTemp << endl;
    }
    //     Compute potential evapotranspiration.
    EvapoTranspiration(sim, u, jtout);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
void EvapoTranspiration(Simulation &sim, uint32_t u, int jtout)
//     Function EvapoTranspiration() computes the rate of reference evapotranspiration
//  and related variables. The following subroutines and functions are called for each
//  hour: sunangle, cloudcov(), clcor(), refalbed(), VaporPressure(), clearskyemiss(), del(),
//  gam().
//     Global variables used:
//        declination, Elevation, Latitude, pi, Radiation, RelativeHumidity,
//        SitePar, AirTemp.
//     Global variables computed: CloudCoverRatio, CloudTypeCorr, albedo, rnet, ReferenceETP,
//        es1hour, es2hour, ReferenceTransp.
//     argument: jtout - index indicating if output is required.
{
    State &state = sim.states[u];
    const double stefb = 5.77944E-08; // the Stefan-Boltzman constant, in W m-2 K-4 (= 1.38E-12 * 41880)
    const double c12 = 0.125;         // c12 ... c15 are constant parameters.
    const double c13 = 0.0439;
    const double c14 = 0.030;
    const double c15 = 0.0576;
    int iamhr = 0; // earliest time in day for computing cloud cover
    int ipmhr = 0; // latest time in day for computing cloud cover
    double cosz;   // cosine of sun angle from zenith for this hour
    double suna;   // sun angle from horizon, degrees at this hour
    double isr;    // hourly extraterrestrial radiation in W / m**2
                   //      Start hourly loop
    for (int ihr = 0; ihr < 24; ihr++)
    {
        double ti = ihr + 0.5; // middle of the hourly interval
                               //      The following subroutines and functions are called for each
                               //  hour: sunangle, cloudcov, clcor, refalbed .
        sunangle(ti, sim.latitude, declination, SolarNoon, cosz, suna);
        isr = tmpisr * cosz;
        state.hours[ihr].cloud_cov = cloudcov(state.hours[ihr].radiation, isr, cosz);
        //      clcor is called to compute cloud-type correction.
        //      iamhr and ipmhr are set.
        state.hours[ihr].cloud_cor = clcor(ihr, SitePar[7], isr, cosz, state.day_length, state.hours[ihr].radiation, SolarNoon);
        if ((cosz >= 0.1736) && (iamhr == 0))
            iamhr = ihr;
        if ((ihr >= 12) && (cosz <= 0.1736) && (ipmhr == 0))
            ipmhr = ihr - 1;
        //      refalbed is called to compute the reference albedo for each hour.
        albedo[ihr] = refalbed(isr, state.hours[ihr].radiation, cosz, suna);
        if (jtout > 1)
        {
            // write output to file
            ofstream File18(fs::path("output") / (string(sim.profile_name) + ".TM2"), ios::app);
            if (ihr == 0)
            {
                File18 << "    ***********" << endl
                       << "  date =";
                File18.width(12);
                File18 << sim.states[u].date << " Day of Year =";
                File18.width(4);
                File18 << sim.day_start + u << endl;
                File18 << " hour    cosz       isr   sunangle " << endl;
            }
            File18.width(5);
            File18 << ihr;
            File18.precision(4);
            File18.width(10);
            File18 << cosz;
            File18.width(10);
            File18 << isr;
            File18.width(10);
            File18 << suna << endl;
        }
    } //   End of 1st hourly loop
      //     Zero some variables that will later be used for summation.
    ReferenceTransp = 0;
    Rn = 0;                                         // daily net radiation
    double rnet[24];                                // hourly net radiation
    double cltcoram = state.hours[iamhr].cloud_cor; //  cloud type correction during early morning
    double cltcorpm = state.hours[ipmhr].cloud_cor; //  cloud type correction during late afternoon
                                                    //
    for (int ihr = 0; ihr < 24; ihr++)              //  2nd hourly loop
    {
        double ti = ihr + 0.5;                          // middle of the hourly interval
                                                        //      Compute saturated vapor pressure (svp), using function VaporPressure().
                                                        //      The actual vapor pressure (vp) is computed from svp and the
                                                        //  relative humidity. Compute vapor pressure deficit (vpd). This
                                                        //  procedure is based on the CIMIS algorithm.
        double svp = VaporPressure(AirTemp[ihr]);       // saturated vapor pressure, mb
        double vp = 0.01 * RelativeHumidity[ihr] * svp; // vapor pressure, mb
        double vpd = svp - vp;                          // vapor pressure deficit, mb.
                                                        //   Get cloud cover and cloud correction for night hours
        if (ihr < iamhr || ihr > ipmhr)
        {
            state.hours[ihr].cloud_cov = 0;
            state.hours[ihr].cloud_cor = 0;
        }
        //     The hourly net radiation is computed using the CIMIS algorithm (Dong et al., 1988):
        //     rlonin, the hourly incoming long wave radiation, is computed from ea0, cloud cover
        //  (CloudCoverRatio), air temperature (tk),  stefb, and cloud type correction (CloudTypeCorr).
        //     rnet, the hourly net radiation, W m-2, is computed from the global radiation, the albedo,
        //  the incoming long wave radiation, and the outgoing longwave radiation.
        double tk = AirTemp[ihr] + 273.161; // hourly air temperature in Kelvin.
        double ea0 = clearskyemiss(vp, tk); // clear sky emissivity for long wave radiation
                                            //     Compute incoming long wave radiation:
        double rlonin = (ea0 * (1 - state.hours[ihr].cloud_cov) + state.hours[ihr].cloud_cov) * stefb * pow(tk, 4) - state.hours[ihr].cloud_cor;
        rnet[ihr] = (1 - albedo[ihr]) * state.hours[ihr].radiation + rlonin - stefb * pow(tk, 4);
        Rn += rnet[ihr];
        //     The hourly reference evapotranspiration ReferenceETP is computed by the
        //  CIMIS algorithm using the modified Penman equation:
        //     The weighting ratio (w) is computed from the functions del() (the slope of the saturation
        //  vapor pressure versus air temperature) and gam() (the psychometric constant).
        double w = del(tk, svp) / (del(tk, svp) + gam(sim.elevation, AirTemp[ihr])); // coefficient of the Penman equation
                                                                                     //     The wind function (fu2) is computed using different sets of parameters
                                                                                     //  for day-time and night-time. The parameter values are as suggested by CIMIS.
        double fu2;                                                                  // wind function for computing evapotranspiration
        if (state.hours[ihr].radiation <= 0)
            fu2 = c12 + c13 * WindSpeed[ihr];
        else
            fu2 = c14 + c15 * WindSpeed[ihr];
        //     hlathr, the latent heat for evaporation of water (W m-2 per mm at this hour) is
        //  computed as a function of temperature.
        double hlathr = 878.61 - 0.66915 * (AirTemp[ihr] + 273.161);
        //     ReferenceETP, the hourly reference evapotranspiration, is now computed by the
        //  modified Penman equation.
        ReferenceETP[ihr] = w * rnet[ihr] / hlathr + (1 - w) * vpd * fu2;
        if (ReferenceETP[ihr] < 0)
            ReferenceETP[ihr] = 0;
        //     ReferenceTransp is the sum of ReferenceETP
        ReferenceTransp += ReferenceETP[ihr];
        //     es1hour and es2hour are computed as the hourly potential evapotranspiration due to
        //  radiative and aerodynamic factors, respectively.
        //     es1hour and ReferenceTransp are not computed for periods of negative net radiation.
        es2hour[ihr] = (1 - w) * vpd * fu2;
        if (rnet[ihr] > 0)
            es1hour[ihr] = w * rnet[ihr] / hlathr;
        else
            es1hour[ihr] = 0;
    } //   end of 2nd hourly loop
}
