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
    double daywnd(double, double, double, double, double, double);
    double daytmp(Simulation &, uint32_t, double, double, uint32_t, double, double);
    void AverageAirTemperatures(Hour[24], double &, double &, double &);
    void ComputeDayLength(uint32_t, int32_t, double, double, double &, double &, double &, double &, double &, double &);
    double tdewhour(Simulation &, uint32_t, uint32_t, double, double, double, double, double, double, double, double);
    double SimulateRunoff(Simulation &, uint32_t, double, double, uint32_t);
    void EvapoTranspiration(State &, double, double, double, double, double);
}


//     Definition of file scope variables:
double declination, // daily declination angle, in radians.
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
    ComputeDayLength(sim.day_start + u, sim.year, sim.latitude, sim.longitude, declination, tmpisr, state.solar_noon, state.day_length, sunr, suns);
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
    double t1 = sunr + SitePar[1];             // the hour at which wind begins to blow (SitePar(1) hours after sunrise).
    double t2 = state.solar_noon + SitePar[2]; // the hour at which wind speed is maximum (SitePar(2) hours after solar noon).
    double t3 = suns + SitePar[3];             // the hour at which wind stops to blow (SitePar(3) hours after sunset).
    double wnytf = SitePar[4];                 // used for estimating night time wind (from time t3 to time t1 next day).
                                               //
    for (int ihr = 0; ihr < 24; ihr++)         //  Start hourly loop.
    {
        Hour &hour = state.hours[ihr];
        double ti = ihr + 0.5;                                          // time in the middle of each hourly interval.
        double sinb = sd + cd * cos(pi * (ti - state.solar_noon) / 12); // sine of the solar elevation.
                                                                        //     Compute hourly global radiation, using function dayrad.
        hour.radiation = dayrad(ti, radsum, sinb, c11);
        //     Compute hourly temperature, using function daytmp.
        hour.temperature = daytmp(sim, u, ti, SitePar[8], LastDayWeatherData, sunr, suns);
        //     Compute hourly dew point temperature, using function tdewhour.
        hour.dew_point = tdewhour(sim, u, LastDayWeatherData, ti, hour.temperature, sunr, state.solar_noon, SitePar[8], SitePar[12], SitePar[13], SitePar[14]);
        //     Compute hourly relative humidity, using function dayrh.
        hour.humidity = dayrh(state.hours[ihr].temperature, hour.dew_point);
        //     Compute hourly wind speed, using function daywnd, and daily sum of wind.
        hour.wind_speed = daywnd(ti, sim.climate[u].Wind, t1, t2, t3, wnytf);
    }
    //     Compute average daily temperature, using function AverageAirTemperatures.
    AverageAirTemperatures(state.hours, AvrgDailyTemp, DayTimeTemp, NightTimeTemp);
    //     Compute potential evapotranspiration.
    EvapoTranspiration(sim.states[u], sim.latitude, sim.elevation, declination, tmpisr, SitePar[7]);
}
