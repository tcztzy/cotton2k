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

tuple<double> ComputeDayLength(const int &, const double &, const double &);

extern "C"
{
    double dayrad(double, double, double, double);
    double dayrh(double, double);
}

double daytmp(double, const int &, const double &, const Climstruct[400]);

double tdewhour(double, double, const int &, const Climstruct[400]);

double daywnd(double, double, double, double, double, double);

void AverageAirTemperatures();

void EvapoTranspiration(int, const string &, const string &, const int &, const double &, const double &);

double cloudcov(double, double, double);

double clcor(int, double, double, double, const double &);

double del(double, double);

double gam(double, double);

double refalbed(double, double, double, double);

void sunangle(double, double &, double &, const double &);

double SimulateRunoff(double, const int &, const int &, const Climstruct[400]);

//     Definition of file scope variables:
double declination,      // daily declination angle, in radians.
SolarNoon,        // time of solar noon, hours.
sunr,             // time of sunrise, hours.
suns,             // time of sunset, hours.
tmpisr;           // extraterrestrial radiation, W / m2.
//
//     All daily weather data, read from input file, are stored in the structure:
//         Clim[400];    --  defined in global
//     The values are extracted from this structure by function GetFromClim(), see 
//  file "GeneralFunctions.cpp"
//////////////////////////////////////////////////////////////////////////////
tuple<double> DayClim(const string &ProfileName, const string &Date, const int &Daynum, const int &DayOfSimulation,
                      const int &DayStart, const int &DayFinish, const double &Latitude, const double &Longitude,
                      Climstruct Clim[400])
//     The function DayClim() is called daily from SimulateThisDay(). It calls the
//  the following functions:
//     ComputeDayLength(), GetFromClim(), SimulateRunoff(), AverageAirTemperatures(), dayrad(),
//     daytemp(), EvapoTranspiration(), tdewhour(), dayrh(), daywnd()
//         Global variables referenced:
//    DayStart, declination, LastDayWeatherData, Latitude, 
//    OutIndex, pi, SitePar
//         Global variables set:
//    AirTemp, bPollinSwitch, DewPointTemp, Radiation, RelativeHumidity, WindSpeed
//
{
//     Compute day length and related variables:
    double DayLength;
    tie(DayLength) = ComputeDayLength(Daynum, Latitude, Longitude);
//
    double xlat = Latitude * pi / 180;         // latitude converted to radians.
    double cd = cos(xlat) * cos(declination);  // amplitude of the sine of the solar height.
    double sd = sin(xlat) * sin(declination);  // seasonal offset of the sine of the solar height.
//     The computation of the daily integral of global radiation (from
//  sunrise to sunset) is based on Spitters et al. (1986). 
    const double c11 = 0.4;    // constant parameter.
    double radsum;       // daily radiation integral.
    if (fabs(sd / cd) >= 1)
        radsum = 0;    //  arctic circle
    else {
//     dsbe is the integral of sinb * (1 + c11 * sinb) from sunrise to sunset,
        double dsbe = acos(-sd / cd) * 24 / pi * (sd + c11 * sd * sd + 0.5 * c11 * cd * cd)
                      + 12 * (cd * (2 + 3 * c11 * sd)) * sqrt(1 - (sd / cd) * (sd / cd)) / pi;
//     The daily radiation integral is computed for later use in function Radiation.
//  Daily radiation intedral is converted from langleys to Watt m-2, and divided by dsbe.
//      11.630287 = 1000000 / 3600 / 23.884
        radsum = GetFromClim(Clim, "rad", Daynum) * 11.630287 / dsbe;
    }
//     Set 'pollination switch' for rainy days (as in GOSSYM).
    double rainToday; // The amount of rain today, mm
    rainToday = GetFromClim(Clim, "rain", Daynum);
    if (rainToday >= 2.5)
        bPollinSwitch = false;
    else
        bPollinSwitch = true;
//     Call SimulateRunoff() only if the daily rainfall is more than 2 mm.
//     Note: this is modified from the original GOSSYM - RRUNOFF routine. It is called here
//  for rainfall only, but it is not activated when irrigation is applied.
    double runoffToday = 0; // amount of runoff today, mm
    if (rainToday >= 2.0) {
        runoffToday = SimulateRunoff(rainToday, Daynum, DayStart, Clim);
        if (runoffToday < rainToday)
            rainToday -= runoffToday;
        else
            rainToday = 0;
        int j = Daynum - DayStart;  // days from start of simulation
        Clim[j].Rain = rainToday;
    }
    Scratch21[DayOfSimulation - 1].runoff = runoffToday;
//     Set period for detailed output of weather variables, if requested.
    static int j1 = 0;
    static int j2 = 0;
    if (OutIndex[15] > 0 && Daynum <= DayStart) {
        // NOTE: Used to dialog input DayStart DayFinish
        j1 = DayStart;
        j2 = DayFinish;
    }
    int jtout;  // jtout > 0 if output is required
    if (OutIndex[15] > 0 && Daynum >= j1 && Daynum <= j2)
        jtout = OutIndex[15];
    else
        jtout = 0;
//     Parameters for the daily wind function are now computed:
//     Note:  SitePar[] are site specific parameters.
    double t1 = sunr + SitePar[1];      // the hour at which wind begins to blow (SitePar(1) hours after sunrise).
    double t2 = SolarNoon + SitePar[2]; // the hour at which wind speed is maximum (SitePar(2) hours after solar noon).
    double t3 = suns + SitePar[3];      // the hour at which wind stops to blow (SitePar(3) hours after sunset).
    double wnytf = SitePar[4]; // used for estimating night time wind (from time t3 to time t1 next day).
//
    for (int ihr = 0; ihr < 24; ihr++)  //  Start hourly loop.
    {
        double ti = ihr + 0.5;   // time in the middle of each hourly interval.
        double sinb = sd + cd * cos(pi * (ti - SolarNoon) / 12);  // sine of the solar elevation.
//     Compute hourly global radiation, using function dayrad.
        Radiation[ihr] = dayrad(ti, radsum, sinb, c11);
//     Compute hourly temperature, using function daytmp.
        AirTemp[ihr] = daytmp(ti, Daynum, DayLength, Clim);
//     Compute hourly dew point temperature, using function tdewhour.
        DewPointTemp[ihr] = tdewhour(ti, AirTemp[ihr], Daynum, Clim);
//     Compute hourly relative humidity, using function dayrh.
        RelativeHumidity[ihr] = dayrh(AirTemp[ihr], DewPointTemp[ihr]);
//     Compute hourly wind speed, using function daywnd, and daily sum of wind.
        WindSpeed[ihr] = daywnd(ti, GetFromClim(Clim, "wind", Daynum), t1, t2, t3, wnytf);
    }
//     Write output file if requested.
    if (jtout > 1) {
        ofstream File18(fs::path("output") / (ProfileName + ".TM2"), ios::app);
        File18 << endl;
        File18 << " hour     radiation    temperature             RelativeHumidity           wind" << endl;
        for (int ihr = 0; ihr < 24; ihr++)  //  Start hourly loop.
        {
            File18.width(5);
            File18 << ihr;
            File18.precision(5);
            File18.width(15);
            File18 << Radiation[ihr];
            File18.width(15);
            File18 << AirTemp[ihr];
            File18.width(15);
            File18 << RelativeHumidity[ihr];
            File18.width(15);
            File18 << WindSpeed[ihr] << endl;
        }
    }
//     Compute average daily temperature, using function AverageAirTemperatures.
    AverageAirTemperatures();
//     Write output file if requested.
    if (jtout > 1) {
        ofstream File18(fs::path("output") / (ProfileName + ".TM2"), ios::app);
        File18 << endl << " DayTimeTemp, NightTimeTemp, AvrgDailyTemp =   ";
        File18.precision(2);
        File18.width(10);
        File18 << DayTimeTemp;
        File18.width(10);
        File18 << NightTimeTemp;
        File18.width(10);
        File18 << AvrgDailyTemp << endl;
    }
//     Compute potential evapotranspiration.
    EvapoTranspiration(jtout, ProfileName, Date, Daynum, Latitude, DayLength);
    return make_tuple(DayLength);
}

//////////////////////////////////////////////////////////////////////////////////////////
tuple<double> ComputeDayLength(const int &Daynum, const double &Latitude, const double &Longitude)
//     Function ComputeDayLength() computes day length, declination, time of
//  solar noon, and extra-terrestrial radiation for this day. The CIMIS
//  (California Irrigation Management Information System) algorithms are used. 
//     Global variables referenced here:  
//  iyear, Latitude, Longitude, pi, 
//     Global variables set here:   
//  DayLength, declination
//
{
//     Convert day of year to corresponding angle in radians (xday). It uses function
//  LeapYear() (see file GeneralFunctions.cpp)
    double xday = 2 * pi * (Daynum - 1) / (365 + LeapYear(iyear));
//     Compute declination angle for this day. The equation used here for computing it 
//  is taken from the CIMIS algorithm. 
    declination = .006918 - .399912 * cos(xday) + .070257 * sin(xday)
                  - .006758 * cos(2 * xday) + .000907 * sin(2 * xday)
                  - .002697 * cos(3 * xday) + .001480 * sin(3 * xday);
//     Compute extraterrestrial radiation in W m-2. The 'solar constant' (average
//  value = 1367 W m-2) is corrected for this day's distance between earth and the sun. The 
//  equation used here is from the CIMIS algorithm, which is based on the work of Iqbal (1983).
    tmpisr = 1367 * (1.00011 + .034221 * cos(xday)
                     + .00128 * sin(xday) + .000719 * cos(2 * xday)
                     + .000077 * sin(2 * xday));
//     Time of solar noon (SolarNoon) is computed by the CIMIS algorithm,
//  using a correction for longitude (f), and the date correction (exday).
//     It is assumed that the time zone is "geographically correct". For
//  example, longitude between 22.5 and 37.5 East is in time zone GMT+2,
//  and longitude between 112.5 and 127.5 West is in time zone GMT-8.
//     All daily times in the model are computed by this method.
    double exday = (.000075 + .001868 * cos(xday) - .032077 * sin(xday)
                    - .014615 * cos(2 * xday) - .04089 * sin(2 * xday)) * 12. / pi;
    double st = 15 * (int) (Longitude / 15);
    double f = (Longitude - st) / 15;
    if (Longitude > 0)        // east of Greenwich
    {
        if (f > 0.5)
            f -= 1;
    } else                      // west  of Greenwich
    {
        if (f < -0.5)
            f += 1;
    }
    SolarNoon = 12 - f - exday;
//     Compute day length, by commonly used equations, from latitude and declination of
//  this day. Times of sunrise and of sunset are computed from solar noon and day length.
    double xlat = Latitude * pi / 180;
    double ht = -tan(xlat) * tan(declination);
    if (ht > 1)
        ht = 1;         //  arctic circle
    else if (ht < -1)
        ht = -1;
    double DayLength = 2 * acos(ht) * 12 / pi;
    sunr = SolarNoon - DayLength / 2;
    suns = sunr + DayLength;
    return make_tuple(DayLength);
}

///////////////////////////////////////////////////////////////////////////////////////////
double daytmp(double ti, const int &Daynum, const double &DayLength, const Climstruct Clim[400])
//     Function daytmp() computes and returns the hourly values of air temperature,
//  using the measured daily maximum and minimum.
//     The algorithm is described in Ephrath et al. (1996). It is based on 
//  the following assumptions:
//     The time of minimum daily temperature is at sunrise.
//     The time of maximum daily temperature is SitePar[8] hours after solar noon.
//     Many models assume a sinusoidal curve of the temperature during the day,
//  but actual data deviate from the sinusoidal curve in the following characteristic 
//  way: a faster increase right after sunrise, a near plateau maximum during several 
//  hours in the middle of the day, and a rather fast decrease by sunset. The physical
//  reason for this is a more efficient mixing of heated air from ground level into the 
//  atmospheric boundary layer, driven by strong lapse temperature gradients buoyancy.
//     Note: ** will be used for "power" as in Fortran notation.
//     A first order approximation is
//       daytmp = tmin + (tmax-tmin) * st * tkk / (tkk + daytmp - tmin)
//  where
//       st = sin(pi * (ti - SolarNoon + dayl / 2) / (dayl + 2 * SitePar[8]))
//     Since daytmp appears on both sides of the first equation, it
//  can be solved and written explicitly as:
//      daytmp = tmin - tkk/2 + 0.5 * sqrt(tkk**2 + 4 * amp * tkk * st)
//  where the amplitude of tmin and tmax is calculated as
//      amp = (tmax - tmin) * (1 + (tmax - tmin) / tkk)
//     This ensures that temperature still passes through tmin and tmax values.
//     The value of tkk was determined by calibration as 15.
//     This algorithm is used for the period from sunrise to the time of maximum temperature, 
//  hmax. A similar algorithm is used for the time from hmax to sunset, but the value of the 
//  minimum temperature of the next day (mint_tomorrow) is used instead of mint_today.
//     Night air temperature is described by an exponentially declining curve.
//     For the time from sunset to mid-night:
//        daytmp = (mint_tomorrow - sst * exp((dayl - 24) / tcoef)
//               + (sst - mint_tomorrow) * exp((suns - ti) / tcoef))
//               / (1 - exp((dayl - 24) / tcoef))
//  where
//        tcoef is a time coefficient, determined by calibration as 4
//        sst is the sunset temperature, determined by the daytime equation as:
//        sst = mint_tomorrow - tkk / 2 + 0.5 * sqrt(tkk**2 + 4 * amp * tkk * sts)
//  where
//        sts  = sin(pi * dayl / (dayl + 2 * SitePar[8]))
//        amp = (tmax - mint_tomorrow) * (1 + (tmax - mint_tomorrow) / tkk)
//      For the time from midnight to sunrise, similar equations are used, but the minimum 
//  temperature of this day (mint_today) is used instead of mint_tomorrow, and the maximum 
//  temperature of the previous day (maxt_yesterday) is used instead of maxt_today. Also, 
//  (suns-ti-24) is used for the time variable instead of (suns-ti).
//      These exponential equations for night-time temperature ensure that the curve will 
//  be continuous with the daytime equation at sunset, and will pass through the minimum 
//  temperature at sunrise.
//
//  Input argument:
//     ti - time of day (hours).
//  Global variables used: 
//     DayLength, LastDayWeatherData, pi, SitePar, SolarNoon, sunr, suns
//
{
    const double tkk = 15;     // The temperature increase at which the sensible heat flux is
    //  doubled, in comparison with the situation without buoyancy.
    const double tcoef = 4;    // time coefficient for the exponential part of the equation.
    double hmax = SolarNoon + SitePar[8]; // hour of maximum temperature
    int im1 = Daynum - 1;     // day of year yesterday
    int ip1 = Daynum + 1;     // day of year tomorrow
    if (ip1 > LastDayWeatherData)
        ip1 = Daynum;
//
    double amp;          // amplitude of temperatures for a period.
    double sst;          // the temperature at sunset.
    double st;           // computed from time of day, used for daytime temperature.
    double sts;          // intermediate variable for computing sst.
    double HourlyTemperature;   // computed temperature at time ti.
//
    if (ti <= sunr)        //  from midnight to sunrise
    {
        amp = (GetFromClim(Clim, "tmax", im1) - GetFromClim(Clim, "tmin", Daynum))
              * (1 + (GetFromClim(Clim, "tmax", im1) - GetFromClim(Clim, "tmin", Daynum)) / tkk);
        sts = sin(pi * DayLength / (DayLength + 2 * SitePar[8]));
//  compute temperature at sunset:
        sst = GetFromClim(Clim, "tmin", Daynum) - tkk / 2 + 0.5 * sqrt(tkk * tkk + 4 * amp * tkk * sts);
        HourlyTemperature = (GetFromClim(Clim, "tmin", Daynum) - sst * exp((DayLength - 24) / tcoef)
                             + (sst - GetFromClim(Clim, "tmin", Daynum)) * exp((suns - ti - 24) / tcoef))
                            / (1 - exp((DayLength - 24) / tcoef));
    } else if (ti <= hmax)    //  from sunrise to hmax
    {
        amp = (GetFromClim(Clim, "tmax", Daynum) - GetFromClim(Clim, "tmin", Daynum))
              * (1 + (GetFromClim(Clim, "tmax", Daynum) - GetFromClim(Clim, "tmin", Daynum)) / tkk);
        st = sin(pi * (ti - SolarNoon + DayLength / 2.) / (DayLength + 2 * SitePar[8]));
        HourlyTemperature = GetFromClim(Clim, "tmin", Daynum) - tkk / 2
                            + 0.5 * sqrt(tkk * tkk + 4 * amp * tkk * st);
    } else if (ti <= suns)    //  from hmax to sunset
    {
        amp = (GetFromClim(Clim, "tmax", Daynum) - GetFromClim(Clim, "tmin", ip1))
              * (1 + (GetFromClim(Clim, "tmax", Daynum) - GetFromClim(Clim, "tmin", ip1)) / tkk);
        st = sin(pi * (ti - SolarNoon + DayLength / 2) / (DayLength + 2 * SitePar[8]));
        HourlyTemperature = GetFromClim(Clim, "tmin", ip1) - tkk / 2
                            + 0.5 * sqrt(tkk * tkk + 4 * amp * tkk * st);
    } else                   //  from sunset to midnight
    {
        amp = (GetFromClim(Clim, "tmax", Daynum) - GetFromClim(Clim, "tmin", ip1))
              * (1 + (GetFromClim(Clim, "tmax", Daynum) - GetFromClim(Clim, "tmin", ip1)) / tkk);
        sts = sin(pi * DayLength / (DayLength + 2 * SitePar[8]));
        sst = GetFromClim(Clim, "tmin", ip1) - tkk / 2 + 0.5 * sqrt(tkk * tkk + 4 * amp * tkk * sts);
        HourlyTemperature = (GetFromClim(Clim, "tmin", ip1) - sst * exp((DayLength - 24) / tcoef)
                             + (sst - GetFromClim(Clim, "tmin", ip1)) * exp((suns - ti) / tcoef))
                            / (1. - exp((DayLength - 24) / tcoef));
    }
    return HourlyTemperature;
//     Reference:
//     Ephrath, J.E., Goudriaan, J. and Marani, A. 1996. Modelling
//  diurnal patterns of air temperature, radiation, wind speed and
//  relative humidity by equations from daily characteristics.
//  Agricultural Systems 51:377-393.
}

///////////////////////////////////////////////////////////////////////////////////////////
double tdewhour(double ti, double tt, const int &Daynum, const Climstruct Clim[400])
//     Function tdewhour() computes the hourly values of dew point 
//  temperature from average dew-point and the daily estimated range. 
//  This range is computed as a regression on maximum and minimum temperatures.
//  Input arguments:
//     ti - time of day (hours).
//     tt - air temperature C at this time of day.
//  Global variables used: 
//    LastDayWeatherData, SitePar, SolarNoon, sunr, suns.
//
{
    int im1 = Daynum - 1;     // day of year yeaterday
    int ip1 = Daynum + 1;     // day of year tomorrow
    if (ip1 > LastDayWeatherData)
        ip1 = Daynum;
    double tdewhr;  // the dew point temperature (c) of this hour.
    double tdmin;   // minimum of dew point temperature.
    double tdrange; // range of dew point temperature.
    double hmax = SolarNoon + SitePar[8]; // time of maximum air temperature
//
    if (ti <= sunr)       // from midnight to sunrise
    {
        tdrange = SitePar[12] + SitePar[13] * GetFromClim(Clim, "tmax", im1)
                  + SitePar[14] * GetFromClim(Clim, "tmin", Daynum);
        if (tdrange < 0)
            tdrange = 0;
        tdmin = GetFromClim(Clim, "tdew", im1) - tdrange / 2;
        tdewhr = tdmin + tdrange * (tt - GetFromClim(Clim, "tmin", Daynum))
                         / (GetFromClim(Clim, "tmax", im1) - GetFromClim(Clim, "tmin", Daynum));
    } else if (ti <= hmax)     // from sunrise to hmax
    {
        tdrange = SitePar[12] + SitePar[13] * GetFromClim(Clim, "tmax", Daynum)
                  + SitePar[14] * GetFromClim(Clim, "tmin", Daynum);
        if (tdrange < 0) tdrange = 0;
        tdmin = GetFromClim(Clim, "tdew", Daynum) - tdrange / 2;
        tdewhr = tdmin + tdrange * (tt - GetFromClim(Clim, "tmin", Daynum))
                         / (GetFromClim(Clim, "tmax", Daynum) - GetFromClim(Clim, "tmin", Daynum));
    } else     //  from hmax to midnight
    {
        tdrange = SitePar[12] + SitePar[13] * GetFromClim(Clim, "tmax", Daynum)
                  + SitePar[14] * GetFromClim(Clim, "tmin", ip1);
        if (tdrange < 0) tdrange = 0;
        tdmin = GetFromClim(Clim, "tdew", ip1) - tdrange / 2;
        tdewhr = tdmin + tdrange * (tt - GetFromClim(Clim, "tmin", ip1))
                         / (GetFromClim(Clim, "tmax", Daynum) - GetFromClim(Clim, "tmin", ip1));
    }
    return tdewhr;
}

///////////////////////////////////////////////////////////////////////////////////////////
double daywnd(double ti, double wind, double t1, double t2, double t3, double wnytf)
//     The daywnd function computes the hourly values of wind speed
//  (m / sec), estimated from the measured total daily wind run.
//     Input arguments:
//       t1 = the hour at which day-time wind begins to blow.
//       t2 = the hour at which day-time wind speed is maximum.
//       t3 = the hour at which day-time wind ceases to blow.
//       ti = the hour of the day.
//       wind = the daily total wind run (km per day).
//       wnytf = Factor for estimating night-time wind (from time
//               t3 to time t1 next day).
//
//     The algorithm is described by Ephrath et al. (1996). It is based on the 
//	following assumptions:
//     Although the variability of wind speed during any day is very
//  large, the diurnal wind speed curves appear to be characterized by
//  the following repetitive pattern: increase in wind speed from time
//  t1 in the morning to time t2 in the afternoon, decrease from t2 to t3
//  in the evening, and a low constant wind speed at night, from t3 to t1
//  in the next day.
//     The values of t1, t2, and t3 have been determined in the calling routine:
//  t1 is SitePar(1) hours after sunrise, t2 is SitePar(2) hours after solar
//  noon, and t3 is SitePar(3) hours after sunset. These parameters are site-
//  specific. They are 1, 3, and 0, respectively, for the San Joaquin valley of 
//  California and for Arizona, and 1, 4, and 2, respectively, for the coastal 
//  plain of israel.
//     The wind speed during the night, from t3 to t1 next day (wmin)
//  is assumed to be proportional to the daily total wind run. The
//  ratio wnytf is also site-specific, SitePar(4), ( 0.008 for San Joaquin 
//  and Arizona, 0.0025 for the coastal plain of Israel). wmin is the minimum 
//  wind speed from t1 to t3.
//      wtday is computed by subtracting the daily integral of wmin, after 
//  converting it from m/sec to km/day, from the total daily wind run (wndt).
//      wmax, the maximum wind speed at time t2 (minus wmin), is
//  computed from wtday and converted to m/sec.
//      daywnd from t1 to t2 is now computed as an increasing
//  sinusoidal function from wmin to wmin + wmax, and it is computed from
//  t2 to t3 as a decreasing sinusoidal function from wmin + wmax to  wmin.
//
{
    double HourlyWind;
//   constants related to t1, t2, t3 :
    double sf1 = 4 * (t2 - t1);
    double sf2 = 4 * (t3 - t2);
    double wmin = wind * wnytf;  //  the constant minimum wind speed during the night (m/sec).
    double wtday = wind - wmin * 3.6 * 24;  //  integral of wind run from t1 to t3, minus wmin (km).
    double wmax = wtday * 2 * pi / 3.6 / (sf1 + sf2);  //  the maximum wind speed (m per sec), above wmin.
//
    if (ti >= t1 && ti < t2)
        HourlyWind = wmin + wmax * sin(2 * pi * (ti - t1) / sf1);
    else if (ti >= t2 && ti < t3)
        HourlyWind = wmin + wmax * sin(2 * pi * (ti - (2 * t2 - t3)) / sf2);
    else if (ti >= t3 || ti < t1)
        HourlyWind = wmin;
    return HourlyWind;
//
//     Reference:
//     Ephrath, J.E., Goudriaan, J. and Marani, A. 1996. Modelling
//  diurnal patterns of air temperature, radiation, wind speed and
//  relative humidity by equations from daily characteristics.
//  Agricultural Systems 51:377-393.
}

///////////////////////////////////////////////////////////////////////////////////////////
void AverageAirTemperatures()
//     Function AverageAirTemperatures() calculates daily average temperatures, daytime
//  average and night time average. 
//     Global variables referenced: 
//        AirTemp[], Radiation[].
//     Global variables computed: 
//        AvrgDailyTemp, NightTimeTemp, DayTimeTemp.
//
{
    int nn1 = 0;  // counter of night hours
    int nn2 = 0;  // counter of daytime hours
    AvrgDailyTemp = 0;
    NightTimeTemp = 0;
    DayTimeTemp = 0;
//
    for (int ihr = 0; ihr < 24; ihr++)   // hourly loop
    {
        AvrgDailyTemp += AirTemp[ihr];
        if (Radiation[ihr] <= 0) {
            NightTimeTemp += AirTemp[ihr];
            nn1++;
        } else {
            DayTimeTemp += AirTemp[ihr];
            nn2++;
        }
    }
//
    AvrgDailyTemp = AvrgDailyTemp / 24;
    NightTimeTemp = NightTimeTemp / nn1;
    DayTimeTemp = DayTimeTemp / nn2;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
void
EvapoTranspiration(int jtout, const string &ProfileName, const string &Date, const int &Daynum, const double &Latitude,
                   const double &DayLength)
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
    const double stefb = 5.77944E-08; // the Stefan-Boltzman constant, in W m-2 K-4 (= 1.38E-12 * 41880)
    const double c12 = 0.125;         // c12 ... c15 are constant parameters.
    const double c13 = 0.0439;
    const double c14 = 0.030;
    const double c15 = 0.0576;
    int iamhr = 0;  // earliest time in day for computing cloud cover
    int ipmhr = 0;  // latest time in day for computing cloud cover
    double cosz;    // cosine of sun angle from zenith for this hour
    double suna;    // sun angle from horizon, degrees at this hour
    double isr;     // hourly extraterrestrial radiation in W / m**2
//      Start hourly loop
    for (int ihr = 0; ihr < 24; ihr++) {
        double ti = ihr + 0.5;           // middle of the hourly interval
//      The following subroutines and functions are called for each
//  hour: sunangle, cloudcov, clcor, refalbed . 
        sunangle(ti, cosz, suna, Latitude);
        isr = tmpisr * cosz;
        CloudCoverRatio[ihr] = cloudcov(Radiation[ihr], isr, cosz);
//      clcor is called to compute cloud-type correction.
//      iamhr and ipmhr are set.
        CloudTypeCorr[ihr] = clcor(ihr, SitePar[7], isr, cosz, DayLength);
        if ((cosz >= 0.1736) && (iamhr == 0))
            iamhr = ihr;
        if ((ihr >= 12) && (cosz <= 0.1736) && (ipmhr == 0))
            ipmhr = ihr - 1;
//      refalbed is called to compute the reference albedo for each hour.
        albedo[ihr] = refalbed(isr, Radiation[ihr], cosz, suna);
        if (jtout > 1) {
// write output to file
            ofstream File18(fs::path("output") / (ProfileName + ".TM2"), ios::app);
            if (ihr == 0) {
                File18 << "    ***********" << endl << "  date =";
                File18.width(12);
                File18 << Date << " Day of Year =";
                File18.width(4);
                File18 << Daynum << endl;
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
    }       //   End of 1st hourly loop
//     Zero some variables that will later be used for summation.
    ReferenceTransp = 0;
    Rn = 0; // daily net radiation
    double rnet[24]; // hourly net radiation
    double cltcoram = CloudTypeCorr[iamhr];  //  cloud type correction during early morning
    double cltcorpm = CloudTypeCorr[ipmhr];  //  cloud type correction during late afternoon
//
    for (int ihr = 0; ihr < 24; ihr++)  //  2nd hourly loop
    {
        double ti = ihr + 0.5;           // middle of the hourly interval
//      Compute saturated vapor pressure (svp), using function VaporPressure().
//      The actual vapor pressure (vp) is computed from svp and the
//  relative humidity. Compute vapor pressure deficit (vpd). This
//  procedure is based on the CIMIS algorithm.
        double svp = VaporPressure(AirTemp[ihr]);  // saturated vapor pressure, mb
        double vp = 0.01 * RelativeHumidity[ihr] * svp;  // vapor pressure, mb
        double vpd = svp - vp;             // vapor pressure deficit, mb.
//   Get cloud cover and cloud correction for night hours
        if (ihr < iamhr || ihr > ipmhr) {
            CloudCoverRatio[ihr] = 0;
            CloudTypeCorr[ihr] = 0;
        }
//     The hourly net radiation is computed using the CIMIS algorithm (Dong et al., 1988):
//     rlonin, the hourly incoming long wave radiation, is computed from ea0, cloud cover 
//  (CloudCoverRatio), air temperature (tk),  stefb, and cloud type correction (CloudTypeCorr).
//     rnet, the hourly net radiation, W m-2, is computed from the global radiation, the albedo,
//  the incoming long wave radiation, and the outgoing longwave radiation.
        double tk = AirTemp[ihr] + 273.161;  // hourly air temperature in Kelvin.
        double ea0 = clearskyemiss(vp, tk);  // clear sky emissivity for long wave radiation
//     Compute incoming long wave radiation:
        double rlonin = (ea0 * (1 - CloudCoverRatio[ihr]) + CloudCoverRatio[ihr])
                        * stefb * pow(tk, 4) - CloudTypeCorr[ihr];
        rnet[ihr] = (1 - albedo[ihr]) * Radiation[ihr] + rlonin - stefb * pow(tk, 4);
        Rn += rnet[ihr];
//     The hourly reference evapotranspiration ReferenceETP is computed by the
//  CIMIS algorithm using the modified Penman equation:
//     The weighting ratio (w) is computed from the functions del() (the slope of the saturation
//  vapor pressure versus air temperature) and gam() (the psychometric constant).
        double w = del(tk, svp) / (del(tk, svp) + gam(Elevation, AirTemp[ihr])); // coefficient of the Penman equation
//     The wind function (fu2) is computed using different sets of parameters
//  for day-time and night-time. The parameter values are as suggested by CIMIS.
        double fu2;  // wind function for computing evapotranspiration
        if (Radiation[ihr] <= 0)
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
    }  //   end of 2nd hourly loop
}

/////////////////////////////////////////////////////////////////////////////////
double clearskyemiss(double vp, double tk)
//     Function clearskyemiss() estimates clear sky emissivity for long wave radiation.
//  Input arguments:
//     vp - vapor pressure of the air in KPa
//     tk - air temperature in K.
//
{
//     Convert vp to mbars
    double vp1 = vp * 10; // vapor pressure of the air in mbars.
//     Compute clear sky emissivity by the method of Idso (1981)
    double ea0 = 0.70 + 5.95e-05 * vp1 * exp(1500 / tk);
    if (ea0 > 1)
        ea0 = 1;
    return ea0;
}

//    Reference:
//      Idso, S.B. 1981. A set of equations for full spectrum and 8- to 14-um
// and 10.5- to 12.5- um thermal radiation from cloudless skies. Water
// Resources Res. 17:295.
/////////////////////////////////////////////////////////////////////////////////
double cloudcov(double radihr, double isr, double cosz)
//     Function cloudcov() computes cloud cover for this hour from radiation data, using
//  the CIMIS algorithm. The return value is cloud cover ratio ( 0 to 1 )
//  Input arguments:
//       radihr = hourly global radiation in W m-2 .
//       isr = hourly extraterrestrial radiation in W m-2 . 
//       cosz = cosine of sun angle from zenith.
//     This algorithm is described by Dong et al. (1988). Cloud cover fraction is estimated 
//  as a function of the ratio of actual solar radiation to extraterrestrial radiation. 
//  The parameters of this function have been based on California data.
//      The equation is for daylight hours, when the sun is not less
//  than 10 degrees above the horizon (coszhr > 0.1736).
//
{
    double p1 = 1.333;  //    p1, p2, p3 are constant parameters.
    double p2 = 1.7778;
    double p3 = 0.294118;
    double rasi = 0;   // ratio of radihr to isr.
    double clcov = 0;  // computed cloud cover.
//
    if (isr > 0)
        rasi = radihr / isr;
    if (cosz > 0.1736 && rasi <= p1 / p2) {
        if (rasi >= 0.375)
            clcov = pow((p1 - p2 * rasi), p3);
        else
            clcov = pow((p1 - p2 * 0.375), p3);
        if (clcov < 0)
            clcov = 0;
    }
    return clcov;
}

//      Reference:
//      Dong, A., Prashar, C.K. and Grattan, S.R. 1988. Estimation of
// daily and hourly net radiation. CIMIS Final Report June 1988, pp.
// 58-79.
/////////////////////////////////////////////////////////////////////////////////
double clcor(int ihr, double ck, double isrhr, double coszhr, const double &DayLength)
//     Function clcor() computes cloud type correction, using the CIMIS algorithm.
//     Input arguments:
//            ck = cloud type correction factor (data for this location).
//            coszhr = cosine of sun angle from zenith.
//            ihr = time of day, hours.
//            isrhr = hourly extraterrestrial radiation in W m-2 .
//     Global variables used: DayLength, Radiation[], SolarNoon, pi
//     Note: This algorithm is described by Dong et al. (1988). ck is the
//  cloud-type correction used in the Monteith equation for estimating
//  net radiation. The value of this correction depends on site and
//  time of year. Regional ck values for California are given by Dong
//  et al. (1988). In the San Joaquin valley of California ck is almost
//  constant from April to October, with an average value of 60. The value
//  of ck is site-dependant, assumed to be constant during the growing season.
//     The daily ck is converted to an hourly value for
//  clear or partly cloudy sky (rasi >= 0.375) and when the sun is at
//  least 10 degrees above the horizon.
//     Evening, night and early morning cloud type correction is temporarily
//  assigned 0. It is later assigned the values of first or 
// last non-zero values (in the calling routine). 
//
{
    double rasi = 0;  //  ratio of Radiation to isrhr.
    if (isrhr > 0)
        rasi = Radiation[ihr] / isrhr;
    if (coszhr >= 0.1736 && rasi >= 0.375) {
        double angle = pi * (ihr - SolarNoon + 0.5) / DayLength; // hour angle (from solar noon) in radians.
        return ck * pi / 2 * cos(angle);
    } else return 0;
}

//      Reference:
//      Dong, A., Prashar, C.K. and Grattan, S.R. 1988. Estimation of
// daily and hourly net radiation. CIMIS Final Report June 1988, pp.
// 58-79.
/////////////////////////////////////////////////////////////////////////////////
double del(double tk, double svp)
//     Function del() computes the slope of the saturation vapor
//  pressure (svp, in mb) versus air temperature (tk, in K).
//     This algorithm is the same as used by CIMIS.
//
{
    double a = pow((double) 10, (-0.0304 * tk));
    double b = pow(tk, 2);
    double c = pow((double) 10, (-1302.88 / tk));
    return (6790.5 - 5.02808 * tk + 4916.8 * a * b + 174209 * c) * svp / b;
}

/////////////////////////////////////////////////////////////////////////////////
double gam(double elev, double tt)
//     Function gam() computes the psychometric constant at elevation (elev), m above
//  sea level, and air temperature, C (tt).
//     This algorithm is the same as used by CIMIS.
//
{
    double bp = 101.3 - .01152 * elev + 5.44e-07 * pow(elev, 2); //  barometric pressure, KPa, at this elevation.
    return 0.000646 * bp * (1 + 0.000946 * tt);
}

/////////////////////////////////////////////////////////////////////////////////
double refalbed(double isrhr, double rad, double coszhr, double sunahr)
//     Function refalbed() computes the reference crop albedo, using the
//  CIMIS algorithm.
//     This algorithm is described by Dong et al. (1988). Albedo is
//  estimated as a function of sun elevation above the horizon (suna)
//  for clear or partly cloudy sky (rasi >= 0.375) and when the sun is
//  at least 10 degrees above the horizon. 
//     For very cloudy sky, or when solar altitude is below 10
//  degrees, the following albedo value is assumed: (p4)+ 0.26
//     Input arguments:
//        isrhr = hourly extraterrestrial radiation in W m-2 .
//        rad = hourly global radiation in W / m-2 .
//        coszhr = cosine of sun angle from zenith.
//        sunahr = sun angle from horizon, degrees.
//
{
    const double p1 = 0.00158; //  p1 ... p4 are constant parameters.
    const double p2 = 0.386;
    const double p3 = 0.0188;
    const double p4 = 0.26;
    double rasi = 0;     //   ratio of rad to isrhr
    double refalb = p4;  // the reference albedo
//
    if (isrhr > 0)
        rasi = rad / isrhr;
    if (coszhr > 0.1736 && rasi >= 0.375) {
        refalb = p1 * sunahr + p2 * exp(-p3 * sunahr);
        if (refalb > p4)
            refalb = p4;
    }
    return refalb;
}

//      Reference:
//      Dong, A., Prashar, C.K. and Grattan, S.R. 1988. Estimation of
// daily and hourly net radiation. CIMIS Final Report June 1988, pp. 58-79.
/////////////////////////////////////////////////////////////////////////////////
void sunangle(double ti, double &coszhr, double &sunahr, const double &Latitude)
//     sunangle.cpp : computes sun angle for any time of day. 
//     Input argument:
//        ti = time of day, hours.
//      Output arguments:
//        coszhr = cosine of sun angle from zenith for this hour.
//        sunahr = sun angle from horizon, degrees.
//
{
//      The latitude is converted to radians (xlat). 
    double xlat = Latitude * pi / 180;           // latitude in radians.
    double cd;  // amplitude of the sine of the solar height, computed as
    // the product of cosines of latitude and declination angles.
    cd = cos(xlat) * cos(declination);
    double sd;  // seasonal offset of the sine of the solar height, computed
    // as the product of sines of latitude and declination angles.
    sd = sin(xlat) * sin(declination);
    double hrangle = 15 * (ti - SolarNoon) * pi / 180; // hourly angle converted to radians
    coszhr = sd + cd * cos(hrangle);
    if (coszhr <= 0) {
        coszhr = 0;
        sunahr = 0;
    } else if (coszhr >= 1) {
        coszhr = 1;
        sunahr = 90;
    } else
        sunahr = fabs(acos(coszhr) * 180 / pi - 90);
}

///////////////////////////////////////////////////////////////////////////////
double SimulateRunoff(double rain, const int &Daynum, const int &DayStart, const Climstruct Clim[400])
//     This function is called from DayClim() and is executed on each day with raifall more 
//  than 2 mm. It computes the runoff and the retained portion of the rainfall. Note: This
//  function is based on the code of GOSSYM. No changes have been made from the original GOSSYM 
//  code (except translation to C++). It has not been validated by actual field measurement.
//     It calculates the portion of rainfall that is lost to runoff, and reduces rainfall to the
//  amount which is actually infiltrated into the soil. It uses the soil conservation service 
//  method of estimating runoff.                                               
//     References:
//  - Brady, Nyle C. 1984. The nature and properties of soils, 9th ed. Macmillan Publishing Co.
//  - Schwab, Frevert, Edminster, and Barnes. 1981. Soil and water conservation engineering, 
//  3rd ed. John Wiley & Sons, Inc.                               
//
//     The following global variables are referenced here:
//  ClayVolumeFraction, Irrig (structure), NumIrrigations, SandVolumeFraction.
//     The argument used here:  rain = today,s rainfall.
//     The return value:  the amount of water (mm) lost by runoff.
{
    static bool bFirst = true; // if this is the first time the function is called.
    static int iGroup;   // soil group number (by clay and sand in upper soil layer)
    static double d01;    // Adjustment of curve number for soil groups A,B,C.
//     The following is computed only the first time the function is called.
//     Infiltration rate is estimated from the percent sand and percent clay in the Ap layer.
//  If clay content is greater than 35%, the soil is assumed to have a higher runoff potential,
//  if clay content is less than 15% and sand is greater than 70%, a lower runoff potential is 
//  assumed. Other soils (loams) assumed moderate runoff potential. No 'impermeable' (group D)
//  soils are assumed.  References: Schwab, Brady.
    if (bFirst) {
        if (SandVolumeFraction[0] > 0.70 && ClayVolumeFraction[0] < 0.15) {
//     Soil group A = 1, low runoff potential
            iGroup = 1;
            d01 = 1.0;
        } else if (ClayVolumeFraction[0] > 0.35) {
//     Soil group C = 3, high runoff potential
            iGroup = 3;
            d01 = 1.14;
        } else {
//     Soil group B = 2, moderate runoff potential
            iGroup = 2;
            d01 = 1.09;
        }
        bFirst = false;
    }
//     Loop to accumulate 5-day antecedent rainfall (mm) which will affect the soil's ability
//  to accept new rainfall. This also includes all irrigations.
    int i01 = Daynum - 5;
    if (i01 < DayStart)
        i01 = DayStart;
    int i02 = Daynum;
    double PreviousWetting = 0; // five day total (before this day) of rain and irrigation, mm
    for (int Dayn = i01; Dayn < i02; Dayn++) {
        double amtirr = 0; // mm water applied on this day by irrigation
        for (int i = 0; i < NumIrrigations; i++) {
            if (Dayn == Irrig[i].day)
                amtirr = Irrig[i].amount;
        }
        PreviousWetting += amtirr + GetFromClim(Clim, "rain", Dayn);
    }
//
    double d02; // Adjusting curve number for antecedent rainfall conditions.
    if (PreviousWetting < 36) {
//  low moisture, low runoff potential.
        if (iGroup == 1)
            d02 = 0.71;
        else if (iGroup == 2)
            d02 = 0.78;
        else if (iGroup == 3)
            d02 = 0.83;
    } else if (PreviousWetting > 53) {
//  wet conditions, high runoff potential.
        if (iGroup == 1)
            d02 = 1.24;
        else if (iGroup == 2)
            d02 = 1.15;
        else if (iGroup == 3)
            d02 = 1.10;
    } else {
//  moderate conditions
        d02 = 1.00;
    }
//  Assuming straight rows, and good cropping practice:
    double crvnum = 78.0; // Runoff curve number, unadjusted for moisture and soil type.
    crvnum = crvnum * d01 * d02; // adjusted curve number
    double d03; // maximum potential difference between rainfall and runoff.
    d03 = 25400 / crvnum - 254;
//
    double runoff; // simulated runoff on this day
    if (rain <= 0.2 * d03)
        runoff = 0; // no runoff
    else
        runoff = pow((rain - 0.2 * d03), 2) / (rain + 0.8 * d03);
    return runoff;
}