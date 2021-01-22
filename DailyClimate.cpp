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
    double refalbed(double, double, double, double);
    double del(double, double);
    double gam(double, double);
    double cloudcov(double, double, double);
    double daywnd(double, double, double, double, double, double);
    double clcor(uint8_t, double, double, double, double, double, double);
}

double daytmp(Simulation &, uint32_t, double, const double &);

double tdewhour(Simulation &, uint32_t, double, double);

void AverageAirTemperatures();

void EvapoTranspiration(int, const string &, const string &, const int &, const double &, const double &);

void sunangle(double, double &, double &, const double &);

double SimulateRunoff(Simulation &, uint32_t, double);

//     Definition of file scope variables:
double declination, // daily declination angle, in radians.
    SolarNoon,      // time of solar noon, hours.
    sunr,           // time of sunrise, hours.
    suns,           // time of sunset, hours.
    tmpisr;         // extraterrestrial radiation, W / m2.
//////////////////////////////////////////////////////////////////////////////
tuple<double> DayClim(Simulation &sim, uint32_t u)
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
    //     Compute day length and related variables:
    double DayLength;
    tie(DayLength) = ComputeDayLength(sim.day_start + u, sim.latitude, sim.longitude);
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
        runoffToday = SimulateRunoff(sim, u, rainToday);
        if (runoffToday < rainToday)
            rainToday -= runoffToday;
        else
            rainToday = 0;
        sim.climate[u].Rain = rainToday;
    }
    Scratch21[u].runoff = runoffToday;
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
        Radiation[ihr] = dayrad(ti, radsum, sinb, c11);
        //     Compute hourly temperature, using function daytmp.
        AirTemp[ihr] = daytmp(sim, u, ti, DayLength);
        //     Compute hourly dew point temperature, using function tdewhour.
        DewPointTemp[ihr] = tdewhour(sim, u, ti, AirTemp[ihr]);
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
    EvapoTranspiration(jtout, sim.profile_name, sim.states[u].date, sim.day_start + u, sim.latitude, DayLength);
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
    declination = .006918 - .399912 * cos(xday) + .070257 * sin(xday) - .006758 * cos(2 * xday) + .000907 * sin(2 * xday) - .002697 * cos(3 * xday) + .001480 * sin(3 * xday);
    //     Compute extraterrestrial radiation in W m-2. The 'solar constant' (average
    //  value = 1367 W m-2) is corrected for this day's distance between earth and the sun. The
    //  equation used here is from the CIMIS algorithm, which is based on the work of Iqbal (1983).
    tmpisr = 1367 * (1.00011 + .034221 * cos(xday) + .00128 * sin(xday) + .000719 * cos(2 * xday) + .000077 * sin(2 * xday));
    //     Time of solar noon (SolarNoon) is computed by the CIMIS algorithm,
    //  using a correction for longitude (f), and the date correction (exday).
    //     It is assumed that the time zone is "geographically correct". For
    //  example, longitude between 22.5 and 37.5 East is in time zone GMT+2,
    //  and longitude between 112.5 and 127.5 West is in time zone GMT-8.
    //     All daily times in the model are computed by this method.
    double exday = (.000075 + .001868 * cos(xday) - .032077 * sin(xday) - .014615 * cos(2 * xday) - .04089 * sin(2 * xday)) * 12. / pi;
    double st = 15 * (int)(Longitude / 15);
    double f = (Longitude - st) / 15;
    if (Longitude > 0) // east of Greenwich
    {
        if (f > 0.5)
            f -= 1;
    }
    else // west  of Greenwich
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
        ht = 1; //  arctic circle
    else if (ht < -1)
        ht = -1;
    double DayLength = 2 * acos(ht) * 12 / pi;
    sunr = SolarNoon - DayLength / 2;
    suns = sunr + DayLength;
    return make_tuple(DayLength);
}

///////////////////////////////////////////////////////////////////////////////////////////
double daytmp(Simulation &sim, uint32_t u, double ti, const double &DayLength)
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
    const double tkk = 15; // The temperature increase at which the sensible heat flux is
    //  doubled, in comparison with the situation without buoyancy.
    const double tcoef = 4;               // time coefficient for the exponential part of the equation.
    double hmax = SolarNoon + SitePar[8]; // hour of maximum temperature
    int im1 = u > 1 ? u - 1 : 0;      // day of year yesterday
    ClimateStruct &yesterday = sim.climate[im1];
    ClimateStruct &today = sim.climate[u];
    int ip1 = u + 1; // day of year tomorrow
    if (ip1 > LastDayWeatherData)
        ip1 = u;
    ClimateStruct &tomorrow = sim.climate[ip1];
    //
    double amp;               // amplitude of temperatures for a period.
    double sst;               // the temperature at sunset.
    double st;                // computed from time of day, used for daytime temperature.
    double sts;               // intermediate variable for computing sst.
    double HourlyTemperature; // computed temperature at time ti.
                              //
    if (ti <= sunr)           //  from midnight to sunrise
    {
        amp = (yesterday.Tmax - today.Tmin) * (1 + (yesterday.Tmax - today.Tmin) / tkk);
        sts = sin(pi * DayLength / (DayLength + 2 * SitePar[8]));
        //  compute temperature at sunset:
        sst = today.Tmin - tkk / 2 + 0.5 * sqrt(tkk * tkk + 4 * amp * tkk * sts);
        HourlyTemperature = (today.Tmin - sst * exp((DayLength - 24) / tcoef) + (sst - today.Tmin) * exp((suns - ti - 24) / tcoef)) / (1 - exp((DayLength - 24) / tcoef));
    }
    else if (ti <= hmax) //  from sunrise to hmax
    {
        amp = (today.Tmax - today.Tmin) * (1 + (today.Tmax - today.Tmin) / tkk);
        st = sin(pi * (ti - SolarNoon + DayLength / 2.) / (DayLength + 2 * SitePar[8]));
        HourlyTemperature = today.Tmin - tkk / 2 + 0.5 * sqrt(tkk * tkk + 4 * amp * tkk * st);
    }
    else if (ti <= suns) //  from hmax to sunset
    {
        amp = (today.Tmax - tomorrow.Tmin) * (1 + (today.Tmax - tomorrow.Tmin) / tkk);
        st = sin(pi * (ti - SolarNoon + DayLength / 2) / (DayLength + 2 * SitePar[8]));
        HourlyTemperature = tomorrow.Tmin - tkk / 2 + 0.5 * sqrt(tkk * tkk + 4 * amp * tkk * st);
    }
    else //  from sunset to midnight
    {
        amp = (today.Tmax - tomorrow.Tmin) * (1 + (today.Tmax - tomorrow.Tmin) / tkk);
        sts = sin(pi * DayLength / (DayLength + 2 * SitePar[8]));
        sst = tomorrow.Tmin - tkk / 2 + 0.5 * sqrt(tkk * tkk + 4 * amp * tkk * sts);
        HourlyTemperature = (tomorrow.Tmin - sst * exp((DayLength - 24) / tcoef) + (sst - tomorrow.Tmin) * exp((suns - ti) / tcoef)) / (1. - exp((DayLength - 24) / tcoef));
    }
    return HourlyTemperature;
    //     Reference:
    //     Ephrath, J.E., Goudriaan, J. and Marani, A. 1996. Modelling
    //  diurnal patterns of air temperature, radiation, wind speed and
    //  relative humidity by equations from daily characteristics.
    //  Agricultural Systems 51:377-393.
}

///////////////////////////////////////////////////////////////////////////////////////////
double tdewhour(Simulation &sim, uint32_t u, double ti, double tt)
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
    int im1 = u > 1 ? u - 1 : 0;      // day of year yesterday
    ClimateStruct &yesterday = sim.climate[im1];
    ClimateStruct &today = sim.climate[u];
    int ip1 = u + 1; // day of year tomorrow
    if (ip1 > LastDayWeatherData)
        ip1 = u;
    ClimateStruct &tomorrow = sim.climate[ip1];
    double tdewhr;                        // the dew point temperature (c) of this hour.
    double tdmin;                         // minimum of dew point temperature.
    double tdrange;                       // range of dew point temperature.
    double hmax = SolarNoon + SitePar[8]; // time of maximum air temperature
                                          //
    if (ti <= sunr)                       // from midnight to sunrise
    {
        tdrange = SitePar[12] + SitePar[13] * yesterday.Tmax + SitePar[14] * today.Tmin;
        if (tdrange < 0)
            tdrange = 0;
        tdmin = yesterday.Tdew - tdrange / 2;
        tdewhr = tdmin + tdrange * (tt - today.Tmin) / (yesterday.Tmax - today.Tmin);
    }
    else if (ti <= hmax) // from sunrise to hmax
    {
        tdrange = SitePar[12] + SitePar[13] * today.Tmax + SitePar[14] * today.Tmin;
        if (tdrange < 0)
            tdrange = 0;
        tdmin = today.Tdew - tdrange / 2;
        tdewhr = tdmin + tdrange * (tt - today.Tmin) / (today.Tmax - today.Tmin);
    }
    else //  from hmax to midnight
    {
        tdrange = SitePar[12] + SitePar[13] * today.Tmax + SitePar[14] * tomorrow.Tmin;
        if (tdrange < 0)
            tdrange = 0;
        tdmin = tomorrow.Tdew - tdrange / 2;
        tdewhr = tdmin + tdrange * (tt - tomorrow.Tmin) / (today.Tmax - tomorrow.Tmin);
    }
    return tdewhr;
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
    int nn1 = 0; // counter of night hours
    int nn2 = 0; // counter of daytime hours
    AvrgDailyTemp = 0;
    NightTimeTemp = 0;
    DayTimeTemp = 0;
    //
    for (int ihr = 0; ihr < 24; ihr++) // hourly loop
    {
        AvrgDailyTemp += AirTemp[ihr];
        if (Radiation[ihr] <= 0)
        {
            NightTimeTemp += AirTemp[ihr];
            nn1++;
        }
        else
        {
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
void EvapoTranspiration(int jtout, const string &ProfileName, const string &Date, const int &Daynum, const double &Latitude,
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
        sunangle(ti, cosz, suna, Latitude);
        isr = tmpisr * cosz;
        CloudCoverRatio[ihr] = cloudcov(Radiation[ihr], isr, cosz);
        //      clcor is called to compute cloud-type correction.
        //      iamhr and ipmhr are set.
        CloudTypeCorr[ihr] = clcor(ihr, SitePar[7], isr, cosz, DayLength, Radiation[ihr], SolarNoon);
        if ((cosz >= 0.1736) && (iamhr == 0))
            iamhr = ihr;
        if ((ihr >= 12) && (cosz <= 0.1736) && (ipmhr == 0))
            ipmhr = ihr - 1;
        //      refalbed is called to compute the reference albedo for each hour.
        albedo[ihr] = refalbed(isr, Radiation[ihr], cosz, suna);
        if (jtout > 1)
        {
            // write output to file
            ofstream File18(fs::path("output") / (ProfileName + ".TM2"), ios::app);
            if (ihr == 0)
            {
                File18 << "    ***********" << endl
                       << "  date =";
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
    } //   End of 1st hourly loop
      //     Zero some variables that will later be used for summation.
    ReferenceTransp = 0;
    Rn = 0;                                 // daily net radiation
    double rnet[24];                        // hourly net radiation
    double cltcoram = CloudTypeCorr[iamhr]; //  cloud type correction during early morning
    double cltcorpm = CloudTypeCorr[ipmhr]; //  cloud type correction during late afternoon
                                            //
    for (int ihr = 0; ihr < 24; ihr++)      //  2nd hourly loop
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
            CloudCoverRatio[ihr] = 0;
            CloudTypeCorr[ihr] = 0;
        }
        //     The hourly net radiation is computed using the CIMIS algorithm (Dong et al., 1988):
        //     rlonin, the hourly incoming long wave radiation, is computed from ea0, cloud cover
        //  (CloudCoverRatio), air temperature (tk),  stefb, and cloud type correction (CloudTypeCorr).
        //     rnet, the hourly net radiation, W m-2, is computed from the global radiation, the albedo,
        //  the incoming long wave radiation, and the outgoing longwave radiation.
        double tk = AirTemp[ihr] + 273.161; // hourly air temperature in Kelvin.
        double ea0 = clearskyemiss(vp, tk); // clear sky emissivity for long wave radiation
                                            //     Compute incoming long wave radiation:
        double rlonin = (ea0 * (1 - CloudCoverRatio[ihr]) + CloudCoverRatio[ihr]) * stefb * pow(tk, 4) - CloudTypeCorr[ihr];
        rnet[ihr] = (1 - albedo[ihr]) * Radiation[ihr] + rlonin - stefb * pow(tk, 4);
        Rn += rnet[ihr];
        //     The hourly reference evapotranspiration ReferenceETP is computed by the
        //  CIMIS algorithm using the modified Penman equation:
        //     The weighting ratio (w) is computed from the functions del() (the slope of the saturation
        //  vapor pressure versus air temperature) and gam() (the psychometric constant).
        double w = del(tk, svp) / (del(tk, svp) + gam(Elevation, AirTemp[ihr])); // coefficient of the Penman equation
                                                                                 //     The wind function (fu2) is computed using different sets of parameters
                                                                                 //  for day-time and night-time. The parameter values are as suggested by CIMIS.
        double fu2;                                                              // wind function for computing evapotranspiration
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
    } //   end of 2nd hourly loop
}

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
    double xlat = Latitude * pi / 180; // latitude in radians.
    double cd;                         // amplitude of the sine of the solar height, computed as
    // the product of cosines of latitude and declination angles.
    cd = cos(xlat) * cos(declination);
    double sd; // seasonal offset of the sine of the solar height, computed
    // as the product of sines of latitude and declination angles.
    sd = sin(xlat) * sin(declination);
    double hrangle = 15 * (ti - SolarNoon) * pi / 180; // hourly angle converted to radians
    coszhr = sd + cd * cos(hrangle);
    if (coszhr <= 0)
    {
        coszhr = 0;
        sunahr = 0;
    }
    else if (coszhr >= 1)
    {
        coszhr = 1;
        sunahr = 90;
    }
    else
        sunahr = fabs(acos(coszhr) * 180 / pi - 90);
}

///////////////////////////////////////////////////////////////////////////////
double SimulateRunoff(Simulation &sim, uint32_t u, double rain)
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
    static int iGroup;         // soil group number (by clay and sand in upper soil layer)
    static double d01;         // Adjustment of curve number for soil groups A,B,C.
                               //     The following is computed only the first time the function is called.
                               //     Infiltration rate is estimated from the percent sand and percent clay in the Ap layer.
                               //  If clay content is greater than 35%, the soil is assumed to have a higher runoff potential,
                               //  if clay content is less than 15% and sand is greater than 70%, a lower runoff potential is
                               //  assumed. Other soils (loams) assumed moderate runoff potential. No 'impermeable' (group D)
                               //  soils are assumed.  References: Schwab, Brady.
    if (bFirst)
    {
        if (SandVolumeFraction[0] > 0.70 && ClayVolumeFraction[0] < 0.15)
        {
            //     Soil group A = 1, low runoff potential
            iGroup = 1;
            d01 = 1.0;
        }
        else if (ClayVolumeFraction[0] > 0.35)
        {
            //     Soil group C = 3, high runoff potential
            iGroup = 3;
            d01 = 1.14;
        }
        else
        {
            //     Soil group B = 2, moderate runoff potential
            iGroup = 2;
            d01 = 1.09;
        }
        bFirst = false;
    }
    //     Loop to accumulate 5-day antecedent rainfall (mm) which will affect the soil's ability
    //  to accept new rainfall. This also includes all irrigations.
    int i01 = u - 5;
    if (i01 < 0)
        i01 = 0;
    double PreviousWetting = 0; // five day total (before this day) of rain and irrigation, mm
    for (int Dayn = i01; Dayn < u; Dayn++)
    {
        double amtirr = 0; // mm water applied on this day by irrigation
        for (int i = 0; i < NumIrrigations; i++)
        {
            if (Dayn == Irrig[i].day)
                amtirr = Irrig[i].amount;
        }
        PreviousWetting += amtirr + sim.climate[Dayn].Rain;
    }
    //
    double d02; // Adjusting curve number for antecedent rainfall conditions.
    if (PreviousWetting < 36)
    {
        //  low moisture, low runoff potential.
        if (iGroup == 1)
            d02 = 0.71;
        else if (iGroup == 2)
            d02 = 0.78;
        else if (iGroup == 3)
            d02 = 0.83;
    }
    else if (PreviousWetting > 53)
    {
        //  wet conditions, high runoff potential.
        if (iGroup == 1)
            d02 = 1.24;
        else if (iGroup == 2)
            d02 = 1.15;
        else if (iGroup == 3)
            d02 = 1.10;
    }
    else
    {
        //  moderate conditions
        d02 = 1.00;
    }
    //  Assuming straight rows, and good cropping practice:
    double crvnum = 78.0;        // Runoff curve number, unadjusted for moisture and soil type.
    crvnum = crvnum * d01 * d02; // adjusted curve number
    double d03;                  // maximum potential difference between rainfall and runoff.
    d03 = 25400 / crvnum - 254;
    //
    double runoff; // simulated runoff on this day
    if (rain <= 0.2 * d03)
        runoff = 0; // no runoff
    else
        runoff = pow((rain - 0.2 * d03), 2) / (rain + 0.8 * d03);
    return runoff;
}