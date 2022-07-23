use crate::utils::fmin;
use crate::{
    albedo, es1hour, es2hour, AirTemp, AvrgDailyTemp, CloudCoverRatio, CloudTypeCorr, DayTimeTemp,
    Daynum, DewPointTemp, Elevation, GetFromClim, NightTimeTemp, Profile, Radiation, ReferenceETP,
    ReferenceTransp, RelativeHumidity, Rn, SitePar, WindSpeed, CLIMATE_METRIC_IRRD,
    CLIMATE_METRIC_TDEW, CLIMATE_METRIC_TMAX, CLIMATE_METRIC_TMIN, CLIMATE_METRIC_WIND,
};
use chrono::{DateTime, Datelike, Duration, FixedOffset, NaiveDate, TimeZone, Timelike};

#[derive(Debug, Clone, Copy)]
pub struct Atmosphere {
    date: NaiveDate,
    /// daily declination angle, in radians.
    pub declination: f64,
    timezone: FixedOffset,
    sunrise: DateTime<FixedOffset>,
    solar_noon: DateTime<FixedOffset>,
    sunset: DateTime<FixedOffset>,
    pub daylength: Duration,
    /// extraterrestrial radiation, W / m2.
    tmpisr: f64,
}

pub fn num_hours(duration: Duration) -> f64 {
    duration.num_seconds() as f64 / Duration::hours(1).num_seconds() as f64
}

fn hour(datetime: DateTime<FixedOffset>, tz: FixedOffset) -> f64 {
    let localtime = datetime.with_timezone(&tz);
    (localtime.hour() * 3600 + localtime.minute() * 60 + localtime.second()) as f64 / 3600.
}

fn hours(hours: f64) -> Duration {
    Duration::seconds((hours * 3600.).round() as i64)
}

impl Atmosphere {
    /// Computes day length, declination, time of solar noon, and extra-terrestrial radiation for this day.
    /// The CIMIS (California Irrigation Management Information System) algorithms are used.
    pub fn new(date: NaiveDate, longitude: f64, latitude: f64) -> Self {
        let xday = 2. * std::f64::consts::PI * date.ordinal0() as f64
            / chrono::NaiveDate::from_ymd(date.year(), 12, 31).ordinal() as f64;
        // Compute declination angle for this day. The equation used here for computing it is taken from the CIMIS algorithm.
        let declination = 0.006918 - 0.399912 * xday.cos() + 0.070257 * xday.sin()
            - 0.006758 * (2. * xday).cos()
            + 0.000907 * (2. * xday).sin()
            - 0.002697 * (3. * xday).cos()
            + 0.001480 * (3. * xday).sin();
        // Compute extraterrestrial radiation in W m-2. The 'solar constant' (average value = 1367 W m-2) is corrected for this day's distance between earth and the sun. The equation used here is from the CIMIS algorithm, which is based on the work of Iqbal (1983).
        let tmpisr = 1367.
            * (1.00011
                + 0.034221 * xday.cos()
                + 0.00128 * xday.sin()
                + 0.000719 * (2. * xday).cos()
                + 0.000077 * (2. * xday).sin());
        // Time of solar noon is computed by the CIMIS algorithm, using a correction for longitude (f), and the date correction (exday).
        //
        // It is assumed that the time zone is "geographically correct".
        //
        // For example, longitude between 22.5 and 37.5 East is in time zone GMT+2, and longitude between 112.5 and 127.5 West is in time zone GMT-8.
        //
        // All daily times in the model are computed by this method.
        let exday = (0.000075 + 0.001868 * xday.cos()
            - 0.032077 * xday.sin()
            - 0.014615 * (2. * xday).cos()
            - 0.04089 * (2. * xday).sin())
        .to_degrees();
        let timezone = FixedOffset::east(
            ((longitude + 7.5) / 15.).floor() as i32 * Duration::hours(1).num_seconds() as i32,
        );
        let solar_noon = FixedOffset::east(((longitude + exday) * 240.).round() as i32)
            .from_local_date(&date)
            .unwrap()
            .and_hms(12, 0, 0);
        // Compute day length, by commonly used equations, from latitude and declination of this day.
        // Times of sunrise and of sunset are computed from solar noon and day length.
        let xlat = latitude.to_radians();
        let ht = -xlat.tan() * declination.tan();
        let daylength = Duration::seconds(
            (2. * ((if ht > 1. {
                1.
            } else if ht < -1. {
                -1.
            } else {
                ht
            })
            .acos())
            .to_degrees()
                * 240.) as i64,
        );
        let sunrise = solar_noon - daylength / 2;
        let sunset = sunrise + daylength;
        Atmosphere {
            date,
            declination,
            timezone,
            sunrise,
            solar_noon,
            sunset,
            daylength,
            tmpisr,
        }
    }
    /// All daily weather data, read from input file, are stored in the structure:
    /// Clim[400];    --  defined in global
    /// The values are extracted from this structure by function GetFromClim()
    ///
    /// The function DayClim() is called daily from SimulateThisDay().
    /// It calls the the following functions:
    /// GetFromClim(), SimulateRunoff(),
    /// AverageAirTemperatures(), dayrad(), daytemp(), EvapoTranspiration(),
    /// tdewhour(), dayrh(), daywnd()
    ///
    /// Global variables referenced:
    ///
    /// Daynum, DayFinish, DayStart, declination, SitePar
    ///
    /// Global variables set:
    /// AirTemp, bPollinSwitch, DewPointTemp, Radiation, RelativeHumidity, WindSpeed
    pub unsafe fn meteorology(&self, profile: &Profile) {
        // latitude converted to radians.
        let xlat = profile.latitude.to_radians();
        // amplitude of the sine of the solar height.
        let cd = xlat.cos() * self.declination.cos();
        // seasonal offset of the sine of the solar height.
        let sd = xlat.sin() * self.declination.sin();
        // The computation of the daily integral of global radiation (from sunrise to sunset) is based on Spitters et al. (1986).
        // constant parameter.
        const c11: f64 = 0.4;
        // daily radiation integral.
        let radsum = if (sd / cd).abs() >= 1. {
            //  arctic circle
            0.
        } else {
            // dsbe is the integral of sinb * (1 + c11 * sinb) from sunrise to sunset,
            let dsbe = (-sd / cd).acos() * 24. / std::f64::consts::PI
                * (sd + c11 * sd * sd + 0.5 * c11 * cd * cd)
                + 12. * (cd * (2. + 3. * c11 * sd)) * (1. - (sd / cd) * (sd / cd)).sqrt()
                    / std::f64::consts::PI;
            // The daily radiation integral is computed for later use in function Radiation.
            // Daily radiation intedral is converted from langleys to Watt m-2, and divided by dsbe.
            //   11.630287 = 1000000 / 3600 / 23.884
            GetFromClim(CLIMATE_METRIC_IRRD, self.date.ordinal() as i32) * 11.630287 / dsbe
        };
        // Parameters for the daily wind function are now computed:
        // Note:  SitePar[] are site specific parameters.
        let t1 = self.sunrise + hours(SitePar[1]); // the hour at which wind begins to blow (SitePar(1) hours after sunrise).
        let t2 = self.solar_noon + hours(SitePar[2]); // the hour at which wind speed is maximum (SitePar(2) hours after solar noon).
        let t3 = self.sunset + hours(SitePar[3]); // the hour at which wind stops to blow (SitePar(3) hours after sunset).
        let wnytf = SitePar[4]; // used for estimating night time wind (from time t3 to time t1 next day).
        for ihr in 0..24 as usize {
            // time in the middle of each hourly interval.
            let ti = self
                .timezone
                .from_local_date(&self.date)
                .and_hms_opt(ihr as u32, 30, 0)
                .unwrap();
            // sine of the solar elevation.
            let sinb = sd + cd * ((num_hours(ti - self.solar_noon) * 15.).to_radians()).cos();
            // Compute hourly global radiation, using function dayrad.
            Radiation[ihr] = dayrad(radsum, sinb, c11);
            // Compute hourly temperature, using function daytmp.
            AirTemp[ihr] = self.daytmp(profile, ti);
            // Compute hourly dew point temperature, using function tdewhour.
            DewPointTemp[ihr] = self.tdewhour(profile, ti, AirTemp[ihr]);
            // Compute hourly relative humidity, using function dayrh.
            RelativeHumidity[ihr] = dayrh(AirTemp[ihr], DewPointTemp[ihr]);
            // Compute hourly wind speed, using function daywnd, and daily sum of wind.
            WindSpeed[ihr] = daywnd(
                ti,
                GetFromClim(CLIMATE_METRIC_WIND, Daynum),
                t1,
                t2,
                t3,
                wnytf,
                self.timezone,
            );
        }
        //     Compute average daily temperature, using function
        //     AverageAirTemperatures.
        AverageAirTemperatures();
        //     Compute potential evapotranspiration.
        EvapoTranspiration(profile, self.date, self);
    }

    /// Computes and returns the hourly values of air temperature, using the measured daily maximum and minimum.
    ///
    /// The algorithm is described in Ephrath et al. (1996).
    /// It is based on the following assumptions:
    /// * The time of minimum daily temperature is at sunrise.
    /// * The time of maximum daily temperature is SitePar[8] hours after solar noon.
    ///
    /// Many models assume a sinusoidal curve of the temperature during the day, but actual data deviate from the sinusoidal curve in the following characteristic way:
    /// a faster increase right after sunrise, a near plateau maximum during several hours in the middle of the day, and a rather fast decrease by sunset.
    /// The physical reason for this is a more efficient mixing of heated air from ground level into the atmospheric boundary layer, driven by strong lapse temperature gradients buoyancy.
    ///
    /// Note: ** will be used for "power" as in Fortran notation.
    ///
    /// A first order approximation is
    ///     daytmp = tmin + (tmax-tmin) * st * tkk / (tkk + daytmp - tmin)
    /// where
    ///     st = sin(pi * (ti - SolarNoon + dayl / 2) / (dayl + 2 * SitePar[8]))
    ///
    /// Since daytmp appears on both sides of the first equation, it can be solved and written explicitly as:
    ///     daytmp = tmin - tkk/2 + 0.5 * sqrt(tkk**2 + 4 * amp * tkk * st)
    /// where the amplitude of tmin and tmax is calculated as
    ///     amp = (tmax - tmin) * (1 + (tmax - tmin) / tkk)
    ///
    /// This ensures that temperature still passes through tmin and tmax values.
    ///
    /// The value of tkk was determined by calibration as 15.
    ///
    /// This algorithm is used for the period from sunrise to the time of maximum temperature, hmax.
    /// A similar algorithm is used for the time from hmax to sunset, but the value of the minimum temperature of the next day (mint_tomorrow) is used instead of mint_today.
    ///
    /// Night air temperature is described by an exponentially declining curve.
    ///
    /// For the time from sunset to mid-night:
    ///     daytmp = (mint_tomorrow - sst * exp((dayl - 24) / tcoef)
    ///              + (sst - mint_tomorrow) * exp((suns - ti) / tcoef))
    ///              / (1 - exp((dayl - 24) / tcoef))
    /// where
    ///       tcoef is a time coefficient, determined by calibration as 4
    ///       sst is the sunset temperature, determined by the daytime equation as:
    ///       sst = mint_tomorrow - tkk / 2 + 0.5 * sqrt(tkk**2 + 4 * amp * tkk * sts)
    /// where
    ///       sts  = sin(pi * dayl / (dayl + 2 * SitePar[8]))
    ///       amp = (tmax - mint_tomorrow) * (1 + (tmax - mint_tomorrow) / tkk)
    ///
    /// For the time from midnight to sunrise, similar equations are used,
    /// but the minimum temperature of this day (mint_today) is used instead of mint_tomorrow,
    /// and the maximum temperature of the previous day (maxt_yesterday) is used instead of maxt_today.
    /// Also, (suns-ti-24) is used for the time variable instead of (suns-ti).
    ///
    /// These exponential equations for night-time temperature ensure that the curve will be continuous with the daytime equation at sunset, and will pass through the minimum temperature at sunrise.
    ///
    /// Input argument:
    /// * `ti` - time of day (hours).
    /// Global variables used:
    ///
    /// Daynum, pi, SitePar, SolarNoon, sunr, suns
    unsafe fn daytmp(&self, profile: &Profile, ti: DateTime<FixedOffset>) -> f64 {
        // The temperature increase at which the sensible heat flux is doubled, in comparison with the situation without buoyancy.
        const tkk: f64 = 15.;
        // time coefficient for the exponential part of the equation.
        const tcoef: f64 = 4.;
        // hour of maximum temperature
        let hmax = self.solar_noon + hours(SitePar[8]);
        // day of year yesterday
        let im1 = self.date.ordinal0();
        // day of year tomorrow
        let ip1 = if self.date.ordinal() + 1 > profile.last_day_weather_data.ordinal() {
            self.date.ordinal()
        } else {
            self.date.ordinal() + 1
        };
        let amp; // amplitude of temperatures for a period.
        let sst; // the temperature at sunset.
        let st; // computed from time of day, used for daytime temperature.
        let sts; // intermediate variable for computing sst.
        let HourlyTemperature; // computed temperature at time ti.
        let one_day = Duration::days(1);
        //
        if ti <= self.sunrise {
            //  from midnight to sunrise
            amp = (GetFromClim(CLIMATE_METRIC_TMAX, im1 as i32)
                - GetFromClim(CLIMATE_METRIC_TMIN, self.date.ordinal() as i32))
                * (1.
                    + (GetFromClim(CLIMATE_METRIC_TMAX, im1 as i32)
                        - GetFromClim(CLIMATE_METRIC_TMIN, self.date.ordinal() as i32))
                        / tkk);
            sts = (std::f64::consts::PI * self.daylength.num_seconds() as f64
                / (self.daylength + hours(2. * SitePar[8])).num_seconds() as f64)
                .sin();
            //  compute temperature at sunset:
            sst = GetFromClim(CLIMATE_METRIC_TMIN, self.date.ordinal() as i32) - tkk / 2.
                + 0.5 * (tkk * tkk + 4. * amp * tkk * sts).sqrt();
            HourlyTemperature = (GetFromClim(CLIMATE_METRIC_TMIN, self.date.ordinal() as i32)
                - sst * (num_hours(self.daylength - one_day) / tcoef).exp()
                + (sst - GetFromClim(CLIMATE_METRIC_TMIN, self.date.ordinal() as i32))
                    * (num_hours(self.sunset - ti - one_day) / tcoef).exp())
                / (1. - (num_hours(self.daylength - one_day) / tcoef).exp());
        } else if ti <= hmax {
            //  from sunrise to hmax
            amp = (GetFromClim(CLIMATE_METRIC_TMAX, self.date.ordinal() as i32)
                - GetFromClim(CLIMATE_METRIC_TMIN, self.date.ordinal() as i32))
                * (1.
                    + (GetFromClim(CLIMATE_METRIC_TMAX, self.date.ordinal() as i32)
                        - GetFromClim(CLIMATE_METRIC_TMIN, self.date.ordinal() as i32))
                        / tkk);
            st = (std::f64::consts::PI * num_hours(ti - self.solar_noon + self.daylength / 2)
                / num_hours(self.daylength + hours(2. * SitePar[8])))
            .sin();
            HourlyTemperature = GetFromClim(CLIMATE_METRIC_TMIN, self.date.ordinal() as i32)
                - tkk / 2.
                + 0.5 * (tkk * tkk + 4. * amp * tkk * st).sqrt();
        } else if ti <= self.sunset {
            //  from hmax to sunset
            amp = (GetFromClim(CLIMATE_METRIC_TMAX, self.date.ordinal() as i32)
                - GetFromClim(CLIMATE_METRIC_TMIN, ip1 as i32))
                * (1.
                    + (GetFromClim(CLIMATE_METRIC_TMAX, self.date.ordinal() as i32)
                        - GetFromClim(CLIMATE_METRIC_TMIN, ip1 as i32))
                        / tkk);
            st = (std::f64::consts::PI * num_hours(ti - self.solar_noon + self.daylength / 2)
                / num_hours(self.daylength + hours(2. * SitePar[8])))
            .sin();
            HourlyTemperature = GetFromClim(CLIMATE_METRIC_TMIN, ip1 as i32) - tkk / 2.
                + 0.5 * (tkk * tkk + 4. * amp * tkk * st).sqrt();
        } else
        //  from sunset to midnight
        {
            amp = (GetFromClim(CLIMATE_METRIC_TMAX, self.date.ordinal() as i32)
                - GetFromClim(CLIMATE_METRIC_TMIN, ip1 as i32))
                * (1.
                    + (GetFromClim(CLIMATE_METRIC_TMAX, self.date.ordinal() as i32)
                        - GetFromClim(CLIMATE_METRIC_TMIN, ip1 as i32))
                        / tkk);
            sts = (std::f64::consts::PI * self.daylength.num_seconds() as f64
                / (self.daylength + hours(2. * SitePar[8])).num_seconds() as f64)
                .sin();
            sst = GetFromClim(CLIMATE_METRIC_TMIN, ip1 as i32) - tkk / 2.
                + 0.5 * (tkk * tkk + 4. * amp * tkk * sts).sqrt();
            HourlyTemperature = (GetFromClim(CLIMATE_METRIC_TMIN, ip1 as i32)
                - sst * (num_hours(self.daylength - one_day) / tcoef).exp()
                + (sst - GetFromClim(CLIMATE_METRIC_TMIN, ip1 as i32))
                    * (num_hours(self.sunset - ti) / tcoef).exp())
                / (1. - (num_hours(self.daylength - one_day) / tcoef).exp());
        }
        return HourlyTemperature;
        //     Reference:
        //     Ephrath, J.E., Goudriaan, J. and Marani, A. 1996. Modelling
        //  diurnal patterns of air temperature, radiation, wind speed and
        //  relative humidity by equations from daily characteristics.
        //  Agricultural Systems 51:377-393.
    }
    /// Computes the hourly values of dew point temperature from average dew-point and the daily estimated range.
    /// This range is computed as a regression on maximum and minimum temperatures.
    /// Input arguments:
    /// * `ti` - time of day (hours).
    /// * `tt` - air temperature C at this time of day.
    /// Global variables used:
    ///
    /// Daynum SitePar, SolarNoon, sunr, suns.
    unsafe fn tdewhour(&self, profile: &Profile, ti: DateTime<FixedOffset>, tt: f64) -> f64 {
        let im1 = Daynum - 1; // day of year yeaterday
        let mut ip1 = Daynum + 1; // day of year tomorrow
        if ip1 > profile.last_day_weather_data.ordinal() as i32 {
            ip1 = Daynum;
        }
        let tdewhr; // the dew point temperature (c) of this hour.
        let tdmin; // minimum of dew point temperature.
        let mut tdrange; // range of dew point temperature.
        let hmax = self.solar_noon + hours(SitePar[8]); // time of maximum air temperature

        if ti <= self.sunrise {
            // from midnight to sunrise
            tdrange = SitePar[12]
                + SitePar[13] * GetFromClim(CLIMATE_METRIC_TMAX, im1)
                + SitePar[14] * GetFromClim(CLIMATE_METRIC_TMIN, Daynum);
            if tdrange < 0. {
                tdrange = 0.;
            }
            tdmin = GetFromClim(CLIMATE_METRIC_TDEW, im1) - tdrange / 2.;
            tdewhr = tdmin
                + tdrange * (tt - GetFromClim(CLIMATE_METRIC_TMIN, Daynum))
                    / (GetFromClim(CLIMATE_METRIC_TMAX, im1)
                        - GetFromClim(CLIMATE_METRIC_TMIN, Daynum));
        } else if ti <= hmax {
            // from sunrise to hmax
            tdrange = SitePar[12]
                + SitePar[13] * GetFromClim(CLIMATE_METRIC_TMAX, Daynum)
                + SitePar[14] * GetFromClim(CLIMATE_METRIC_TMIN, Daynum);
            if tdrange < 0. {
                tdrange = 0.;
            }
            tdmin = GetFromClim(CLIMATE_METRIC_TDEW, Daynum) - tdrange / 2.;
            tdewhr = tdmin
                + tdrange * (tt - GetFromClim(CLIMATE_METRIC_TMIN, Daynum))
                    / (GetFromClim(CLIMATE_METRIC_TMAX, Daynum)
                        - GetFromClim(CLIMATE_METRIC_TMIN, Daynum));
        } else {
            // from hmax to midnight
            tdrange = SitePar[12]
                + SitePar[13] * GetFromClim(CLIMATE_METRIC_TMAX, Daynum)
                + SitePar[14] * GetFromClim(CLIMATE_METRIC_TMIN, ip1);
            if tdrange < 0. {
                tdrange = 0.;
            }
            tdmin = GetFromClim(CLIMATE_METRIC_TDEW, ip1) - tdrange / 2.;
            tdewhr = tdmin
                + tdrange * (tt - GetFromClim(CLIMATE_METRIC_TMIN, ip1))
                    / (GetFromClim(CLIMATE_METRIC_TMAX, Daynum)
                        - GetFromClim(CLIMATE_METRIC_TMIN, ip1));
        }
        return tdewhr;
    }
}
/// Computes the hourly values of global radiation, in W m-2, using the measured daily total global radiation.
/// The algorithm follows the paper of Spitters et al. (1986).
/// It assumes that atmospheric transmission of radiation is lower near the margins of the daylight period,
/// because of an increase in the path length through the atmosphere at lower solar heights.
/// Radiation is therefore assumed to be proportional to sinb * (1 + c11 * sinb), where the value of c11 is set as 0.4 .
///
/// Input arguments:
/// radsum - daily radiation integral.
/// sinb - sine of the solar elevation.
/// c11 - constant parameter (0.4).
fn dayrad(radsum: f64, sinb: f64, c11: f64) -> f64 {
    let HourlyRadiation = radsum * sinb * (1. + c11 * sinb);
    if HourlyRadiation < 0. {
        0.
    } else {
        HourlyRadiation
    }
    //      References:
    //      Spitters, C.J.T., Toussaint, H.A.J.M. and Goudriaan, J. 1986.
    // Separating the diffuse and direct component of global radiation and
    // its implications for modeling canopy photosynthesis. Part I.
    // Components of incoming radiation. Agric. For. Meteorol. 38:217-229.
    //     Ephrath, J.E., Goudriaan, J. and Marani, A. 1996. Modelling
    //  diurnal patterns of air temperature, radiation, wind speed and
    //  relative humidity by equations from daily characteristics.
    //  Agricultural Systems 51:377-393.
}

/// Computes the hourly values of relative humidity, using the hourly air and dew point temperatures.
/// It calls function [vapor_pressure()].
/// If the estimated dew point is higher than the actual air temperature, its value is taken as the air temperature (relative humidity 100%).
///
/// The relative humidity is calculated as the percentage ratio of the saturated vapor pressure at dew point temperature and the saturated vapor pressure at actual air temperature.
///
/// Input arguments:
/// tt - air temperature C at this time of day.
/// tdew - dew point temperature C at this time of day.
fn dayrh(tt: f64, tdew: f64) -> f64 {
    let td = fmin(tt, tdew); // the dew point temperature (C), is assumed to
                             // be tt if tt < tdew.
    let esvp = vapor_pressure(tt); // the saturated vapor pressure in the air (mbar).
    let vpa = vapor_pressure(td); // the actual vapor pressure in the air (mbar).
    let mut RelHumHour = 100. * vpa / esvp; // relative humidity at this time of day, %.
    if RelHumHour < 1. {
        RelHumHour = 1.;
    }
    if RelHumHour > 100. {
        RelHumHour = 100.;
    }
    return RelHumHour;
    //     Reference:
    //     Ephrath, J.E., Goudriaan, J. and Marani, A. 1996. Modelling
    //  diurnal patterns of air temperature, radiation, wind speed and
    //  relative humidity by equations from daily characteristics.
    //  Agricultural Systems 51:377-393.
}

/// Computes the hourly values of wind speed (m / sec), estimated from the measured total daily wind run.
///
/// Input arguments:
///
/// t1 = the hour at which day-time wind begins to blow.
/// t2 = the hour at which day-time wind speed is maximum.
/// t3 = the hour at which day-time wind ceases to blow.
/// ti = the hour of the day.
/// wind = the daily total wind run (km per day).
/// wnytf = Factor for estimating night-time wind (from time t3 to time t1 next day).
///
/// The algorithm is described by Ephrath et al. (1996). It is based on the following assumptions:
///
/// Although the variability of wind speed during any day is very large, the diurnal wind speed curves appear to be characterized by
/// the following repetitive pattern: increase in wind speed from time
/// t1 in the morning to time t2 in the afternoon, decrease from t2 to t3
/// in the evening, and a low constant wind speed at night, from t3 to t1
/// in the next day.
///
/// The values of t1, t2, and t3 have been determined in the calling routine:
/// t1 is SitePar(1) hours after sunrise, t2 is SitePar(2) hours after solar
/// noon, and t3 is SitePar(3) hours after sunset. These parameters are site-
/// specific. They are 1, 3, and 0, respectively, for the San Joaquin valley of
/// California and for Arizona, and 1, 4, and 2, respectively, for the coastal
/// plain of israel.
///
/// The wind speed during the night, from t3 to t1 next day (wmin)
/// is assumed to be proportional to the daily total wind run. The
/// ratio wnytf is also site-specific, SitePar(4), ( 0.008 for San Joaquin
/// and Arizona, 0.0025 for the coastal plain of Israel). wmin is the minimum
/// wind speed from t1 to t3.
///
/// wtday is computed by subtracting the daily integral of wmin, after
/// converting it from m/sec to km/day, from the total daily wind run (wndt).
///
/// wmax, the maximum wind speed at time t2 (minus wmin), is
/// computed from wtday and converted to m/sec.

/// daywnd from t1 to t2 is now computed as an increasing
/// sinusoidal function from wmin to wmin + wmax, and it is computed from
/// t2 to t3 as a decreasing sinusoidal function from wmin + wmax to  wmin.
fn daywnd(
    ti: DateTime<FixedOffset>,
    wind: f64,
    t1: DateTime<FixedOffset>,
    t2: DateTime<FixedOffset>,
    t3: DateTime<FixedOffset>,
    wnytf: f64,
    timezone: FixedOffset,
) -> f64 {
    let HourlyWind;
    // constants related to t1, t2, t3 :
    let sf1 = (t2 - t1) * 4;
    let sf2 = (t3 - t2) * 4;
    // the constant minimum wind speed during the night (m/sec).
    let wmin = wind * wnytf;
    // integral of wind run from t1 to t3, minus wmin (km).
    let wtday = wind - wmin * 3.6 * 24.;
    // the maximum wind speed (m per sec), above wmin.
    let wmax = wtday * 2. * std::f64::consts::PI / 3.6 / num_hours(sf1 + sf2);

    if ti >= t1 && ti < t2 {
        HourlyWind =
            wmin + wmax * (2. * std::f64::consts::PI * num_hours(ti - t1) / num_hours(sf1)).sin();
    } else if ti >= t2 && ti < t3 {
        HourlyWind = wmin
            + wmax
                * (2. * std::f64::consts::PI * hour(ti - (t2 - t3) * 2, timezone) / num_hours(sf2))
                    .sin();
    } else {
        HourlyWind = wmin;
    }
    return HourlyWind;
    // Reference:
    //     Ephrath, J.E., Goudriaan, J. and Marani, A. 1996. Modelling
    //  diurnal patterns of air temperature, radiation, wind speed and
    //  relative humidity by equations from daily characteristics.
    //  Agricultural Systems 51:377-393.
}

/// Calculates daily average temperatures, daytime average and night time average.
unsafe fn AverageAirTemperatures() {
    let mut nn1 = 0; // counter of night hours
    let mut nn2 = 0; // counter of daytime hours
    AvrgDailyTemp = 0.;
    NightTimeTemp = 0.;
    DayTimeTemp = 0.;

    for ihr in 0..24 {
        AvrgDailyTemp += AirTemp[ihr];
        if Radiation[ihr] <= 0. {
            NightTimeTemp += AirTemp[ihr];
            nn1 += 1;
        } else {
            DayTimeTemp += AirTemp[ihr];
            nn2 += 1;
        }
    }
    AvrgDailyTemp = AvrgDailyTemp / 24.;
    NightTimeTemp /= nn1 as f64;
    DayTimeTemp /= nn2 as f64;
}
/// Computes the water vapor pressure in the air (in KPa units) as a function of the air at temperature tt (C). This equation is widely used.
pub fn vapor_pressure(tt: f64) -> f64 {
    0.61078 * (17.269 * tt / (tt + 237.3)).exp()
}
/// Computes the rate of reference evapotranspiration and related variables.
/// The following subroutines and functions are called for each hour: sunangle, cloudcov(), clcor(), refalbed(), [vapor_pressure()],
/// [clearskyemiss()], [del()], [gam()].
unsafe fn EvapoTranspiration(profile: &Profile, date: NaiveDate, atmosphere: &Atmosphere) {
    const stefb: f64 = 5.77944E-08; // the Stefan-Boltzman constant, in W m-2 K-4 (= 1.38E-12 * 41880)
    const c12: f64 = 0.125; // c12 ... c15 are constant parameters.
    const c13: f64 = 0.0439;
    const c14: f64 = 0.030;
    const c15: f64 = 0.0576;
    let mut iamhr: usize = 0; // earliest time in day for computing cloud cover
    let mut ipmhr: usize = 0; // latest time in day for computing cloud cover
    let mut isr: f64; // hourly extraterrestrial radiation in W / m**2

    for ihr in 0..24 as usize {
        // The following subroutines and functions are called for each hour: sunangle, cloudcov, clcor, refalbed .
        // cosine of sun angle from zenith for this hour
        let cosz = sunangle(
            profile.latitude,
            atmosphere
                .timezone
                .from_local_date(&date)
                .and_hms_opt(ihr as u32, 30, 0)
                .unwrap(),
            atmosphere.declination,
            atmosphere.solar_noon,
        );
        // sun angle from horizon, degrees at this hour
        let suna = (cosz.acos().to_degrees() - 90.).abs();
        isr = atmosphere.tmpisr * cosz;
        CloudCoverRatio[ihr] = cloudcov(Radiation[ihr], isr, cosz);
        //      clcor is called to compute cloud-type correction.
        //      iamhr and ipmhr are set.
        CloudTypeCorr[ihr] = clcor(
            ihr,
            SitePar[7],
            isr,
            cosz,
            atmosphere.solar_noon,
            atmosphere.daylength,
        );
        if (cosz >= 0.1736) && (iamhr == 0) {
            iamhr = ihr;
        }
        if ihr >= 12 && cosz <= 0.1736 && ipmhr == 0 {
            ipmhr = ihr - 1;
        }
        //      refalbed is called to compute the reference albedo for each
        //      hour.
        albedo[ihr] = refalbed(isr, Radiation[ihr], cosz, suna);
    }
    // Zero some variables that will later be used for summation.
    ReferenceTransp = 0.;
    Rn = 0.; // daily net radiation
    let mut rnet: [f64; 24] = [0.; 24]; // hourly net radiation
    for ihr in 0..24 as usize {
        // Compute saturated vapor pressure (svp), using function vapor_pressure().
        // The actual vapor pressure (vp) is computed from svp and the relative humidity.
        // Compute vapor pressure deficit (vpd). This procedure is based on the CIMIS algorithm.
        //
        // saturated vapor pressure, mb
        let svp = vapor_pressure(AirTemp[ihr]);
        // vapor pressure, mb
        let vp = 0.01 * RelativeHumidity[ihr] * svp;
        // vapor pressure deficit, mb.
        let vpd = svp - vp;
        // Get cloud cover and cloud correction for night hours
        if ihr < iamhr || ihr > ipmhr {
            CloudCoverRatio[ihr] = 0.;
            CloudTypeCorr[ihr] = 0.;
        }
        // The hourly net radiation is computed using the CIMIS algorithm (Dong et al., 1988): rlonin, the hourly incoming long wave
        // radiation, is computed from ea0, cloud cover
        // (CloudCoverRatio), air temperature (tk),  stefb, and cloud type
        // correction (CloudTypeCorr).
        //
        // rnet, the hourly net radiation, W m-2, is computed from the  global radiation, the albedo, the incoming long wave radiation, and the outgoing longwave radiation.
        // hourly air temperature in Kelvin.
        let tk = AirTemp[ihr] + 273.161;
        // clear sky emissivity for long wave radiation
        let ea0 = clearskyemiss(vp, tk);
        // Compute incoming long wave radiation:
        let rlonin =
            (ea0 * (1. - CloudCoverRatio[ihr]) + CloudCoverRatio[ihr]) * stefb * tk.powi(4)
                - CloudTypeCorr[ihr];
        rnet[ihr] = (1. - albedo[ihr]) * Radiation[ihr] + rlonin - stefb * tk.powi(4);
        Rn += rnet[ihr];
        // The hourly reference evapotranspiration ReferenceETP is computed by the CIMIS algorithm using the modified Penman equation:
        // The weighting ratio (w) is computed from the functions del() (the slope of the saturation vapor pressure versus air temperature) and gam() (the psychometric constant).
        //
        // coefficient of the Penman equation
        let w = del(tk, svp) / (del(tk, svp) + gam(Elevation, AirTemp[ihr]));
        // The wind function (fu2) is computed using different sets of parameters for day-time and night-time.
        // The parameter values are as suggested by CIMIS.
        //
        // wind function for computing evapotranspiration
        let fu2 = if Radiation[ihr] <= 0. {
            c12 + c13 * WindSpeed[ihr]
        } else {
            c14 + c15 * WindSpeed[ihr]
        };
        // hlathr, the latent heat for evaporation of water (W m-2 per mm at this hour) is computed as a function of temperature.
        let hlathr = 878.61 - 0.66915 * (AirTemp[ihr] + 273.161);
        // ReferenceETP, the hourly reference evapotranspiration, is now computed by the modified Penman equation.
        ReferenceETP[ihr] = w * rnet[ihr] / hlathr + (1. - w) * vpd * fu2;
        if ReferenceETP[ihr] < 0. {
            ReferenceETP[ihr] = 0.;
        }
        // ReferenceTransp is the sum of ReferenceETP
        ReferenceTransp += ReferenceETP[ihr];
        // es1hour and es2hour are computed as the hourly potential evapotranspiration due to radiative and aerodynamic factors, respectively.
        // es1hour and ReferenceTransp are not computed for periods of negative net radiation.
        es2hour[ihr] = (1. - w) * vpd * fu2;
        if rnet[ihr] > 0. {
            es1hour[ihr] = w * rnet[ihr] / hlathr;
        } else {
            es1hour[ihr] = 0.;
        }
    }
}

/// Function clearskyemiss() estimates clear sky emissivity for long wave radiation.
///
/// Input arguments:
/// * `vp` - vapor pressure of the air in KPa
/// * `tk` - air temperature in K.
///
/// References:
/// * https://doi.org/10.1029/WR017i002p00295
pub fn clearskyemiss(vp: f64, tk: f64) -> f64 {
    // Convert vp to mbars
    // vapor pressure of the air in mbars.
    let vp1 = vp * 10.;
    // Compute clear sky emissivity by the method of Idso (1981)
    let ea0 = 0.70 + 5.95e-05 * vp1 * (1500. / tk).exp();
    if ea0 > 1. {
        1.
    } else {
        ea0
    }
}

/// Computes cloud cover for this hour from radiation data, using the CIMIS algorithm.
/// The return value is cloud cover ratio ( 0 to 1 )
///
/// Input arguments:
/// * `radihr` - hourly global radiation in W m<sup>-2</sup>.
/// * `isr` - hourly extraterrestrial radiation in W m<sup>-2</sup>.
/// * `cosz` - cosine of sun angle from zenith.
///
/// This algorithm is described by Dong et al. (1988).
/// Cloud cover fraction is estimated as a function of the ratio of actual solar radiation to extraterrestrial radiation.
/// The parameters of this function have been based on California data.
///
/// The equation is for daylight hours, when the sun is not less than 10 degrees above the horizon (coszhr > 0.1736).
fn cloudcov(radihr: f64, isr: f64, cosz: f64) -> f64 {
    const p1: f64 = 1.333; //    p1, p2, p3 are constant parameters.
    const p2: f64 = 1.7778;
    const p3: f64 = 0.294118;
    // ratio of radihr to isr.
    let rasi = if isr > 0. { radihr / isr } else { 0. };
    if cosz > 0.1736 && rasi <= p1 / p2 {
        // computed cloud cover.
        let clcov = (p1 - p2 * if rasi >= 0.375 { rasi } else { 0.375 }).powf(p3);
        if clcov < 0. {
            0.
        } else {
            clcov
        }
    } else {
        0.
    }
    // Reference:
    // Dong, A., Prashar, C.K. and Grattan, S.R. 1988. Estimation of daily and hourly net radiation.
    // CIMIS Final Report June 1988, pp.58-79.
}

/// Computes cloud type correction, using the CIMIS algorithm.
///
/// Input arguments:
/// * `ck` - cloud type correction factor (data for this location).
/// * `coszhr` - cosine of sun angle from zenith.
/// * `ihr` - time of day, hours.
/// * `isrhr` - hourly extraterrestrial radiation in W m-2 .
/// Global variables used: Radiation[], SolarNoon, pi
///
/// NOTE: This algorithm is described by Dong et al. (1988).
/// ck is the cloud-type correction used in the Monteith equation for estimating net radiation.
/// The value of this correction depends on site and time of year.
/// Regional ck values for California are given by Dong et al. (1988).
/// In the San Joaquin valley of California ck is almost constant from April to October, with an average value of 60.
/// The value of ck is site-dependant, assumed to be constant during the growing season.
///
/// The daily ck is converted to an hourly value for clear or partly cloudy sky (rasi >= 0.375) and when the sun is at least 10 degrees above the horizon.
///
/// Evening, night and early morning cloud type correction is temporarily assigned 0. It is later assigned the values of first or last non-zero values (in the calling routine).
unsafe fn clcor(
    ihr: usize,
    ck: f64,
    isrhr: f64,
    coszhr: f64,
    solar_noon: DateTime<FixedOffset>,
    daylength: Duration,
) -> f64 {
    let mut rasi = 0.; //  ratio of Radiation to isrhr.
    if isrhr > 0. {
        rasi = Radiation[ihr] / isrhr;
    }
    if coszhr >= 0.1736 && rasi >= 0.375 {
        let angle = std::f64::consts::PI
            * (solar_noon.with_hour(ihr as u32).unwrap() - solar_noon + Duration::minutes(30))
                .num_seconds() as f64
            / daylength.num_seconds() as f64; // hour angle (from solar noon) in radians.
        ck * std::f64::consts::PI / 2. * angle.cos()
    } else {
        0.
    }
    // Reference:
    // Dong, A., Prashar, C.K. and Grattan, S.R. 1988. Estimation of daily and hourly net radiation.
    // CIMIS Final Report June 1988, pp.58-79.
}

/// Computes the slope of the saturation vapor pressure (svp, in mb) versus air temperature (tk, in K).
/// This algorithm is the same as used by CIMIS.
fn del(tk: f64, svp: f64) -> f64 {
    let a = 10f64.powf(-0.0304 * tk);
    let b = tk.powi(2);
    let c = 10f64.powf(-1302.88 / tk);
    return (6790.5 - 5.02808 * tk + 4916.8 * a * b + 174209. * c) * svp / b;
}

/// Function gam() computes the psychometric constant at elevation (elev), m above sea level, and air temperature, C (tt).
/// This algorithm is the same as used by CIMIS.
fn gam(elev: f64, tt: f64) -> f64 {
    let bp = 101.3 - 0.01152 * elev + 5.44e-07 * elev.powi(2); //  barometric pressure, KPa, at this elevation.
    return 0.000646 * bp * (1. + 0.000946 * tt);
}
/// Computes the reference crop albedo, using the CIMIS algorithm.
///
/// This algorithm is described by Dong et al. (1988). Albedo is estimated as a function of sun elevation above the horizon (suna) for clear or partly cloudy sky (rasi >= 0.375) and when the sun is at least 10 degrees above the horizon.
///
/// For very cloudy sky, or when solar altitude is below 10 degrees, the following albedo value is assumed: (p4)+ 0.26
///
/// Input arguments:
/// * `isrhr` - hourly extraterrestrial radiation in W m-2 .
/// * `rad` - hourly global radiation in W / m-2 .
/// * `coszhr` - cosine of sun angle from zenith.
/// * `sunahr` - sun angle from horizon, degrees.
fn refalbed(isrhr: f64, rad: f64, coszhr: f64, sunahr: f64) -> f64 {
    const p1: f64 = 0.00158; //  p1 ... p4 are constant parameters.
    const p2: f64 = 0.386;
    const p3: f64 = 0.0188;
    const p4: f64 = 0.26;
    let mut rasi = 0.; //   ratio of rad to isrhr

    if isrhr > 0. {
        rasi = rad / isrhr;
    }
    if coszhr > 0.1736 && rasi >= 0.375 {
        let refalb = p1 * sunahr + p2 * (-p3 * sunahr).exp();
        if refalb > p4 {
            p4
        } else {
            refalb
        }
    } else {
        p4
    }
    // Reference:
    // Dong, A., Prashar, C.K. and Grattan, S.R. 1988. Estimation of daily and hourly net radiation.
    // CIMIS Final Report June 1988, pp. 58-79.
}

/// Computes sun angle for any time of day.
///
/// Input argument:
/// * `latitude` - latitude in degree.
/// * `ti` - time of day, hours.
///
/// Output arguments:
///
/// * `coszhr` - cosine of sun angle from zenith for this hour.
/// * `sunahr` - sun angle from horizon, degrees.
fn sunangle(
    latitude: f64,
    ti: DateTime<FixedOffset>,
    declination: f64,
    solar_noon: DateTime<FixedOffset>,
) -> f64 {
    // The latitude is converted to radians (xlat).
    let xlat = latitude.to_radians();
    // amplitude of the sine of the solar height, computed as the product of cosines of latitude and declination angles.
    let cd = xlat.cos() * declination.cos();
    // seasonal offset of the sine of the solar height, computed as the product of sines of latitude and declination angles.
    let sd = xlat.sin() * declination.sin();
    // hourly angle converted to radians
    let hrangle = 15. * num_hours(ti - solar_noon).to_radians();
    let result = sd + cd * hrangle.cos();
    if result <= 0. {
        0.
    } else if result >= 1. {
        1.
    } else {
        result
    }
}
