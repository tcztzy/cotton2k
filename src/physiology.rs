/// computes physiological age
///
/// This function returns the daily 'physiological age' increment,
/// based on hourly temperatures. It is called each day by `SimulateThisDay`.
#[no_mangle]
extern "C" fn PhysiologicalAge(air_temp: *const f64) -> f64 {
    let air_temperature = unsafe { std::slice::from_raw_parts(air_temp, 24) };
    //     The threshold value is assumed to be 12 C (p1). One physiological day is
    //  equivalent to a day with an average temperature of 26 C, and therefore the
    //  heat units are divided by 14 (p2).
    //     A linear relationship is assumed between temperature and heat unit
    //  accumulation in the range of 12 C (p1) to 33 C (p2*p3+p1). the effect of
    //  temperatures higher than 33 C is assumed to be equivalent to that of 33 C.
    //     The following constant Parameters are used in this function:
    let p1 = 12.; // threshold temperature, C
    let p2 = 14.; // temperature, C, above p1, for one physiological day.
    let p3 = 1.5; // maximum value of a physiological day.

    let mut dayfd = 0.; // the daily contribution to physiological age (return value).
    for t in air_temperature {
        let mut tfd = (t - p1) / p2; // the hourly contribution to physiological age.
        if tfd < 0. {
            tfd = 0.;
        }
        if tfd > p3 {
            tfd = p3;
        }
        dayfd += tfd;
    }
    dayfd / 24.
}

/// TODO: effects of pix applied.
#[no_mangle]
extern "C" fn Pix() {}
