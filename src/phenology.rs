static mut avtemp: f64 = 0f64; // average temperature from day of emergence.
static mut sumstrs: f64 = 0f64; // cumulative effect of water and N stresses on date of first square.
#[no_mangle]
extern "C" fn DaysToFirstSquare(
    Daynum: i32,
    DayEmerge: i32,
    AvrgDailyTemp: f64,
    WaterStress: f64,
    NStressVeg: f64,
    var30: f64,
) -> f64
//     This function computes and returns tsq1, the number of days from emergence to first
//  square. It is called from CottonPhenology().
//     The following global variables are referenced here:
//        AvrgDailyTemp, Kday, NStressVeg, VarPar, WaterStress.
//
{
    //     The following constant parameters are used:
    const p1: f64 = 34.;
    const p2: f64 = 132.2;
    const p3: f64 = -7.;
    const p4: f64 = 0.125;
    const p5: f64 = 0.08;
    const p6: f64 = 0.30;
    //     On day of emergence assign initial values to some local variables.
    if Daynum <= DayEmerge {
        unsafe {
            avtemp = AvrgDailyTemp;
            sumstrs = 0f64;
        }
    }
    //      Compute the average temperature from the day of emergence to
    //  now. Compute the number of days from emergence to first
    //  square, as a function of this average temperature. This
    //  relationships is derived from data of K. R. Reddy et al.
    //  (unpublished), CSRU, for Delta cultivars.
    unsafe {
        avtemp = ((Daynum - DayEmerge) as f64 * avtemp + AvrgDailyTemp)
            / (Daynum - DayEmerge + 1) as f64;
        if avtemp > p1 {
            avtemp = p1;
        }
        //      The cumulative effect of water stress and vegetative N stress
        //  is computed. This effect is substacted from TSQ, assuming that the
        //  first square appears earlier under water or N stressed conditions.
        sumstrs += p5 * (1f64 - WaterStress) + p6 * (1f64 - NStressVeg);
        (p2 + avtemp * (p3 + avtemp * p4)) * var30 - sumstrs
    }
}
