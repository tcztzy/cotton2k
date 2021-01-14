/// The function `PsiOnTranspiration` computes and returns the effect of the average soil
/// matrix water potential on transpiration rate. It is called by `WaterUptake`.
/// The argument PsiAverage is the average soil water matrix potential, bars.
#[no_mangle]
extern "C" fn PsiOnTranspiration(psi_average: f64) -> f64 {
    // This is a third degree function with two parameters (a, b). It has the
    // value of 1 when PsiAverage = b - a, and the value of 0 when PsiAverage = - a.

    // The minimum value, however, is set to d, and the maximum value to c.
    let a = 20.;
    let b = 14.;
    let c = 1.00;
    let d = 0.05;
    let rfep = ((a + psi_average) / b).powi(3);
    if rfep > c {
        c
    } else if rfep < d {
        d
    } else {
        rfep
    }
}

/// This function computes soil water hydraulic conductivity for a given value of soil water content, using the Van-Genuchten equation. The units of the computed conductivity are the same as the given saturated conductivity (`SaturatedHydCond`).
#[no_mangle]
extern "C" fn wcond(
    q: f64,
    qr: f64,
    qsat: f64,
    beta: f64,
    saturated_hyd_cond: f64,
    pore_space: f64,
) -> f64
// The following arguments are used:
//   beta  - parameter of the van-genuchten equation.
//   saturated_hyd_cond - saturated hydraulic conductivity (at qsat).
//   pore_space - pore space volume.
//   q - soil water content, cm3 cm-3.
//   qr - residual water content, cm3 cm-3.
//   qsat - saturated water content, cm3 cm-3.
{
    // For very low values of water content (near the residual water content) wcond is 0.
    if (q - qr) < 0.0001 {
        return 0.;
    }
    // Water content for saturated conductivity is minimum of PoreSpace and qsat.

    // For very high values of water content (exceeding the saturated
    // water content or pore space) conductivity is SaturatedHydCond.
    let xsat = if qsat < pore_space { qsat } else { pore_space };
    if q >= xsat {
        return saturated_hyd_cond;
    }
    // The following equation is used (in FORTRAN notation):
    //   WCOND = CONDSAT * ((Q-QR)/(XSAT-QR))**0.5
    //           * (1-(1-((Q-QR)/(XSAT-QR))**(1/GAMA))**GAMA)**2
    let gama = 1. - 1. / beta;
    let gaminv = 1. / gama;
    let sweff = (q - qr) / (xsat - qr); // intermediate variable (effective water content).
    let acoeff = (1. - sweff.powf(gaminv)).powf(gama); // intermediate variable
    let bcoeff = (1. - acoeff).powi(2); // intermediate variable
    sweff.powf(0.5) * bcoeff * saturated_hyd_cond
}

/// This function computes soil water content (cm3 cm-3) for
/// a given value of matrix potential, using the Van-Genuchten equation.
#[no_mangle]
extern "C" fn qpsi(psi: f64, qr: f64, qsat: f64, alpha: f64, beta: f64) -> f64
// The following arguments are used:
//   alpha, beta  - parameters of the van-genuchten equation.
//   psi - soil water matrix potential (bars).
//   qr - residual water content, cm3 cm-3.
//   qsat - saturated water content, cm3 cm-3.
{
    // For very high values of PSI, saturated water content is assumed.
    // For very low values of PSI, air-dry water content is assumed.
    if psi >= -0.00001 {
        qsat
    } else if psi <= -500000f64 {
        qr
    } else {
        // The soil water matric potential is transformed from bars (psi)
        // to cm in positive value (psix).
        let psix = 1000. * (psi + 0.00001).abs();
        // The following equation is used (in FORTRAN notation):
        //   QPSI = QR + (QSAT-QR) / (1 + (ALPHA*PSIX)**BETA)**(1-1/BETA)
        let gama = 1. - 1. / beta;
        let term = 1. + (alpha * psix).powf(beta); //  intermediate variable
        let swfun = qr + (qsat - qr) / term.powf(gama); //  computed water content
        if swfun < (qr + 0.0001) {
            qr + 0.0001
        } else {
            swfun
        }
    }
}

/// This function computes soil water matric potential (in bars) for a given value of soil water content, using the Van-Genuchten equation.
#[no_mangle]
extern "C" fn psiq(q: f64, qr: f64, qsat: f64, alpha: f64, beta: f64) -> f64
// The following arguments are used:
//   alpha, beta  - parameters of the van-genuchten equation.
//   q - soil water content, cm3 cm-3.
//   qr - residual water content, cm3 cm-3.
//   qsat - saturated water content, cm3 cm-3.
{
    // For very low values of water content (near the residual water
    // content) psiq is -500000 bars, and for saturated or higher water
    // content psiq is -0.00001 bars.
    if (q - qr) < 0.00001 {
        return -500000.;
    } else if q >= qsat {
        return -0.00001;
    }
    // The following equation is used (FORTRAN notation):
    // PSIX = (((QSAT-QR) / (Q-QR))**(1/GAMA) - 1) **(1/BETA) / ALPHA
    let gama = 1. - 1. / beta;
    let gaminv = 1. / gama;
    let term = ((qsat - qr) / (q - qr)).powf(gaminv); //  intermediate variable
    let mut psix = (term - 1.).powf(1. / beta) / alpha;
    if psix < 0.01 {
        psix = 0.01;
    }
    // psix (in cm) is converted to bars (negative value).
    psix = (0.01 - psix) * 0.001;
    if psix < -500000. {
        psix = -500000.;
    }
    if psix > -0.00001 {
        psix = -0.00001;
    }
    psix
}

/// This function computes soil water osmotic potential (in bars, positive value).
#[no_mangle]
extern "C" fn PsiOsmotic(q: f64, qsat: f64, ec: f64) -> f64
// The following arguments are used:
//   q - soil water content, cm3 cm-3.
//   qsat - saturated water content, cm3 cm-3.
//   ec - electrical conductivity of saturated extract (mmho/cm)
{
    if ec > 0f64 {
        let result = 0.36 * ec * qsat / q;
        if result > 6f64 {
            6f64
        } else {
            result
        }
    } else {
        0f64
    }
}

/// This function computes the effect of temperature on the rate of mineralization of organic mineralizable nitrogen. It is based on GODWIN and JONES (1991).
#[no_mangle]
extern "C" fn SoilTemperatureEffect(tt: f64) -> f64
// The following argument is used:  tt - soil temperature (C).
{
    // The following constant parameters are used:
    let tfpar1 = 0.010645;
    let tfpar2 = 0.12979;
    // The temperature function of CERES is replaced by the function
    // suggested by Vigil and Kissel (1995):
    //   tfm = 0.010645 * exp(0.12979 * tt)
    // Note: tfm = 0.5 for 29.66 C, tfm = 1 for 35 C, tfm = 2 for 40.34 C.
    let tfm = tfpar1 * std::f64::consts::E.powf(tfpar2 * tt);
    if tfm < 0f64 {
        0f64
    } else if tfm > 2f64 {
        2f64
    } else {
        tfm
    }
}

#[no_mangle]
extern "C" fn SoilWaterEffect(
    volumetric_water_content: f64,
    field_capacity: f64,
    volumetric_water_content_at_permanent_wilting_point: f64,
    volumetric_water_content_saturated: f64,
    xx: f64,
) -> f64
//     This function computes the effect of soil moisture on the rate of mineralization of 
//  organic mineralizable nitrogen, and on the rates of urea hydrolysis and nitrification.
//     It is based on Godwin and Jones (1991).
//     The following global variables are referenced:
//       FieldCapacity, thetar, thts, VolWaterContent.
//     The argument xx is 0.5 when used for mineralization and urea hydrolysis,
//  or 1.0 when used for nitrification.
//     l, k are layer and column of this cell.
//
{
    let wf = // the effect of soil moisture on process rate.
    if volumetric_water_content <= field_capacity {
        // Soil water content less than field capacity:
        (volumetric_water_content - volumetric_water_content_at_permanent_wilting_point) / (field_capacity - volumetric_water_content_at_permanent_wilting_point)}
    else{
        // Soil water content more than field capacity:
        1f64 - xx * (volumetric_water_content - field_capacity) / (volumetric_water_content_saturated - field_capacity)
    };

    if wf < 0f64 {
        0f64
    } else {
        wf
    }
}

/// This function calculates soil mechanical resistance of cell l,k. It is computed
/// on the basis of parameters read from the input and calculated in RootImpedance().
///
/// It is called from PotentialRootGrowth().
///
///  The function has been adapted, without change, from the code of GOSSYM. Soil mechanical
/// resistance is computed as an empirical function of bulk density and water content.
/// It should be noted, however, that this empirical function is based on data for one type
/// of soil only, and its applicability for other soil types is questionable. The effect of soil
/// moisture is only indirectly reflected in this function. A new module (SoilWaterOnRootGrowth)
/// has therefore been added in COTTON2K to simulate an additional direct effect of soil
/// moisture on root growth.
///
/// The minimum value of rtimpd of this and neighboring soil cells is used to compute
/// rtpct. The code is based on a segment of RUTGRO in GOSSYM, and the values of the p1 to p3
/// parameters are based on GOSSYM usage:
#[no_mangle]
extern "C" fn SoilMechanicResistance(rtimpdmin: f64) -> f64 {
    let p1 = 1.046;
    let p2 = 0.034554;
    let p3 = 0.5;

    // effect of soil mechanical resistance on root growth (the return value).
    let rtpct = p1 - p2 * rtimpdmin;
    if rtpct > 1f64 {
        1f64
    } else if rtpct < p3 {
        p3
    } else {
        rtpct
    }
}

/// This function computes the aggregation factor for 2 mixed soil materials.
#[no_mangle]
extern "C" fn form(c0: f64, d0: f64, g0: f64) -> f64
// Arguments referenced:
//   c0 - heat conductivity of first material
//   d0 - heat conductivity of second material	
//   g0 - shape factor for these materials
{
    (2f64 / (1f64 + (c0 / d0 - 1f64) * g0) + 1f64 / (1f64 + (c0 / d0 - 1f64) * (1f64 - 2f64 * g0)))
        / 3f64
}
