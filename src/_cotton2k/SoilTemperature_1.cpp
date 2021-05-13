// File SoilTemperature_1.cpp
//   List of functions in this file:
//       ColumnShading()
//       SoilTemperature()
//       SoilTemperatureInit()
//
#include <cmath>
#include "global.h"
#include "exceptions.h"
#include "Simulation.hpp"

using namespace std;

void SoilTemperatureInit(Simulation &);

// SoilTemperature_2
void EnergyBalance(Simulation &, uint32_t, int, int, double, double, const double &, double[20]);

// SoilTemperature_3
void SoilHeatFlux(double, int, int, int, int, double);

void PredictEmergence(Simulation &, unsigned int, int);

//////////////////////////
void ColumnShading(State &state, double rracol[20], double day_emerge, double row_space, unsigned int plant_row_column)
//     This function computes light interception by crop canopy and shading
//  of soil columns by the plants. It is called from SimulateThisDay().
//
//     The following global variables are referenced here:
//       isw, LeafAreaIndex, PlantHeight, PlantRowColumn, nk.
//     The following global variables are set here:    LightIntercept, rracol.
{
    //     Before emergence: no light interception and no shading. LightIntercept is
    //  assigned zero, and the rracol array is assigned 1.
    if (state.daynum < day_emerge || isw <= 0 || day_emerge <= 0)
    {
        LightIntercept = 0;
        for (int k = 0; k < nk; k++)
            rracol[k] = 1;
        return;
    }
    //     Compute the maximum leaf area index until this day (lmax).
    static double lmax; // maximum leaf area index.
    if (state.daynum <= day_emerge)
        lmax = 0;
    else if (state.leaf_area_index > lmax)
        lmax = state.leaf_area_index;
    //     Light interception is computed by two methods:
    //     (1) It is assumed to be proportional to the ratio of plant height to row spacing.
    double zint; // light interception computed from plant height.
    zint = 1.0756 * state.plant_height / row_space;
    //     (2) It is computed as a function of leaf area index. If LeafAreaIndex is not greater
    //  than 0.5 lfint is a linear function of it.
    double lfint; // light interception computed from leaf area index.
    if (state.leaf_area_index <= 0.5)
        lfint = 0.80 * state.leaf_area_index;
    else
        //     If the leaf area index is greater than 0.5, lfint is computed as an exponential
        //  function of LeafAreaIndex.
        lfint = 1 - exp(0.07 - 1.16 * state.leaf_area_index);
    //     If lfint is greater then zint, LightIntercept is their average value.
    //  Otherwise, if the LeafAreaIndex is decreasing, it is lfint. Else it is zint.
    if (lfint > zint)
        LightIntercept = 0.5 * (zint + lfint);
    else if (state.leaf_area_index < lmax)
        LightIntercept = lfint;
    else
        LightIntercept = zint;
    //     The value of LightIntercept is between zero and one.
    if (LightIntercept < 0)
        LightIntercept = 0;
    if (LightIntercept > 1)
        LightIntercept = 1;
    //     Loop of soil columns.
    double sw = 0; // sum of column widths
    double sw0;    // sum of column widths up to location of plant row.
    double sw1;    // distance from middle of a column to the plant row, cm.
    int j, k0;     // number of columns from plant row location.
    for (int k = 0; k < nk; k++)
    {
        if (k <= plant_row_column)
        {
            //     When the column is on the left of the plant row.
            j = plant_row_column - k;
            sw += wk(j, row_space);
            sw0 = sw;
            sw1 = sw - wk(j, row_space) / 2;
            k0 = j;
        }
        else
        {
            //     When the column is on the right of the plant row.
            sw += wk(k, row_space);
            sw1 = sw - sw0 - wk(k, row_space) / 2;
            k0 = k;
        }
        //     Relative shading is computed as a function of sw1 and plant height, modified by LightIntercept.
        //  rracol is the fraction of radiation received by a soil column. Iys minimum value is 0.05 .
        double shade; // relative shading of a soil column.
        if (sw1 >= state.plant_height)
            shade = 0;
        else
        {
            shade = 1 - (sw1 / state.plant_height) * (sw1 / state.plant_height);
            if (LightIntercept < zint && state.leaf_area_index < lmax)
                shade = shade * LightIntercept / zint;
        }
        rracol[k0] = 1 - shade;
        if (rracol[k0] < 0.05)
            rracol[k0] = 0.05;
    } // end for k
}

////////////////////////////////////////////////////////////////////////////////////////////////
void SoilTemperature(Simulation &sim, uint32_t u, double rracol[20])
//     This is the main part of the soil temperature sub-model. It is called daily from
//  SimulateThisDay(). It calls the following functions:
//  EnergyBalance(), PredictEmergence(), SoilHeatFlux(), SoilTemperatureInit().
//
//     The following global variables are referenced here:
//       ActualTranspiration, AirTemp,
//       DewPointTemp, dl, FoliageTemp, isw, nk, nl,
//       PlantRowColumn, ReferenceETP, ReferenceTransp, RelativeHumidity,
//       rracol, SitePar, thad, wk
//     The following global variables are set here:
//       ActualSoilEvaporation, CumEvaporation, DeepSoilTemperature, es, SoilTemp,
//       SoilTempDailyAvrg, VolWaterContent,
{
    State &state = sim.states[u];
    if (u == 0)
        SoilTemperatureInit(sim);
    //     Compute dts, the daily change in deep soil temperature (C), as
    //  a site-dependent function of Daynum.
    double dts = 2 * pi * SitePar[10] / 365 * cos(2 * pi * (sim.day_start + u - SitePar[11]) / 365);
    //     Define iter1 and dlt for hourly time step.
    int iter1 = 24;     // number of iterations per day.
    double dlt = 3600;  // time (seconds) of one iteration.
    int kk = 1;         // number of soil columns for executing computations.
                        //     If there is no canopy cover, no horizontal heat flux is assumed, kk = 1.
                        //  Otherwise it is equal to the number of columns in the slab.
    double shadeav = 0; // average shaded area in all shaded soil columns.
                        //     isw defines the type of soil temperature computation.
    if (isw > 1)
    {
        double shadetot = 0; // sum of shaded area in all shaded soil columns.
        int nshadedcol = 0;  // number of at least partially shaded soil columns.
        kk = nk;
        for (int k = 0; k < nk; k++)
            if (rracol[k] <= 0.99)
            {
                shadetot += 1 - rracol[k];
                nshadedcol++;
            }
        //
        if (nshadedcol > 0)
            shadeav = shadetot / nshadedcol;
    }
    //     Set daily averages of soil temperature to zero.
    for (int l = 0; l < nl; l++)
        for (int k = 0; k < nk; k++)
            SoilTempDailyAvrg[l][k] = 0;
    //     es and ActualSoilEvaporation are computed as the average for the whole soil
    //  slab, weighted by column widths.
    double es = 0; // potential evaporation rate, mm day-1
    state.actual_soil_evaporation = 0;
    //     Start hourly loop of iterations.
    for (int ihr = 0; ihr < iter1; ihr++)
    {
        Hour &hour = state.hours[ihr];
        //     Update the temperature of the last soil layer (lower boundary conditions).
        DeepSoilTemperature += dts * dlt / 86400;
        double etp0 = 0; // actual transpiration (mm s-1) for this hour
        if (state.evapotranspiration > 0.000001)
            etp0 = state.actual_transpiration * hour.ref_et / state.evapotranspiration / dlt;
        //
        //     Compute vertical transport for each column
        //
        for (int k = 0; k < kk; k++)
        {
            //     Set SoilTemp for the lowest soil layer.
            SoilTemp[nl - 1][k] = DeepSoilTemperature;
            //     Compute transpiration from each column, weighted by its relative shading.
            double etp1 = 0; // actual hourly transpiration (mm s-1) for a column.
            if (shadeav > 0.000001)
                etp1 = etp0 * (1 - rracol[k]) / shadeav;
            double ess = 0; //   evaporation rate from surface of a soil column (mm / sec).
            //     The potential evaporation rate (escol1k) from a column is the sum of the radiation
            //  component of the Penman equation(es1hour), multiplied by the relative radiation reaching
            //  this column, and the wind and vapor deficit component of the Penman equation (es2hour).
            double escol1k; // potential evaporation fron soil surface of a column, mm per hour.
            escol1k = hour.et1 * rracol[k] + hour.et2;
            es += escol1k * wk(k, sim.row_space);
            //     Compute actual evaporation from soil surface. update VolWaterContent of
            //  the soil soil cell, and add to daily sum of actual evaporation.
            double evapmax = 0.9 * (VolWaterContent[0][k] - thad[0]) * 10 *
                             dl(0); // maximum possible evaporatio from a soil cell near the surface.
            if (escol1k > evapmax)
                escol1k = evapmax;
            VolWaterContent[0][k] -= 0.1 * escol1k / dl(0);
            state.actual_soil_evaporation += escol1k * wk(k, sim.row_space);
            ess = escol1k / dlt;
            //     Call EnergyBalance to compute soil surface and canopy temperature.
            EnergyBalance(sim, u, ihr, k, ess, etp1, state.plant_height, rracol);

        } // End loop of k
        //     Compute soil temperature flux in the vertical direction.
        //     Assign iv = 1, layer = 0, nn = nl.
        int iv = 1;          // indicates vertical (=1) or horizontal (=0) flux.
        int nn = nl;         // number of array members for heat flux.
        int layer = 0;       // soil layer number
        double tsolav[maxl]; // hourly average soil temperature C, of a soil layer.
                             //     Loop over kk columns, and call SoilHeatFlux().
        for (int k = 0; k < kk; k++)
            SoilHeatFlux(dlt, iv, nn, layer, k, sim.row_space);
        //     If no horizontal heat flux is assumed, make all array members of SoilTemp equal to the
        //  value computed for the first column. Also, do the same for array memebers of VolWaterContent.
        if (isw <= 1)
            for (int l = 0; l < nl; l++)
                for (int k = 1; k < nk; k++)
                {
                    SoilTemp[l][k] = SoilTemp[l][0];
                    if (l == 0)
                        VolWaterContent[l][k] = VolWaterContent[l][0];
                }
        //
        //     Compute horizontal transport for each layer
        //
        //     Compute soil temperature flux in the horizontal direction, when isw = 2.
        //     Assign iv = 0 and nn = nk. Start loop for soil layers, and call SoilHeatFlux.
        if (isw > 1)
        {
            iv = 0;
            nn = nk;
            for (int l = 0; l < nl; l++)
            {
                layer = l;
                SoilHeatFlux(dlt, iv, nn, layer, l, sim.row_space);
            }
        }
        //     Compute average temperature of soil layers, in degrees C.
        for (int l = 0; l < nl; l++)
        {
            tsolav[l] = 0;
            for (int k = 0; k < nk; k++)
            {
                SoilTempDailyAvrg[l][k] += SoilTemp[l][k];
                tsolav[l] += SoilTemp[l][k] - 273.161;
            }
            tsolav[l] = tsolav[l] / nk;
        }
        //     Compute average temperature of foliage, in degrees C. The average is weighted by the canopy
        //  shading of each column, only columns which are shaded 5% or more by canopy are used.
        double tfc = 0;     // average foliage temperature, weighted by shading in each column
        double shading = 0; // sum of shaded area in all shaded columns, used to compute TFC
        for (int k = 0; k < nk; k++)
        {
            if (rracol[k] <= 0.95)
            {
                tfc += (FoliageTemp[k] - 273.161) * (1 - rracol[k]);
                shading += 1 - rracol[k];
            }
        }
        if (shading >= 0.01)
            tfc = tfc / shading;
        //     If emergence date is to be simulated, call PredictEmergence().
        if (isw == 0 && sim.day_start + u >= sim.day_plant)
            PredictEmergence(sim, u, ihr);
    } // end of hourly loop
      //      At the end of the day compute actual daily evaporation and its cumulative sum.
    if (kk == 1)
    {
        es = es / wk(1, sim.row_space);
        state.actual_soil_evaporation /= wk(1, sim.row_space);
    }
    else
    {
        es = es / sim.row_space;
        state.actual_soil_evaporation /= sim.row_space;
    }
    sim.states[u].cumulative_evaporation += state.actual_soil_evaporation;
    if (state.kday > 0)
    {
        Scratch21[u].es = es;
    }
    //  compute daily averages.
    for (int l = 0; l < nl; l++)
        for (int k = 0; k < nk; k++)
            SoilTempDailyAvrg[l][k] = SoilTempDailyAvrg[l][k] / iter1;
}
/*
                           References.

     Benjamin, J.G., Ghaffarzadeh, M.R. and Cruse, R.M. 1990.
  Coupled water and heat transport in ridged soils. Soil Sci. Soc.
  Am. J. 54:963-969.

     Chen, J. 1984. Uncoupled multi-layer model for the transfer of
  sensible and latent heat flux densities from vegetation. Boundary-
  Layer Meteorology 28:213-225.

     Chen, J. 1985. A graphical extrapolation method to determine
  canopy resistance from measured temperature and humidity profiles
  above a crop canopy. Agric. For. Meteorol. 37:75-88.

     Clothier, B.E., Clawson, K.L., Pinter, P.J.Jr., Moran, M.S.,
  Reginato, R.J. and Jackson, R.D. 1986. Estimation of soil heat flux
  from net radiation during the growth of alfalfa. Agric. For.
  Meteorol. 37:319-329.

     Costello, T.A. and Braud, H.J. Jr. 1989. Thermal diffusivity
  of soil by nonlinear regression analysis of soil temperature data.
  Trans. ASAE 32:1281-1286.

     De Vries, D.A. 1963. Thermal properties of soils. In: W.R. Van
  Wijk (ed) Physics of plant environment, North Holland, Amsterdam,
  pp 210-235.

     Deardorff, J.W. 1978. Efficient prediction of ground surface
  temperature and moisture with inclusion of a layer of vegetation.
  J. Geophys. Res. 83 (C4):1889-1903.

     Dong, A., Prashar, C.K. and Grattan, S.R. 1988. Estimation of
  daily and hourly net radiation. CIMIS Final Report June 1988, pp.
  58-79.

     Ephrath, J.E., Goudriaan, J. and Marani, A. 1996. Modelling
  diurnal patterns of air temperature, radiation, wind speed and
  relative humidity by equations from daily characteristics.
  Agricultural Systems 51:377-393.

     Hadas, A. 1974. Problem involved in measuring the soil thermal
  conductivity and diffusivity in a moist soil. Agric. Meteorol.
  13:105-113.

     Hadas, A. 1977. Evaluation of theoretically predicted thermal
  conductivities of soils under field and laboratory conditions. Soil
  Sci. Soc. Am. J. 41:460-466.

     Hanks, R.J., Austin, D.D. and Ondrechen, W.T. 1971. Soil
  temperature estimation by a numerical method. Soil Sci. Soc. Am.
  Proc. 35:665-667.

     Hares, M.A. and Novak, M.D. 1992. Simulation of surface energy
  balance and soil temperature under strip tillage: I. Model
  description. Soil Sci. Soc. Am. J. 56:22-29.

     Hares, M.A. and Novak, M.D. 1992. Simulation of surface energy
  balance and soil temperature under strip tillage: II. Field test.
  Soil Sci. Soc. Am. J. 56:29-36.

     Horton, E. and Wierenga, P.J. 1983. Estimating the soil heat
  flux from observations of soil temperature near the surface. Soil
  Sci. Soc. Am. J. 47:14-20.

     Horton, E., Wierenga, P.J. and Nielsen, D.R. 1983. Evaluation
  of methods for determining apparent thermal diffusivity of soil
  near the surface. Soil Sci. Soc. Am. J. 47:25-32.

     Horton, R. 1989. Canopy shading effects on soil heat and water
  flow. Soil Sci. Soc. Am. J. 53:669-679.

     Horton, R., and Chung, S-O, 1991. Soil Heat Flow. Ch. 17 in: Hanks,
  J., and Ritchie, J.T., (Eds.) Modeling Plant and Soil Systems. Am. Soc.
  Agron., Madison, WI, pp 397-438.

     Iqbal, M. 1983. An Introduction to Solar Radiation. Academic
  Press.

     Kimball, B.A., Jackson, R.D., Reginato, R.J., Nakayama, F.S.
  and Idso, S.B. 1976. Comparison of field-measured and calculated
  soil heat fluxes. Soil Sci. Soc. Am. J. 40:18-28.

     Lettau, B. 1971. Determination of the thermal diffusivity in
  the upper layers of a natural ground cover. Soil Sci. 112:173-177.

     Monin, A.S. 1973. Boundary layers in planetary atmospheres. In:
  P. Morrel (ed.), Dynamic meteorology, D. Reidel Publishing Company,
  Boston, pp. 419-458.

     Spitters, C.J.T., Toussaint, H.A.J.M. and Goudriaan, J. 1986.
  Separating the diffuse and direct component of global radiation and
  its implications for modeling canopy photosynthesis. Part I.
  Components of incoming radiation. Agric. For. Meteorol. 38:217-229.

     Wierenga, P.J. and de Wit, C.T. 1970. Simulation of heat flow
  in soils. Soil Sci. Soc. Am. Proc. 34:845-848.

     Wierenga, P.J., Hagan, R.M. and Nielsen, D.R. 1970. Soil
  temperature profiles during infiltration and redistribution of cool
  and warm irrigation water. Water Resour. Res. 6:230-238.

     Wierenga, P.J., Nielsen, D.R. and Hagan, R.M. 1969. Thermal
  properties of soil based upon field and laboratory measurements.
  Soil Sci. Soc. Am. Proc. 33:354-360.
*/
////////////////////////////////////////////////////////////////////////
void SoilTemperatureInit(Simulation &sim)
//     This function is called from SoilTemperature() at the start of the simulation. It sets
//  initial values to soil and canopy temperatures.
//     The following arguments are set in this function:
//  jt1, jt2 - input of start and stop of output of soil temperatures.
//
//     The following global variables are referenced here:
//  Clim (structure), nl, SitePar.
//     The following global variables are set here:
//  DeepSoilTemperature, SoilTemp.
{
    //     Compute initial values of soil temperature: It is assumed that at the start of simulation
    //  the temperature of the first soil layer (upper boundary) is equal to the average air temperature
    //  of the previous five days (if climate data not available - start from first climate data).
    //     NOTE: For a good simulation of soil temperature, it is recommended to start simulation at
    //  least 10 days before planting date. This means that climate data should be available for
    //  this period. This is especially important if emergence date has to be simulated.
    int idd = 0;     // number of days minus 4 from start of simulation.
    double tsi1 = 0; // Upper boundary (surface layer) initial soil temperature, C.
    for (int i = 0; i < 5; i++)
        tsi1 += sim.climate[i].Tmax + sim.climate[i].Tmin;
    tsi1 = tsi1 / 10;
    //     The temperature of the last soil layer (lower boundary) is computed as a sinusoidal function
    //  of day of year, with site-specific parameters.
    DeepSoilTemperature = SitePar[9] + SitePar[10] * sin(2 * pi * (sim.day_start - SitePar[11]) / 365);
    //     SoilTemp is assigned to all columns, converted to degrees K.
    tsi1 += 273.161;
    DeepSoilTemperature += 273.161;
    for (int l = 0; l < nl; l++)
    {
        //     The temperatures of the other soil layers are linearly interpolated.
        //  tsi = computed initial soil temperature, C, for each layer
        double tsi = ((nl - l - 1) * tsi1 + l * DeepSoilTemperature) / (nl - 1);
        for (int k = 0; k < nk; k++)
            SoilTemp[l][k] = tsi;
    }
}
