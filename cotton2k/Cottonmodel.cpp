// CottonModel.cpp  "CottonModel" is the simulation module of the Cotton2K
// model. Version 4.0 was written by A. Marani June 2004. Compiled by Microsoft
// Visual C++.Net 2003. This file defines the following functions in this file:
//       Class C2KApp:
//    Message map and constructor
//    InitInstance()
//    ExitInstance()
//    GetProfilesList()
//    DailySimulation()
//    DoAdjustments()
//    SimulateThisDay()
//
#include "Cottonmodel.h"

#include "CottonSimulation.h"
#include "GeneralFunctions.h"

//
#ifdef _DEBUG
#define new DEBUG_NEW
#endif
///////////////////////////////////////////////////////////////////////////////
//    Class C2KApp
///////////////////////////////////////////////////////////////////////////////
// C2KApp construction
C2KApp::C2KApp() {}
///////////////////////////////////////////////////////////////////////////////
bool C2KApp::DoAdjustments()
//     This function is called from DailySimulation(). It checks if plant
//     adjustment data
//  are available for this day and calls the necessary functions to compute
//  adjustment. It calls PlantAdjustments(), SimulateThisDay(),
//  WriteStateVariables()
//
//     The following global variable are referenced:  DayEmerge, Daynum, Kday,
//     The following global variable are set:   MapDataDate, nadj,
//     NumAdjustDays.
//
{
    //     Check if plant map data are available for this day. If there are no
    //     more adjustments, return.
    static int kprevadj =
        0;  // day after emergence of the previous plant map adjustment.
    int sumsad = 0;  // sum for checking if any adjustments are due
    for (int i = 0; i < 30; i++) sumsad += MapDataDate[i];
    if (sumsad <= 0) return false;
    //     Loop for all adjustment data, and check if there is an adjustment for
    //     this day.
    for (int i = 0; i < 30; i++) {
        if (Daynum == MapDataDate[i]) {
            //     Compute NumAdjustDays, the number of days for retroactive
            //     adjustment. This can not be more
            //  than 12 days, limited by the date of the previous adjustment.
            NumAdjustDays = Kday - kprevadj;
            if (NumAdjustDays > 12) NumAdjustDays = 12;
            //     Loop for six possible adjustments. On each iteration call
            //     first PlantAdjustments(), which
            //  will assign true to nadj(jj) if adjustment is necessary, and
            //  compute the necessary parameters.
            for (int jj = 0; jj < 5; jj++) {
                PlantAdjustments(i, jj);
                //     If adjustment is necessary, rerun the simulation for the
                //     previous NumAdjustDays (number
                //  of days) and call WriteStateVariables() to write state
                //  variables in scratch structure.
                if (nadj[jj])
                    for (int j1 = 0; j1 < NumAdjustDays; j1++) {
                        SimulateThisDay();
                        if (Kday > 0) WriteStateVariables(true);
                    }  // end for j1, and if nadj
            }          // end for jj
            //     After finishing this adjustment date, set kprevadj (date of
            //     previous adjustment, to
            //  be used for next adjustment), and assign zero to the present
            //  msadte, and to array nadj[].
            kprevadj = MapDataDate[i] - DayEmerge + 1;
            MapDataDate[i] = 0;
            for (int jj = 0; jj < 5; jj++) nadj[jj] = false;
            continue;
        }  // end if Daynum
    }      // end do i
    return true;
}
//////////////////////////////////////////////////
void C2KApp::SimulateThisDay()
//     This function executes all the simulation computations in a day. It is
//     called from
//  DailySimulation(), and DoAdjustments().   It calls the following functions:
//     ColumnShading(), DayClim(), SoilTemperature(),
//     SoilProcedures(), SoilNitrogen(), SoilSum(), PhysiologicalAge(), Pix(),
//     Defoliate(), Stress(), GetNetPhotosynthesis(), PlantGrowth(),
//     CottonPhenology(), PlantNitrogen(), CheckDryMatterBal(),
//     PlantNitrogenBal(), SoilNitrogenBal(), SoilNitrogenAverage(),
//
//     The following global variables are referenced here:  DayEmerge,
//     DayFinish,
//  DayStart, iyear, Kday, LastDayWeatherData, LeafAreaIndex, pixday.
//
//     The following global variables are set here:
//  bEnd, Date, DayInc, Daynum, DayOfSimulation, isw, Kday.
//
{
    //    Compute Daynum (day of year), Date, and DayOfSimulation (days from
    //    start of simulation).
    Daynum++;
    DayOfSimulation = Daynum - DayStart + 1;
    //    Compute Kday (days from emergence).
    if (DayEmerge <= 0)
        Kday = 0;
    else
        Kday = Daynum - DayEmerge + 1;
    if (Kday < 0) Kday = 0;
    //     The following functions are executed each day (also before
    //     emergence).
    ColumnShading();    // computes light interception and soil shading.
    DayClim();          // computes climate variables for today.
    SoilTemperature();  // executes all modules of soil and canopy temperature.
    SoilProcedures();   // executes all other soil processes.
    SoilNitrogen();     // computes nitrogen transformations in the soil.
    SoilSum();          // computes totals of water and N in the soil.
    //     The following is executed each day after plant emergence:
    if (Daynum >= DayEmerge && isw > 0) {
        //     If this day is after emergence, assign to isw the value of 2.
        isw = 2;
        DayInc = PhysiologicalAge();  // computes physiological age
        if (pixday[0] > 0) Pix();     // effects of pix applied.
        Defoliate();                  // effects of defoliants applied.
        Stress();                     // computes water stress factors.
        GetNetPhotosynthesis();       // computes net photosynthesis.
        PlantGrowth();                // executes all modules of plant growth.
        CottonPhenology();    // executes all modules of plant phenology.
        PlantNitrogen();      // computes plant nitrogen allocation.
        CheckDryMatterBal();  // checks plant dry matter balance.
        //     If the relevant output flag is not zero, compute soil nitrogen
        //     balance and soil
        //  nitrogen averages by layer, and write this information to files.
        if (false) {
            PlantNitrogenBal();     // checks plant nitrogen balance.
            SoilNitrogenBal();      // checks soil nitrogen balance.
            SoilNitrogenAverage();  // computes average soil nitrogen by layers.
        }
    }
    //     Check if the date to stop simulation has been reached, or if this is
    //     the last day
    //  with available weather data. Simulation will also stop when no leaves
    //  remain on the plant. bEnd = true  indicates stopping this simulation.
    //
    if (Daynum >= DayFinish || Daynum >= LastDayWeatherData ||
        (Kday > 10 && LeafAreaIndex < 0.0002))
        bEnd = true;
}
