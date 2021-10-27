// CottonModel.cpp  "CottonModel" is the simulation module of the Cotton2K
// model. Version 4.0 was written by A. Marani June 2004. Compiled by Microsoft
// Visual C++.Net 2003. This file defines the following functions in this file:
//       Class C2KApp:
//    Message map and constructor
//    InitInstance()
//    ExitInstance()
//    GetProfilesList()
//    RunTheModel()
//    DailySimulation()
//    DoAdjustments()
//    SimulateThisDay()
//
#include "Cottonmodel.h"

#include <boost/algorithm/string.hpp>
#include <iostream>

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
/////////////////////////////////////////////////////////////////////////////
bool C2KApp::InitInstance(std::string JobFileName)
//     InitInstance() starts the application. it calls the functions:
//  GetProfilesList() and RunTheModel().
//     Global variables set:  FrameTitle
//
{
    //     Build and set title of the main frame.
    FrameTitle = "";
    char Slash = '\\';
    int nLen = JobFileName.length();
    //
    for (int i = nLen - 1; i >= 0; i--) {
        if (JobFileName[i] == Slash) break;
        FrameTitle = JobFileName[i] + FrameTitle;
    }
    //     Get the list of the "Profiles" to run.
    GetProfilesList(JobFileName.c_str());
    //     Now run the model.
    RunTheModel();
    return true;
}
/////////////////////////////////////////////////////////////////////////
void C2KApp::GetProfilesList(std::string JobFileName)
//     Function GetProfilesList() opens the "JOB" file, reads it, and gets
//  from it the profile list. It is called from InitInstance().
//     Input argument: name of the JOB file (with path)
//
{
    //     Open the selected Job file for input.
    ifstream DataFile(JobFileName, ios::in);
    if (DataFile.fail()) {
        //     When the file can not be opened:
        DataFile.close();
        std::cerr << JobFileName << " cannot be open!" << std::endl;
        return;
    }
    char m_TempString[100];
    //     Skip the first line of the file.
    DataFile.getline(m_TempString, 99);
    //     Read the simulation profiles from the next lines of the file, and
    //     store them
    //  in the StringArray ProfileArray.
    ProfileArray.clear();
    int readlength = 99;
    while (readlength > 0) {
        //     Create array of Profile names
        DataFile.get(m_TempString, 21);
        if (DataFile.eof() == 1) break;
        std::string m_String = m_TempString;
        boost::algorithm::erase_all(m_String, " ");
        readlength = m_String.length();
        if (readlength > 4)
            ProfileArray.push_back(m_String.substr(0, readlength - 4));
        DataFile.getline(m_TempString, 99);
    }
    DataFile.close();
}
/////////////////////////////////////////////////////////////////////////////
void C2KApp::RunTheModel()
//     This function calls the following functions for each profile:
//          ReadInput(), DailySimulation() and  DataOutput()
//     Global variables set: ProfileName
//     Global variables referenced:  DayFinish, DayStart
//
{
    for (int i = 0; i < ProfileArray.size(); i++) {
        ProfileName = ProfileArray[i];
        //     Read the input data for this simulation
        ReadInput();
        //     Do daily simulations
        DailySimulation();
        //     Write output data
        DataOutput();
    }
    //     End of simulation of all profiles
    std::cout << " Simulation Ended." << std::endl;
}
////////////////////////////////////////////////////////////////////////////////
void C2KApp::DailySimulation()
//     This function controls the dynamic phase of the simulation, allowing
//  for in-run adjustments when there is an input of plant map adjustments.
//     It calls the functions:
//        DoAdjustments(), SimulateThisDay(), WriteStateVariables()
//
//     The following global variable are referenced:   DayStart, Kday.
//     The following global variable are set:   bEnd, Daynum.
//
{
    Daynum = DayStart - 1;
    bEnd = false;
    //     Start the daily loop. If variable bEnd has been assigned a value
    //  of true end simulation.
    while (!bEnd) {
        bool bAdjustToDo = DoAdjustments();
        std::cout << "\rProcess (" << Daynum - DayStart + 2 << "/"
                  << DayFinish - DayStart + 1 << ")" << std::flush;
        //     Execute simulation for this day.
        SimulateThisDay();
        //     If there are pending plant adjustments, call
        //     WriteStateVariables() to write
        //  state variables of this day in a scratch file.
        if (bAdjustToDo) WriteStateVariables(true);
    }  // end while
}
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
//     DoyToDate(), ColumnShading(), DayClim(), SoilTemperature(),
//     SoilProcedures(), SoilNitrogen(), SoilSum(), PhysiologicalAge(), Pix(),
//     Defoliate(), Stress(), GetNetPhotosynthesis(), PlantGrowth(),
//     CottonPhenology(), PlantNitrogen(), CheckDryMatterBal(),
//     PlantNitrogenBal(), SoilNitrogenBal(), SoilNitrogenAverage(),
//     DailyOutput();
//
//     The following global variables are referenced here:  DayEmerge,
//     DayFinish,
//  DayStart, iyear, Kday, LastDayWeatherData, LeafAreaIndex, OutIndex, pixday.
//
//     The following global variables are set here:
//  bEnd, Date, DayInc, Daynum, DayOfSimulation, isw, Kday.
//
{
    //    Compute Daynum (day of year), Date, and DayOfSimulation (days from
    //    start of simulation).
    Daynum++;
    Date = DoyToDate(Daynum, iyear);
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
        if (OutIndex[20] > 0) {
            PlantNitrogenBal();     // checks plant nitrogen balance.
            SoilNitrogenBal();      // checks soil nitrogen balance.
            SoilNitrogenAverage();  // computes average soil nitrogen by layers.
        }
    }
    //     Call DailyOutput for reporting some simulation data for this day.
    DailyOutput();
    //     Check if the date to stop simulation has been reached, or if this is
    //     the last day
    //  with available weather data. Simulation will also stop when no leaves
    //  remain on the plant. bEnd = true  indicates stopping this simulation.
    //
    if (Daynum >= DayFinish || Daynum >= LastDayWeatherData ||
        (Kday > 10 && LeafAreaIndex < 0.0002))
        bEnd = true;
}

int main(int argc, char **argv) {
    if (argc < 2) {
        std::cerr << "Job file path should be provided!" << std::endl;
        return 1;
    }
    std::string job_filename = argv[1];
    C2KApp app;
    app.InitInstance(job_filename);
    return 0;
}
