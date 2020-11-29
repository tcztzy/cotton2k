//  File WriteOutput_1.cpp
//
//   functions in this file:
// OpenOutputFiles()
// DailyOutput{}
// output1()
// DataOutput()
// WriteLine22()
//
#include <ctime>
#include "CottonSimulation.h"
#include "GeneralFunctions.h"

/////////////////////////////////////////////////////////////
void OpenOutputFiles(const string& m_fileDesc, const string& ProfileName)
//     This function opens the output files that will be used by
//  this simulation. It is called from function ReadProfileFile().
//
//     The following global variables are referenced here:     OutIndex[]
//     The following argument is used:
//       m_fileDesc = the description of this profile file.
//
//     The file extensions are as follows:
//                                         Extension  OutIndex 
//                                         ---------   ------
//  Output files which are always opened:
//     Summary of the input data               B01
//     Summary output file                     S01
//     Output data for plotting charts         PLT
//     Detailed output file                    F01
//
//  Optional output files opened if the relevant OutIndex parameter is not zero:
//     Soil map data for plotting              SMP     8,9,10,11 or 12
//     Plant map data for plotting             PLM     7 or 13
//
//  For advanced users only: Use a text editor for editing the profile file:
//     Computed hourly weather data            TM1     15
//     Detailed computed hourly weather data   TM2     15 (value=2)
//     Computed hourly soil temperature data   TMS     16
//     Simulated soil water data               WAT     17
//     Carbon balance of plants data           CHB     18
//     Leaf water potential and related data   LWP     19
//     Simulated soil nitrogen data            NB0     20
//     Simulated soil nitrogen data            NB1     20
//     Simulated soil nitrogen data            NB2     20
//     Simulated soil nitrogen data            NB3     20
//     Simulated soil nitrogen data            NB4     20
//     Simulated plant root data               RUT     22
//
{
//     Get the date of running this simulation.  
      time_t rawtime;
      struct tm * timeinfo;
      char buffer [80];
      time(&rawtime);
      timeinfo = localtime(&rawtime);
      strftime(buffer, 80, "%A, %B %d, %Y", timeinfo);
      string Rundate = buffer;
//     Open file F01 and write header to it.
      ofstream File46(fs::path("output") / (ProfileName + ".F01"), ios::out);
      File46.unsetf(ios::left);
	  File46.width(50);
	  File46 << "COTTON2K Version 4.0 (2004)" << endl;
	  File46.width(62);
	  File46 << "A simulation model for irrigated cotton in arid regions" << endl;
	  File46.width(50);
	  File46 << "Written by Avishalom Marani" << endl << endl;
      File46.setf(ios::left);
	  File46 << "Profile Name:    ";
	  File46.width(20);
	  File46 << ProfileName << endl;
	  File46 << "Simulation Date: ";
	  File46.width(30);
	  File46 << Rundate << endl;
	  File46 << "Description:     ";
	  File46.width(55);
	  File46 << m_fileDesc << endl << endl;
      if (DayEmerge > 0)
		    File46 << " Defined Germination on " << DoyToDate(DayEmerge, iyear) << 
			" (Day of Year = " << DayEmerge << " )" << endl << endl;
//     The following lines in the header depend on whether metric or 
//  English units are used, or if results are per plant or per unit area.
      if (OutIndex[2] == 0)       // per plant
	  {
         if (OutIndex[1] == 0)    // metric units
		 {
            File46 << "                       Per plant                      Per Plant" << endl;
			File46 << "                  ------------------  -----------------------------------" << endl;
			File46 << " DAE     Date     Height  LAI  Nodes  Sites  Squares  Green  Open   Absc  Yield" << endl;
            File46 << "                   (cm)                               bolls  bolls  sites kg/ha" << endl;
	 	 }
         else                     // English units
		 {
            File46 << "                       Per plant                      Per Plant" << endl;
			File46 << "                  ------------------  -----------------------------------" << endl;
			File46 << " DAE     Date     Height  LAI  Nodes  Sites  Squares  Green  Open   Absc  Yield" << endl;
            File46 << "                   (in)                               bolls  bolls  sites lbs/a" << endl;
 		 }
	  }
      else                        // per unit area
	  {
         if (OutIndex[1] == 0)    // metric units 
		 {
            File46 << "                       Per plant                      Per sq m" << endl;
			File46 << "                  ------------------  -----------------------------------" << endl;
			File46 << " DAE     Date     Height  LAI  Nodes  Sites  Squares  Green  Open   Absc  Yield" << endl;
            File46 << "                   (cm)                               bolls  bolls  sites kg/ha" << endl;
		 }
         else                     // English units
		 {
            File46 << "                       Per plant                    1000s per acre" << endl;
			File46 << "                  ------------------  -----------------------------------" << endl;
			File46 << " DAE     Date     Height  LAI  Nodes  Sites  Squares  Green  Open   Absc  Yield" << endl;
            File46 << "                   (in)                               bolls  bolls  sites lbs/a" << endl;
		 }
	  }
	  File46 << endl;
//
// Open the input information file B01, and write its header. 
      ofstream File20(fs::path("output") / (ProfileName + ".B01"), ios::out);
      File20.unsetf(ios::left);
	  File20.width(50);
	  File20 << "COTTON2K Version 4.0 (2003)" << endl;
	  File20.width(62);
	  File20 << "A simulation model for irrigated cotton in arid regions" << endl;
	  File20.width(50);
	  File20 << "Written by Avishalom Marani" << endl << endl;
      File20.setf(ios::left);
	  File20 << "Profile Name:    ";
	  File20.width(20);
	  File20 << ProfileName << endl;
	  File20 << "Simulation Date: ";
	  File20.width(30);
	  File20 << Rundate << endl;
	  File20 << "Description:     ";
	  File20.width(55);
	  File20 << m_fileDesc << endl << endl;
//
// Open the summary output file S01, and write its header. 
      ofstream File22(fs::path("output") / (ProfileName + ".S01"), ios::out);
	  File22 << endl << "                         SUMMARY OF SIMULATION RUN" << endl;
      File22.setf(ios::left);
	  File22 << "Profile Name:    ";
	  File22.width(20);
	  File22 << ProfileName << endl;
	  File22 << "Simulation Date: ";
	  File22.width(30);
	  File22 << Rundate << endl;
	  File22 << "Description:     ";
	  File22.width(55);
	  File22 << m_fileDesc << endl << endl;
      if (DayEmerge > 0)
		    File22 << " Defined Germination on " << DoyToDate(DayEmerge, iyear) << 
			" (Day of Year = " << DayEmerge << " )" << endl << endl;
      File22 << "                                   PLANT                   GREEN  OPEN" << endl;
      File22 << "     EVENT         DATE       DAE  HEIGHT NODES  LAI  SQARS BOLLS  BOLLS YIELD" << endl;
      if (OutIndex[1] == 0 && OutIndex[2] == 0)
	     File22 << "                                    (cm)               (- per plant -)  (kg/ha)" << endl;
      else if (OutIndex[1] == 1 && OutIndex[2] == 0) 
         File22 << "                                    (in)               (- per plant -)  (lb/ac)" << endl;
      else if (OutIndex[1] == 0 && OutIndex[2] == 1) 
         File22 << "                                    (cm)               (- per sq.m  -)  (kg/ha)" << endl;
      else if (OutIndex[1] == 1 && OutIndex[2] == 1) 
         File22 << "                                    (in)               (1000 per acre)  (lb/ac)" << endl;
//
//      The following files are optional output. Opening these files
//  is determined by the output flags in the Profile file.
//      When any of the output flags 8 to 12 is non-zero, open output
//  file for soil map plotting data, as unit 23.
      int imap; // indicator for soil-map or plant-map output.
      imap = OutIndex[8] + OutIndex[9] + OutIndex[10] + OutIndex[11] + OutIndex[12];
      if (imap > 0) 
         ofstream File23(fs::path("output") / (ProfileName + ".SMP"), ios::out);
//
//      When any of the output flags 7 or 13 is non-zero, open output
//  file for plant map data, as unit 24.
      imap = OutIndex[7] + OutIndex[13];
      if (imap > 0) 
         ofstream File24(fs::path("output") / (ProfileName + ".PLM"), ios::out);
//
//      The following files are used for checking during model
// development. They are not intended for general use.
//      For the soil water data output, two files are opened. The extensions are '.WAT'
//  and ".WA2". These files are opened when the flag OutIndex[17] is non-zero.
      if (OutIndex[17] > 0) 
      {
         ofstream File42(fs::path("output") / (ProfileName + ".WAT"), ios::out);
         ofstream File56(fs::path("output") / (ProfileName + ".WA2"), ios::out);
         File56 << "    DATE    vertical water content  for columns 0 and 9";
		 File56 << endl;
      }
//
//      For the root data output,  the extension is '.RUT'. The file
//  is opened when the flag OutIndex[22] is non-zero.
      if (OutIndex[22] > 0) 
         ofstream File34(fs::path("output") / (ProfileName + ".RUT"), ios::out);
//
//      When the output flag 15 is non-zero: file '.TM1' which is used for output
//  of computed 24-hour weather data for plotting is opened, and its header line is written.
//      When the output flag 15 is greater than 1, file '.TM2' which is used for detailed output 
//  of computed 24-hour weather data (for checking) is also opened.
//      When the output flag 16 is non-zero, file '.TMS', which
//  outputs simulated soil temperatue data, is opened as unit 19.
      if (OutIndex[15] > 0) 
	  {
         ofstream File17(fs::path("output") / (ProfileName + ".TM1"), ios::out);
         File17 <<  "     TIME         RAD         TMP          RH      WND        ETREF";
		 File17 << endl;
	  }
      if (OutIndex[15] > 1) 
         ofstream File18(fs::path("output") / (ProfileName + ".TM2"), ios::out);
//
//      When the output flag 18 is non-zero, file '.CHB', which outputs plant carbon balance data, 
//  is opened, and the heading for this file is written.
      if (OutIndex[18] > 0) 
	  {
         ofstream File36(fs::path("output") / (ProfileName + ".CHB"), ios::out);
         File36 << "    DATE       CDSTEM   CDLEAF     CDPET    CDROOT    CSTRES    STEMWT    LEAFWT     PETWT    ROOTWT    SUMPDR      NR";
		 File36 << endl;
	  }
//
//      When the output flag 19 is non-zero, file '.LWP', which outputs simulated leaf water
//  potential and related data, is opened, and the heading for this file is written.
      if (OutIndex[19] > 0) 
	  {
         ofstream File44(fs::path("output") / (ProfileName + ".LWP"), ios::out);
         File44 << " KDAY  PSISOIL   COND      RSOIL     RROOT     RSHOOT    RLEAF     PSILN     PSILD ";
		 File44 << endl;
	  }
//
//      When the output flag 20 is non-zero, files '.NB0' to 'NB4', which output nitrogen data, 
//  are opened and the headings for these files are written.
      if (OutIndex[20] > 0) 
	  {
         ofstream File35(fs::path("output") / (ProfileName + ".NB0"), ios::out);
         File35 << " KDAY  VNO3C                              VNH4C";
		 File35 << endl;
         File35 << "       0 - 30  30 - 60  60 - 90 90 - 120  0 - 30  30 - 60  60 - 90 90 - 120";
		 File35 << endl;
         File35 << "                   in parts per million (per volume)" << endl;
         ofstream File47(fs::path("output") / (ProfileName + ".NB1"), ios::out);
         File47 << " KDAY  BALSN    FERN     ORGN   UPTAKEN    SOILN   SNLOSS    BALPN    ADDN   PLANTN    NLOSS";
		 File47 << endl;
         ofstream File37(fs::path("output") / (ProfileName + ".NB2"), ios::out);
         File37 << "    CHECKING  NITROGEN IN PLANT PARTS  ";
		 File37 << endl;
         File37 << " KDAY SLEAFN   PETN  STEMN  ROOTN  BURRN  SEEDN SQRN   LEAFCN  PETCN STEMCN ROOTCN  BURCN SEEDCN PETCNO3";
		 File37 << endl;
         ofstream File38(fs::path("output") / (ProfileName + ".NB3"), ios::out);
         File38 << "    CHECKING  NITROGEN REQUIREMENTS";
		 File38 << endl;
         File38 << " KDAY RQNLEF RQNPET RQNSTM RQNRUT   REQV RQNSQR RQNSED RQNBUR   REQF   REQTOT";
		 File38 << endl;
         ofstream File39(fs::path("output") / (ProfileName + ".NB4"), ios::out);
         File39 << "    CHECKING  NITROGEN SUPPLY POOLS";
		 File39 << endl;
         File39 << " KDAY SUPNO3 SUPNH4   UPTN LEAFRS PETRS  STEMRS ROOTRS BURRES   RESN   PETN  NPOOL";
		 File39 << endl;
	  }
}
////////////////////////////////////////////////////////////////////////////
void DailyOutput(const string& ProfileName, const string& Date)
//     DailyOutput() writes output at the end of each day. It is called from SimulateThisDay().
//  This function calls WriteStateVariables(), cotplt(), and output1().
//
//     The following global variables are referenced here:
//       Daynum, DayOfSimulation, DayStartPlantMaps, DayStopPlantMaps, dl, NumFruitBranches, 
//       NumPreFruNodes, OutIndex, PlantMapFreq, RowSpace, VolNo3NContent, wk,    
//     The following global variables are set here:      MainStemNodes, SumNO3N90.
{
//     1. Compute some variables needed for output:
//     Main stem node count (MainStemNodes) starts from the cotyledonary node ( # = 0 ).
//  After first square, last prefruiting node becomes first fruiting branch node.
      if ( NumFruitBranches[0] <= 0 ) 
         MainStemNodes = NumPreFruNodes - 1;
      else
         MainStemNodes = NumPreFruNodes + NumFruitBranches[0] - 2;
//     Compute SumNO3N90 as the total nitrate N, up to a depth of 90 cm of soil, in kg / ha
      SumNO3N90 = 0; 
      double sumdl = 0; // depth to the end of a layer
      for (int l = 0; l < nl; l++)
	  {
	     sumdl += dl[l];
	     if (sumdl > 90) 
		    break;
         for (int k = 0; k < nk; k++)
            SumNO3N90 += VolNo3NContent[l][k] * dl[l] * wk[k];
      }
      SumNO3N90 = SumNO3N90 * 100 / RowSpace;
//
//     2. Call WriteStateVariables() which saves values of all important state variables for
//  this day in structure Scratch21, which will be used for output at the end of the simulation.
      if (DayOfSimulation > 0)
           WriteStateVariables(false, Date);
//
//     3. Write plant map output by calling cotplt():
//     If the output flags (OutIndex(7) or OutIndex(13)) indicate that plant map output is  
//  requested, and if it is after first square - call function cotplt() if day of year is between 
//  DayStartPlantMaps and DayStopPlantMaps. 
      if ( FruitingCode[0][0][0] > 0 && (OutIndex[7] + OutIndex[13]) > 0) 
         if ( Daynum >= DayStartPlantMaps && Daynum <= DayStopPlantMaps ) 
		 {
            int idum; // days from start of output period.
            idum = Daynum - DayStartPlantMaps;
//     Plant maps will be executed at PlantMapFreq day intervals:
            if ( (idum % PlantMapFreq) == 0 ) 
			{
               if ( OutIndex[7] > 0 ) 
                  cotplt(1, ProfileName, Date);
               if ( OutIndex[13] > 0 ) 
                  cotplt(4, ProfileName, Date);
            }
         }
//
//     4. Call output1() to write output to F01 and S01 files:
      output1(ProfileName, Date);
//     Check for end of the simulation and assign true to bEnd.
      if ( Daynum >= DayFinish || LeafAreaIndex < 0.0002 || Daynum >= LastDayWeatherData ) 
          bEnd = true;
}
//////////////////////////
void output1(const string& ProfileName, const string& Date)
//     This function is a collection of write statements for output of model results. 
//  It writes daily data to files F01 and S01. It is called each day from DailyOutput().
//
//     The following global variables are referenced here:
//       AbscisedFruitSites, Date, Daynum, DayEmerge, FirstBloom, FirstSquare, isw, Kday, 
//       LeafAreaIndex, LastDayWeatherData, LintYield, MainStemNodes, NumFruitSites, NumGreenBolls, 
//       NumOpenBolls, NumSquares, OutIndex, PlantHeight, PlantPopulation
{
       ofstream File46(fs::path("output") / (ProfileName + ".F01"), ios::app);
//     Date of first square and date of first bloom are written to file F01.
      if ( Daynum == FirstSquare  )
			 File46 << " ** First square on " << Date << " **" << endl;
      if ( Daynum == FirstBloom ) 
			 File46 << " ** First bloom  on " << Date << " **" << endl;
//     If it is after emergence, the following output will be written:
      if ( Daynum >= DayEmerge && isw > 0 ) 
	  {
//     Convertion for site numbers etc. to per area basis (metric or English units)
          double conversion = 1; // per plant
          if (OutIndex[2] == 1) 
		  {
             if (OutIndex[1] == 0)          // metric units
                conversion = PlantPopulation / 10000;   // per sq m
             else                           // English units
                conversion = PlantPopulation / 2469;    // 1000s per acre
          }
//     Write Kday, Date, plant height, LAI, main stem nodes, numbers (per plant, or converted
//  to per m2, or to 1000s per acre) of fruiting sites, squares, green bolls, open bolls, 
//  and cumulated abscised sites, and LintYield (in kg per ha, or lbs per acre).
          File46.unsetf(ios::left);
          File46.width(4);
          File46 << Kday;
          File46.width(12);
          File46 << Date;
		  File46.setf(ios::fixed);
          File46.width(7);
          File46.precision(1);
          if (OutIndex[1] == 0)     // metric units
			 File46 << PlantHeight;
		  else
			 File46 << PlantHeight / 2.54;   // English
          File46.width(6);
          File46.precision(2);
	      File46 << LeafAreaIndex;
          File46.width(6);
          File46.precision(0);
		  File46 << MainStemNodes;
          File46.width(9);
          if (OutIndex[1] == 0)    
             File46.precision(2);
		  else
             File46.precision(1);
		  File46 << NumFruitSites * conversion;
          File46.width(8);
		  File46 << NumSquares * conversion;
          File46.width(7);
		  File46 << NumGreenBolls * conversion;
          File46.width(7);
		  File46 << NumOpenBolls * conversion;
          File46.width(7);
		  File46 << AbscisedFruitSites * conversion;
          File46.width(6);
          File46.precision(0);
          if (OutIndex[1] == 0)    
			 File46 << LintYield;
		  else
			 File46 << LintYield * 0.893;
		  File46 << endl;
//
//    Also, write information to the 'summary file' (S01) if this is
//  the day of first square, first bloom, or day of last actual weather data.
          ofstream File22(fs::path("output") / (ProfileName + ".S01"), ios::app);
          if ( Daynum == LastDayWeatherData )
		  {
			 File22 << " Last act. wthr.";
             File22.width(12);
             File22 << Date;
		     File22.setf(ios::fixed);
             File22.width(5);
             File22 << Kday;
             File22.width(7);
             File22.precision(1);
             if (OutIndex[1] == 0)     // metric units
			    File22 << PlantHeight;
		     else
			    File22 << PlantHeight / 2.54;   // English
             File22.width(6);
             File22.precision(0);
             File22 << MainStemNodes;
             File22.width(6);
             File22.precision(3);
             File22 << LeafAreaIndex;
             File22.width(6);
             if (OutIndex[1] == 0)     // metric units
                 File22.precision(1);
			 else
                 File22.precision(0);
             File22 << NumSquares * conversion;
             File22.width(6);
             File22 << NumGreenBolls * conversion;
             File22.width(6);
             File22 << NumOpenBolls * conversion;
             File22.width(7);
             File22.precision(0);
             if (OutIndex[1] == 0)    
			    File22 << LintYield;
		     else
			    File22 << LintYield * 0.893;
             File22 << endl;
		  }
      }
}
/////////////////////////////////////////////////////////////////////////////////////
tuple<string> DataOutput(const string& ProfileName, const string& Date, const int& DayStart)
//     This function is called from RunSimulation() at the end of the simulation. It
//  gets the data from structure Scratch21 and writes summary data in file *.S01.
//     It calls the functions WriteLine22(), outputplt(), output2(), output3(), output4(),
//  output5(), output6(), and output7() for executing further output.
//
//     The following global variables are referenced here:
//       DefoliationDate, OutIndex, PlantPopulation, NumOpenBolls.
//
//     The following global variables are set here (extracted from the structure Scratch21):
//       Date, Kday, LeafAreaIndex, LintYield, MainStemNodes, PlantHeight.
//
{
      bool ifg1 = false;  // flags to check if output has been written.
      bool ifg2 = false;  // flags to check if output has been written.
      bool ifg3 = false;  // flags to check if output has been written.
      bool ifg4 = false;  // flags to check if output has been written.
      double i00; // number of squares, per unit area
      double i01; // number of green bolls, per unit area.
      double i02; // number of open bolls, per unitn area.
      double sixpct = NumOpenBolls * 0.6; // 60 percent of the final number of open bolls.
      string date = Date;
      ofstream File22(fs::path("output") / (ProfileName + ".S01"), ios::app);
//     Start reading data from struct Scratch21.
//     Read data for each day, and compute the local variables i00,
//  i01, i02. Convert from per plant to per sq m, or to English units.
      for (int irec = 0; irec < 400; irec++)
	  {		  
		 if (Scratch21[irec].daynum <= 0)
			 break;
         date = Scratch21[irec].date;
         i01 = Scratch21[irec].numGreenBolls;
         Kday = Scratch21[irec].kday;
         LeafAreaIndex = Scratch21[irec].leafAreaIndex;
         MainStemNodes = Scratch21[irec].mainStemNodes;
         i02 = Scratch21[irec].numOpenBolls;
         i00 = Scratch21[irec].numSquares;
         LintYield = Scratch21[irec].lintYield;
         PlantHeight = Scratch21[irec].plantHeight;
//     Convert height and LintYield to English units if necessary.
         if (OutIndex[1] == 1) 
		 {
            PlantHeight = PlantHeight / 2.54;
            LintYield = LintYield * 0.893;
		 }
         if (OutIndex[2] == 1) 
//     Convert to per area basis, metric or English units.
		 {
            if (OutIndex[1] == 0) // metric units: per sq m
			{
               i00 = i00 * PlantPopulation / 10000;
               i01 = i01 * PlantPopulation / 10000;
               i02 = i02 * PlantPopulation / 10000;
			}
            else                  // English units: 1000s per acre
			{
               i00 = i00 * PlantPopulation / 2469;
               i01 = i01 * PlantPopulation / 2469;
               i02 = i02 * PlantPopulation / 2469;
			}
		 }
//     Report first day of open bolls, date of 60% open bolls, and date
//  of defoliation.
         if ( ! ifg1 && Scratch21[irec].daynum == FirstSquare)
		 {
            File22 << " First square    " << date;
			WriteLine22(File22, i00, i01, i02);
            ifg1 = true;
		 }
         if ( ! ifg2 && Scratch21[irec].daynum == FirstBloom)
		 {
            File22 << " First Bloom     " << date;
			WriteLine22(File22, i00, i01, i02);
            ifg2 = true;
		 }
         if ( ! ifg3 && i02 > 0.1 ) 
		 {
            File22 << " 1st open boll   " << date;
			WriteLine22(File22, i00, i01, i02);
            ifg3 = true;
		 }
         if ( ! ifg4 && sixpct <= i02 && i02 > 0. ) 
		 {
            File22 << " 60% open bolls  " << date;
			WriteLine22(File22, i00, i01, i02);
            ifg4 = true;
		 }
         for (int i = 0; i < 5; i++)
		 {
            if (Scratch21[irec].daynum == DefoliationDate[i]) 
			{
               File22 << " Defoliation     " << date;
			   WriteLine22(File22, i00, i01, i02);
			}
		 } 
	  }
//     When end of file is reached report final LintYield in file 22.
      File22 << " Max yield on    " << date;
	  WriteLine22(File22, i00, i01, i02);
//     Call procedures for printing output
      outputplt(ProfileName);
      if ( OutIndex[6] > 0 ) 
		   output2(ProfileName);
      if ( OutIndex[3] > 0 )
		  output3(ProfileName);
      if ( OutIndex[4] > 0 )
		  output4(ProfileName);
      if ( OutIndex[5] > 0 )
		   output5(ProfileName);
      output6(ProfileName);
      if ( OutIndex[8] + OutIndex[9] + OutIndex[10] + OutIndex[11] + OutIndex[12] > 0)
           output7(ProfileName, DayStart);
      return make_tuple(date);
}
////////////////////////
void WriteLine22(ofstream &File22, double i00, double i01, double i02)
//     This function writes a formatted line in file *.S01. It is called from DataOutput().
//     Global variables referenced: 
//       Kday, LeafAreaIndex, LintYield, MainStemNodes, PlantHeight. 
{
		    File22.setf(ios::fixed); 
			File22.width(5);
			File22 << Kday;
			File22.width(7);
			File22.precision(1);
			File22 << PlantHeight;
			File22.width(6);
			File22 << MainStemNodes;
			File22.width(6);
			File22.precision(2);
			File22 << LeafAreaIndex;
            if (OutIndex[1] == 0 || OutIndex[2] == 0) 
			   File22.precision(1);
			else
			   File22.precision(0);
			File22.width(6);
			File22 << i00;
			File22.width(6);
			File22 << i01;
			File22.width(6);
			File22 << i02;
		    File22.precision(0);
			File22.width(6);
			File22 << LintYield << endl;
}

