//  GettingInput_1.cpp
//
//   functions in this file:
// ReadInput()
// ReadProfileFile()
// ReadCalibrationData()
// InitializeGrid()
// WriteInitialInputData()
//
#include "CottonSimulation.h"
#include "GeneralFunctions.h"
#include "resource.h"
#include <math.h>

//
// Definitions of File scope variables:
   bool bLat,                // true if latitude is south, false if north
        bLong;               // true if longitude is west, false if east
   int  nSiteNum,            // index number for site. 
        LastDayOfActualWeather,   // last day of actual weather. 
        nVarNum;             // index number for cultivar. 
   double SkipRowWidth,      // the smaller distance between skip rows, cm
        PlantsPerM;          // average number of plants pre meter of row.
   string m_mulchdata,      // string containing input data of mulching
        VarName,              // name of the cultivar
        SiteName;             // name of the site
/////////////////////////////////////////////////////////////
void ReadInput(const string& ProfileName)
//     This is the main function for reading input. It is called from RunTheModel().
//     The following global variables are set here:
//        PlantWeightAtStart , SoilNitrogenAtStart
//     The following global variables are referenced here:
//        ReserveC, TotalLeafWeight, TotalRootWeight, TotalSoilNh4N, TotalSoilNo3N, 
//        TotalSoilUreaN, TotalStemWeight.
//
{
//     The following functions are called to read initial values of some variables from 
//  input files, or initialize them otherwise.
	InitializeGlobal();
	ReadProfileFile(ProfileName);
	ReadCalibrationData();
	LastDayOfActualWeather = OpenClimateFile();
	InitializeGrid();
	ReadSoilImpedance();
    WriteInitialInputData(ProfileName);
	InitSoil();
	ReadAgriculturalInput(ProfileName);
	ReadPlantMapInput();
	InitializeSoilData();
	InitializeSoilTemperature();
	InitializeRootData();
//     initialize some variables at the start of simulation.
    SoilNitrogenAtStart = TotalSoilNo3N + TotalSoilNh4N + TotalSoilUreaN;
    PlantWeightAtStart = TotalRootWeight + TotalStemWeight + TotalLeafWeight + ReserveC;
}
/////////////////////////////////////////////////////////////////////////////
void ReadProfileFile(const string& ProfileName)
//     This function opens and reads the profile file. It is called from ReadInput().
//  It calls GetLineData(), DateToDoy() and OpenOutputFiles().
//     The following global or file-scope variables are set here:
//  ActWthFileName, AgrInputFileName, bEnd, bLat, bLong, CO2EnrichmentFactor, DayEmerge,
//  DayEndCO2, DayFinish, DayPlant, DayStart, DayStartCO2, DayStartPlantMaps, DayStartSoilMaps,
//  DayStopPlantMaps, DayStopSoilMaps, Elevation, isw, iyear, Latitude, Longitude, m_mulchdata, 
//  MulchIndicator, nSiteNum, nVarNum, OutIndex, PlantmapFileName, PlantMapFreq, PlantsPerM, 
//  PrdWthFileName, RowSpace, SkipRowWidth, SoilHydFileName, SoilInitFileName, SoilMapFreq.
//
{
    fs::path strFileName = fs::path("profiles") / (ProfileName + ".pro"); // file name with path
    CString strMessage;
//     If file does not exist, or can not be opened, display message
    if (!fs::exists(strFileName))
    {
         AfxFormatString1(strMessage, IDS_FILE_NOT_EXISTS, strFileName.string().c_str());
         AfxMessageBox(strMessage);
         return;
    }
    ifstream DataFile(strFileName, ios::in);
    if ( DataFile.fail() )
    {
          AfxMessageBox(("Error opening " + strFileName.string() + ".").c_str());
          DataFile.close();
          return;
    }
//     Line #1: Read file description.
	string Dummy = GetLineData(DataFile);
    string m_fileDesc; // Description of the Profile file
    if (Dummy.length() > 20)
    {
        m_fileDesc = Dummy.substr(20);
        m_fileDesc.erase(m_fileDesc.find_last_not_of(" \r\n\t\f\v") + 1);
    }
    else   
        m_fileDesc = "";
//     Line #2: Read dates of emergence, start and end of simulation, and planting date. 
	string DateEmerge, DateSimStart, DateSimEnd, DatePlant;
    Dummy = GetLineData(DataFile);
    int nLength = Dummy.length();
    if (nLength >= 14)
    {
       DateEmerge =  Dummy.substr(0,14);
       DateEmerge.erase(remove(DateEmerge.begin(), DateEmerge.end(), ' '), DateEmerge.end());
    }
    if (nLength >= 23)
    {
       DateSimStart = Dummy.substr(15,14);
       DateSimStart.erase(remove(DateSimStart.begin(), DateSimStart.end(), ' '), DateSimStart.end());
    }
    if (nLength >= 38)
    {
       DateSimEnd = Dummy.substr(30,14);
       DateSimEnd.erase(remove(DateSimEnd.begin(), DateSimEnd.end(), ' '), DateSimEnd.end());
    }
    if (nLength >= 53)
    {
       DatePlant = Dummy.substr(45,14);
       DatePlant.erase(remove(DatePlant.begin(), DatePlant.end(), ' '), DatePlant.end());
    }
    else 
       DatePlant = "";
//     For advanced users only: if there is CO2 enrichment, read also CO2 factor, DOY dates 
//	for start and stop of enrichment (these are left blank if there is no CO2 enrichment).
    if (nLength > 76)
	{
       string ttt = Dummy.substr(60,10);
       ttt.erase(remove(ttt.begin(), ttt.end(), ' '), ttt.end());
       CO2EnrichmentFactor = atof(ttt.c_str());
       ttt = Dummy.substr(70,5);
       ttt.erase(remove(ttt.begin(), ttt.end(), ' '), ttt.end());
       DayStartCO2 = atoi(ttt.c_str());
       ttt = Dummy.substr(75);
       ttt.erase(remove(ttt.begin(), ttt.end(), ' '), ttt.end());
       DayEndCO2 = atoi(ttt.c_str());
	}
    else  
       CO2EnrichmentFactor = 0;
//     Line #3: Names of weather files: actual and predicted. 
    Dummy = GetLineData(DataFile);
    nLength = Dummy.length();
    if (nLength > 1)
    {
        ActWthFileName = Dummy.substr(0,20).c_str();
        ActWthFileName.Remove(' ');
    }
    if (nLength > 21)
    {
        PrdWthFileName = Dummy.substr(20,20).c_str();
        PrdWthFileName.Remove(' ');
    }
    else 
        PrdWthFileName = "";
//     For advanced users only: If soil mulch is used, read relevant parameters.
    if (nLength > 41)
	{
           m_mulchdata = Dummy.substr(40);
	       MulchIndicator = atoi (m_mulchdata.substr(0,10).c_str() );
           if (MulchIndicator > 0)
           {
              MulchTranSW = atof (m_mulchdata.substr(10, 10).c_str());
              MulchTranLW = atof (m_mulchdata.substr(20, 10).c_str());
			  DayStartMulch = atoi (m_mulchdata.substr(30, 5).c_str());
			  DayEndMulch = atoi (m_mulchdata.substr(35, 5).c_str());
              if (DayEndMulch <= 0)
                  DayEndMulch = DateToDoy(DateSimEnd.c_str(), iyear);
           }
	}
    else
	{
		  m_mulchdata = "";
	      MulchIndicator = 0;
	}
//     Line #4: Names of files for soil hydraulic data, soil initial
//  conditions, agricultural input, and plant map adjustment.
    Dummy = GetLineData(DataFile);
    nLength = Dummy.length();
    if (nLength > 1)
    {
       SoilHydFileName = Dummy.substr(0,20).c_str();
       SoilHydFileName.Remove(' ');
    }
    if (nLength > 20)
    {
       SoilInitFileName = Dummy.substr(20,20).c_str();
       SoilInitFileName.Remove(' ');
    }
    if (nLength > 40)
    {
       AgrInputFileName = Dummy.substr(40,20).c_str();
       AgrInputFileName.Remove(' ');
    }
    if (nLength > 60)
    {
       PlantmapFileName = Dummy.substr(60).c_str();
       PlantmapFileName.Remove(' ');
    }
    else 
       PlantmapFileName = "";
//     Line #5: Latitude and longitude of this site, elevation (in m
//  above sea level), and the index number for this geographic site.
    Dummy = GetLineData(DataFile);
    nLength = Dummy.length();
	bLat = false;
	bLong = false;
    if (nLength > 1)
	{
          Latitude = atof(Dummy.substr(0,10).c_str());
		  if (Latitude < 0)
		  {
			  bLat = true;
			  Latitude = -Latitude;
		  }
	}
    if (nLength >= 20)
	{
          Longitude = atof(Dummy.substr(10,10).c_str());
		  if (Longitude < 0)
		  {
			  bLong = true;
			  Longitude = -Longitude;
		  }
	}
    if (nLength >= 30)
          Elevation = atof(Dummy.substr(20,10).c_str());
    if (nLength > 30)
          nSiteNum = atoi(Dummy.substr(30).c_str());
//     Line #6: Row spacing in cm, skip-row spacing in cm (blank or 0 
//  for no skip rows), number of plants per meter of row, and index
//  number for the cultivar.
    Dummy = GetLineData(DataFile);
    nLength = Dummy.length();
    if (nLength > 1)
          RowSpace =  atof(Dummy.substr(0,10).c_str());
    if (nLength >= 20)
          SkipRowWidth = atof(Dummy.substr(10,10).c_str());
    if (nLength >= 30)
          PlantsPerM = atof(Dummy.substr(20,10).c_str());
    if (nLength > 30)
          nVarNum = atoi(Dummy.substr(30).c_str());
//     Line #7: Frequency in days for output of soil maps, and dates 
//  for start and stop of this output (blank or 0 if no such output is
//  required. Same is repeated for output of plant maps.
    string SoilMapStartDate, SoilMapStopDate, PlantMapStartDate, PlantMapStopDate; 
    Dummy = GetLineData(DataFile);
    nLength = Dummy.length();
    if (nLength > 9)
          SoilMapFreq = atoi(Dummy.substr(0,10).c_str());
    if (nLength >= 16)
    {
          SoilMapStartDate =  Dummy.substr(14,11);
          SoilMapStartDate.erase(remove(SoilMapStartDate.begin(), SoilMapStartDate.end(), ' '), SoilMapStartDate.end());
    }
    if (nLength >= 31)
    {
          SoilMapStopDate =  Dummy.substr(29,11);
          SoilMapStopDate.erase(remove(SoilMapStopDate.begin(), SoilMapStopDate.end(), ' '), SoilMapStopDate.end());
    }
    if (nLength > 41)
          PlantMapFreq = atoi(Dummy.substr(40,10).c_str());
    if (nLength >= 56)
    {
          PlantMapStartDate =  Dummy.substr(54,11);
          PlantMapStartDate.erase(remove(PlantMapStartDate.begin(), PlantMapStartDate.end(), ' '), PlantMapStartDate.end());
    }
    if (nLength >= 71)
    {
          PlantMapStopDate =  Dummy.substr(69,11);
          PlantMapStopDate.erase(remove(PlantMapStopDate.begin(), PlantMapStopDate.end(), ' '), PlantMapStopDate.end());
    }
//     Line #8: 23 output flags.
// - Line 8 consists of zeros and ones.  "1" tells the simulator to
//   produce a particular report and "0" indicates that no report
//   should be produced.  
    for (int n = 0; n < 24; n++)
        OutIndex[n] = 0;
    Dummy = GetLineData(DataFile);
    nLength = Dummy.length();
    for (int n = 0; n < 23 && 3*n < nLength ; n++)
    {
          int n1 = 3 * n;
          OutIndex[n+1] = atoi(Dummy.substr(n1, 3).c_str()); 
    }
    DataFile.close();
//     Calendar dates of emergence, planting, start and stop of simulation, start and stop of 
// output of soil slab and plant maps are converted to DOY dates by calling function DateToDoy.
    iyear = atoi(DateSimStart.substr(7,4).c_str());
    DayStart = DateToDoy(DateSimStart.c_str(), iyear);
    DayEmerge = DateToDoy(DateEmerge.c_str(), iyear);
    DayFinish = DateToDoy(DateSimEnd.c_str(), iyear);
    DayPlant = DateToDoy(DatePlant.c_str(), iyear);
    DayStartSoilMaps = DateToDoy(SoilMapStartDate.c_str(), iyear);
    DayStopSoilMaps = DateToDoy(SoilMapStopDate.c_str(), iyear);
    DayStartPlantMaps = DateToDoy(PlantMapStartDate.c_str(), iyear);
    DayStopPlantMaps = DateToDoy(PlantMapStopDate.c_str(), iyear);
//     If the output frequency indicators are zero, they are set to 999.
      if ( SoilMapFreq <= 0 ) 
		   SoilMapFreq = 999;
      if ( PlantMapFreq <= 0 )
		   PlantMapFreq = 999;
//     If the date of emergence has not been given, emergence will be
//  simulated by the model. In this case, isw = 0, and a check is
//  performed to make sure that the date of planting has been given.
	  if ( DayEmerge <= 0 ) 
	  {
         isw = 0;
         if (DayPlant <= 0) 
		 {
            string msg = " planting date or emergence date must";
            msg += " be given in the profile file !!";
			AfxMessageBox(msg.c_str());
            bEnd = true;
         }
	  }
//     If the date of emergence has been given in the input: isw = 1 if simulation
//  starts before emergence, or isw = 2 if simulation starts at emergence.
      else if ( DayEmerge > DayStart ) 
         isw = 1;
      else
      {
         isw = 2;
         Kday = 1;
      }
//     Call function OpenOutputFiles() to open the output files.
      OpenOutputFiles(m_fileDesc, ProfileName);
}
//////////////////////////////////////////////////////////
void ReadCalibrationData()
//     This function reads the values of the calibration parameters
//  from input files. It is called from ReadInput(). It calls GetLineData().
//
//     The following global or file-scope variables are set here:
//  SiteName, SitePar, VarName, VarPar
{
//     Open file of variety file list. 
    fs::path strFileName = fs::path("data") / "vars" / "varlist.dat";
    CString strMessage;
//     If file does not exist, display message and and open a new file
    if (!fs::exists(strFileName))
    {
         AfxFormatString1(strMessage, IDS_FILE_NOT_EXISTS, strFileName.string().c_str());
         AfxMessageBox(strMessage);
         return;
    }
    ifstream DataFile(strFileName, ios::in);
    if ( DataFile.fail() )
    {
          AfxMessageBox(("Error opening " + strFileName.string() + ".").c_str());
          DataFile.close();
          return;
    }
//
    string Dummy, VarFile;
    // FIXME: if it go through to the last blank line, it wouldn't break because it doesn't reach the eof
    for (int m_idx = 0; m_idx < 1000; m_idx++)
    {
        if (DataFile.eof() == 1)
            break;
        Dummy = GetLineData(DataFile);
        int nLength = Dummy.length();
	    int num;
	    string Name, FileName;
        if (nLength >= 4)
		{
           num = atoi(Dummy.substr(0,4).c_str());
		}
        if (nLength >= 25)
		{
           Name = Dummy.substr(5,20);
           Name.erase(remove(Name.begin(), Name.end(), ' '), Name.end());
		}
        if (nLength >= 45)
		{
           FileName = Dummy.substr(40,20);
           FileName.erase(remove(FileName.begin(), FileName.end(), ' '), FileName.end());
		}
	    if (num == nVarNum)
		{
           VarFile = FileName;
		   VarName = Name;
		   break;
		}
	}
    DataFile.close();
//
    strFileName = fs::path("data") / "vars" / VarFile;
//  If file does not exist, or can not be opened, display message
    if (!fs::exists(strFileName))
    {
         AfxFormatString1(strMessage, IDS_FILE_NOT_EXISTS, strFileName.string().c_str());
         AfxMessageBox(strMessage);
         return;
    }
    ifstream DataFile1(strFileName, ios::in);
    if ( DataFile1.fail() )
    {
          AfxMessageBox(("Error opening " + strFileName.string() + ".").c_str());
          DataFile1.close();
          return;
    }
//     Read values of variety related parameters
	GetLineData(DataFile1);  // skip 1st line
	for (int i = 1; i <= 60; i++)
	{
	    Dummy = GetLineData(DataFile1); 
		VarPar[i] = atof (Dummy.substr(0,20).c_str());
	}
    DataFile1.close();
//     Open file of site file list. 
    strFileName = fs::path("data") / "site" / "sitelist.dat";
//     If file does not exist, or can not be opened, display message 
    if (!fs::exists(strFileName))
    {
         AfxFormatString1(strMessage, IDS_FILE_NOT_EXISTS, strFileName.string().c_str());
         AfxMessageBox(strMessage);
         return;
    }
    string SiteFile;
    ifstream DataFile2(strFileName, ios::in);
    if ( DataFile2.fail() )
    {
          AfxMessageBox(("Error opening " + strFileName.string() + ".").c_str());
          DataFile.close();
          return;
    }
//
    for (int m_idx = 0; m_idx < 1000; m_idx++)
    {
        if (DataFile2.eof() == 1)
            break;

        Dummy = GetLineData(DataFile2);
        int nLength = Dummy.length();
	    int num;
	    string Name, FileName;
        if (nLength >= 4)
		{
           num =  atoi(Dummy.substr(0,4).c_str());
		}
        if (nLength >= 25)
		{
           Name = Dummy.substr(5,20);
           Name.erase(remove(Name.begin(), Name.end(), ' '), Name.end());
		}
        if (nLength >= 45)
		{
           FileName = Dummy.substr(40,20);
           FileName.erase(remove(FileName.begin(), FileName.end(), ' '), FileName.end());
		}
	    if (num == nSiteNum)
		{
           SiteFile = FileName;
		   SiteName = Name;
		   break;
		}
	}
    DataFile2.close();
//
    strFileName = fs::path("data") / "site" / SiteFile;
//     If file does not exist, or can not be opened, display message 
    if (!fs::exists(strFileName))
    {
         AfxFormatString1(strMessage, IDS_FILE_NOT_EXISTS, strFileName.string().c_str());
         AfxMessageBox(strMessage);
         return;
    }
    ifstream DataFile3(strFileName, ios::in);
    if ( DataFile3.fail() )
    {
          AfxMessageBox(("Error opening " + strFileName.string() + ".").c_str());
          DataFile3.close();
          return;
    }
//     Read values of site related parameters
	GetLineData(DataFile3);  // skip 1st line
	for (int i = 1; i <= 20; i++)
	{
	    Dummy = GetLineData(DataFile3); 
		SitePar[i] = atof (Dummy.substr(0,20).c_str());
	}
    DataFile3.close();
}
//////////////////////////////////////////////////////////
void InitializeGrid()
//     This function initializes the soil grid variables. It is executed once 
//  at the beginning of the simulation. It is called from ReadInput().
//
//     The following global or file-scope variables are set here:
//  DensityFactor, dl, nk, nl, PerPlantArea, PlantPopulation, 
//  PlantRowColumn, PlantRowLocation, RowSpace, wk. 
//     The following global variables are referenced here:
//  PlantsPerM, SkipRowWidth, VarPar, maxk, maxl.
{
//     PlantRowLocation is the distance from edge of slab, cm, of the plant row.
    PlantRowLocation = 0.5 * RowSpace;
    if (SkipRowWidth > 1) 
	{
//     If there is a skiprow arrangement, RowSpace and PlantRowLocation are redefined.
       RowSpace = 0.5 * (RowSpace + SkipRowWidth); // actual width of the soil slab (cm)
       PlantRowLocation = 0.5 * SkipRowWidth;
    }
//     Compute PlantPopulation - number of plants per hectar, and PerPlantArea - the average 
//  surface area per plant, in dm2, and the empirical plant density factor (DensityFactor). 
//  This factor will be used to express the effect of plant density on some plant 
//  growth rate functions.  Note that DensityFactor = 1 for 5 plants per sq m (or 50000 per ha).
	PlantPopulation = PlantsPerM / RowSpace * 1000000;
    PerPlantArea = 1000000 / PlantPopulation;
    DensityFactor = exp(VarPar[1] * (5 - PlantPopulation / 10000));
//     Define the numbers of rows and columns in the soil slab (nl, nk).
//  Define the depth, in cm, of consecutive nl layers. 
//     Note that maxl and maxk are defined as constants in file "global.h".
    nl = maxl;
	nk = maxk;
	dl[0] = 2;
	dl[1] = 2;
	dl[2] = 2;
	dl[3] = 4;
	for (int l = 4; l < maxl-2; l++) 
        dl[l] = 5;
	dl[maxl-2] = 10;
	dl[maxl-1] = 10;
//      The width of the slab columns is computed by dividing the row
//  spacing by the number of columns. It is assumed that slab width is
//  equal to the average row spacing, and column widths are uniform.
//      Note: wk is an array - to enable the option of non-uniform 
//  column widths in the future.
//      PlantRowColumn (the column including the plant row) is now computed from 
//  PlantRowLocation (the distance of the plant row from the edge of the slab).
      double sumwk = 0; // sum of column widths
      PlantRowColumn = 0;
      for ( int k = 0; k < nk; k++)
	  {
         wk[k] = RowSpace / nk;
         sumwk = sumwk + wk[k];
         if ( PlantRowColumn == 0 && sumwk > PlantRowLocation ) 
		 {
            if ((sumwk-PlantRowLocation) >  (0.5*wk[k])) 
                PlantRowColumn = k - 1;
            else                                        
                PlantRowColumn = k;
		 }
	  }
}
//////////////////////////////////////////////////////////
void WriteInitialInputData(const string& ProfileName)
//     This function writes the input data to File20 (*.B01). It is executed once 
//  at the beginning of the simulation. It is called by ReadInput().
//
//     The following global or file-scope variables are set here:
//  DayEndMulch, DayStartMulch, MulchTranLW, MulchTranSW.
//     The following global or file-scope variables are referenced here:
//  ActWthFileName, AgrInputFileName, CO2EnrichmentFactor, DayEmerge, DayEndCO2, DayFinish, 
//  DayPlant, DayStart, DayStartCO2, Elevation, iyear, Latitude, Longitude, m_mulchdata, 
//  maxk, MulchIndicator, OutIndex, PerPlantArea, PlantmapFileName, PlantsPerM, 
//  PrdWthFileName, RowSpace, SiteName, SkipRowWidth, SoilHydFileName,
//  SoilInitFileName, VarName.
{
      ofstream File20(fs::path("output") / (ProfileName + ".B01"), ios::app);
      File20 << "    Latitude.. "; 
      File20.setf(ios::fixed);
	  File20.width(8);
	  File20.precision(2);
	  File20 << fabs(Latitude);
	  File20.width(7);
	  if (bLat)
		  File20 << "  South";
      else
		  File20 << "  North";
      File20 << "         Longitude.. "; 
	  File20.width(8);
	  File20.precision(2);
	  File20 << fabs(Longitude);
	  File20.width(7);
	  if (bLong)
		  File20 << "   West";
      else
		  File20 << "   East";
	  File20 << endl;
//
      if (OutIndex[1] == 0)     // write profile data in metric units
	  {
         File20 << "    Elevation (m).... "; 
	     File20.width(8);
	     File20.precision(2);
	     File20 << Elevation;
	  }
	  else
	  {
         File20 << "    Elevation (ft).... "; 
	     File20.width(8);
	     File20.precision(2);
	     File20 << Elevation * 3.28;
	  }
	  File20 << endl;
//
	  File20 << "    Start Simulation " << DoyToDate(DayStart, iyear);
	  File20 << "       Stop Simulation.... " << DoyToDate(DayFinish, iyear) << endl;
//
	  File20 << "    Planting date...." << DoyToDate(DayPlant, iyear);
	  if ( DayEmerge <= 0 ) 
	     File20 << "       Emergence date is simulated   " << endl;
      else
         File20 << "       Emergence date...   " << DoyToDate(DayEmerge, iyear) << endl;
//
      if (OutIndex[1] == 0)     // write profile data in metric units
	  {
         File20 << "    Row Spacing (cm) "; 
	     File20.width(8);
	     File20.precision(2);
	     File20 << RowSpace;
         File20 << "          Plants Per Row-m... "; 
	     File20.width(8);
	     File20 << PlantsPerM << endl;
         File20 << "    Skip Width (cm). "; 
	     File20.width(8);
	     File20 << SkipRowWidth;
         File20 << "          Plants Per Ha...... "; 
	     File20.width(8);
	     File20.precision(1);
	     File20 << PlantPopulation << endl;
	  }
	  else
	  {
         File20 << "    Row Spacing (in) "; 
	     File20.width(8);
	     File20.precision(2);
	     File20 << RowSpace / 2.54;
         File20 << "          Plants Per Row-ft.. "; 
	     File20.width(8);
	     File20 << PlantsPerM * 0.305 << endl;
         File20 << "    Skip Width (in). "; 
	     File20.width(8);
	     File20 << SkipRowWidth / 2.54;
         File20 << "          Plants Per Acre.... "; 
	     File20.width(8);
	     File20.precision(1);
	     File20 << PlantPopulation * 0.405 << endl;
	  }
//
     if ( CO2EnrichmentFactor > 1 )     //  write CO2 enrichment input data.
	  {
		 File20 << "          CO2 enrichment factor              ...  ";
	     File20.width(8);
	     File20.precision(4);
	     File20 << CO2EnrichmentFactor << endl;
		 File20 << "    from ...... " << DoyToDate(DayStartCO2, iyear);
		 File20 << "    to ...... " << DoyToDate(DayEndCO2, iyear) << endl;
	  }
//
	  if (m_mulchdata.length() > 0 && MulchIndicator > 0)
	  {
			  File20 << "   Polyethylene mulch cover. Transmissivity values are: " << endl;
			  File20 << " For short waves:  ";
	  	      File20.width(8);
	          File20.precision(3);
			  File20 << MulchTranSW;
			  File20 << " For long waves:  ";
	  	      File20.width(8);
	          File20.precision(3);
			  File20 << MulchTranLW << endl;
			  File20 << " From Day of Year  ";
	  	      File20.width(4);
			  File20 << DayStartMulch;
			  File20 << " to Day of Year  ";
	  	      File20.width(4);
			  File20 << DayEndMulch << endl;
			  if (MulchIndicator == 1)
				  File20 << " All soil surface covered by mulch." << endl;
			  else  
			  {
	  	          File20.width(6);
	              File20.precision(2);
				  if (MulchIndicator == 2)
				      File20 << RowSpace / maxk;
				  else if (MulchIndicator == 3)
				      File20 << 2 * RowSpace / maxk;
				  File20 << " cm on each side of planr rows not covered by mulch." << endl;
			  }
	  }
//     Write names of the other input files
      if (ActWthFileName.GetLength() > 0)
      {
         File20 << "    Actual Weather Input File:     " << ActWthFileName << endl; 
         File20 << "    Last date read from Actual Weather File: " 
                << DoyToDate(LastDayOfActualWeather, iyear) << endl; 
      } 
      if (PrdWthFileName.GetLength() > 0)
         File20 << "    Predicted Weather Input File:  " << PrdWthFileName << endl; 
      File20 << "    Cultural Input File:           " << AgrInputFileName << endl; 
      File20 << "    Initial Soil Data Input File:  " << SoilInitFileName << endl; 
      File20 << "    Soil Hydrology Input File:     " << SoilHydFileName << endl; 
      if (PlantmapFileName.GetLength() > 0)
         File20 << "    Plant Map Adjustment File:     " << PlantmapFileName << endl << endl; 
//   Write names of the site and the variety
      File20 << "    Site...     " << SiteName << endl; 
      File20 << "    Variety...  " << VarName << endl << endl; 
}