//     File GeneralFunctions.cpp contains the following functions, which are used in
//  several places in the model:
//  String handling -
//      GetLineData()
//  Date conversions - 
//      DateToDoy()
//      DoyToDate()
//      LeapYear()
//  Soil functions -
//      psiq()
//      qpsi()
//      wcond()
//      PsiOsmotic()
//  Extracting climate data -
//      GetFromClim()
//
#include "global.h"
#include "GeneralFunctions.h"

///////////////////////////////////////////////////////////////////////
string GetLineData(ifstream &DataFile)
//     This function safely reads a line from input. All data are read into a temporary character
//  string with the "getline" command. It stops at end of line, and converts it to string.
//     Input        : the ifstream file which is being read.
//     return value : string with the trailing blanks removed.
//
{
     char tmp[91];   // Temporarily stores a character array.
     DataFile.getline(tmp, 90);
     string LineString = tmp;
     LineString.erase(LineString.find_last_not_of(" \r\n\t\f\v") + 1);
     return LineString;
}
/////////////////////////////////////////////////////////////////////////
// Date conversion functions:
//
int DateToDoy(string Date, int m_YearStart) 
//     This function converts calendar date string to day of year, and
//  allows for leap years and dates in the following year.
//     Day of year, sometimes called "Julian date", counts the days from 
//  the beginning of the calendar year. If the simulation continues to the
//  next year, count continues.
//
//     Arguments input:
//            Date = calendar date string as 'dd-MON-yyyy' (11-character string).
//            m_YearStart = year of start of simulation (4 digit integer).
//     Return value:
//            jday =  Day of year (Julian date).
//
//     NOTE: This function does not check for correct input. Make sure that
//  Date is an 11-character string with dd as a valid 2-digit day of month,
//  MON a valid 3-letter month name, and yyyy a 4-digit year number equal to
//  m_YearStart or m_YearStart+1
//
{
	static string MonthName[] =
	{ "JAN","FEB","MAR","APR","MAY","JUN", "JUL","AUG","SEP","OCT","NOV","DEC" };
    static int i0[] =
	{  0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30 };

      Date.erase(remove(Date.begin(), Date.end(), ' '), Date.end());
//   If the date string is blank, or old style date format, return 0.
      if ( Date.empty() )
          return 0;
      else if (Date.substr(2,1) == "/") 
          throw Cotton2KException(Date + " == ERROR: Old date format used " );

	  int day = atoi (Date.substr(0,2).c_str());  //  Convert characters to integers for day.
      int iy = atoi (Date.substr(7).c_str());    //  Convert characters to integers for year.
      string stmon = Date.substr(3,3);  //  Get string month

      int month = 0 ;
	  for (int i = 0; i < 12; i++)
		  if (stmon == MonthName[i])
		  {
			  month = i + 1;
			  break;
		  }

      if (month == 0)
          throw Cotton2KException(" Error in month definition:  " + stmon);
//     Adjust number of days in February for leap years.
      i0[2] = 28 + LeapYear(iy);
//     Compute jday.
      int jday = 0;
      for (int i = 0; i < month; i++)
			jday += i0[i];
      jday += day;
//   Add correction if this is the next calendar year of simulation.
      int nexty = m_YearStart + 1; //    next year
      if (iy == nexty)
           jday += 365 + LeapYear(m_YearStart);
	  return jday;
}
///////////////////////////////////////////////////////////////////////////
string DoyToDate(int Doy, int m_YearStart) 
//     This function converts day of year (sometimes called 'Julian date')
// to calendar date, allowing for leap years and for days in the following year.
//     Arguments input:
//            Doy = Day of year.
//            m_YearStart = year number (4 digit) at start of simulation.
//     Return value:
//         Calender date as an 11-character string 'dd-MON-yyyy'
//
//     NOTE: This function does not check for correct input. Make sure that Doy 
// is an integer between 1 and 730 (or 731 if m_YearStart or m_YearStart+1 is a leap year),
// and m_YearStart is a 4-digit number).
//
{
      static string MonthName[] =
	  { "JAN","FEB","MAR","APR","MAY","JUN", "JUL","AUG","SEP","OCT","NOV","DEC" };
      static int mday[] =
	  { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };

      if (Doy <= 0) 
          return "           ";

//  Adjust number of days in February for leap years.
      mday[1] = 28 + LeapYear(m_YearStart);

	  int iday = Doy;
	  int iy = m_YearStart;
//  If this is the following year
      int jadd = 365 + LeapYear(m_YearStart);
	  if (Doy > jadd) 
	  {  
	     iday = Doy - jadd;
	     iy ++;
         mday[1] = 28 + LeapYear(iy);
	  }
//   Compute month and day from Julian date.
      int month = 1;
      for (int i = 0; i < 12; i++)
	  {
         if ( iday <= mday[i] )
             break;
         month ++;
         iday = iday - mday[i];
	  }
//   Construct date string.
	  string Date1 = "  ";
      char Cday[3];
      sprintf(Cday, "%2d", iday);
	  if ( Cday[0] == ' ')
               Cday[0] = '0';
      Date1.replace(0, 1, 1, Cday[0]);
      Date1.replace(1, 1, 1, Cday[1]);

	  string Date2 = "    ";
      char Cyear[5];
      sprintf(Cyear, "%4d", iy);
      Date2.replace(0, 1, 1, Cyear[0]);
      Date2.replace(1, 1, 1, Cyear[1]);
      Date2.replace(2, 1, 1, Cyear[2]);
      Date2.replace(3, 1, 1, Cyear[3]);

	  string DateOut;
	  DateOut = Date1 + "-" + MonthName[month - 1] + "-" + Date2;
      return DateOut;
}
/////////////////////////////////////////////////////////////////
int LeapYear(int nYear) 
//    This function returns 1 if this is a leap year, or 0 if not.
//    Argument input:
//           nYear = the year (4 digit integer).
//    Return value:
//           nLeap =  1 if this is a leap year, 0 if not.
//
{
      int nLeap = 0;
      if ( nYear % 4 == 0 )
          nLeap = 1;
      if ( nYear % 100 == 0 ) 
          nLeap = 0; 
      if ( nYear % 400 == 0 )
          nLeap = 1;

	  return nLeap;
}
//////////////////////////////////////////////////////////////////
double psiq ( double q, double qr, double qsat, double alpha, double beta )
//     This function computes soil water matric potential (in bars)
//  for a given value of soil water content, using the Van-Genuchten equation.
//
//     The following arguments are used:
//        alpha, beta  - parameters of the van-genuchten equation.
//        q - soil water content, cm3 cm-3.
//        qr - residual water content, cm3 cm-3.
//        qsat - saturated water content, cm3 cm-3.
//
{
//      For very low values of water content (near the residual water
//  content) psiq is -500000 bars, and for saturated or higher water
//  content psiq is -0.00001 bars.
      if ( (q-qr) < 0.00001 )
          return -500000;
      else if ( q >= qsat ) 
          return -0.00001;
//     The following equation is used (FORTRAN notation):
//      PSIX = (((QSAT-QR) / (Q-QR))**(1/GAMA) - 1) **(1/BETA) / ALPHA
      double gama = 1 - 1 / beta;
      double gaminv = 1 / gama;
      double term = (qsat - qr) / ( q - qr);  //  intermediate variable
      term = pow (term, gaminv);
      double psix = pow ((term - 1), (1 / beta)) / alpha;
      if ( psix < 0.01 ) 
           psix = 0.01;
//      psix (in cm) is converted to bars (negative value).
      psix = ( 0.01 - psix ) * 0.001;
      if ( psix < -500000 ) 
          psix = - 500000;
      if ( psix > -0.00001 ) 
          psix = - 0.00001;
	  return psix;
}
////////////////////////////////////////////////////////////////////
double qpsi ( double psi, double qr, double qsat, double alpha, double beta )
//     This function computes soil water content (cm3 cm-3) for
//  a given value of matric potential, using the Van-Genuchten equation.
//
//     The following arguments are used:
//        alpha, beta  - parameters of the van-genuchten equation.
//        psi - soil water matric potential (bars).
//        qr - residual water content, cm3 cm-3.
//        qsat - saturated water content, cm3 cm-3.
//
{
//     For very high values of PSI, saturated water content is assumed.
//     For very low values of PSI, air-dry water content is assumed.
      if ( psi >= -0.00001 )     
          return qsat;
      else if ( psi <= -500000 ) 
          return qr;
//     The soil water matric potential is transformed from bars (psi)
//  to cm in positive value (psix).
      double psix = 1000 * fabs(psi + 0.00001);
//     The following equation is used (in FORTRAN notation):
//      QPSI = QR + (QSAT-QR) / (1 + (ALPHA*PSIX)**BETA)**(1-1/BETA)
      double gama = 1 - 1 / beta;
      double term = 1 + pow( (alpha * psix), beta);  //  intermediate variable
      double swfun = qr + (qsat - qr) / pow(term, gama);  //  computed water content
      if (swfun < (qr + 0.0001)) 
          swfun = qr + 0.0001;
      return swfun;
}////////////////////////////////////////////////////////////////////////
double wcond ( double q, double qr, double qsat, double beta, double SaturatedHydCond, double PoreSpace )
//     This function computes soil water hydraulic conductivity
//  for a given value of soil water content, using the Van-Genuchten
//  equation. The units of the computed conductivity are the same as the given
//  saturated conductivity (SaturatedHydCond).
//
//     The following arguments are used:
//        beta  - parameter of the van-genuchten equation.
//        SaturatedHydCond - saturated hydraulic conductivity (at qsat).
//        PoreSpace - pore space volume.
//        q - soil water content, cm3 cm-3.
//        qr - residual water content, cm3 cm-3.
//        qsat - saturated water content, cm3 cm-3.
//
{
//
//     For very low values of water content (near the residual water
//  content) wcond is 0.
      if ( (q-qr) < 0.0001 )
          return 0;
//     Water content for saturated conductivity is minimum of PoreSpace and qsat.
//     For very high values of water content (exceeding the saturated
//  water content or pore space) conductivity is SaturatedHydCond.
      double xsat = min( qsat, PoreSpace );
      if ( q >= xsat) 
          return SaturatedHydCond;
//      The following equation is used (in FORTRAN notation):
//      WCOND = CONDSAT * ((Q-QR)/(XSAT-QR))**0.5
//             * (1-(1-((Q-QR)/(XSAT-QR))**(1/GAMA))**GAMA)**2
      double gama   = 1 - 1 / beta;
      double gaminv = 1 / gama;
      double sweff  = (q - qr) / (xsat - qr);  // intermediate variable (effective water content).
      double acoeff = pow ( (1 - pow(sweff, gaminv)),  gama);  // intermediate variable
      double bcoeff = pow ( (1 - acoeff),  2);  // intermediate variable
      double conductivity = pow ( sweff, 0.5) * bcoeff * SaturatedHydCond;
	  return conductivity;
}
///////////////////////////////////////////////////////////////////////////////
double PsiOsmotic ( double q, double qsat, double ec)
//      This function computes soil water osmotic potential (in bars, positive value).
//
//     The following arguments are used:
//        q - soil water content, cm3 cm-3.
//        qsat - saturated water content, cm3 cm-3.
//        ec - electrical conductivity of saturated extract (mmho/cm)
//
{
    double ReturnValue;
	if (ec > 0) 
	{
	     ReturnValue = 0.36 * ec * qsat / q;
         if (ReturnValue > 6)
             ReturnValue = 6;
		 return ReturnValue;
	}
    else   
	     return 0;
}
///////////////////////////////////////////////////////////////////////////////
double GetFromClim(string item, int Doy)
//     This function extracts daily climate values for day of year Doy
//  from the structure Clim.
//     Input arguments:
//       item - string defining which item to extract.
//       Doy -  defines day of year to extract.
//
{
      int i;
	  for (i = 0; i < 400; i++)
	  {
		if (Clim[i].nDay == Doy)
	            		break;
	  }
//
	  if (i > 399)
         i = 399;
	  if (Doy < Clim[0].nDay)
         i = 0;
	  if (item == "tmin")
		  return Clim[i].Tmin;
	  else if (item == "tmax")
		  return Clim[i].Tmax;
	  else if (item == "rad")
		  return Clim[i].Rad;
	  else if (item == "rain")
		  return Clim[i].Rain;
	  else if (item == "wind")
		  return Clim[i].Wind;
	  else if (item == "tdew")
		  return Clim[i].Tdew;
	  else
		  return -99;
}
