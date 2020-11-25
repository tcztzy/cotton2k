//  File WriteOutput_2.cpp
//
//   functions in this file:
// cotplt()
// sitecode()
// bollsize()
//
#include "CottonSimulation.h"
#include "GeneralFunctions.h"
#include "resource.h"

// File scope variables:
double bsize[3][30][5]; // average boll size (g seeds + lint per boll)
int    mcode[3][30][5]; // code indicating state of fruit.
/////////////////////////////////////////////////////////////
void cotplt(int mtype, const string& ProfileName)
//     This function executes output of plant maps. It plots the cotton plant with the 
//  corresponding symbols at the correct locations. It is called from function DailyOutput().
//     It calls bollsize() and sitecode().
//
//     The following global and file scope variables are referenced here:
//       bsize, Date, FruitingCode, mcode, NumFruitBranches, NumNodes, NumVegBranches.
//     The following argument is referenced here:
//       mtype - defines the type of plant map requested:
//        (1) - Average plant map (mcode); function sitecode() is called.
//        (4) - Average boll size (bsize); function bollsize() is called.
{
//     The following constants are used:
	  CString char1[7] = { "-X", "-*", "-$", "-A", "-A", "-A", "-B" };
	  CString char2[7] = { "X-", "*-", "$-", "A-", "A-", "A-", "B-" };
	  CString char3[13] = { "-1", "-2", "-3", "-4", "-5", "-6", "-7", "-8",
                         "-9", "-0", "-A", "-$", "-X" };
	  CString char4[13] = { "1-", "2-", "3-", "4-", "5-", "6-", "7-", "8-",
                         "9-", "0-", "A-", "$-", "X-" };
	  char chari = 'I'; // symbols used in plots.
      ofstream File24("Output\\" + ProfileName + ".PLM", ios::app);
//     Print two blank lines to begin a new plot.
 	  File24 << endl << endl;
      if ( mtype == 1 ) 
	  {
         sitecode();  // compute mcode
//     Write heading and legends for the plant map.
	     File24 << "            Average Plant Fruit Code on " << Date << endl;
		 File24 << "  X = Square          B = Bloom or young green boll" << endl;
		 File24 << "  * = Set green boll  $ = mature boll" << endl;
		 File24 << "  A = Abscised fruit  I = Monopodial node" << endl << endl;
	  }
      else if ( mtype == 4 ) 
	  {
         bollsize(); // compute bsize
//     Write heading and legends for the plant map.
	     File24 << "            Boll Weight on Location on " << Date << endl;
		 File24 << "  A ==> Abscised fruit  X ==> Square" << endl;
		 File24 << "  0 ==> 0   to  0.5 g   1 ==> 0.5 to 1.5 g" << endl;
		 File24 << "  2 ==> 1.5 to  2.5 g   3 ==> 2.5 to 3.5 g" << endl;
		 File24 << "  4 ==> 3.5 to  4.5 g   5 ==> 4.5 to 5.5 g" << endl;
		 File24 << "  6 ==> 5.5 to  5.5 g   7 ==> 6.5 to 7.5 g" << endl;
		 File24 << "  8 ==> 7.5 to  8.5 g   9 ==> 8.5 or more" << endl << endl;
      }
//     Blank pri and prt.
      char pri[3][30]; // character defined as 'I' for existing stem nodes.
      CString prt[3][30][5]; // array of character symbols used in plot.
	  for (int k = 0; k < 3; k++)
		  for (int l = 0; l < 30; l++)
		  {
			  pri[k][l] = ' ';
			  for (int m = 0; m < 5; m++)
				  prt[k][l][m] = "  ";
		  }
//     Loop over all vegetative stems. Separate loops over fruiting branches and nodes 
//  for odd numbered fruiting branches (right side of the plotted plant) and even numbered 
//  fruiting branches (left side).
      int icode; // code used to define the character symbol to plot.
      int m0; // number of node on fruiting branch in reverse order.
      int nbrch; // the number of fruiting branches on a vegetative stem.
      int nnid; // the number of nodes on a fruiting branch.
	  for (int k = 0; k < NumVegBranches; k++)
	  {
         nbrch = NumFruitBranches[k];
//     Right side of the graph.
		 for (int l = 0; l < nbrch; l+=2)
		 {
            nnid  = NumNodes[k][l];
			for (int m = 0; m < nnid; m++)
			{
               if ( FruitingCode[k][l][m] > 0 ) 
			   {
//     Array prt for all existing sites is set with symbols for plotting. mcode is assigned 
//  symbol sets char1 and char2, whereas bsize is assigned symbol sets char3 and char4, for 
//  right and left sides, respectively.
//     'X' is a square, 'A' is an absciced site, '$' is an open boll. For plotting mcode, 
//  symbols 'B' and '*' are young and full sized green bolls, respectively. For boll size 
//  plots, symbols '0' to '9' are used for size groups.
//     Array pri for all existing monopodial stem nodes is assigned symbol chari ('I').
                  if ( mtype == 1 ) 
				  {
                     icode = mcode[k][l][m];
                     if ( icode > 0 )
					    prt[k][l][m] = char1[icode-1];
				  }
                  else if ( mtype == 4 ) 
				  {
                     icode = (int) (bsize[k][l][m] + 0.5);
	                 if (icode > 9)
					    icode = 9;
                     if ( bsize[k][l][m] < 0.5 )
					    icode  = 10;
                     if ( FruitingCode[k][l][m] == 4 || FruitingCode[k][l][m] == 5 )
                        icode = 11;
                     if ( FruitingCode[k][l][m] == 1 )
					    icode  = 13;
                     if ( icode > 0 )
					    prt[k][l][m] = char3[icode-1];
				  }
                  pri[k][l] = chari;
			   } // if FruitingCode
			} // for m
		 } // for l
//     Repeat for the left side of the graph.
		 for (int l = 1; l < nbrch; l+=2)
		 {
            nnid  = NumNodes[k][l];
			for (int m = 0; m < nnid; m++)
			{
               m0 = 4 - m;
               if ( FruitingCode[k][l][m] > 0 ) 
			   {
                  if ( mtype == 1 ) 
				  {
                     icode = mcode[k][l][m];
                     if ( icode > 0 )
					    prt[k][l][m0] = char2[icode-1];
				  }
                  else if ( mtype == 4 ) 
				  {
                     icode = (int) (bsize[k][l][m] + 0.5);
	                 if (icode > 9)
					    icode = 9;
                     if ( bsize[k][l][m] < 0.5 )
					    icode  = 10;
                     if ( FruitingCode[k][l][m] == 4 || FruitingCode[k][l][m] == 5 )
                        icode = 11;
                     if ( FruitingCode[k][l][m] == 1 )
					    icode  = 13;
                     if ( icode > 0 )
					    prt[k][l][m0] = char4[icode-1];
				  }
                  pri[k][l] = chari;
			   }  // if FruitingCode
			} // for m
		 } // for l
	  } // for k
//     Define nbrch as the number of fruiting branches on the main
//  stem. Plotting sequence depends on whether this is odd or even.
      int lx; // number of fruiting branch in reverse order.
      int lx1; // lx minus 1.
      nbrch = NumFruitBranches[0];
      if ( nbrch % 2 != 0 ) 
	  {
         for (int l = 0; l < nbrch; l+= 2)
		 {
            lx = nbrch - 1 - l;
            lx1 = lx - 1;
			File24 << "           " << pri[0][lx];
			for (int m = 0; m < 5; m++)
				File24 << prt[0][lx][m];
			File24 << "            " << pri[1][lx];
			for (int m = 0; m < 5; m++)
				File24 << prt[1][lx][m];
			File24 << "            " << pri[2][lx];
			for (int m = 0; m < 5; m++)
				File24 << prt[2][lx][m];
			File24 << endl;
            if ( lx1 >= 0 )
			{
				File24 << " ";
  			    for (int m = 0; m < 5; m++)
				   File24 << prt[0][lx1][m];
			    File24 << pri[0][lx1] << "            ";
  			    for (int m = 0; m < 5; m++)
				   File24 << prt[1][lx1][m];
			    File24 << pri[1][lx1] << "            " ;
  			    for (int m = 0; m < 5; m++)
				   File24 << prt[2][lx1][m];
			    File24 << pri[2][lx1] << endl;
			}
         }
	  }
      else
	  {
         for (int l = 0; l < nbrch; l+= 2)
		 {
            lx = nbrch - 1 - l;
            lx1 = lx - 1;
			File24 << " ";
  			for (int m = 0; m < 5; m++)
			   File24 << prt[0][lx][m];
			File24 << pri[0][lx] << "            ";
  			for (int m = 0; m < 5; m++)
			   File24 << prt[1][lx][m];
			File24 << pri[1][lx] << "            " ;
  			for (int m = 0; m < 5; m++)
			   File24 << prt[2][lx][m];
			File24 << pri[2][lx] << endl;
            if ( lx1 >= 0 )
			{
			   File24 << "           " << pri[0][lx1];
			   for (int m = 0; m < 5; m++)
			    	File24 << prt[0][lx1][m];
			   File24 << "            " << pri[1][lx1];
			   for (int m = 0; m < 5; m++)
			    	File24 << prt[1][lx1][m];
			   File24 << "            " << pri[2][lx1];
			   for (int m = 0; m < 5; m++)
			    	File24 << prt[2][lx1][m];
			   File24 << endl;
			}
		 }
	  }
//
      File24 << "           " << pri[0][0];
      File24 << "                      " << pri[1][0];
      File24 << "                      " << pri[2][0] << endl;
      File24 << "           " << pri[0][0];
      File24 << "                      " << pri[1][0];
      File24 << "                      " << pri[2][0] << endl;
      File24 << "           " << pri[0][0];
      File24 << "                      " << pri[1][0];
      File24 << "                      " << pri[2][0] << endl;
      File24 << endl << endl;
}
////////////////////////
void sitecode()
//     This function is called from function cotplt(). It sets mcode to draw the plant map. 
//  A fruit will be included in the average plant map if the fraction of fruit remaining at 
//  each site is equal to or larger then the average retention for each fruit class.
//
//     The following global variables are referenced here:
//  FruitingCode, FruitFraction, NumFruitBranches, NumNodes, NumVegBranches.
//     The following file scope variable is set here:  mcode.
//
//     The following values are used for mcode: 1 = square;
//  2 = full sized green boll; 3 = mature boll; 4 = fruit abscised as
//  a square; 5 = fruit abscised as a boll; 7 = young green boll.
//     The average retention is computed from the data for these groups.
{
//     Numbers (per plant) of the following:
      double sumfru1 = 0; // retained squares
	  double sumfru2 = 0; // old green bolls
	  double sumfru3 = 0; // open bolls
	  double sumfru7 = 0; // young green bolls
//     Numbers (per plant) of sites with the following:
      double fru1no = 0; // squares
	  double fru2no = 0; // old green bolls
	  double fru3no = 0; // open bolls
	  double fru7no = 0; // young green bolls
//     Loop over all possible fruiting sites. The number of retained squares or bolls, 
//  as well as the number of sites, are computed for each group of fruits (1, 2, 3, and 7, 
//  as classified for FruitingCode).
	  for (int k = 0; k < NumVegBranches; k++)
	  {
         int nbrch; // the number of fruiting branches on a vegetative stem.
         nbrch = NumFruitBranches[k];
		 for (int l = 0; l < nbrch; l++)
		 {
            int nnid;  // the number of nodes on a fruiting branch.
            nnid  = NumNodes[k][l];
			for (int m = 0; m < nnid; m++)
			{
               if ( FruitingCode[k][l][m] == 1 ) 
			   {
                  sumfru1 += FruitFraction[k][l][m];
                  fru1no++;
			   }
               else if ( FruitingCode[k][l][m] == 2 ) 
			   {
                  sumfru2 += FruitFraction[k][l][m];
                  fru2no++;
			   }
               else if ( FruitingCode[k][l][m] == 3) 
			   {
                  sumfru3 += FruitFraction[k][l][m];
                  fru3no++;
			   }
               else if ( FruitingCode[k][l][m] == 7 ) 
			   {
                  sumfru7 += FruitFraction[k][l][m];
                  fru7no++;
               } // if FruitingCode
			} // for m
         } // for l
      } // for k
//     Compute average retention per site for each group.
      double avgfru1 = 0, avgfru2 = 0, avgfru3 = 0, avgfru7 = 0; // average retention (per site)
//       of squares, old green bolls, open bolls, and young green bolls, respectively.
      if ( fru1no > 0 ) 
		  avgfru1 = sumfru1 / fru1no;
      if ( fru2no > 0 )
		  avgfru2 = sumfru2 / fru2no;
      if ( fru3no > 0 ) 
		  avgfru3 = sumfru3 / fru3no;
      if ( fru7no > 0 ) 
		  avgfru7 = sumfru7 / fru7no;
//     Loop over all possible fruiting sites, and check for existence of square or boll.
//  If it is a square, define mcode = 1 if retention is above the average, mcode = 4 if below.
//  If it is a full sized green boll, define mcode = 2 if  retention is above the average, 
//  mcode = 5 if below. If it is an open boll, define mcode = 3 if retention is above
//  the average, mcode = 5 if below. If it is a young green boll, define mcode = 7 if retention
//  is above the average, mcode = 5 if below. For abscised squares or bolls (fcode is 4 or 5) 
//  define mcode as 4 or 5, respectively.
	  for (int k = 0; k < NumVegBranches; k++)
	  {
         int nbrch = NumFruitBranches[k];
		 for (int l = 0; l < nbrch; l++)
		 {
            int nnid  = NumNodes[k][l];
			for (int m = 0; m < nnid; m++)
			{
               if ( FruitingCode[k][l][m] == 1 ) 
			   {
                  if ( FruitFraction[k][l][m] >= avgfru1 ) 
                     mcode[k][l][m] = 1;
                  else
                     mcode[k][l][m] = 4;
			   }
               else if ( FruitingCode[k][l][m] == 2 ) 
			   {
                  if ( FruitFraction[k][l][m] >= avgfru2 ) 
                     mcode[k][l][m] = 2;
                  else
                     mcode[k][l][m] = 5;
			   }
               else if ( FruitingCode[k][l][m] == 3 ) 
			   {
                  if ( FruitFraction[k][l][m] >= avgfru3 ) 
                     mcode[k][l][m] = 3;
                  else
                     mcode[k][l][m] = 5;
			   }
               else if ( FruitingCode[k][l][m] == 4 )
				   mcode[k][l][m] = 4;
               else if ( FruitingCode[k][l][m] == 5 )
				   mcode[k][l][m] = 5;
               else if ( FruitingCode[k][l][m] == 7) 
			   {
                  if ( FruitFraction[k][l][m] > avgfru7 ) 
                     mcode[k][l][m] = 7;
                  else
                     mcode[k][l][m] = 5;
			   }
			}
         }
      }
}
////////////////////////////
void bollsize()
//     This function is called from function cotplt().
//
//     The following global variables are referenced here:
//       BollWeight, FruitingCode, FruitFraction, NumFruitBranches, NumNodes, 
//       NumVegBranches, VarPar.
//     The following fie scope variable is set here:    bsize
{
//     Loop over all vegetative stems, fruiting branches and nodes:
//     Compute boll size, bsize, as the average weight of seed cotton per boll, for each site 
//  with a boll. bsize will be zero when no bolls are retained. It can not be larger than a 
//  maximum value defined by VarPar(11).
	  for (int k = 0; k < NumVegBranches; k++)
	  {
         int nbrch; // the number of fruiting branches on a vegetative stem.
         nbrch = NumFruitBranches[k];
		 for (int l = 0; l < nbrch; l++)
		 {
            int nnid; // the number of nodes on a fruiting branch.
            nnid  = NumNodes[k][l];
			for (int m = 0; m < nnid; m++)
			{
               if ( FruitingCode[k][l][m] == 2 || FruitingCode[k][l][m] == 3 || FruitingCode[k][l][m] == 7)  
			   {
                  if ( FruitFraction[k][l][m] > 0.001 ) 
                     bsize[k][l][m] = BollWeight[k][l][m] / FruitFraction[k][l][m];
                  else
                     bsize[k][l][m] = 0;
                  if ( bsize[k][l][m] > VarPar[11] )
                       bsize[k][l][m] = VarPar[11];
			   }
			} // for m
		 } // for l
      } // for k
}
