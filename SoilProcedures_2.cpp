// File SoilProcedures_2.cpp
//
//   functions in this file:
// CapillaryFlow() 
// Drain()
// DripFlow()    
// CellDistance()
// 
#include "global.h"
#include "GeneralFunctions.h"

double Drain();
double CellDistance(int, int, int, int);
// SoilProcedures_3
void WaterFlux(double[], double[], double[], double[], double[], double[], int, int, int, long);
void NitrogenFlow(int, double[], double[], double[], double[], double[]);

//////////////////////////
void CapillaryFlow(const int& DayStart)
//     This function computes the capillary water flow between soil cells. It is called by 
//  SoilProcedures(), noitr times per day.  The number of iterations (noitr) has been computed 
//  in SoilProcedures() as a function of the amount of water applied. It is executed only once 
//  per day if no water is applied by rain or irrigation.
//     It calls functions:   Drain(), NitrogenFlow(), psiq(), PsiOsmotic(), WaterFlux().
//
//     The following global variables are referenced:
//       alpha, ,vanGenuchtenBeta, Daynum, DayStart, dl, ElCondSatSoilToday, nk, nl, PoreSpace,
//       RowSpace, SoilHorizonNum, thad ,thts, WaterTableLayer, wk.
//     The following global variables are set:
//       CumWaterDrained, SoilPsi, VolNo3NContent, VolUreaNContent, VolWaterContent.
{
      static long numiter; // counter used for WaterFlux() calls.
      double wk1[40];      // dummy array for passing values of array wk.
//     Set initial values in first day.
      if (Daynum <= DayStart) 
	  {
         numiter = 0;
         for (int l = 0; l < nl; l++)
                wk1[l] = 0;
      }
//     Increase the counter numiter, and compute the updated values of SoilPsi in each
//  soil cell by calling functions psiq() and PsiOsmotic().
      numiter++;
      for (int l = 0; l < nl; l++)
	  {
         int j = SoilHorizonNum[l];  //  the soil horizon number
         for (int k = 0; k < nk; k++)
            SoilPsi[l][k] = psiq(VolWaterContent[l][k],thad[l],thts[l],alpha[j],vanGenuchtenBeta[j])
                          - PsiOsmotic ( VolWaterContent[l][k], thts[l], ElCondSatSoilToday);
      }
//    
      int nlx = nl; // The last layer without a water table.
      if (WaterTableLayer < nlx)
          nlx = WaterTableLayer - 1;
      int iv; //  direction indicator: iv = 1 for vertical flow in each column;
                                 //    iv = 0 for horizontal flow in each layer.
      double q01[40]; // one dimensional array of a layer or a column of previous
                      // values of VolWaterContent.
      double q1[40];  // one dimensional array of a layer or a column of VolWaterContent.
      double psi1[40];// one dimensional array of a layer or a column of SoilPsi.
      double nit[40]; // one dimensional array of a layer or a column of VolNo3NContent.
      double nur[40]; // one dimensional array of a layer or a column of VolUreaNContent.
//
//     VERTICAL FLOW in each column. the direction indicator iv is set to 1.
      iv = 1;  
//     Loop over all columns. Temporary one-dimensional arrays are defined for each column: 
//  assign the VolWaterContent[] values to temporary one-dimensional arrays q1 and q01. Assign 
//  SoilPsi, VolNo3NContent and VolUreaNContent values to arrays psi1, nit and nur, respectively.
      for (int k = 0; k < nk; k++)
	  {
         for (int l = 0; l < nlx; l++)
		 {
            q1[l] = VolWaterContent[l][k];
            q01[l] = VolWaterContent[l][k];
            psi1[l] = SoilPsi[l][k] + PsiOsmotic ( VolWaterContent[l][k], thts[l], ElCondSatSoilToday);
            nit[l] = VolNo3NContent[l][k];
            nur[l] = VolUreaNContent[l][k];
         } // end loop l
//     Call the following functions: WaterFlux() calculates the water flow caused by potential 
//  gradients; NitrogenFlow() computes the movement of nitrates caused by the flow of water. 
         WaterFlux( q1, psi1, dl, thad, thts, PoreSpace, nlx, iv, 0, numiter );
         NitrogenFlow(nl, q01, q1, dl, nit, nur);
//     Reassign the updated values of q1, nit, nur and psi1 back to 
//  VolWaterContent, VolNo3NContent, VolUreaNContent and SoilPsi.
         for (int l = 0; l < nlx; l++)
		 {
            VolWaterContent[l][k] = q1[l];
            VolNo3NContent[l][k] = nit[l];
            VolUreaNContent[l][k] = nur[l];
            SoilPsi[l][k] = psi1[l] - PsiOsmotic ( VolWaterContent[l][k], thts[l], ElCondSatSoilToday);
         } // end loop l
      } // end loop k
      double pp1[40]; // one dimensional array of a layer or a column of PP.
      double qr1[40]; // one dimensional array of a layer or a column of THAD.
      double qs1[40]; // one dimensional array of a layer or a column of THTS.
//
//     HORIZONTAL FLUX in each layer. The direction indicator iv is set to 0.
      iv = 0;
//     Loop over all layers. Define the horizon number j for this layer. Temporary 
//  one-dimensional arrays are defined for each layer: assign the VolWaterContent values 
//  to  q1 and q01. Assign SoilPsi, VolNo3NContent, VolUreaNContent, thad and thts values 
//  of the soil cells to arrays psi1, nit, nur, qr1 and qs1, respectively.
      for (int l = 0; l < nlx; l++)
	  {
         for (int k = 0; k < nk; k++)
		 {
            q1[k] = VolWaterContent[l][k];
            q01[k] = VolWaterContent[l][k];
            psi1[k] = SoilPsi[l][k] + PsiOsmotic ( VolWaterContent[l][k], thts[l], ElCondSatSoilToday);
            qr1[k] = thad[l];
            qs1[k] = thts[l];
            pp1[k] = PoreSpace[l];
            nit[k] = VolNo3NContent[l][k];
            nur[k] = VolUreaNContent[l][k];
            wk1[k] = wk[k];
         }
//     Call subroutines WaterFlux(), and NitrogenFlow() to compute water nitrate and 
//  urea transport in the layer. 
         WaterFlux( q1, psi1, wk1, qr1, qs1, pp1, nk, iv, l, numiter );
         NitrogenFlow(nk, q01, q1, wk1, nit, nur);
//     Reassign the updated values of q1, nit, nur and psi1 back to 
//  VolWaterContent, VolNo3NContent, VolUreaNContent and SoilPsi.
         for (int k = 0; k < nk; k++)
		 {
            VolWaterContent[l][k] = q1[k];
            SoilPsi[l][k] = psi1[k] - PsiOsmotic ( VolWaterContent[l][k], thts[l], ElCondSatSoilToday);
            VolNo3NContent[l][k] = nit[k];
            VolUreaNContent[l][k] = nur[k];
         } // snd loop k
      } // end loop l
//     Call Drain() to move excess water down in the column and compute drainage out
//  of the column. Update cumulative drainage.
      double WaterDrainedOut  = 0;    // water drained out of the slab, mm.
      WaterDrainedOut += Drain();
      if ( WaterDrainedOut > 0 ) 
         CumWaterDrained += 10 * WaterDrainedOut / RowSpace;
//  Compute the soil water potential for all soil cells.
      for (int l = 0; l < nl; l++)
	  {
         int j = SoilHorizonNum[l];
         for (int k = 0; k < nk; k++)
		 {
            SoilPsi[l][k] = psiq(VolWaterContent[l][k],thad[l],thts[l],alpha[j],vanGenuchtenBeta[j])
                          - PsiOsmotic ( VolWaterContent[l][k], thts[l], ElCondSatSoilToday);
         }
      }
}
////////////////////////////////////////////////////////////////////////////////////
double Drain()
//     This function computes the gravity flow of water in the slab, and returns the 
//  drainage of water out of the slab. It is called from GravityFlow() and CapillaryFlow().
//
//     The following global variables are referenced:
//       dl, FieldCapacity, MaxWaterCapacity, nk, nl, NO3FlowFraction, PoreSpace, RowSpace, 
//       WaterTableLayer, wk
//     The following global variables are set:
//       SoilNitrogenLoss, VolNo3NContent, VolUreaNContent, VolWaterContent, 
{
     int nlx = nl; // last soil layer for computing drainage.
     if (WaterTableLayer < nlx) 
		 nlx = WaterTableLayer;
     double oldvh2oc[20]; // stores previous values of VolWaterContent.
     double nitconc; // nitrate N concentration in the soil solution.
     double nurconc; // urea N concentration in the soil solution.
//     The following is executed if this is not the bottom layer.
     for (int l = 0; l < nlx-1; l++)
	 {
//     Compute the average water content (avwl) of layer l. Store the
//  water content in array oldvh2oc.
         double avwl = 0;     // average water content in a soil layer
         for (int k = 0; k < nk; k++)
		 {
            avwl += VolWaterContent[l][k] * wk[k] / RowSpace;
            oldvh2oc[k] = VolWaterContent[l][k];
         }
//     Upper limit of water content in free drainage..
         double uplimit = MaxWaterCapacity[l]; 
//
//     Check if the average water content exceeds uplimit for this layer, and if it does, 
//  compute amount (wmov) to be moved to the next layer from each cell.
         double wmov; // amount of water moving out of a cell.
         if ( avwl > uplimit ) 
		 {
             wmov = avwl - uplimit;
             wmov = wmov * dl[l] / dl[l+1];
             for (int k = 0; k < nk; k++)
			 {
//     Water content of all soil cells in this layer will be uplimit. the amount (qmv) 
//  to be added to each cell of the next layer is computed (corrected for non uniform 
//  column widths). The water content in the next layer is computed.
                 VolWaterContent[l][k] = uplimit;
                 VolWaterContent[l+1][k] += wmov * wk[k] * nk / RowSpace;
//     The concentrations of nitrate and urea N in the soil solution are
//  computed and their amounts in this layer and in the next one are updated.
                 double qvout; // amount of water moving out of a cell.
                 qvout = (oldvh2oc[k] - uplimit);
                 if ( qvout > 0 ) 
			     {
                     nitconc = VolNo3NContent[l][k] / oldvh2oc[k];
	                 if (nitconc < 1.e-30)
						 nitconc = 0;
                     nurconc = VolUreaNContent[l][k] / oldvh2oc[k];
	                 if (nurconc < 1.e-30)
						 nurconc = 0;
                     VolNo3NContent[l][k] = VolWaterContent[l][k] * nitconc;
                     VolUreaNContent[l][k] = VolWaterContent[l][k] * nurconc;
//     Only a part ( NO3FlowFraction ) of N is moved with water draining.
                     double vno3mov; // amount of nitrate N moving out of a cell.
                     double vnurmov; // amount of urea N moving out of a cell.
                     vno3mov = qvout * nitconc;
                     VolNo3NContent[l+1][k] += NO3FlowFraction[l] * vno3mov  * dl[l] / dl[l+1];
                     VolNo3NContent[l][k] += (1-NO3FlowFraction[l]) * vno3mov;
                     vnurmov = qvout * nurconc;
                     VolUreaNContent[l+1][k] += NO3FlowFraction[l] * vnurmov * dl[l] / dl[l+1];
                     VolUreaNContent[l][k] += (1 - NO3FlowFraction[l]) * vnurmov;
                 }
             }
         }
//     If the average water content is not higher than uplimit, start another loop over columns.  
		 else
		 {
             for (int k = 0; k < nk; k++)
			 {
//  Check each soil cell if water content exceeds uplimit,
                 if ( VolWaterContent[l][k] > uplimit ) 
			     {
                    wmov = VolWaterContent[l][k] - uplimit;
                    VolWaterContent[l][k] = uplimit;
                    VolWaterContent[l+1][k] += wmov * dl[l] / dl[l+1];
                    nitconc = VolNo3NContent[l][k] / oldvh2oc[k];
	                if (nitconc < 1.e-30) 
					    nitconc = 0;
                    nurconc = VolUreaNContent[l][k] / oldvh2oc[k];
	                if (nurconc < 1.e-30)
				        nurconc = 0;
                    VolNo3NContent[l][k] = VolWaterContent[l][k] * nitconc;
                    VolUreaNContent[l][k] = VolWaterContent[l][k] * nurconc;
//
                    VolNo3NContent[l+1][k] += NO3FlowFraction[l] * wmov * nitconc * dl[l] / dl[l+1];
                    VolUreaNContent[l+1][k] += NO3FlowFraction[l] * wmov * nurconc * dl[l] / dl[l+1];
                    VolNo3NContent[l][k] += (1- NO3FlowFraction[l]) * wmov * nitconc;
                    VolUreaNContent[l][k] += (1- NO3FlowFraction[l]) * wmov * nurconc;
                 } // end if Vol...
			 }// end loop k
         } // end if avwl...
	 } // end loop l
//     For the lowermost soil layer, loop over columns:
//     It is assumed that the maximum amount of water held at the lowest soil layer (nlx-1)
//  of the slab is equal to FieldCapacity. If water content exceeds MaxWaterCapacity, compute
//  the water drained out (Drainage), update water, nitrate and urea, compute nitrogen lost
//  by drainage, and add it to the cumulative N loss SoilNitrogenLoss.
	 double Drainage; // drainage of water out of the slab, cm3 (return value)
	 Drainage = 0; 
     for (int k = 0; k < nk; k++)
     {		
         if ( VolWaterContent[nlx-1][k] > MaxWaterCapacity[nlx-1] ) 
		 {
             Drainage += (VolWaterContent[nlx-1][k] - MaxWaterCapacity[nlx-1]) * dl[nlx-1] * wk[k];
             nitconc = VolNo3NContent[nlx-1][k] / oldvh2oc[k];
	         if (nitconc < 1.e-30)
			     nitconc = 0;
             nurconc = VolUreaNContent[nlx-1][k] / oldvh2oc[k];
	         if (nurconc < 1.e-30)
			     nurconc = 0;
             double saven; //  intermediate variable for computing N loss.
             saven = (VolNo3NContent[nlx-1][k] + VolUreaNContent[nlx-1][k]) * dl[nlx-1] * wk[k];
             VolWaterContent[nlx-1][k] = MaxWaterCapacity[nlx-1];
             VolNo3NContent[nlx-1][k] = nitconc * MaxWaterCapacity[nlx-1];
             VolUreaNContent[nlx-1][k] = nurconc * MaxWaterCapacity[nlx-1];
             SoilNitrogenLoss += saven - (VolNo3NContent[nlx-1][k] + VolUreaNContent[nlx-1][k]) * dl[nlx-1] * wk[k];
         } // end if Vol...
     } // end loop k
	 return Drainage;  
}
/////////////////////////////////////////////////////////////////////////////////////
void DripFlow(double Drip)
//     This function computes the water redistribution in the soil after irrigation 
//  by a drip system. It also computes the resulting redistribution of nitrate and urea N. 
//  It is called by SoilProcedures() noitr times per day. It calls function CellDistrib().
//     The following argument is used:
//  Drip - amount of irrigation applied by the drip method, mm.
//
//     The following global variables are referenced:
//       dl, LocationColumnDrip, LocationLayerDrip, MaxWaterCapacity, 
//       nk, nl, NO3FlowFraction, PoreSpace, RowSpace, wk
//
//     The following global variables are set:
//       CumWaterDrained, SoilNitrogenLoss, VolWaterContent, VolNo3NContent, VolUreaNContent
//
{
	  double dripw[40]; // amount of water applied, or going from one ring of
                        // soil cells to the next one, cm3. (array)
	  double dripn[40]; // amount of nitrate N applied, or going from one ring of soil
                        // soil cells to the next one, mg. (array)
	  double dripu[40]; // amount of urea N applied, or going from one ring of soil
                        // soil cells to the next one, mg. (array)
	  for (int i = 0; i < 40; i++)
	  {
		  dripw[i] = 0;
		  dripn[i] = 0;
		  dripu[i] = 0;
	  }
//     Incoming flow of water (Drip, in mm) is converted to dripw(0), in cm3 per slab.
      dripw[0] = Drip * RowSpace * .10;
//     Wetting the cell in which the emitter is located.
      double h2odef; // the difference between the maximum water capacity (at a water content
                     // of uplimit) of this ring of soil cell, and the actual water content, cm3.
	  int l0 = LocationLayerDrip;  //  layer where the drip emitter is situated
      int k0 = LocationColumnDrip; //  column where the drip emitter is situated
//     It is assumed that wetting cannot exceed MaxWaterCapacity of this cell. Compute
//  h2odef, the amount of water needed to saturate this cell.
      h2odef = ( MaxWaterCapacity[l0] - VolWaterContent[l0][k0] ) * dl[l0] * wk[k0];
//      If maximum water capacity is not exceeded - update VolWaterContent of
//  this cell and exit the function.
      if (dripw[0] <= h2odef) 
	  {
         VolWaterContent[l0][k0] += dripw[0] / ( dl[l0] * wk[k0] );
         return;
      }
//      If maximum water capacity is exceeded - calculate the excess of water flowing out of 
//  this cell (in cm3 per slab) as dripw[1]. The next ring of cells (kr=1) will receive it
//  as incoming water flow.
      dripw[1] = dripw[0] - h2odef;
//      Compute the movement of nitrate N to the next ring
      double cnw = 0;  //  concentration of nitrate N in the outflowing water
      if ( VolNo3NContent[l0][k0] > 1.e-30) 
	  {
         cnw = VolNo3NContent[l0][k0] / ( VolWaterContent[l0][k0] + dripw[0] / ( dl[l0] * wk[k0] ) );
//     cnw is multiplied by dripw[1] to get dripn[1], the amount of nitrate N going out
//  to the next ring of cells. It is assumed, however, that not more than a proportion  
//  (NO3FlowFraction) of the nitrate N in this cell can be removed in one iteration.
         if ( ( cnw * MaxWaterCapacity[l0] ) < ( NO3FlowFraction[l0] * VolNo3NContent[l0][k0] ) ) 
		 {
            dripn[1] = NO3FlowFraction[l0] * VolNo3NContent[l0][k0] * dl[l0] * wk[k0];
            VolNo3NContent[l0][k0] = (1 - NO3FlowFraction[l0]) * VolNo3NContent[l0][k0];
         }
		 else
         {
			dripn[1] = dripw[1] * cnw;
            VolNo3NContent[l0][k0] = MaxWaterCapacity[l0] * cnw;
         }
      }
//     The movement of urea N to the next ring is computed similarly.
	  double cuw = 0;  //  concentration of urea N in the outflowing water
      if ( VolUreaNContent[l0][k0] > 1.e-30) 
	  {
         cuw = VolUreaNContent[l0][k0]  / ( VolWaterContent[l0][k0] + dripw[0] / ( dl[l0] * wk[k0] ) );
         if ( ( cuw*MaxWaterCapacity[l0] ) < ( NO3FlowFraction[l0] * VolUreaNContent[l0][k0] ) )
		 {
            dripu[1] = NO3FlowFraction[l0] * VolUreaNContent[l0][k0] * dl[l0] * wk[k0];
            VolUreaNContent[l0][k0] = (1 - NO3FlowFraction[l0]) * VolUreaNContent[l0][k0];
         }
		 else
		 {
			dripu[1] = dripw[1] * cuw;
            VolUreaNContent[l0][k0] = MaxWaterCapacity[l0] * cuw;
		 }
      }
      double defcit[40][20]; // array of the difference between water capacity and
                             // actual water content in each cell of the ring
//     Set VolWaterContent of the cell in which the drip is located to MaxWaterCapacity.
      VolWaterContent[l0][k0] = MaxWaterCapacity[l0];
//     Loop of concentric rings of cells, starting from ring 1.
//     Assign zero to the sums sv, st, sn, sn1, su and su1.
      for (int kr = 1; kr < maxl; kr++)
	  {
         double uplimit; //  upper limit of soil water content in a soil cell
         double sv = 0; // sum of actual water content in a ring of cells, cm3
         double st = 0; // sum of total water capacity in a ring of cells, cm3
         double sn = 0; // sum of nitrate N content in a ring of cells, mg.
         double sn1 = 0;// sum of movable nitrate N content in a ring of cells, mg
         double su = 0; // sum of urea N content in a ring of cells, mg
         double su1 = 0;// sum of movable urea N content in a ring of cells, mg
	     double radius = 6 * kr;// radius (cm) of the wetting ring
         double dist; // distance (cm) of a cell center from drip location
//     Loop over all soil cells
         for(int l = 1; l < nl; l++)
		 {
//     Upper limit of water content is the porespace volume in layers below the water table, 
//  MaxWaterCapacity in other layers.
            if (l >= WaterTableLayer) 
                   uplimit = PoreSpace[l];
            else
                   uplimit = MaxWaterCapacity[l];
            for (int k = 0; k < nk; k++)
			{
//     Compute the sums sv, st, sn, sn1, su and su1 within the radius limits of this ring. The 
//  array defcit is the sum of difference between uplimit and VolWaterContent of each cell.
               dist = CellDistance(l, k, l0, k0); 
               if (dist <= radius && dist > (radius - 6) )
			   {
                  sv += VolWaterContent[l][k] * dl[l] * wk[k];
                  st += uplimit * dl[l] * wk[k];
                  sn += VolNo3NContent[l][k] * dl[l] * wk[k];
                  sn1 += VolNo3NContent[l][k] * dl[l] * wk[k] * NO3FlowFraction[l];
                  su += VolUreaNContent[l][k] * dl[l] * wk[k];
                  su1 +=  VolUreaNContent[l][k] * dl[l] * wk[k] * NO3FlowFraction[l];
                  defcit[l][k] = uplimit - VolWaterContent[l][k];
			   }
               else
                  defcit[l][k] = 0;
			} // end loop k
		 } // end loop l
//     Compute the amount of water needed to saturate all the cells in this ring (h2odef). 
         h2odef = st - sv;
//     Test if the amount of incoming flow, dripw(kr), is greater than  h2odef.
         if ( dripw[kr] <= h2odef ) 
		 {
//     In this case, this will be the last wetted ring.
//     Update VolWaterContent in this ring, by wetting each cell in proportion
//  to its defcit. Update VolNo3NContent and VolUreaNContent of the cells in this ring
//  by the same proportion. this is executed for all the cells in the ring.
            for (int l = 1; l < nl; l++)
			{
               for (int k = 0; k < nk; k++)
			   {
                  dist = CellDistance(l, k, l0, k0);
                  if (dist <= radius && dist > (radius - 6)) 
				  {
                     VolWaterContent[l][k] += dripw[kr] * defcit[l][k] / h2odef;
                     VolNo3NContent[l][k] += dripn[kr] * defcit[l][k] / h2odef;
                     VolUreaNContent[l][k] += dripu[kr] * defcit[l][k] / h2odef;
                  }
               } // end loop k
			} // end loop l
            return;
         } // end if dripw
//     If dripw(kr) is greater than h2odef, calculate cnw and cuw as the concentration of nitrate
//  and urea N in the total water of this ring after receiving the incoming water and nitrogen. 
         cnw = (sn + dripn[kr]) / (sv + dripw[kr]);
         cuw = (su + dripu[kr]) / (sv + dripw[kr]);
         double drwout = dripw[kr] - h2odef;  //  the amount of water going out of a ring, cm3.
//     Compute the nitrate and urea N going out of this ring, and their amount lost from this 
//  ring. It is assumed that not more than a certain part of the total nitrate or urea N 
//  (previously computed as sn1 an su1) can be lost from a ring in one iteration. drnout and 
//  xnloss are adjusted accordingly. druout and xuloss are computed similarly for urea N.
         double drnout = drwout * cnw;  //  the amount of nitrate N going out of a ring, mg
         double xnloss = 0; // the amount of nitrate N lost from a ring, mg
         if ( drnout > dripn[kr] ) 
		 {
            xnloss = drnout - dripn[kr];
            if ( xnloss > sn1 ) 
			{
               xnloss = sn1;
               drnout = dripn[kr] + xnloss;
            }
         }
         double druout = drwout * cuw;  //  the amount of urea N going out of a ring, mg
         double xuloss = 0;             // the amount of urea N lost from a ring, mg
         if ( druout > dripu[kr] ) 
		 {
            xuloss = druout - dripu[kr];
            if ( xuloss > su1 ) 
			{
               xuloss = su1;
               druout = dripu[kr] + xuloss;
            }
         }
//     For all the cells in the ring, as in the 1st cell, saturate VolWaterContent to uplimit, 
//  and update VolNo3NContent and VolUreaNContent.
         for (int l = 1; l < nl; l++)
		 {
            if (l >= WaterTableLayer) 
                  uplimit = PoreSpace[l];
            else
                  uplimit = MaxWaterCapacity[l];
//
            for (int k = 0; k < nk; k++)
			{
               dist = CellDistance(l, k, l0, k0);
               if (dist <= radius && dist > (radius - 6)) 
			   {
                  VolWaterContent[l][k] = uplimit;
                  if ( xnloss <= 0 ) 
                     VolNo3NContent[l][k] = uplimit * cnw;
                  else
                     VolNo3NContent[l][k] = VolNo3NContent[l][k] * (1.- xnloss / sn);
                  if ( xuloss <= 0 ) 
                     VolUreaNContent[l][k] = uplimit * cuw;
                  else
                     VolUreaNContent[l][k] = VolUreaNContent[l][k] * (1.- xuloss / su);
               } // end if dist
            } // end loop k
         } // end loop l
//     The outflow of water, nitrate and urea from this ring will be the inflow into the next ring. 
         if ( kr < (nl-l0-1) && kr < maxl-1) 
		 {
            dripw[kr+1] = drwout;
            dripn[kr+1] = drnout;
            dripu[kr+1] = druout;
		 }
         else
//     If this is the last ring, the outflowing water will be added to the drainage, 
//  CumWaterDrained, the outflowing nitrogen to SoilNitrogenLoss.
		 {
            CumWaterDrained += 10 * drwout / RowSpace;
            SoilNitrogenLoss += drnout + druout;
            return;
		 } // end if kr...
//     Repeat all these procedures for the next ring.
	  } // end loop kr
}
/////////////////////////
double CellDistance(int l, int k, int l0, int k0)
//     This function computes the distance between the centers of cells l,k an l0,k0
//  It is called from DripFlow().
{
//     Compute vertical distance between centers of l and l0
      double xl = 0;  // vertical distance (cm) between cells 
      if (l > l0) 
	  {
         for (int il = l0; il <= l; il++)
               xl += dl[il];
         xl -= (dl[l] + dl[l0]) * 0.5;
	  }
      else if (l < l0) 
	  {
         for (int il = l0; il >= l; il--)
               xl += dl[il];
         xl -= (dl[l] + dl[l0]) * 0.5;
      }
//     Compute horizontal distance between centers of k and k0
      double xk = 0;  // horizontal distance (cm) between cells
      if (k > k0) 
	  {
         for (int ik = k0; ik <= k; ik++)
               xk += wk[ik];
         xk -= (wk[k] + wk[k0]) * 0.5;
      }
	  else if (k < k0) 
	  {
         for (int ik = k0; ik >= k; ik--)
               xk += wk[ik];
         xk -= (wk[k] + wk[k0]) * 0.5;
      }
//     Compute diagonal distance between centers of cells
      return  sqrt(xl*xl + xk*xk);
}
