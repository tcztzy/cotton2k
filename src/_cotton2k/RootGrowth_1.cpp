//  RootGrowth_1.cpp
//
//   functions in this file:
// PotentialRootGrowth()
// RootImpedance()
// ComputeActualRootGrowth()
//
#include "global.h"
#include "exceptions.h"
#include "Simulation.hpp"

using namespace std;
//////////////////////////////////////////////////
//                   THE COTTON ROOT SUB-MODEL.
//     The following is a documentation of the root sub-model used in COTTON2K. It is
//  derived from the principles of RHIZOS, as implemented in GOSSYM and in GLYCIM,
//  and from some principles of ROOTSIMU (Hoogenboom and Huck, 1986). It is devised
//  to be generally applicable, and may be used with root systems of different crops by
//  redefining the parameters, which are set here as constants, and some of them are
//  set in function InitializeRootData(). These parameters are of course specific for the
//  crop species, and perhaps also for cultivars or cultivar groups.
//
//     This is a two-dimensional model and it may be used with soil cells of different
//  sizes. The grid can be defined by the modeler. The maximum numbers of layers
//  and columns are given by the parameters maxl and maxk, respectively. These are set
//  to 40 and 20, in this version of COTTON2K. The grid is set in function InitializeGrid().
//
//     The whole slab is being simulated. Thus, non-symmetrical processes (such as
//  side-dressing of fertilizers or drip-irrigation) can be handled. The plant is assumed to
//  be situated at the center of the soil slab, or off-center for skip-rows. Adjoining soil
//  slabs are considered as mirror-images of each other. Alternate-row drip systems (or any
//  other agricultural input similarly situated) are located at one edge of the slab.
//
//     The root mass in each cell is made up of NumRootAgeGroups classes, whose number is to be
//  defined by the modeler. The maximum number of classes is 3 in this version of COTTON2K.
//
//     The following functions account for root growth morphology: TapRootGrowth() describes
//  growth of the taproot, and LateralRootGrowth() describes growth of the lateral roots.
//
//     The calling sequence of the root submodel modules is as follows:
//     InitializeRootData() is called from ReadInput()
//  at the start of the simulation (see their code in file gettinginput_2.cpp) .
//     PotentialRootGrowth() and ActualRootGrowth() are called each day from PlantGrowth().
//     PotentialRootGrowth() calls RootImpedance(), SoilMechanicResistance(), SoilAirOnRootGrowth(),
//  SoilNitrateOnRootGrowth(), SoilTemOnRootGrowth().
//     ActualRootGrowth() calls RedistRootNewGrowth(), TapRootGrowth(), LateralRootGrowth(),
//  RootAging(), RootDeath(), RootCultivation(), RootSummation().
//////////////////////////
void RootImpedance(SoilCell soil_cells[40][20])
//     This function calculates soil mechanical impedance to root growth, rtimpd(l,k),
//  for all soil cells. It is called from PotentialRootGrowth(). The impedance is a function
//  of bulk density and water content in each soil soil cell. No changes have been made
//  in the original GOSSYM code.
//
//     The following global variables are referenced here:
//  BulkDensity, gh2oc, impede, inrim, ncurve, nk, nl,
//  SoilHorizonNum, tstbd.
//     The following global variables are set here:    RootImpede.
{
    for (int l = 0; l < nl; l++)
    {
        int j = SoilHorizonNum[l];
        double Bd = BulkDensity[j]; // bulk density for this layer
                                    //
        int jj;
        for (jj = 0; jj < inrim; jj++)
        {
            if (Bd <= tstbd[jj][0])
                break;
        }
        int j1 = jj;
        if (j1 > inrim - 1)
            j1 = inrim - 1;
        int j0 = jj - 1;
        if (j0 < 0)
            j0 = 0;
        //
        for (int k = 0; k < nk; k++)
        {
            double Vh2o = soil_cells[l][k].water_content / Bd;
            int ik;
            for (ik = 0; ik < ncurve; ik++)
            {
                if (Vh2o <= gh2oc[ik])
                    break;
            }
            int i1 = ik;
            if (i1 > ncurve - 1)
                i1 = ncurve - 1;
            int i0 = ik - 1;
            if (i0 < 0)
                i0 = 0;
            //
            if (j1 == 0)
            {
                if (i1 == 0 || Vh2o <= gh2oc[i1])
                    RootImpede[l][k] = impede[j1][i1];
                else
                    RootImpede[l][k] = impede[j1][i0] - (impede[j1][i0] - impede[j1][i1]) * (Vh2o - gh2oc[i0]) / (gh2oc[i1] - gh2oc[i0]);
            }
            else
            {
                if (i1 == 0 || Vh2o <= gh2oc[i1])
                    RootImpede[l][k] = impede[j0][i1] - (impede[j0][i1] - impede[j1][i1]) * (tstbd[j0][i1] - Bd) / (tstbd[j0][i1] - tstbd[j1][i1]);
                else
                {
                    double temp1 = impede[j0][i1] - (impede[j0][i1] - impede[j1][i1]) * (tstbd[j0][i1] - Bd) / (tstbd[j0][i1] - tstbd[j1][i1]);
                    double temp2 = impede[j0][i0] - (impede[j0][i0] - impede[j1][i1]) * (tstbd[j0][i0] - Bd) / (tstbd[j0][i0] - tstbd[j1][i0]);
                    RootImpede[l][k] = temp2 + (temp1 - temp2) * (Vh2o - gh2oc[i0]) / (gh2oc[i1] - gh2oc[i0]);
                }
            }
        }
    }
    //
}
