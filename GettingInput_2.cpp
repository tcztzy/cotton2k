//  File GettingInput_2.cpp
//     Consists of the following functions:
// InitializeSoilData()
// ReadSoilHydraulicData()
// InitializeRootData()
// InitializeSoilTemperature()
// form()
//
#include <boost/algorithm/string.hpp>
#include <iostream>

#include "CottonSimulation.h"
#include "GeneralFunctions.h"

//
#include <math.h>
#ifdef _DEBUG
#define new DEBUG_NEW
#endif
//////////////////////////////////////////////////////////
void InitializeRootData()
//     This function initializes the root submodel parameters and variables. It
//     is called
//  by ReadInput(). it is executed once at the beginning of the simulation.
//
//     Global or file scope variables referenced:
//        dl, PlantRowColumn, nk, nl, PerPlantArea, RowSpace.
//     Global or file scope variables set:
//        ActualRootGrowth[maxl][maxk], cgind[3], DepthLastRootLayer,
//	      LastTaprootLayer, LateralRootFlag[maxl], NumLayersWithRoots,
// NumRootAgeGroups,
//        PotGroRoots[maxl][maxk], RootAge[maxl][maxk], RootColNumLeft[maxl],
//        RootColNumRight[maxl], RootGroFactor[maxl][maxk],
//        RootWeight[maxl][maxk][3], TapRootLength, TotalRootWeight.
//
{
    //     The parameters of the root model are defined for each root class:
    //       grind(i), cuind(i), thtrn(i), trn(i), thdth(i), dth(i).
    NumRootAgeGroups = 3;
    cgind[0] = 1;
    cgind[1] = 1;
    cgind[2] = 0.10;
    double rlint = 10;  // Vertical interval, in cm, along the taproot, for
                        // initiating lateral roots.
    int ll = 1;         // Counter for layers with lateral roots.
    double sumdl =
        0;  // Distance from soil surface to the middle of a soil layer.
    for (int l = 0; l < nl; l++) {
        //     Using the value of rlint (interval between lateral roots), the
        //  layers from which lateral roots may be initiated are now computed.
        //  LateralRootFlag[l] is assigned a value of 1 for these layers.
        LateralRootFlag[l] = 0;
        if (l > 0) sumdl += 0.5 * dl[l - 1];
        sumdl += 0.5 * dl[l];
        if (sumdl >= ll * rlint) {
            LateralRootFlag[l] = 1;
            ll++;
        }
    }
    //     All the state variables of the root system are initialized to zero.
    for (int l = 0; l < nl; l++) {
        if (l < 3) {
            RootColNumLeft[l] = PlantRowColumn - 1;
            RootColNumRight[l] = PlantRowColumn + 2;
        } else if (l < 7) {
            RootColNumLeft[l] = PlantRowColumn;
            RootColNumRight[l] = PlantRowColumn + 1;
        } else {
            RootColNumLeft[l] = 0;
            RootColNumRight[l] = 0;
        }
        //
        for (int k = 0; k < nk; k++) {
            PotGroRoots[l][k] = 0;
            RootGroFactor[l][k] = 1;
            ActualRootGrowth[l][k] = 0;
            RootAge[l][k] = 0;
            for (int i = 0; i < 3; i++) RootWeight[l][k][i] = 0;
        }
    }
    //
    RootWeight[0][PlantRowColumn - 1][0] = 0.0020;
    RootWeight[0][PlantRowColumn][0] = 0.0070;
    RootWeight[0][PlantRowColumn + 1][0] = 0.0070;
    RootWeight[0][PlantRowColumn + 2][0] = 0.0020;
    RootWeight[1][PlantRowColumn - 1][0] = 0.0040;
    RootWeight[1][PlantRowColumn][0] = 0.0140;
    RootWeight[1][PlantRowColumn + 1][0] = 0.0140;
    RootWeight[1][PlantRowColumn + 2][0] = 0.0040;
    RootWeight[2][PlantRowColumn - 1][0] = 0.0060;
    RootWeight[2][PlantRowColumn][0] = 0.0210;
    RootWeight[2][PlantRowColumn + 1][0] = 0.0210;
    RootWeight[2][PlantRowColumn + 2][0] = 0.0060;
    RootWeight[3][PlantRowColumn][0] = 0.0200;
    RootWeight[3][PlantRowColumn + 1][0] = 0.0200;
    RootWeight[4][PlantRowColumn][0] = 0.0150;
    RootWeight[4][PlantRowColumn + 1][0] = 0.0150;
    RootWeight[5][PlantRowColumn][0] = 0.0100;
    RootWeight[5][PlantRowColumn + 1][0] = 0.0100;
    RootWeight[6][PlantRowColumn][0] = 0.0050;
    RootWeight[6][PlantRowColumn + 1][0] = 0.0050;
    //     Start loop for all soil layers containing roots.
    DepthLastRootLayer = 0;
    TotalRootWeight = 0;
    for (int l = 0; l < 7; l++) {
        DepthLastRootLayer += dl[l];  // compute total depth to the last layer
                                      // with roots (DepthLastRootLayer).
        //     For each soil soil cell with roots, compute total root weight
        //  per plant (TotalRootWeight), and convert RootWeight from g per plant
        //  to g per cell.
        for (int k = 0; k < nk; k++) {
            for (int i = 0; i < 3; i++) {
                TotalRootWeight += RootWeight[l][k][i];
                RootWeight[l][k][i] =
                    RootWeight[l][k][i] * 0.01 * RowSpace / PerPlantArea;
            }
            //     initialize RootAge to a non-zero value for each cell
            //     containing roots.
            if (RootWeight[l][k][0] > 0) RootAge[l][k] = 0.01;
        }
    }
    //     Initial value of taproot length, TapRootLength, is computed to the
    // middle of the last layer with roots. The last soil layer with
    // taproot, LastTaprootLayer, is defined.
    NumLayersWithRoots = 7;
    TapRootLength = (DepthLastRootLayer - 0.5 * dl[NumLayersWithRoots - 1]);
    LastTaprootLayer = 6;
}
//////////////////////////////////////////////////////////////////////
double form(double c0, double d0, double g0)
//     This function computes the aggregation factor for 2 mixed soil materials.
//  Arguments referenced:
//     c0	- heat conductivity of first material
//     d0	- heat conductivity of second material
//     g0	- shape factor for these materials
//
{
    double form0 = (2 / (1 + (c0 / d0 - 1) * g0) +
                    1 / (1 + (c0 / d0 - 1) * (1 - 2 * g0))) /
                   3;
    return form0;
}
