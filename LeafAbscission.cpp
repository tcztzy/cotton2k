//  LeafAbscission.cpp
//
//   functions in this file:
// LeafAbscission()
// PreFruitLeafAbscission()
// MainStemLeafAbscission()
// FruitNodeLeafAbscission()
// DefoliationLeafAbscission()
// SortArray()
//
#include "CottonSimulation.h"
#include "GeneralFunctions.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif
//
//////////////////////////////////////////////////
void LeafAbscission()
//     This function simulates leaf abscission. It is called from
//  CottonPhenology(). It calls the following functions:
//        FruitNodeLeafAbscission(), MainStemLeafAbscission(),
//        PreFruitLeafAbscission(), DefoliationLeafAbscission().
//
//     The following global variables are referenced here:
//        Daynum, DayFirstDef, NumFruitBranches, NumVegBranches,
//        PerPlantArea, TotalLeafArea, TotalLeafWeight.
//
//     The following global variables are set here:
//        AbscisedLeafWeight, LeafAreaIndex, ReserveC.
//
{
    //     If there are almost no leaves, this routine is not executed.
    if (LeafAreaIndex <= 0.0001) return;
    //     Compute droplf as a function of LeafAreaIndex.
    double droplf = drop_leaf_age(LeafAreaIndex);  // leaf age until its abscission.
    //     Call PreFruitLeafAbscission() to simulate the physiological
    //     abscission of
    //  prefruiting node leaves.
    PreFruitLeafAbscission(droplf);
    //     Loop for all vegetative branches and fruiting branches, and call
    //     MainStemLeafAbscission()
    //  for each fruiting branch to simulate the physiological abscission of the
    //  other leaves.
    for (int k = 0; k < NumVegBranches; k++) {
        int nbrch;  // number of fruiting branches on a vegetative branch.
        nbrch = NumFruitBranches[k];
        for (int l = 0; l < nbrch; l++) MainStemLeafAbscission(k, l, droplf);
    }
    //     Call DefoliationLeafAbscission() to simulate leaf abscission caused
    //     by defoliants.
    if (DayFirstDef > 0 && Daynum >= DayFirstDef) DefoliationLeafAbscission();
    //     If the reserves in the leaf are too high, add the lost reserves
    //  to AbscisedLeafWeight and adjust ReserveC.
    if (ReserveC > 0) {
        double resmax;  // maximum possible amount of reserve C in the leaves.
        double resratio =
            0.20;  // maximum possible ratio of reserve C to leaf dry weight.
        resmax = resratio * TotalLeafWeight();
        if (ReserveC > resmax) {
            AbscisedLeafWeight += ReserveC - resmax;
            ReserveC = resmax;
        }
    }
    //     Compute the resulting LeafAreaIndex but do not let it get too small.
    LeafAreaIndex = TotalLeafArea() / PerPlantArea;
    if (LeafAreaIndex < 0.0001) LeafAreaIndex = 0.0001;
}
/////////////////////////
void PreFruitLeafAbscission(double droplf)
//     This function simulates the abscission of prefruiting node
//  leaves. It is called from function LeafAbscission().
//
//     The following global variables are referenced here:
//        DayInc, Daynum, DayFirstDef, FirstSquare, LeafAreaIndex, LeafNConc,
//        NumPreFruNodes, PetioleNConc, pixcon.
//
//     The following global variable are set here:
//        AbscisedLeafWeight, AgeOfPreFruNode, CumPlantNLoss, LeafAreaPreFru,
//        LeafNitrogen, LeafWeightPreFru, NumAbscisedLeaves, PetioleNitrogen,
//        PetioleWeightPreFru, PixInPlants, TotalLeafArea, TotalLeafWeight,
//        TotalPetioleWeight.
//
//     The following  argument is used in this function:
//        droplf - leaf age until it is abscised.
//
{
    //     Loop over all prefruiting nodes. If it is after first square,
    //  node age is updated here.
    for (int j = 0; j < NumPreFruNodes; j++) {
        if (FirstSquare > 0) AgeOfPreFruNode[j] += DayInc;
        //     The leaf on this node is abscised if its age has reached
        //  droplf, and if there is a leaf here, and if LeafAreaIndex is not too
        //  small:
        //     Update TotalLeafArea, AbscisedLeafWeight, TotalLeafWeight,
        //     TotalPetioleWeight,
        //  PixInPlants, LeafNitrogen, CumPlantNLoss.
        //     Assign zero to LeafAreaPreFru, PetioleWeightPreFru and
        //     LeafWeightPreFru of this leaf. If a defoliation was applied, add
        //     1 to the counter NumAbscisedLeaves.
        if (AgeOfPreFruNode[j] >= droplf && LeafAreaPreFru[j] > 0 &&
            LeafAreaIndex > 0.1) {
            LeafArea[NodeLayerPreFru[j]] -= LeafAreaPreFru[j];
            AbscisedLeafWeight += LeafWeightPreFru[j] + PetioleWeightPreFru[j];
            TotalPetioleWeight -= PetioleWeightPreFru[j];
            PixInPlants -= LeafWeightPreFru[j] * pixcon;
            LeafNitrogen -= LeafWeightPreFru[j] * LeafNConc;
            PetioleNitrogen -= PetioleWeightPreFru[j] * PetioleNConc;
            CumPlantNLoss += LeafWeightPreFru[j] * LeafNConc +
                             PetioleWeightPreFru[j] * PetioleNConc;
            LeafAreaPreFru[j] = 0;
            LeafWeightPreFru[j] = 0;
            PetioleWeightPreFru[j] = 0;
            if (DayFirstDef > 0 &&
                Daynum > DayFirstDef)  // if defoliation has been applied
                NumAbscisedLeaves++;
        }
    }
}
/////////////////////////
void MainStemLeafAbscission(int k, int l, double droplf)
//     This function simulates the abscission of main stem leaves
//  on node l of vegetative branch k. It is called from function
//  LeafAbscission(). It calls function FruitNodeLeafAbscission().
//
//     The following global variables are referenced here:
//        Daynum, DayFirstDef, LeafAge, LeafAreaIndex, LeafNConc,
//        NumNodes, PetioleNConc, pixcon.
//
//     The following global variable are set here:
//        AbscisedLeafWeight, CumPlantNLoss, LeafAreaMainStem, LeafNitrogen,
//        LeafWeightMainStem, NumAbscisedLeaves, PetioleWeightMainStem,
//        PetioleNitrogen, PixInPlants, TotalLeafArea, TotalLeafWeight,
//        TotalPetioleWeight.
//
//     The following arguments are used in this function:
//        droplf - leaf age until it is abscised.
//        k, l - numbers of this vegetative branch and fruiting branch.
//
{
    //     The leaf on this main stem node is abscised if its age has reached
    //     droplf,
    //  and if there is a leaf here, and if LeafAreaIndex is not too small:
    //     Update TotalLeafArea, AbscisedLeafWeight, TotalLeafWeight,
    //     TotalPetioleWeight,
    //  PixInPlants, LeafNitrogen, CumPlantNLoss.
    //     Assign zero to LeafAreaMainStem, PetioleWeightMainStem and
    //     LeafWeightMainStem of this leaf. If this is after defoliation, add 1
    //     to the counter NumAbscisedLeaves.
    if (LeafAge[k][l][0] > droplf && LeafAreaMainStem[k][l] > 0 &&
        LeafAreaIndex > 0.1) {
        AbscisedLeafWeight +=
            LeafWeightMainStem[k][l] + PetioleWeightMainStem[k][l];
        TotalPetioleWeight -= PetioleWeightMainStem[k][l];
        PixInPlants -= LeafWeightMainStem[k][l] * pixcon;
        LeafNitrogen -= LeafWeightMainStem[k][l] * LeafNConc;
        PetioleNitrogen -= PetioleWeightMainStem[k][l] * PetioleNConc;
        CumPlantNLoss += LeafWeightMainStem[k][l] * LeafNConc +
                         PetioleWeightMainStem[k][l] * PetioleNConc;
        LeafArea[NodeLayer[k][l]] -= LeafAreaMainStem[k][l];
        LeafAreaMainStem[k][l] = 0;
        LeafWeightMainStem[k][l] = 0;
        PetioleWeightMainStem[k][l] = 0;
        if (DayFirstDef > 0 &&
            Daynum > DayFirstDef)  // if defoliation has been applied
            NumAbscisedLeaves++;
    }
    //     Loop over all nodes on this fruiting branch and call
    //     FruitNodeLeafAbscission().
    int nnid = NumNodes[k][l];  // node number on this fruiting branch.
    for (int m = 0; m < nnid; m++) FruitNodeLeafAbscission(k, l, m, droplf);
}
/////////////////////////
void FruitNodeLeafAbscission(int k, int l, int m, double droplf)
//     This function simulates the abscission of fruiting node
//  leaves on node m of fruiting branch l of vegetative branch k. It is
//  called from function MainStemLeafAbscission().
//
//     The following global variables are referenced here:
//        Daynum, DayFirstDef, LeafAge, LeafAreaIndex, LeafNConc, PetioleNConc,
//        pixcon.
//
//     The following global variables are set here:
//        AbscisedLeafWeight, CumPlantNLoss, LeafAreaNodes, LeafNitrogen,
//        LeafWeightNodes, NumAbscisedLeaves, PetioleNitrogen,
//        PetioleWeightNodes, TotalLeafArea, TotalLeafWeight,
//        TotalPetioleWeight, PixInPlants.
//
//     The following arguments are used in this function:
//        droplf - leaf age until it is abscised.
//        k, l - numbers of this vegetative branch and fruiting branch.
//        m - node number on this fruiting branch.
//
{
    //     The leaf on this fruiting node is abscised if its age has reached
    //     droplf, and
    //  if there is a leaf here, and if LeafAreaIndex is not too small:
    //     Update TotalLeafArea, AbscisedLeafWeight, TotalLeafWeight,
    //     TotalPetioleWeight,
    //  PixInPlants, LeafNitrogen, CumPlantNLoss,
    //     Assign zero to LeafAreaNodes, PetioleWeightNodes and LeafWeightNodes
    //     of this leaf. If this is after defoliation, add 1 to the
    //     NumAbscisedLeaves.
    if (LeafAge[k][l][m] >= droplf && LeafAreaNodes[k][l][m] > 0 &&
        LeafAreaIndex > 0.1) {
        AbscisedLeafWeight +=
            LeafWeightNodes[k][l][m] + PetioleWeightNodes[k][l][m];
        TotalPetioleWeight -= PetioleWeightNodes[k][l][m];
        PixInPlants -= LeafWeightNodes[k][l][m] * pixcon;
        LeafNitrogen -= LeafWeightNodes[k][l][m] * LeafNConc;
        PetioleNitrogen -= PetioleWeightNodes[k][l][m] * PetioleNConc;
        CumPlantNLoss += LeafWeightNodes[k][l][m] * LeafNConc +
                         PetioleWeightNodes[k][l][m] * PetioleNConc;
        LeafArea[NodeLayer[k][l]] -= LeafAreaNodes[k][l][m];
        LeafAreaNodes[k][l][m] = 0;
        LeafWeightNodes[k][l][m] = 0;
        PetioleWeightNodes[k][l][m] = 0;
        if (DayFirstDef > 0 &&
            Daynum > DayFirstDef)  // if defoliation has been applied
            NumAbscisedLeaves++;
    }
}
/////////////////////////
void DefoliationLeafAbscission()
//     This function simulates leaf abscission caused by defoliants.
//     It is called from function LeafAbscission().
//
//     The following global variables are referenced here:
//        Daynum, DayFirstDef, LeafNConc, NumPreFruNodes, PercentDefoliation,
//        PetioleNConc, pixcon.
//
//     The following global variable are set here:
//        AbscisedLeafWeight, CumPlantNLoss, LeafAreaMainStem, LeafAreaNodes,
//        LeafAreaPreFru, LeafNitrogen, LeafWeightMainStem, LeafWeightNodes,
//        LeafWeightPreFru, NumAbscisedLeaves, PetioleNitrogen,
//        PetioleWeightMainStem, PetioleWeightNodes, PetioleWeightPreFru,
//        PixInPlants, TotalLeafArea, TotalLeafWeight, TotalPetioleWeight.
//
{
    //     When this is the first day of defoliation - if there are any leaves
    //     left on the
    //  prefruiting nodes, they will be shed at this stage.
    if (Daynum == DayFirstDef) {
        for (int j = 0; j < NumPreFruNodes; j++) {
            if (LeafAreaPreFru[j] > 0) {
                LeafArea[NodeLayerPreFru[j]] -= LeafAreaPreFru[j];
                LeafAreaPreFru[j] = 0;
                AbscisedLeafWeight +=
                    LeafWeightPreFru[j] + PetioleWeightPreFru[j];
                TotalPetioleWeight -= PetioleWeightPreFru[j];
                PixInPlants -= LeafWeightPreFru[j] * pixcon;
                LeafNitrogen -= LeafWeightPreFru[j] * LeafNConc;
                PetioleNitrogen -= PetioleWeightPreFru[j] * PetioleNConc;
                CumPlantNLoss += LeafWeightPreFru[j] * LeafNConc +
                                 PetioleWeightPreFru[j] * PetioleNConc;
                LeafWeightPreFru[j] = 0;
                PetioleWeightPreFru[j] = 0;
            }  //  if LeafAreaPreFru
        }      //  for j
    }          //  if Daynum
    //     When this is after the first day of defoliation - count the
    //  number of existing leaves (lefcnt) and sort them by age
    if (Daynum <= DayFirstDef) return;
    //
    double SortByAge[450];  // array of site ages to be sorted
    int indexk[450];        // index associated with k
    int indexl[450];        // index associated with l
    int indexm[450];        // index associated with m
    int lefcnt = 0;         // counter of existing leaves
    for (int k = 0; k < NumVegBranches; k++) {
        int nbrch = NumFruitBranches[k];  // number of fruiting branches on a
                                          // vegetative branch.
        for (int l = 0; l < nbrch; l++) {
            if (LeafWeightMainStem[k][l] > 0) {
                SortByAge[lefcnt] = AgeOfSite[k][l][0];
                indexk[lefcnt] = k;
                indexl[lefcnt] = l;
                indexm[lefcnt] = 66;  // 66 indicates this leaf is at the base
                                      // of the fruiting branch
                lefcnt++;
            }
            int nnid = NumNodes[k][l];  // total existing node number on this
                                        // fruiting branch.
            for (int m = 0; m < nnid; m++) {
                if (LeafWeightNodes[k][l][m] > 0) {
                    SortByAge[lefcnt] = AgeOfSite[k][l][m];
                    indexk[lefcnt] = k;
                    indexl[lefcnt] = l;
                    indexm[lefcnt] = m;
                    lefcnt++;
                }
            }  // for m
        }      // for l
    }          // for k
               //     Call function to sort the leaves by their age
    SortArray(lefcnt, SortByAge, indexk, indexl, indexm);
    //     Compute the number of leaves to be shed on this day
    //     (numLeavesToShed).
    int numLeavesToShed;  // the computed number of leaves to be shed.
    numLeavesToShed = (int)(lefcnt * PercentDefoliation / 100);
    //  Execute leaf shedding according to leaf age (by order of lefcnt).
    for (int i = 0; i < lefcnt; i++) {
        if (numLeavesToShed <= 0) break;
        if (numLeavesToShed > 0 && SortByAge[i] > 0) {
            double pixlos;  // amount of pix lost, g per plant.
            int k = indexk[i];
            int l = indexl[i];
            int m = indexm[i];
            if (m == 66)  // main stem leaves
            {
                AbscisedLeafWeight +=
                    LeafWeightMainStem[k][l] + PetioleWeightMainStem[k][l];
                TotalPetioleWeight -= PetioleWeightMainStem[k][l];
                pixlos = LeafWeightMainStem[k][l] * pixcon;
                LeafNitrogen -= LeafWeightMainStem[k][l] * LeafNConc;
                PetioleNitrogen -= PetioleWeightMainStem[k][l] * PetioleNConc;
                CumPlantNLoss += LeafWeightMainStem[k][l] * LeafNConc +
                                 PetioleWeightMainStem[k][l] * PetioleNConc;
                LeafArea[NodeLayer[k][l]] -= LeafAreaMainStem[k][l];
                LeafAreaMainStem[k][l] = 0;
                LeafWeightMainStem[k][l] = 0;
                PetioleWeightMainStem[k][l] = 0;
            } else  // leaves on fruit nodes
            {
                AbscisedLeafWeight +=
                    LeafWeightNodes[k][l][m] + PetioleWeightNodes[k][l][m];
                TotalPetioleWeight -= PetioleWeightNodes[k][l][m];
                pixlos = LeafWeightNodes[k][l][m] * pixcon;
                LeafNitrogen -= LeafWeightNodes[k][l][m] * LeafNConc;
                PetioleNitrogen -= PetioleWeightNodes[k][l][m] * PetioleNConc;
                CumPlantNLoss += LeafWeightNodes[k][l][m] * LeafNConc +
                                 PetioleWeightNodes[k][l][m] * PetioleNConc;
                LeafArea[NodeLayer[k][l]] -= LeafAreaNodes[k][l][m];
                LeafAreaNodes[k][l][m] = 0;
                LeafWeightNodes[k][l][m] = 0;
                PetioleWeightNodes[k][l][m] = 0;
            }  // if m
            PixInPlants -= pixlos;
            numLeavesToShed--;
            NumAbscisedLeaves++;
        }  // if numLeavesToShed
    }      // for i
}
///////////////////////////////////////////////////////////////////////////////
void SortArray(int size, double InData[], int indexk[], int indexl[],
               int indexm[])
//     This function sorts an array of values by its value (larger is first)
//  together with three indexes associated with each value. It is called from
//  DefoliationLeafAbscission().
//     Arguments:
//        size -    size of the arrays
//        InData - array to be sorted by its values (descending)
//        indexk - first array of associated indexes
//        indexl - second array of associated indexes
//        indexm - third array of associated indexes
//
{
    //     Zeroise temporary output arrays
    int outk[450], outl[450], outm[450];
    double OutData[450];
    for (int i = 0; i < size; i++) {
        outk[i] = 0;
        outl[i] = 0;
        outm[i] = 0;
        OutData[i] = 0;
    }
    //     Loop over arrays, extract values of the ith member
    for (int i = 0; i < size; i++) {
        double value = InData[i];
        int n0 = indexk[i];
        int n1 = indexl[i];
        int n2 = indexm[i];
        for (int j = 0; j < size; j++) {
            if (value > OutData[j]) {
                for (int k = size - 1; k > j; k--) {
                    OutData[k] = OutData[k - 1];
                    outk[k] = outk[k - 1];
                    outl[k] = outl[k - 1];
                    outm[k] = outm[k - 1];
                }
                OutData[j] = value;
                outk[j] = n0;
                outl[j] = n1;
                outm[j] = n2;
                break;
            }  // if value
        }      // for j
    }          // for i
               //     Transfer values back to the input arrays
    for (int i = 0; i < size; i++) {
        indexk[i] = outk[i];
        indexl[i] = outl[i];
        indexm[i] = outm[i];
        InData[i] = OutData[i];
    }
}
