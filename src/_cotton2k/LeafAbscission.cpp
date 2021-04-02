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
#include <cstdint>
#include "global.h"
#include "Simulation.hpp"

using namespace std;

void PreFruitLeafAbscission(State &, double, const int &, const int &, const double &);

void MainStemLeafAbscission(State &, int, int, double, const int &);

void FruitNodeLeafAbscission(State &, int, int, int, double, const int &);

void DefoliationLeafAbscission(State &);

extern "C"
{
    void SortArray(size_t, double[], int32_t *, int32_t *, int32_t *);
}
//////////////////////////////////////////////////
void LeafAbscission(Simulation &sim, uint32_t u)
//     This function simulates leaf abscission. It is called from
//  CottonPhenology(). It calls the following functions:
//        FruitNodeLeafAbscission(), MainStemLeafAbscission(),
//        PreFruitLeafAbscission(), DefoliationLeafAbscission().
//
//     The following global variables are referenced here:
//        DayFirstDef, NumFruitBranches, NumVegBranches,
//        PerPlantArea, TotalLeafArea, TotalLeafWeight.
//
//     The following global variables are set here:
//        AbscisedLeafWeight, LeafAreaIndex, ReserveC.
//
{
    State &state = sim.states[u];
    //     If there are almost no leaves, this routine is not executed.
    if (state.leaf_area_index <= 0.0001)
        return;
    //     Compute droplf as a function of LeafAreaIndex.
    const double vdrop1 = 140, vdrop2 = 1; // constant parameters to compute droplf.
    double droplf;                         // leaf age until its abscission.
    droplf = vdrop1 - vdrop2 * state.leaf_area_index;
    //     Call PreFruitLeafAbscission() to simulate the physiological abscission of
    //  prefruiting node leaves.
    PreFruitLeafAbscission(state, droplf, state.daynum, sim.first_square, state.day_inc);
    //     Loop for all vegetative branches and fruiting branches, and call MainStemLeafAbscission()
    //  for each fruiting branch to simulate the physiological abscission of the other leaves.
    for (int k = 0; k < state.number_of_vegetative_branches; k++)
    {
        for (int l = 0; l < state.vegetative_branches[k].number_of_fruiting_branches; l++)
            MainStemLeafAbscission(state, k, l, droplf, state.daynum);
    }
    //     Call DefoliationLeafAbscission() to simulate leaf abscission caused by defoliants.
    if (DayFirstDef > 0 && state.daynum >= DayFirstDef)
        DefoliationLeafAbscission(state);
    //     If the reserves in the leaf are too high, add the lost reserves
    //  to AbscisedLeafWeight and adjust ReserveC.
    if (ReserveC > 0)
    {
        double resmax;          // maximum possible amount of reserve C in the leaves.
        double resratio = 0.20; // maximum possible ratio of reserve C to leaf dry weight.
        resmax = resratio * TotalLeafWeight;
        if (ReserveC > resmax)
        {
            state.abscised_leaf_weight += ReserveC - resmax;
            ReserveC = resmax;
        }
    }
    //     Compute the resulting LeafAreaIndex but do not let it get too small.
    state.leaf_area_index = TotalLeafArea / PerPlantArea;
    if (state.leaf_area_index < 0.0001)
        state.leaf_area_index = 0.0001;
}

/////////////////////////
void PreFruitLeafAbscission(State &state, double droplf, const int &Daynum, const int &FirstSquare, const double &DayInc)
//     This function simulates the abscission of prefruiting node
//  leaves. It is called from function LeafAbscission().
//
//     The following global variables are referenced here:
//        DayInc, DayFirstDef, FirstSquare, LeafAreaIndex, LeafNConc, NumPreFruNodes,
//        PetioleNConc.
//
//     The following global variable are set here:
//        AbscisedLeafWeight, AgeOfPreFruNode, CumPlantNLoss, LeafAreaPreFru,
//        LeafNitrogen, LeafWeightPreFru, NumAbscisedLeaves, PetioleNitrogen,
//        PetioleWeightPreFru, TotalLeafArea, TotalLeafWeight, TotalPetioleWeight.
//
//     The following  argument is used in this function:
//        droplf - leaf age until it is abscised.
//
{
    //     Loop over all prefruiting nodes. If it is after first square,
    //  node age is updated here.
    for (int j = 0; j < NumPreFruNodes; j++)
    {
        if (FirstSquare > 0)
            AgeOfPreFruNode[j] += DayInc;
        //     The leaf on this node is abscised if its age has reached
        //  droplf, and if there is a leaf here, and if LeafAreaIndex is not too small:
        //     Update TotalLeafArea, AbscisedLeafWeight, TotalLeafWeight, TotalPetioleWeight,
        //  LeafNitrogen, CumPlantNLoss.
        //     Assign zero to LeafAreaPreFru, PetioleWeightPreFru and LeafWeightPreFru of this leaf.
        //     If a defoliation was applied, add 1 to the counter NumAbscisedLeaves.
        if (AgeOfPreFruNode[j] >= droplf && LeafAreaPreFru[j] > 0 && state.leaf_area_index > 0.1)
        {
            TotalLeafArea -= LeafAreaPreFru[j];
            state.abscised_leaf_weight += LeafWeightPreFru[j] + PetioleWeightPreFru[j];
            TotalLeafWeight -= LeafWeightPreFru[j];
            TotalPetioleWeight -= PetioleWeightPreFru[j];
            LeafNitrogen -= LeafWeightPreFru[j] * state.leaf_nitrogen_concentration;
            PetioleNitrogen -= PetioleWeightPreFru[j] * state.petiole_nitrogen_concentration;
            state.cumulative_nitrogen_loss += LeafWeightPreFru[j] * state.leaf_nitrogen_concentration + PetioleWeightPreFru[j] * state.petiole_nitrogen_concentration;
            LeafAreaPreFru[j] = 0;
            LeafWeightPreFru[j] = 0;
            PetioleWeightPreFru[j] = 0;
            if (DayFirstDef > 0 && Daynum > DayFirstDef) // if defoliation has been applied
                NumAbscisedLeaves++;
        }
    }
}

/////////////////////////
void MainStemLeafAbscission(State &state, int k, int l, double droplf, const int &Daynum)
//     This function simulates the abscission of main stem leaves
//  on node l of vegetative branch k. It is called from function
//  LeafAbscission(). It calls function FruitNodeLeafAbscission().
//
//     The following global variables are referenced here:
//        DayFirstDef, LeafAge, LeafAreaIndex, LeafNConc,
//        NumNodes, PetioleNConc.
//
//     The following global variable are set here:
//        AbscisedLeafWeight, CumPlantNLoss, LeafAreaMainStem, LeafNitrogen,
//        LeafWeightMainStem, NumAbscisedLeaves, PetioleWeightMainStem,
//        PetioleNitrogen, TotalLeafArea, TotalLeafWeight, TotalPetioleWeight.
//
//     The following arguments are used in this function:
//        droplf - leaf age until it is abscised.
//        k, l - numbers of this vegetative branch and fruiting branch.
//
{
    //     The leaf on this main stem node is abscised if its age has reached droplf,
    //  and if there is a leaf here, and if LeafAreaIndex is not too small:
    //     Update TotalLeafArea, AbscisedLeafWeight, TotalLeafWeight, TotalPetioleWeight,
    //  LeafNitrogen, CumPlantNLoss.
    //     Assign zero to LeafAreaMainStem, PetioleWeightMainStem and LeafWeightMainStem of this leaf.
    //     If this is after defoliation, add 1 to the counter NumAbscisedLeaves.
    MainStemLeaf &main_stem_leaf = state.vegetative_branches[k].fruiting_branches[l].main_stem_leaf;
    if (state.vegetative_branches[k].fruiting_branches[l].nodes[0].leaf.age > droplf && main_stem_leaf.leaf_area > 0 && state.leaf_area_index > 0.1)
    {
        state.abscised_leaf_weight += main_stem_leaf.leaf_weight + main_stem_leaf.petiole_weight;
        TotalLeafWeight -= main_stem_leaf.leaf_weight;
        TotalPetioleWeight -= main_stem_leaf.petiole_weight;
        LeafNitrogen -= main_stem_leaf.leaf_weight * state.leaf_nitrogen_concentration;
        PetioleNitrogen -= main_stem_leaf.petiole_weight * state.petiole_nitrogen_concentration;
        state.cumulative_nitrogen_loss += main_stem_leaf.leaf_weight * state.leaf_nitrogen_concentration + main_stem_leaf.petiole_weight * state.petiole_nitrogen_concentration;
        TotalLeafArea -= main_stem_leaf.leaf_area;
        main_stem_leaf.leaf_area = 0;
        main_stem_leaf.leaf_weight = 0;
        main_stem_leaf.petiole_weight = 0;
        if (DayFirstDef > 0 && Daynum > DayFirstDef) // if defoliation has been applied
            NumAbscisedLeaves++;
    }
    //     Loop over all nodes on this fruiting branch and call FruitNodeLeafAbscission().
    for (int m = 0; m < state.vegetative_branches[k].fruiting_branches[l].number_of_fruiting_nodes; m++)
        FruitNodeLeafAbscission(state, k, l, m, droplf, Daynum);
}

/////////////////////////
void FruitNodeLeafAbscission(State &state, int k, int l, int m, double droplf, const int &Daynum)
//     This function simulates the abscission of fruiting node
//  leaves on node m of fruiting branch l of vegetative branch k. It is
//  called from function MainStemLeafAbscission().
//
//     The following global variables are referenced here:
//        DayFirstDef, LeafAge, LeafAreaIndex, LeafNConc, PetioleNConc.
//
//     The following global variables are set here:
//        AbscisedLeafWeight, CumPlantNLoss, LeafAreaNodes, LeafNitrogen, LeafWeightNodes,
//        NumAbscisedLeaves, PetioleNitrogen, PetioleWeightNodes,
//        TotalLeafArea, TotalLeafWeight, TotalPetioleWeight.
//
//     The following arguments are used in this function:
//        droplf - leaf age until it is abscised.
//        k, l - numbers of this vegetative branch and fruiting branch.
//        m - node number on this fruiting branch.
//
{
    //     The leaf on this fruiting node is abscised if its age has reached droplf, and
    //  if there is a leaf here, and if LeafAreaIndex is not too small:
    //     Update TotalLeafArea, AbscisedLeafWeight, TotalLeafWeight, TotalPetioleWeight,
    //  LeafNitrogen, CumPlantNLoss,
    //     Assign zero to LeafAreaNodes, PetioleWeightNodes and LeafWeightNodes of this leaf.
    //     If this is after defoliation, add 1 to the NumAbscisedLeaves.
    FruitingSite &site = state.vegetative_branches[k].fruiting_branches[l].nodes[m];
    if (site.leaf.age >= droplf && site.leaf.area > 0 && state.leaf_area_index > 0.1)
    {
        state.abscised_leaf_weight += site.leaf.weight + site.petiole.weight;
        TotalLeafWeight -= site.leaf.weight;
        TotalPetioleWeight -= site.petiole.weight;
        LeafNitrogen -= site.leaf.weight * state.leaf_nitrogen_concentration;
        PetioleNitrogen -= site.petiole.weight * state.petiole_nitrogen_concentration;
        state.cumulative_nitrogen_loss += site.leaf.weight * state.leaf_nitrogen_concentration + site.petiole.weight * state.petiole_nitrogen_concentration;
        TotalLeafArea -= site.leaf.area;
        site.leaf.area = 0;
        site.leaf.weight = 0;
        site.petiole.weight = 0;
        if (DayFirstDef > 0 && Daynum > DayFirstDef) // if defoliation has been applied
            NumAbscisedLeaves++;
    }
}

/////////////////////////
void DefoliationLeafAbscission(State &state)
//     This function simulates leaf abscission caused by defoliants.
//     It is called from function LeafAbscission().
//
//     The following global variables are referenced here:
//        DayFirstDef, LeafNConc, NumPreFruNodes, PercentDefoliation,
//        PetioleNConc.
//
//     The following global variable are set here:
//        AbscisedLeafWeight, CumPlantNLoss, LeafAreaMainStem, LeafAreaNodes, LeafAreaPreFru,
//        LeafNitrogen, LeafWeightMainStem, LeafWeightNodes, LeafWeightPreFru, NumAbscisedLeaves,
//        PetioleNitrogen, PetioleWeightMainStem, PetioleWeightNodes, PetioleWeightPreFru,
//        TotalLeafArea, TotalLeafWeight, TotalPetioleWeight.
//
{
    //     When this is the first day of defoliation - if there are any leaves left on the
    //  prefruiting nodes, they will be shed at this stage.
    if (state.daynum == DayFirstDef)
    {
        for (int j = 0; j < NumPreFruNodes; j++)
        {
            if (LeafAreaPreFru[j] > 0)
            {
                TotalLeafArea -= LeafAreaPreFru[j];
                LeafAreaPreFru[j] = 0;
                state.abscised_leaf_weight += LeafWeightPreFru[j] + PetioleWeightPreFru[j];
                TotalLeafWeight -= LeafWeightPreFru[j];
                TotalPetioleWeight -= PetioleWeightPreFru[j];
                LeafNitrogen -= LeafWeightPreFru[j] * state.leaf_nitrogen_concentration;
                PetioleNitrogen -= PetioleWeightPreFru[j] * state.petiole_nitrogen_concentration;
                state.cumulative_nitrogen_loss += LeafWeightPreFru[j] * state.leaf_nitrogen_concentration + PetioleWeightPreFru[j] * state.petiole_nitrogen_concentration;
                LeafWeightPreFru[j] = 0;
                PetioleWeightPreFru[j] = 0;
            } //  if LeafAreaPreFru
        }     //  for j
    }         //  if Daynum
              //     When this is after the first day of defoliation - count the
              //  number of existing leaves (lefcnt) and sort them by age
    if (state.daynum <= DayFirstDef)
        return;
    //
    double SortByAge[450]; // array of site ages to be sorted
    int indexk[450];       // index associated with k
    int indexl[450];       // index associated with l
    int indexm[450];       // index associated with m
    int lefcnt = 0;        // counter of existing leaves
    for (int k = 0; k < state.number_of_vegetative_branches; k++)
    {
        for (int l = 0; l < state.vegetative_branches[k].number_of_fruiting_branches; l++)
        {
            if (state.vegetative_branches[k].fruiting_branches[l].main_stem_leaf.leaf_weight > 0)
            {
                SortByAge[lefcnt] = state.vegetative_branches[k].fruiting_branches[l].nodes[0].age;
                indexk[lefcnt] = k;
                indexl[lefcnt] = l;
                indexm[lefcnt] = 66; // 66 indicates this leaf is at the base of the fruiting branch
                lefcnt++;
            }
            for (int m = 0; m < state.vegetative_branches[k].fruiting_branches[l].number_of_fruiting_nodes; m++)
            {
                if (state.vegetative_branches[k].fruiting_branches[l].nodes[m].leaf.weight > 0)
                {
                    SortByAge[lefcnt] = state.vegetative_branches[k].fruiting_branches[l].nodes[m].age;
                    indexk[lefcnt] = k;
                    indexl[lefcnt] = l;
                    indexm[lefcnt] = m;
                    lefcnt++;
                }
            } // for m
        }     // for l
    }         // for k
              //     Call function to sort the leaves by their age
    SortArray(lefcnt, SortByAge, indexk, indexl, indexm);
    //     Compute the number of leaves to be shed on this day (numLeavesToShed).
    int numLeavesToShed; // the computed number of leaves to be shed.
    numLeavesToShed = (int)(lefcnt * PercentDefoliation / 100);
    //  Execute leaf shedding according to leaf age (by order of lefcnt).
    for (int i = 0; i < lefcnt; i++)
    {
        if (numLeavesToShed <= 0)
            break;
        if (numLeavesToShed > 0 && SortByAge[i] > 0)
        {
            int k = indexk[i];
            int l = indexl[i];
            int m = indexm[i];
            FruitingSite &site = state.vegetative_branches[k].fruiting_branches[l].nodes[m];
            if (m == 66) // main stem leaves
            {
                MainStemLeaf &main_stem_leaf = state.vegetative_branches[k].fruiting_branches[l].main_stem_leaf;
                state.abscised_leaf_weight += main_stem_leaf.leaf_weight + main_stem_leaf.petiole_weight;
                TotalLeafWeight -= main_stem_leaf.leaf_weight;
                TotalPetioleWeight -= main_stem_leaf.petiole_weight;
                LeafNitrogen -= main_stem_leaf.leaf_weight * state.leaf_nitrogen_concentration;
                PetioleNitrogen -= main_stem_leaf.petiole_weight * state.petiole_nitrogen_concentration;
                state.cumulative_nitrogen_loss += main_stem_leaf.leaf_weight * state.leaf_nitrogen_concentration + main_stem_leaf.petiole_weight * state.petiole_nitrogen_concentration;
                TotalLeafArea -= main_stem_leaf.leaf_area;
                main_stem_leaf.leaf_area = 0;
                main_stem_leaf.leaf_weight = 0;
                main_stem_leaf.petiole_weight = 0;
            }
            else // leaves on fruit nodes
            {
                state.abscised_leaf_weight += site.leaf.weight + site.petiole.weight;
                TotalLeafWeight -= site.leaf.weight;
                TotalPetioleWeight -= site.petiole.weight;
                LeafNitrogen -= site.leaf.weight * state.leaf_nitrogen_concentration;
                PetioleNitrogen -= site.petiole.weight * state.petiole_nitrogen_concentration;
                state.cumulative_nitrogen_loss += site.leaf.weight * state.leaf_nitrogen_concentration + site.petiole.weight * state.petiole_nitrogen_concentration;
                TotalLeafArea -= site.leaf.area;
                site.leaf.area = 0;
                site.leaf.weight = 0;
                site.petiole.weight = 0;
            } // if m
            numLeavesToShed--;
            NumAbscisedLeaves++;
        } // if numLeavesToShed
    }     // for i
}