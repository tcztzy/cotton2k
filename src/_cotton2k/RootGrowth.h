#pragma once
#include <cinttypes>
#include "Simulation.hpp"
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

void RootSummation(State &state, int NumRootAgeGroups, double row_space, double per_plant_area);

double RootCultivation(SoilCell[40][20], int, double, double, double);

void LateralRootGrowthLeft(State &, int, int, unsigned int, double);

void LateralRootGrowthRight(State &, int, int, unsigned int, double);

void InitiateLateralRoots();
