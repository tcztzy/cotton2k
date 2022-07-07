use crate::utils::fmin;
use crate::{
    ActualBollGrowth, ActualBurrGrowth, ActualSquareGrowth, ActualStemGrowth, AgeOfPreFruNode,
    BurrNConc, BurrNitrogen, BurrWeightGreenBolls, BurrWeightOpenBolls,
    CarbonAllocatedForRootGrowth, CottonWeightGreenBolls, CottonWeightOpenBolls, DayEmerge, Daynum,
    ExtraCarbon, Gintot, Kday, LeafAge, LeafNConc, LeafNitrogen, NStressFruiting, NStressRoots,
    NStressVeg, NitrogenStress, NumFruitBranches, NumNodes, NumPreFruNodes, NumVegBranches,
    PetioleNConc, PetioleNO3NConc, PetioleNitrogen, RootNConc, RootNitrogen, SeedNConc,
    SeedNitrogen, SquareNConc, SquareNitrogen, StemNConc, StemNitrogen, StemWeight, SupplyNH4N,
    SupplyNO3N, TotalActualLeafGrowth, TotalActualPetioleGrowth, TotalLeafWeight,
    TotalPetioleWeight, TotalRequiredN, TotalRootWeight, TotalSquareWeight, TotalStemWeight,
};
#[derive(Debug, Clone, Copy)]
pub struct PlantNitrogen {
    // for plant_nitrogen
    /// daily added nitrogen to fruit, g per plant.
    pub addnf: f64,
    /// daily added nitrogen to root, g per plant.
    pub addnr: f64,
    /// daily added nitrogen to vegetative shoot, g per plant.
    pub addnv: f64,
    /// amount of nitrogen not used for growth of plant parts.
    pub xtran: f64,
    /// reserve N in leaves, in g per plant.
    pub leafrs: f64,
    /// reserve N in petioles, in g per plant.
    pub petrs: f64,
    /// reserve N in stems, in g per plant.
    pub stemrs: f64,
    /// reserve N in roots, in g per plant.
    pub rootrs: f64,
    /// reserve N in burrs, in g per plant.
    pub burres: f64,
    /// nitrogen requirement for fruit growth.
    pub reqf: f64,
    /// total nitrogen requirement for plant growth.
    pub reqtot: f64,
    /// nitrogen requirement for vegetative shoot growth.
    pub reqv: f64,
    /// nitrogen requirement for burr growth.
    pub rqnbur: f64,
    /// nitrogen requirement for leaf growth.
    pub rqnlef: f64,
    /// nitrogen requirement for petiole growth.
    pub rqnpet: f64,
    /// nitrogen requirement for root growth.
    pub rqnrut: f64,
    /// nitrogen requirement for seed growth.
    pub rqnsed: f64,
    /// nitrogen requirement for square growth.
    pub rqnsqr: f64,
    /// nitrogen requirement for stem growth.
    pub rqnstm: f64,
    /// total nitrogen available for growth.
    pub npool: f64,
    /// nitrogen uptake from the soil, g per plant.
    pub uptn: f64,
}

impl PlantNitrogen {
    pub fn new() -> Self {
        PlantNitrogen {
            addnf: 0.,
            addnr: 0.,
            addnv: 0.,
            xtran: 0.,
            leafrs: 0.,
            petrs: 0.,
            stemrs: 0.,
            rootrs: 0.,
            burres: 0.,
            reqf: 0.,
            reqtot: 0.,
            reqv: 0.,
            rqnbur: 0.,
            rqnlef: 0.,
            rqnpet: 0.,
            rqnrut: 0.,
            rqnsed: 0.,
            rqnsqr: 0.,
            rqnstm: 0.,
            npool: 0.,
            uptn: 0.,
        }
    }
    /// This function simulates the nitrogen accumulation and distribution in cotton plants, and computes nitrogen stresses. It is called from SimulateThisDay().
    ///
    /// The maximum and minimum N concentrations for the various organs are modified from those reported by: Jones et. al. (1974): Development of a nitrogen balance for cotton growth models: a first approximation. Crop Sci. 14:541-546.
    ///
    /// The following parameters are used in all plant N routines:
    ///
    /// |        |Growth requirement|Uptake requirement| Minimum content |
    /// |--------|------------------|------------------|-----------------|
    /// |leaves  | lefcn0    = .064 | vnreqlef  = .042 | vlfnmin   = .018|
    /// |petioles| petcn0    = .032 | vnreqpet  = .036 | vpetnmin  = .005|
    /// |stems   | stmcn0    = .036 | vnreqstm  = .012 | vstmnmin  = .006|
    /// |roots   | rootcn0   = .018 | vnreqrt   = .010 | vrtnmin   = .010|
    /// |burrs   | burcn0    = .026 | vnreqbur  = .012 | vburnmin  = .006|
    /// |seeds   | seedcn0   = .036 | seedcn1   = .045 |                 |
    /// |squares | sqrcn0    = .024 | vnreqsqr  = .024 |                 |
    ///
    /// The following global variables are referenced here:
    ///
    /// BurrNConc, BurrNitrogen, Kday, LeafNConc, LeafNitrogen,
    /// PetioleNConc, PetioleNConc, PetioleNO3NConc, PetioleNitrogen,
    /// RootNConc, RootNitrogen, SeedNConc, SeedNitrogen, SquareNitrogen,
    /// StemNConc, StemNitrogen.
    ///
    /// The following global and file scope variables are set here:
    ///
    /// addnf, addnr, addnv, burres, leafrs, npool, petrs, reqf, reqtot, reqv,
    /// rootrs, rqnbur, rqnlef, rqnpet, rqnrut, rqnsed, rqnsqr, rqnstm, stemrs,
    /// uptn, xtran.
    pub fn run(&mut self) {
        //     Assign zero to some variables.
        self.leafrs = 0.;
        self.petrs = 0.;
        self.stemrs = 0.;
        self.rootrs = 0.; // reserve N in roots, in g per plant.
        self.reqf = 0.; // nitrogen requirement for fruit growth.
        self.reqtot = 0.; // total nitrogen requirement for plant growth.
        self.reqv = 0.; // nitrogen requirement for vegetative shoot growth.
        self.rqnbur = 0.; // nitrogen requirement for burr growth.
        self.rqnlef = 0.; // nitrogen requirement for leaf growth.
        self.rqnpet = 0.; // nitrogen requirement for petiole growth.
        self.rqnrut = 0.; // nitrogen requirement for root growth.
        self.rqnsed = 0.; // nitrogen requirement for seed growth.
        self.rqnsqr = 0.; // nitrogen requirement for square growth.
        self.rqnstm = 0.; // nitrogen requirement for stem growth.
        self.npool = 0.; // total nitrogen available for growth.
        self.uptn = 0.; // nitrogen uptake from the soil, g per plant.
        self.xtran = 0.; // amount of nitrogen not used for growth of plant parts.
        self.addnf = 0.; // daily added nitrogen to fruit, g per plant.
        self.addnr = 0.; // daily added nitrogen to root, g per plant.
        self.addnv = 0.; // daily added nitrogen to vegetative shoot, g per plant.

        self.nitrogen_requirement(); //  computes the N requirements for growth.
        self.nitrogen_supply(); //  computes the supply of N from uptake and reserves.
        self.nitrogen_allocation(); //  computes the allocation of N in the plant.
        self.nitrogen_content(); // computes the concentrations of N in plant dry matter.
        self.nitrogen_stress(); //  computes nitrogen stress factors.
        self.nitrogen_uptake_requirement(); // computes N requirements for uptake
    }
    /// This function computes the N requirements for growth. It is called from [PlantNitrogen::plant_nitrogen()].
    ///
    /// The following global variables are referenced here:
    ///
    /// * [ActualBollGrowth]
    /// * [ActualBurrGrowth]
    /// * [ActualSquareGrowth]
    /// * [ActualStemGrowth]
    /// * [CarbonAllocatedForRootGrowth]
    /// * [CottonWeightGreenBolls]
    /// * [DayEmerge]
    /// * [Daynum]
    /// * [ExtraCarbon]
    /// * [SeedNitrogen]
    /// * [TotalActualLeafGrowth]
    /// * [TotalActualPetioleGrowth]
    /// * [PetioleNConc]
    /// * [PetioleNO3NConc]
    fn nitrogen_requirement(&mut self) {
        //     The following constant parameters are used:
        const burcn0: f64 = 0.026; //  maximum N content for burrs
        const lefcn0: f64 = 0.064; //  maximum N content for leaves
        const petcn0: f64 = 0.032; //  maximum N content for petioles
        const rootcn0: f64 = 0.018; //  maximum N content for roots
        const seedcn0: f64 = 0.036; //  maximum N content for seeds
        const seedcn1: f64 = 0.045; //  additional requirement of N for existing seed tissue.
        const seedratio: f64 = 0.64; //  ratio of seeds to seedcotton in green bolls.
        const sqrcn0: f64 = 0.024; //  maximum N content for squares
        const stmcn0: f64 = 0.036; //  maximum N content for stems
                                   //     On emergence, assign initial values to petiole N concentrations.
        if unsafe { Daynum <= DayEmerge } {
            unsafe {
                PetioleNConc = petcn0;
                PetioleNO3NConc = petcn0;
            }
        }
        // Compute the nitrogen requirements for growth, by multiplying the daily added dry weight by the maximum N content of each organ.
        // Nitrogen requirements are based on actual growth rates.
        //
        // These N requirements will be used to compute the allocation of N to plant parts and the nitrogen stress factors.
        //
        // All nitrogen requirement variables are in g N per plant.
        // for leaf blade
        self.rqnlef = lefcn0 * unsafe { TotalActualLeafGrowth };
        // for petiole
        self.rqnpet = petcn0 * unsafe { TotalActualPetioleGrowth };
        // for stem
        self.rqnstm = stmcn0 * unsafe { ActualStemGrowth };
        // Add ExtraCarbon to CarbonAllocatedForRootGrowth to compute the total supply of carbohydrates for root growth.
        // for root
        self.rqnrut = rootcn0 * unsafe { CarbonAllocatedForRootGrowth + ExtraCarbon };
        // for squares
        self.rqnsqr = unsafe { ActualSquareGrowth } * sqrcn0;
        // components of seed N requirements.
        // for seed growth
        let rqnsed1 = unsafe { ActualBollGrowth } * seedratio * seedcn0;
        // The N required for replenishing the N content of existing seed tissue (rqnsed2) is added to seed growth requirement.
        let rqnsed2 = if unsafe { CottonWeightGreenBolls > ActualBollGrowth } {
            // existing ratio of N to dry matter in the seeds.
            let rseedn =
                unsafe { SeedNitrogen / ((CottonWeightGreenBolls - ActualBollGrowth) * seedratio) };
            let result = unsafe { CottonWeightGreenBolls - ActualBollGrowth }
                * seedratio
                * (seedcn1 - rseedn);
            if result < 0. {
                0.
            } else {
                result
            }
        } else {
            0.
        };
        self.rqnsed = rqnsed1 + rqnsed2; // total requirement for seeds
        self.rqnbur = unsafe { ActualBurrGrowth } * burcn0; // for burrs
        self.reqf = self.rqnsqr + self.rqnsed + self.rqnbur; // total for fruit
        self.reqv = self.rqnlef + self.rqnpet + self.rqnstm; // total for shoot
        self.reqtot = self.rqnrut + self.reqv + self.reqf; // total N requirement
    }
    /// This function computes the supply of N by uptake from the soil reserves,
    /// it is called from [PlantNitrogen::plant_nitrogen()], it calls [PlantNitrogen::petiole_nitrate()].
    fn nitrogen_supply(&mut self) {
        //     The following constant parameters are used:
        const MobilizNFractionBurrs: f64 = 0.08; //  fraction of N mobilizable for burrs
        const MobilizNFractionLeaves: f64 = 0.09; //  fraction of N mobilizable for leaves and petioles
        const MobilizNFractionStemRoot: f64 = 0.40; //  fraction of N mobilizable for stems and roots
        const vburnmin: f64 = 0.006; //  minimum N contents of burrs
        const vlfnmin: f64 = 0.018; //  minimum N contents of leaves
        const vpetnmin: f64 = 0.005; //  minimum N contents of petioles, non-nitrate fraction.
        const vpno3min: f64 = 0.002; //  minimum N contents of petioles, nitrate fraction.
        const vrtnmin: f64 = 0.010; //  minimum N contents of roots
        const vstmnmin: f64 = 0.006; //  minimum N contents of stems
                                     //     uptn is the total supply of nitrogen to the plant by uptake of
                                     //     nitrate and ammonium.
        self.uptn = unsafe { SupplyNO3N + SupplyNH4N };
        // If total N requirement is less than the supply, define npool as the supply and assign zero to the N reserves in all organs.
        if self.reqtot <= self.uptn {
            self.npool = self.uptn;
            self.leafrs = 0.;
            self.petrs = 0.;
            self.stemrs = 0.;
            self.rootrs = 0.;
            self.burres = 0.;
            self.xtran = 0.;
        } else {
            //     If total N requirement exceeds the supply, compute the nitrogen
            //  reserves in the plant. The reserve N in an organ is defined as a
            //  fraction of the nitrogen content exceeding a minimum N content in
            //  it.
            //     The N reserves in leaves, petioles, stems, roots and burrs of
            //  green bolls are computed, and their N content updated.
            self.leafrs =
                unsafe { (LeafNitrogen - vlfnmin * TotalLeafWeight()) * MobilizNFractionLeaves };
            if self.leafrs < 0. {
                self.leafrs = 0.;
            }
            unsafe {
                LeafNitrogen -= self.leafrs;
            }
            // The petiole N content is subdivided to nitrate and non-nitrate.
            // The nitrate ratio in the petiole N is computed by calling function
            // petiole_nitrate(). Note that the nitrate fraction is more available
            // for redistribution.
            //
            // ratio of NO3 N to total N in petioles.
            let rpetno3 = self.nitrogen_petiole_nitrate();
            // components of reserve N in petioles, for non-NO3 and NO3 origin, respectively.
            let mut petrs1 =
                unsafe { PetioleNitrogen * (1. - rpetno3) - vpetnmin * TotalPetioleWeight }
                    * MobilizNFractionLeaves;
            if petrs1 < 0. {
                petrs1 = 0.;
            }
            let mut petrs2 = unsafe { PetioleNitrogen * rpetno3 - vpno3min * TotalPetioleWeight }
                * MobilizNFractionLeaves;
            if petrs2 < 0. {
                petrs2 = 0.;
            }
            self.petrs = petrs1 + petrs2;
            unsafe {
                PetioleNitrogen -= self.petrs;
            }
            //  Stem N reserves.
            self.stemrs =
                unsafe { StemNitrogen - vstmnmin * TotalStemWeight } * MobilizNFractionStemRoot;
            if self.stemrs < 0. {
                self.stemrs = 0.;
            }
            unsafe {
                StemNitrogen -= self.stemrs;
            }
            //  Root N reserves
            self.rootrs =
                unsafe { RootNitrogen - vrtnmin * TotalRootWeight } * MobilizNFractionStemRoot;
            if self.rootrs < 0. {
                self.rootrs = 0.;
            }
            unsafe {
                RootNitrogen -= self.rootrs;
            }
            //  Burr N reserves
            if unsafe { BurrWeightGreenBolls > 0. } {
                self.burres = unsafe { BurrNitrogen - vburnmin * BurrWeightGreenBolls }
                    * MobilizNFractionBurrs;
                if self.burres < 0. {
                    self.burres = 0.;
                }
                unsafe {
                    BurrNitrogen -= self.burres;
                }
            } else {
                self.burres = 0.;
            }
            // The total reserves, resn, are added to the amount taken up from the soil, for computing npool.
            // Note that N of seeds or squares is not available for redistribution in the plant.
            //
            // total reserve N, in g per plant.
            let resn = self.leafrs + self.petrs + self.stemrs + self.rootrs + self.burres;
            self.npool = self.uptn + resn;
        }
    }
    /// This function computes the allocation of supplied nitrogen to the plant parts.
    fn nitrogen_allocation(&mut self) {
        //     The following constant parameters are used:
        const vseednmax: f64 = 0.70; // maximum proportion of N pool that can be added to seeds
        const vsqrnmax: f64 = 0.65; // maximum proportion of N pool that can be added to squares
        const vburnmax: f64 = 0.65; // maximum proportion of N pool that can be added to burrs
        const vlfnmax: f64 = 0.90; // maximum proportion of N pool that can be added to leaves
        const vstmnmax: f64 = 0.70; // maximum proportion of N pool that can be added to stems
        const vpetnmax: f64 = 0.75; // maximum proportion of N pool that can be added to petioles
                                    //     If total N requirement is less than npool, add N required for growth
                                    //     to the N in
                                    //  each organ, compute added N to vegetative parts, fruiting parts and
                                    //  roots, and compute xtran as the difference between npool and the total N
                                    //  requirements.
        if self.reqtot <= self.npool {
            unsafe {
                LeafNitrogen += self.rqnlef;
                PetioleNitrogen += self.rqnpet;
                StemNitrogen += self.rqnstm;
                RootNitrogen += self.rqnrut;
                SquareNitrogen += self.rqnsqr;
                SeedNitrogen += self.rqnsed;
                BurrNitrogen += self.rqnbur;
            }
            self.addnv = self.rqnlef + self.rqnstm + self.rqnpet;
            self.addnf = self.rqnsqr + self.rqnsed + self.rqnbur;
            self.addnr = self.rqnrut;
            self.xtran = self.npool - self.reqtot;
            return;
        }
        //     If N requirement is greater than npool, execute the following:
        //     First priority is nitrogen supply to the growing seeds. It is assumed
        //     that up to
        //  vseednmax = 0.70 of the supplied N can be used by the seeds. Update seed
        //  N and addnf by the amount of nitrogen used for seed growth, and decrease
        //  npool by this amount. The same procedure is used for each organ,
        //  consecutively.
        let mut useofn; // amount of nitrogen used in growth of a plant organ.
        if self.rqnsed > 0. {
            useofn = fmin(vseednmax * self.npool, self.rqnsed);
            unsafe {
                SeedNitrogen += useofn;
            }
            self.addnf += useofn;
            self.npool -= useofn;
        }
        //     Next priority is for burrs, which can use N up to vburnmax = 0.65 of
        //     the
        //  remaining N pool, and for squares, which can use N up to vsqrnmax = 0.65
        if self.rqnbur > 0. {
            useofn = fmin(vburnmax * self.npool, self.rqnbur);
            unsafe {
                BurrNitrogen += useofn;
            }
            self.addnf += useofn;
            self.npool -= useofn;
        }
        if self.rqnsqr > 0. {
            useofn = fmin(vsqrnmax * self.npool, self.rqnsqr);
            unsafe {
                SquareNitrogen += useofn;
            }
            self.addnf += useofn;
            self.npool -= useofn;
        }
        //     Next priority is for leaves, which can use N up to vlfnmax = 0.90
        //  of the remaining N pool, for stems, up to vstmnmax = 0.70, and for
        //  petioles, up to vpetnmax = 0.75
        if self.rqnlef > 0. {
            useofn = fmin(vlfnmax * self.npool, self.rqnlef);
            unsafe {
                LeafNitrogen += useofn;
            }
            self.addnv += useofn;
            self.npool -= useofn;
        }
        if self.rqnstm > 0. {
            useofn = fmin(vstmnmax * self.npool, self.rqnstm);
            unsafe {
                StemNitrogen += useofn;
            }
            self.addnv += useofn;
            self.npool -= useofn;
        }
        if self.rqnpet > 0. {
            useofn = fmin(vpetnmax * self.npool, self.rqnpet);
            unsafe {
                PetioleNitrogen += useofn;
            }
            self.addnv += useofn;
            self.npool -= useofn;
        }
        //     The remaining npool goes to root growth. If any npool remains
        //  it is defined as xtran.
        if self.rqnrut > 0. {
            useofn = fmin(self.npool, self.rqnrut);
            unsafe {
                RootNitrogen += useofn;
            }
            self.addnr += useofn;
            self.npool -= useofn;
        }
        self.xtran = self.npool;
        if self.xtran > 0. {
            // computes the further allocation of N in the plant
            self.extra_nitrogen_allocation();
        }
    }
    /// This function computes the allocation of extra nitrogen to the plant parts.
    /// It is called from [PlantNitrogen::plant_nitrogen()] if there is a non-zero [xtran].
    fn extra_nitrogen_allocation(&mut self) {
        // If there are any N reserves in the plant, allocate remaining xtran in proportion to the N reserves in each of these organs.
        // Note: all reserves are in g per plant units.
        // reserve N to be added to the burrs.
        let addbur;
        // reserve N to be added to the leaves.
        let addlfn;
        // reserve N to be added to the petioles.
        let addpetn;
        // reserve N to be added to the roots.
        let addrt;
        // reserve N to be added to the stem.
        let addstm;
        // sum of existing reserve N in plant parts.
        let rsum = self.leafrs + self.petrs + self.stemrs + self.rootrs + self.burres;
        if rsum > 0. {
            addlfn = self.xtran * self.leafrs / rsum;
            addpetn = self.xtran * self.petrs / rsum;
            addstm = self.xtran * self.stemrs / rsum;
            addrt = self.xtran * self.rootrs / rsum;
            addbur = self.xtran * self.burres / rsum;
        } else {
            //     If there are no reserves, allocate xtran in proportion to the dry
            //  weights in each of these organs.
            // weight of vegetative plant parts, plus burrs.
            let vegwt = unsafe {
                TotalLeafWeight()
                    + TotalPetioleWeight
                    + TotalStemWeight
                    + TotalRootWeight
                    + BurrWeightGreenBolls
            };
            addlfn = self.xtran * unsafe { TotalLeafWeight() } / vegwt;
            addpetn = self.xtran * unsafe { TotalPetioleWeight } / vegwt;
            addstm = self.xtran * unsafe { TotalStemWeight } / vegwt;
            addrt = self.xtran * unsafe { TotalRootWeight } / vegwt;
            addbur = self.xtran * unsafe { BurrWeightGreenBolls } / vegwt;
        }
        //     Update N content in these plant parts. Note that at this stage of
        //  nitrogen allocation, only vegetative parts and burrs are updated (not
        //  seeds or squares).
        unsafe {
            LeafNitrogen += addlfn;
            PetioleNitrogen += addpetn;
            StemNitrogen += addstm;
            RootNitrogen += addrt;
            BurrNitrogen += addbur;
        }
    }
    // This function computes the concentrations of nitrogen in the dry matter of the plant parts.
    fn nitrogen_content(&mut self) {
        // The following constant parameter is used:
        const seedratio: f64 = 0.64;
        //     Compute N concentration in plant organs as the ratio of N content to
        //     weight of dry matter.
        unsafe {
            if TotalLeafWeight() > 0.00001 {
                LeafNConc = LeafNitrogen / TotalLeafWeight();
            }
            if TotalPetioleWeight > 0.00001 {
                PetioleNConc = PetioleNitrogen / TotalPetioleWeight;
                PetioleNO3NConc = PetioleNConc * self.nitrogen_petiole_nitrate();
            }
            if TotalStemWeight > 0. {
                StemNConc = StemNitrogen / TotalStemWeight;
            }
            if TotalRootWeight > 0. {
                RootNConc = RootNitrogen / TotalRootWeight;
            }
            if TotalSquareWeight > 0. {
                SquareNConc = SquareNitrogen / TotalSquareWeight;
            }
            // weight of seeds in green and mature bolls.
            let xxseed = CottonWeightOpenBolls * (1. - Gintot) + CottonWeightGreenBolls * seedratio;
            if xxseed > 0. {
                SeedNConc = SeedNitrogen / xxseed;
            }
            // weight of burrs in green and mature bolls.
            let xxbur = BurrWeightOpenBolls + BurrWeightGreenBolls;
            if xxbur > 0. {
                BurrNConc = BurrNitrogen / xxbur;
            }
        }
    }
    /// This function computes the nitrogen stress factors.
    fn nitrogen_stress(&mut self) {
        unsafe {
            //     Set the default values for the nitrogen stress coefficients to 1.
            NStressVeg = 1.;
            NStressRoots = 1.;
            NStressFruiting = 1.;
            NitrogenStress = 1.;
            //     Compute the nitrogen stress coefficients. NStressFruiting is the
            //     ratio of
            //  N added actually to the fruits, to their N requirements. NStressVeg is
            //  the same for vegetative shoot growth, and NStressRoots for roots. Also,
            //  an average stress coefficient for vegetative and reproductive organs is
            //  computed as NitrogenStress.
            //     Each stress coefficient has a value between 0 and 1.
            if self.reqf > 0. {
                NStressFruiting = self.addnf / self.reqf;
                if NStressFruiting > 1. {
                    NStressFruiting = 1.;
                }
                if NStressFruiting < 0. {
                    NStressFruiting = 0.;
                }
            }
            if self.reqv > 0. {
                NStressVeg = self.addnv / self.reqv;
                if NStressVeg > 1. {
                    NStressVeg = 1.;
                }
                if NStressVeg < 0. {
                    NStressVeg = 0.;
                }
            }
            if self.rqnrut > 0. {
                NStressRoots = self.addnr / self.rqnrut;
                if NStressRoots > 1. {
                    NStressRoots = 1.;
                }
                if NStressRoots < 0. {
                    NStressRoots = 0.;
                }
            }
            if (self.reqf + self.reqv) > 0. {
                NitrogenStress = (self.addnf + self.addnv) / (self.reqf + self.reqv);
                if NitrogenStress > 1. {
                    NitrogenStress = 1.;
                }
                if NitrogenStress < 0. {
                    NitrogenStress = 0.;
                }
            }
        }
    }
    /// This function computes TotalRequiredN, the nitrogen requirements of the plant - to be used for simulating the N uptake from the soil (in function NitrogenUptake() ) in the next day.
    fn nitrogen_uptake_requirement(&mut self) {
        // The following constant parameters are used:
        const seedcn1: f64 = 0.045; // further requirement for existing seed tissue.
        const seedratio: f64 = 0.64; // the ratio of seeds to seedcotton in green bolls.
        const vnreqlef: f64 = 0.042; // coefficient for computing N uptake requirements of leaves
        const vnreqpet: f64 = 0.036; // coefficient for computing N uptake requirements of petioles
        const vnreqstm: f64 = 0.012; // coefficient for computing N uptake requirements of stems
        const vnreqrt: f64 = 0.010; // coefficient for computing N uptake requirements of roots
        const vnreqsqr: f64 = 0.024; // coefficient for computing N uptake requirements of squares
        const vnreqbur: f64 = 0.012; // coefficient for computing N uptake requirements of burrs
        const voldstm: i32 = 32; // active stem tissue growth period, calendar days.
        unsafe {
            TotalRequiredN = self.reqtot;
            // After the requirements of today's growth are supplied, N is also required for supplying necessary functions in other active plant tissues.
            // Add nitrogen uptake required for leaf and petiole tissue to TotalRequiredN.
            if LeafNConc < vnreqlef {
                TotalRequiredN += TotalLeafWeight() * (vnreqlef - LeafNConc);
            }
            if PetioleNConc < vnreqpet {
                TotalRequiredN += TotalPetioleWeight * (vnreqpet - PetioleNConc);
            }
            // The active stem tissue is the stem formed during the last voldstm days (32 calendar days). add stem requirement to TotalRequiredN.
            //
            // day (from emergence) of oldest actively growing stem tissue.
            let kkday = Kday - voldstm;
            // weight of actively growing stems.
            let grstmwt = if kkday < 1 {
                TotalStemWeight
            } else {
                TotalStemWeight - StemWeight[kkday as usize]
            };
            if StemNConc < vnreqstm {
                TotalRequiredN += grstmwt * (vnreqstm - StemNConc);
            }
            // Compute nitrogen uptake requirement for existing tissues of roots, squares, and seeds and burrs of green bolls. Add it to TotalRequiredN.
            if RootNConc < vnreqrt {
                TotalRequiredN += TotalRootWeight * (vnreqrt - RootNConc);
            }
            if SquareNConc < vnreqsqr {
                TotalRequiredN += TotalSquareWeight * (vnreqsqr - SquareNConc);
            }
            if SeedNConc < seedcn1 {
                TotalRequiredN += CottonWeightGreenBolls * seedratio * (seedcn1 - SeedNConc);
            }
            if BurrNConc < vnreqbur {
                TotalRequiredN += BurrWeightGreenBolls * (vnreqbur - BurrNConc);
            }
        }
    }
    /// This function computes the ratio of NO3 nitrogen to total N in the petioles.
    fn nitrogen_petiole_nitrate(&mut self) -> f64 {
        // The following constant parameters are used:
        //
        // the maximum ratio (of NO3 to total N in petioles).
        const p1: f64 = 0.96;
        // the rate of decline of this ratio with age.
        const p2: f64 = 0.015;
        // the minimum ratio
        const p3: f64 = 0.02;
        // The ratio of NO3 to total N in each individual petiole is computed
        // as a linear function of leaf age. It is assumed that this ratio is
        // maximum for young leaves and is declining with leaf age.
        let mut numl = unsafe { NumPreFruNodes }; // number of petioles computed.
        let mut spetno3 = 0.; // sum of petno3r.
        let mut petno3r; // ratio of NO3 to total N in an individual petiole.
        for j in 0..unsafe { NumPreFruNodes } as usize {
            petno3r = p1 - unsafe { AgeOfPreFruNode[j] } * p2;
            if petno3r < p3 {
                petno3r = p3;
            }
            spetno3 += petno3r;
        }
        // Loop of all the other leaves, with the same computations.
        for k in 0..unsafe { NumVegBranches } as usize {
            let nbrch = unsafe { NumFruitBranches[k] } as usize;
            for l in 0..nbrch {
                let nnid = unsafe { NumNodes[k][l] } as usize;
                numl += nnid as i32;
                for m in 0..nnid {
                    numl += 1;
                    petno3r = p1 - unsafe { LeafAge[k][l][m] } * p2;
                    if petno3r < p3 {
                        petno3r = p3;
                    }
                    spetno3 += petno3r;
                }
            }
        }
        // The return value of the function is the average ratio of NO3 to total N for all the petioles in the plant.
        spetno3 / numl as f64
    }
    /// This function calculates the nitrogen balance in the cotton plant, for diagnostic purposes. It is called from SimulateThisDay().
    ///
    ///  Units are g per plant.
    pub fn balance(&mut self) {
        // addn is the cumulative supply of nitrogen from the soil, including N content of the seedling at emergence.
        /*
        static double addn;
        if (Kday == 1) {addn = RootNitrogen + StemNitrogen + LeafNitrogen;}
        addn += SupplyNO3N + SupplyNH4N;
        double plantn;  // total nitrogen in plant.
        plantn = RootNitrogen + StemNitrogen + LeafNitrogen + PetioleNitrogen +
                SeedNitrogen + BurrNitrogen + SquareNitrogen;
        // CumPlantNLoss is the amount lost in abscised fruit and leaves and dead
        // roots.
        double balpn;  // the plant nitrogen balance, which should be zero.
        balpn = addn - plantn - CumPlantNLoss;*/
    }
}
