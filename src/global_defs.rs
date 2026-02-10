// Auto-generated from bindgen extern statics in global.h.
// Defines all global variables in Rust so C++ global.cpp is no longer required.
#[allow(non_upper_case_globals)]
#[allow(non_snake_case)]
mod global_rust_defs {
    #[allow(unused_imports)]
    use crate::{Climstruct, Irrigation};

    #[no_mangle]
    pub static mut Clim: [Climstruct; 400usize] = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut Irrig: [Irrigation; 150usize] = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut DayEmerge: ::std::os::raw::c_int = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut DayEndMulch: ::std::os::raw::c_int = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut DayFinish: ::std::os::raw::c_int = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut DayFirstDef: ::std::os::raw::c_int = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut DayOfSimulation: ::std::os::raw::c_int = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut Daynum: ::std::os::raw::c_int = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut DayPlant: ::std::os::raw::c_int = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut DayStart: ::std::os::raw::c_int = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut DayStartMulch: ::std::os::raw::c_int = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut DayStartPredIrrig: ::std::os::raw::c_int = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut DayStopPredIrrig: ::std::os::raw::c_int = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut FirstBloom: ::std::os::raw::c_int = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut FirstSquare: ::std::os::raw::c_int = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut IrrigMethod: ::std::os::raw::c_int = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut isw: ::std::os::raw::c_int = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut iyear: ::std::os::raw::c_int = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut Kday: ::std::os::raw::c_int = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut KdayAdjust: ::std::os::raw::c_int = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut LastIrrigation: ::std::os::raw::c_int = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut LastTaprootLayer: ::std::os::raw::c_int = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut LocationColumnDrip: ::std::os::raw::c_int = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut LocationLayerDrip: ::std::os::raw::c_int = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut MainStemNodes: ::std::os::raw::c_int = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut MinDaysBetweenIrrig: ::std::os::raw::c_int = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut MulchIndicator: ::std::os::raw::c_int = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut nk: ::std::os::raw::c_int = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut nl: ::std::os::raw::c_int = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut noitr: ::std::os::raw::c_int = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut NumAbscisedLeaves: ::std::os::raw::c_int = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut NumAdjustDays: ::std::os::raw::c_int = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut NumFruitSites: ::std::os::raw::c_int = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut NumIrrigations: ::std::os::raw::c_int = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut NumLayersWithRoots: ::std::os::raw::c_int = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut NumPreFruNodes: ::std::os::raw::c_int = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut NumRootAgeGroups: ::std::os::raw::c_int = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut NumSheddingTags: ::std::os::raw::c_int = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut NumVegBranches: ::std::os::raw::c_int = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut PlantRowColumn: ::std::os::raw::c_int = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut WaterTableLayer: ::std::os::raw::c_int = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut DefoliationDate: [::std::os::raw::c_int; 5usize] = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut DefoliationMethod: [::std::os::raw::c_int; 5usize] =
        unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut FruitingCode: [[[::std::os::raw::c_int; 5usize]; 30usize]; 3usize] =
        unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut LateralRootFlag: [::std::os::raw::c_int; 40usize] =
        unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut NumFruitBranches: [::std::os::raw::c_int; 3usize] =
        unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut NumNodes: [[::std::os::raw::c_int; 30usize]; 3usize] =
        unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut NodeLayer: [[::std::os::raw::c_int; 30usize]; 3usize] =
        unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut NodeLayerPreFru: [::std::os::raw::c_int; 9usize] = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut RootColNumLeft: [::std::os::raw::c_int; 40usize] = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut RootColNumRight: [::std::os::raw::c_int; 40usize] =
        unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut SoilHorizonNum: [::std::os::raw::c_int; 40usize] = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut bEnd: bool = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut bPollinSwitch: bool = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut nadj: [bool; 5usize] = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut AbscisedFruitSites: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut AbscisedLeafWeight: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut ActualBollGrowth: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut ActualBurrGrowth: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut ActualSquareGrowth: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut ActualSoilEvaporation: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut ActualStemGrowth: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut ActualTranspiration: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut addwtbl: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut AdjAddHeightRate: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut AdjAddMSNodesRate: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut AdjAddSitesRate: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut AdjGreenBollAbsc: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut AdjSquareAbsc: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut AppliedWater: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut AverageLwp: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut AverageLwpMin: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut AverageSoilPsi: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut AvrgDailyTemp: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut BloomWeightLoss: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut BurrNConc: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut BurrNitrogen: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut BurrWeightGreenBolls: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut BurrWeightOpenBolls: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut CarbonAllocatedForRootGrowth: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut CarbonStress: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut conmax: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut CottonWeightGreenBolls: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut CottonWeightOpenBolls: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut CumEvaporation: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut CumFertilizerN: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut CumNetPhotosynth: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut CumNitrogenUptake: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut CumPlantNLoss: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut CumTranspiration: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut CumWaterAdded: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut CumWaterDrained: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut DailyRootLoss: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut DayInc: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut DayTimeTemp: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut dclay: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut DeepSoilTemperature: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut DensityFactor: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut DepthLastRootLayer: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut dsand: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut ElCondSatSoilToday: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut ExtraCarbon: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut FruitGrowthRatio: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut ginp: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut Gintot: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut GreenBollsLost: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut InitialTotalSoilWater: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut IrrigationDepth: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut LeafAreaIndex: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut LeafNConc: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut LeafNitrogen: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut LeafWeightAreaRatio: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut LightIntercept: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut LintYield: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut LwpMax: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut LwpMin: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut MaxIrrigation: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut MineralizedOrganicN: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut MulchTranLW: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut MulchTranSW: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut NetPhotosynthesis: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut NightTimeTemp: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut NitrogenStress: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut NStressFruiting: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut NStressRoots: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut NStressVeg: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut NumGreenBolls: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut NumOpenBolls: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut NumSquares: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut PercentDefoliation: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut PerPlantArea: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut PetioleNConc: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut PetioleNitrogen: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut PetioleNO3NConc: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut pixcon: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut pixda: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut pixdn: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut pixdz: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut PixInPlants: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut PlantHeight: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut PlantPopulation: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut PlantRowLocation: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut PlantWeight: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut PlantWeightAtStart: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut PotGroAllBolls: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut PotGroAllBurrs: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut PotGroAllLeaves: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut PotGroAllPetioles: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut PotGroAllRoots: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut PotGroAllSquares: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut PotGroStem: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut RatioImplicit: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut ReferenceTransp: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut ReserveC: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut Rn: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut RootNConc: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut RootNitrogen: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut RootWeightLoss: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut RowSpace: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut SeedNConc: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut SeedNitrogen: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut SoilNitrogenAtStart: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut SoilNitrogenLoss: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut SquareNConc: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut SquareNitrogen: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut StemNConc: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut StemNitrogen: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut SumNO3N90: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut SupplyNH4N: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut SupplyNO3N: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut TapRootLength: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut TotalActualLeafGrowth: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut TotalActualPetioleGrowth: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut TotalPetioleWeight: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut TotalRequiredN: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut TotalRootWeight: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut TotalSoilNh4N: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut TotalSoilNitrogen: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut TotalSoilNo3N: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut TotalSoilUreaN: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut TotalSoilWater: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut TotalSquareWeight: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut TotalStemWeight: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut WaterStress: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut WaterStressStem: f64 = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut AbscissionLag: [f64; 20usize] = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut ActualRootGrowth: [[f64; 20usize]; 40usize] = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut AgeOfBoll: [[[f64; 5usize]; 30usize]; 3usize] = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut AgeOfPreFruNode: [f64; 9usize] = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut AgeOfSite: [[[f64; 5usize]; 30usize]; 3usize] = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut airdr: [f64; 9usize] = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut AirTemp: [f64; 24usize] = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut albedo: [f64; 24usize] = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut alpha: [f64; 9usize] = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut AvrgNodeTemper: [[[f64; 5usize]; 30usize]; 3usize] =
        unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut AverageLeafAge: [f64; 20usize] = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut beta: [f64; 9usize] = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut BollWeight: [[[f64; 5usize]; 30usize]; 3usize] = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut BulkDensity: [f64; 9usize] = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut BurrWeight: [[[f64; 5usize]; 30usize]; 3usize] = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut cgind: [f64; 3usize] = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut ClayVolumeFraction: [f64; 40usize] = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut CloudCoverRatio: [f64; 24usize] = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut CloudTypeCorr: [f64; 24usize] = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut DefoliantAppRate: [f64; 5usize] = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut DelayNewFruBranch: [f64; 3usize] = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut DelayNewNode: [[f64; 30usize]; 3usize] = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut DewPointTemp: [f64; 24usize] = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut dl: [f64; 40usize] = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut es1hour: [f64; 24usize] = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut es2hour: [f64; 24usize] = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut FieldCapacity: [f64; 40usize] = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut FoliageTemp: [f64; 20usize] = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut FreshOrganicMatter: [[f64; 20usize]; 40usize] = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut FreshOrganicNitrogen: [[f64; 20usize]; 40usize] = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut FruitFraction: [[[f64; 5usize]; 30usize]; 3usize] =
        unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut HeatCapacitySoilSolid: [f64; 40usize] = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut HeatCondDrySoil: [f64; 40usize] = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut HumusNitrogen: [[f64; 20usize]; 40usize] = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut HumusOrganicMatter: [[f64; 20usize]; 40usize] = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut LeafAge: [[[f64; 5usize]; 30usize]; 3usize] = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut LeafAreaMainStem: [[f64; 30usize]; 3usize] = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut LeafArea: [f64; 20usize] = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut LeafAreaIndexes: [f64; 20usize] = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut LeafAreaNodes: [[[f64; 5usize]; 30usize]; 3usize] =
        unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut LeafAreaPreFru: [f64; 9usize] = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut LeafWeightLayer: [f64; 20usize] = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut LeafWeightMainStem: [[f64; 30usize]; 3usize] = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut LeafWeightNodes: [[[f64; 5usize]; 30usize]; 3usize] =
        unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut LeafWeightPreFru: [f64; 9usize] = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut LightInterceptLayer: [f64; 20usize] = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut LwpMinX: [f64; 3usize] = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut LwpX: [f64; 3usize] = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut MarginalWaterContent: [f64; 40usize] = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut MaxWaterCapacity: [f64; 40usize] = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut MulchTemp: [f64; 20usize] = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut NO3FlowFraction: [f64; 40usize] = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut PetioleWeightMainStem: [[f64; 30usize]; 3usize] = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut PetioleWeightNodes: [[[f64; 5usize]; 30usize]; 3usize] =
        unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut PetioleWeightPreFru: [f64; 9usize] = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut PoreSpace: [f64; 40usize] = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut PotGroBolls: [[[f64; 5usize]; 30usize]; 3usize] = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut PotGroBurrs: [[[f64; 5usize]; 30usize]; 3usize] = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut PotGroLeafAreaMainStem: [[f64; 30usize]; 3usize] = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut PotGroLeafAreaNodes: [[[f64; 5usize]; 30usize]; 3usize] =
        unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut PotGroLeafAreaPreFru: [f64; 9usize] = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut PotGroLeafWeightMainStem: [[f64; 30usize]; 3usize] =
        unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut PotGroLeafWeightNodes: [[[f64; 5usize]; 30usize]; 3usize] =
        unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut PotGroLeafWeightPreFru: [f64; 9usize] = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut PotGroPetioleWeightMainStem: [[f64; 30usize]; 3usize] =
        unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut PotGroPetioleWeightNodes: [[[f64; 5usize]; 30usize]; 3usize] =
        unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut PotGroPetioleWeightPreFru: [f64; 9usize] = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut PotGroRoots: [[f64; 20usize]; 40usize] = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut PotGroSquares: [[[f64; 5usize]; 30usize]; 3usize] =
        unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut Radiation: [f64; 24usize] = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut ReferenceETP: [f64; 24usize] = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut RelativeHumidity: [f64; 24usize] = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut rlat1: [f64; 40usize] = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut rlat2: [f64; 40usize] = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut RootAge: [[f64; 20usize]; 40usize] = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut RootGroFactor: [[f64; 20usize]; 40usize] = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut RootImpede: [[f64; 20usize]; 40usize] = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut RootWeight: [[[f64; 3usize]; 20usize]; 40usize] = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut RootWtCapblUptake: [[f64; 20usize]; 40usize] = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut SandVolumeFraction: [f64; 40usize] = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut SaturatedHydCond: [f64; 9usize] = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut ShedByCarbonStress: [f64; 20usize] = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut ShedByNitrogenStress: [f64; 20usize] = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut ShedByWaterStress: [f64; 20usize] = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut SoilPsi: [[f64; 20usize]; 40usize] = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut SoilTemp: [[f64; 20usize]; 40usize] = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut SoilTempDailyAvrg: [[f64; 20usize]; 40usize] = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut SquareWeight: [[[f64; 5usize]; 30usize]; 3usize] = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut StemWeight: [f64; 365usize] = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut thad: [f64; 40usize] = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut thetar: [f64; 40usize] = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut thetas: [f64; 9usize] = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut thts: [f64; 40usize] = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut VarPar: [f64; 61usize] = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut VolNh4NContent: [[f64; 20usize]; 40usize] = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut VolNo3NContent: [[f64; 20usize]; 40usize] = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut VolUreaNContent: [[f64; 20usize]; 40usize] = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut VolWaterContent: [[f64; 20usize]; 40usize] = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut WindSpeed: [f64; 24usize] = unsafe { std::mem::zeroed() };

    #[no_mangle]
    pub static mut wk: [f64; 20usize] = unsafe { std::mem::zeroed() };

    #[export_name = "\u{1}__Z15TotalLeafWeightv"]
    pub unsafe extern "C" fn total_leaf_weight_cpp() -> f64 {
        let mut result = 0.0;
        if crate::FirstSquare <= 0 {
            result += 0.2;
        }

        for i in 0..crate::NumPreFruNodes as usize {
            result += crate::LeafWeightPreFru[i];
        }

        for k in 0..crate::NumVegBranches as usize {
            for l in 0..crate::NumFruitBranches[k] as usize {
                result += crate::LeafWeightMainStem[k][l];
                for m in 0..crate::NumNodes[k][l] as usize {
                    result += crate::LeafWeightNodes[k][l][m];
                }
            }
        }

        result
    }

    #[export_name = "\u{1}__Z13TotalLeafAreav"]
    pub unsafe extern "C" fn total_leaf_area_cpp() -> f64 {
        let mut result = 0.0;
        if crate::FirstSquare <= 0 {
            result += 0.20 * 0.6;
        }

        for i in 0..crate::NumPreFruNodes as usize {
            result += crate::LeafAreaPreFru[i];
        }

        for k in 0..crate::NumVegBranches as usize {
            for l in 0..crate::NumFruitBranches[k] as usize {
                result += crate::LeafAreaMainStem[k][l];
                for m in 0..crate::NumNodes[k][l] as usize {
                    result += crate::LeafAreaNodes[k][l][m];
                }
            }
        }

        result
    }
}
