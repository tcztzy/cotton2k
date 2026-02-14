// Canonical Rust definitions for shared global state used across the model.
// This replaces the previous generated FFI declarations.
#[allow(non_upper_case_globals)]
#[allow(non_snake_case)]
mod global_rust_defs {
    use crate::model_state::for_each_fruiting_site;
    use std::sync::{LazyLock, RwLock};

    #[allow(unused_imports)]
    use crate::{Climstruct, Irrigation};

    pub static Clim: LazyLock<RwLock<[Climstruct; 400usize]>> = LazyLock::new(|| {
        RwLock::new(
            [Climstruct {
                nDay: 0,
                Rad: 0.0,
                Tmax: 0.0,
                Tmin: 0.0,
                Rain: 0.0,
                Wind: 0.0,
                Tdew: 0.0,
            }; 400usize],
        )
    });
    pub static Irrig: LazyLock<RwLock<[Irrigation; 150usize]>> = LazyLock::new(|| {
        RwLock::new(
            [Irrigation {
                day: 0,
                method: 0,
                LocationColumnDrip: 0,
                LocationLayerDrip: 0,
                amount: 0.0,
            }; 150usize],
        )
    });
    pub static mut DayEmerge: ::std::os::raw::c_int = 0;
    pub static mut DayEndMulch: ::std::os::raw::c_int = 0;
    pub static mut DayFinish: ::std::os::raw::c_int = 0;
    pub static mut DayFirstDef: ::std::os::raw::c_int = 0;
    pub static mut DayOfSimulation: ::std::os::raw::c_int = 0;
    pub static mut Daynum: ::std::os::raw::c_int = 0;
    pub static mut DayPlant: ::std::os::raw::c_int = 0;
    pub static mut DayStart: ::std::os::raw::c_int = 0;
    pub static mut DayStartMulch: ::std::os::raw::c_int = 0;
    pub static mut DayStartPredIrrig: ::std::os::raw::c_int = 0;
    pub static mut DayStopPredIrrig: ::std::os::raw::c_int = 0;
    pub static mut FirstBloom: ::std::os::raw::c_int = 0;
    pub static mut FirstSquare: ::std::os::raw::c_int = 0;
    pub static mut IrrigMethod: ::std::os::raw::c_int = 0;
    pub static mut isw: ::std::os::raw::c_int = 0;
    pub static mut iyear: ::std::os::raw::c_int = 0;
    pub static mut Kday: ::std::os::raw::c_int = 0;
    pub static mut KdayAdjust: ::std::os::raw::c_int = 0;
    pub static mut LastIrrigation: ::std::os::raw::c_int = 0;
    pub static mut LastTaprootLayer: ::std::os::raw::c_int = 0;
    pub static mut LocationColumnDrip: ::std::os::raw::c_int = 0;
    pub static mut LocationLayerDrip: ::std::os::raw::c_int = 0;
    pub static mut MainStemNodes: ::std::os::raw::c_int = 0;
    pub static mut MinDaysBetweenIrrig: ::std::os::raw::c_int = 0;
    pub static mut MulchIndicator: ::std::os::raw::c_int = 0;
    pub static mut nk: ::std::os::raw::c_int = 0;
    pub static mut nl: ::std::os::raw::c_int = 0;
    pub static mut noitr: ::std::os::raw::c_int = 0;
    pub static mut NumAbscisedLeaves: ::std::os::raw::c_int = 0;
    pub static mut NumAdjustDays: ::std::os::raw::c_int = 0;
    pub static mut NumFruitSites: ::std::os::raw::c_int = 0;
    pub static mut NumIrrigations: ::std::os::raw::c_int = 0;
    pub static mut NumLayersWithRoots: ::std::os::raw::c_int = 0;
    pub static mut NumPreFruNodes: ::std::os::raw::c_int = 0;
    pub static mut NumRootAgeGroups: ::std::os::raw::c_int = 0;
    pub static mut NumSheddingTags: ::std::os::raw::c_int = 0;
    pub static mut NumVegBranches: ::std::os::raw::c_int = 0;
    pub static mut PlantRowColumn: ::std::os::raw::c_int = 0;
    pub static mut WaterTableLayer: ::std::os::raw::c_int = 0;
    pub static mut DefoliationDate: [::std::os::raw::c_int; 5usize] = [0; 5usize];
    pub static mut DefoliationMethod: [::std::os::raw::c_int; 5usize] = [0; 5usize];
    pub static mut FruitingCode: [[[::std::os::raw::c_int; 5usize]; 30usize]; 3usize] =
        [[[0; 5usize]; 30usize]; 3usize];
    pub static mut LateralRootFlag: [::std::os::raw::c_int; 40usize] = [0; 40usize];
    pub static mut NumFruitBranches: [::std::os::raw::c_int; 3usize] = [0; 3usize];
    pub static mut NumNodes: [[::std::os::raw::c_int; 30usize]; 3usize] = [[0; 30usize]; 3usize];
    pub static mut NodeLayer: [[::std::os::raw::c_int; 30usize]; 3usize] = [[0; 30usize]; 3usize];
    pub static mut NodeLayerPreFru: [::std::os::raw::c_int; 9usize] = [0; 9usize];
    pub static mut RootColNumLeft: [::std::os::raw::c_int; 40usize] = [0; 40usize];
    pub static mut RootColNumRight: [::std::os::raw::c_int; 40usize] = [0; 40usize];
    pub static mut SoilHorizonNum: [::std::os::raw::c_int; 40usize] = [0; 40usize];
    pub static mut bEnd: bool = false;
    pub static mut bPollinSwitch: bool = false;
    pub static mut nadj: [bool; 5usize] = [false; 5usize];
    pub static mut AbscisedFruitSites: f64 = 0.0;
    pub static mut AbscisedLeafWeight: f64 = 0.0;
    pub static mut ActualBollGrowth: f64 = 0.0;
    pub static mut ActualBurrGrowth: f64 = 0.0;
    pub static mut ActualSquareGrowth: f64 = 0.0;
    pub static mut ActualSoilEvaporation: f64 = 0.0;
    pub static mut ActualStemGrowth: f64 = 0.0;
    pub static mut ActualTranspiration: f64 = 0.0;
    pub static mut addwtbl: f64 = 0.0;
    pub static mut AdjAddHeightRate: f64 = 0.0;
    pub static mut AdjAddMSNodesRate: f64 = 0.0;
    pub static mut AdjAddSitesRate: f64 = 0.0;
    pub static mut AdjGreenBollAbsc: f64 = 0.0;
    pub static mut AdjSquareAbsc: f64 = 0.0;
    pub static mut AppliedWater: f64 = 0.0;
    pub static mut AverageLwp: f64 = 0.0;
    pub static mut AverageLwpMin: f64 = 0.0;
    pub static mut AverageSoilPsi: f64 = 0.0;
    pub static mut AvrgDailyTemp: f64 = 0.0;
    pub static mut BloomWeightLoss: f64 = 0.0;
    pub static mut BurrNConc: f64 = 0.0;
    pub static mut BurrNitrogen: f64 = 0.0;
    pub static mut BurrWeightGreenBolls: f64 = 0.0;
    pub static mut BurrWeightOpenBolls: f64 = 0.0;
    pub static mut CarbonAllocatedForRootGrowth: f64 = 0.0;
    pub static mut CarbonStress: f64 = 0.0;
    pub static mut conmax: f64 = 0.0;
    pub static mut CottonWeightGreenBolls: f64 = 0.0;
    pub static mut CottonWeightOpenBolls: f64 = 0.0;
    pub static mut CumEvaporation: f64 = 0.0;
    pub static mut CumFertilizerN: f64 = 0.0;
    pub static mut CumNetPhotosynth: f64 = 0.0;
    pub static mut CumNitrogenUptake: f64 = 0.0;
    pub static mut CumPlantNLoss: f64 = 0.0;
    pub static mut CumTranspiration: f64 = 0.0;
    pub static mut CumWaterAdded: f64 = 0.0;
    pub static mut CumWaterDrained: f64 = 0.0;
    pub static mut DailyRootLoss: f64 = 0.0;
    pub static mut DayInc: f64 = 0.0;
    pub static mut DayTimeTemp: f64 = 0.0;
    pub static mut dclay: f64 = 0.0;
    pub static mut DeepSoilTemperature: f64 = 0.0;
    pub static mut DensityFactor: f64 = 0.0;
    pub static mut DepthLastRootLayer: f64 = 0.0;
    pub static mut dsand: f64 = 0.0;
    pub static mut ElCondSatSoilToday: f64 = 0.0;
    pub static mut ExtraCarbon: f64 = 0.0;
    pub static mut FruitGrowthRatio: f64 = 0.0;
    pub static mut ginp: f64 = 0.0;
    pub static mut Gintot: f64 = 0.0;
    pub static mut GreenBollsLost: f64 = 0.0;
    pub static mut InitialTotalSoilWater: f64 = 0.0;
    pub static mut IrrigationDepth: f64 = 0.0;
    pub static mut LeafAreaIndex: f64 = 0.0;
    pub static mut LeafNConc: f64 = 0.0;
    pub static mut LeafNitrogen: f64 = 0.0;
    pub static mut LeafWeightAreaRatio: f64 = 0.0;
    pub static mut LightIntercept: f64 = 0.0;
    pub static mut LintYield: f64 = 0.0;
    pub static mut LwpMax: f64 = 0.0;
    pub static mut LwpMin: f64 = 0.0;
    pub static mut MaxIrrigation: f64 = 0.0;
    pub static mut MineralizedOrganicN: f64 = 0.0;
    pub static mut MulchTranLW: f64 = 0.0;
    pub static mut MulchTranSW: f64 = 0.0;
    pub static mut NetPhotosynthesis: f64 = 0.0;
    pub static mut NightTimeTemp: f64 = 0.0;
    pub static mut NitrogenStress: f64 = 0.0;
    pub static mut NStressFruiting: f64 = 0.0;
    pub static mut NStressRoots: f64 = 0.0;
    pub static mut NStressVeg: f64 = 0.0;
    pub static mut NumGreenBolls: f64 = 0.0;
    pub static mut NumOpenBolls: f64 = 0.0;
    pub static mut NumSquares: f64 = 0.0;
    pub static mut PercentDefoliation: f64 = 0.0;
    pub static mut PerPlantArea: f64 = 0.0;
    pub static mut PetioleNConc: f64 = 0.0;
    pub static mut PetioleNitrogen: f64 = 0.0;
    pub static mut PetioleNO3NConc: f64 = 0.0;
    pub static mut pixcon: f64 = 0.0;
    pub static mut pixda: f64 = 0.0;
    pub static mut pixdn: f64 = 0.0;
    pub static mut pixdz: f64 = 0.0;
    pub static mut PixInPlants: f64 = 0.0;
    pub static mut PlantHeight: f64 = 0.0;
    pub static mut PlantPopulation: f64 = 0.0;
    pub static mut PlantRowLocation: f64 = 0.0;
    pub static mut PlantWeight: f64 = 0.0;
    pub static mut PlantWeightAtStart: f64 = 0.0;
    pub static mut PotGroAllBolls: f64 = 0.0;
    pub static mut PotGroAllBurrs: f64 = 0.0;
    pub static mut PotGroAllLeaves: f64 = 0.0;
    pub static mut PotGroAllPetioles: f64 = 0.0;
    pub static mut PotGroAllRoots: f64 = 0.0;
    pub static mut PotGroAllSquares: f64 = 0.0;
    pub static mut PotGroStem: f64 = 0.0;
    pub static mut RatioImplicit: f64 = 0.0;
    pub static mut ReferenceTransp: f64 = 0.0;
    pub static mut ReserveC: f64 = 0.0;
    pub static mut Rn: f64 = 0.0;
    pub static mut RootNConc: f64 = 0.0;
    pub static mut RootNitrogen: f64 = 0.0;
    pub static mut RootWeightLoss: f64 = 0.0;
    pub static mut RowSpace: f64 = 0.0;
    pub static mut SeedNConc: f64 = 0.0;
    pub static mut SeedNitrogen: f64 = 0.0;
    pub static mut SoilNitrogenAtStart: f64 = 0.0;
    pub static mut SoilNitrogenLoss: f64 = 0.0;
    pub static mut SquareNConc: f64 = 0.0;
    pub static mut SquareNitrogen: f64 = 0.0;
    pub static mut StemNConc: f64 = 0.0;
    pub static mut StemNitrogen: f64 = 0.0;
    pub static mut SumNO3N90: f64 = 0.0;
    pub static mut SupplyNH4N: f64 = 0.0;
    pub static mut SupplyNO3N: f64 = 0.0;
    pub static mut TapRootLength: f64 = 0.0;
    pub static mut TotalActualLeafGrowth: f64 = 0.0;
    pub static mut TotalActualPetioleGrowth: f64 = 0.0;
    pub static mut TotalPetioleWeight: f64 = 0.0;
    pub static mut TotalRequiredN: f64 = 0.0;
    pub static mut TotalRootWeight: f64 = 0.0;
    pub static mut TotalSoilNh4N: f64 = 0.0;
    pub static mut TotalSoilNitrogen: f64 = 0.0;
    pub static mut TotalSoilNo3N: f64 = 0.0;
    pub static mut TotalSoilUreaN: f64 = 0.0;
    pub static mut TotalSoilWater: f64 = 0.0;
    pub static mut TotalSquareWeight: f64 = 0.0;
    pub static mut TotalStemWeight: f64 = 0.0;
    pub static mut WaterStress: f64 = 0.0;
    pub static mut WaterStressStem: f64 = 0.0;
    pub static mut AbscissionLag: [f64; 20usize] = [0.0; 20usize];
    pub static mut ActualRootGrowth: [[f64; 20usize]; 40usize] = [[0.0; 20usize]; 40usize];
    pub static mut AgeOfBoll: [[[f64; 5usize]; 30usize]; 3usize] =
        [[[0.0; 5usize]; 30usize]; 3usize];
    pub static mut AgeOfPreFruNode: [f64; 9usize] = [0.0; 9usize];
    pub static mut AgeOfSite: [[[f64; 5usize]; 30usize]; 3usize] =
        [[[0.0; 5usize]; 30usize]; 3usize];
    pub static mut airdr: [f64; 9usize] = [0.0; 9usize];
    pub static mut AirTemp: [f64; 24usize] = [0.0; 24usize];
    pub static mut albedo: [f64; 24usize] = [0.0; 24usize];
    pub static mut alpha: [f64; 9usize] = [0.0; 9usize];
    pub static mut AvrgNodeTemper: [[[f64; 5usize]; 30usize]; 3usize] =
        [[[0.0; 5usize]; 30usize]; 3usize];
    pub static mut AverageLeafAge: [f64; 20usize] = [0.0; 20usize];
    pub static mut beta: [f64; 9usize] = [0.0; 9usize];
    pub static mut BollWeight: [[[f64; 5usize]; 30usize]; 3usize] =
        [[[0.0; 5usize]; 30usize]; 3usize];
    pub static mut BulkDensity: [f64; 9usize] = [0.0; 9usize];
    pub static mut BurrWeight: [[[f64; 5usize]; 30usize]; 3usize] =
        [[[0.0; 5usize]; 30usize]; 3usize];
    pub static mut cgind: [f64; 3usize] = [0.0; 3usize];
    pub static mut ClayVolumeFraction: [f64; 40usize] = [0.0; 40usize];
    pub static mut CloudCoverRatio: [f64; 24usize] = [0.0; 24usize];
    pub static mut CloudTypeCorr: [f64; 24usize] = [0.0; 24usize];
    pub static mut DefoliantAppRate: [f64; 5usize] = [0.0; 5usize];
    pub static mut DelayNewFruBranch: [f64; 3usize] = [0.0; 3usize];
    pub static mut DelayNewNode: [[f64; 30usize]; 3usize] = [[0.0; 30usize]; 3usize];
    pub static mut DewPointTemp: [f64; 24usize] = [0.0; 24usize];
    pub static mut dl: [f64; 40usize] = [0.0; 40usize];
    pub static mut es1hour: [f64; 24usize] = [0.0; 24usize];
    pub static mut es2hour: [f64; 24usize] = [0.0; 24usize];
    pub static mut FieldCapacity: [f64; 40usize] = [0.0; 40usize];
    pub static mut FoliageTemp: [f64; 20usize] = [0.0; 20usize];
    pub static mut FreshOrganicMatter: [[f64; 20usize]; 40usize] = [[0.0; 20usize]; 40usize];
    pub static mut FreshOrganicNitrogen: [[f64; 20usize]; 40usize] = [[0.0; 20usize]; 40usize];
    pub static mut FruitFraction: [[[f64; 5usize]; 30usize]; 3usize] =
        [[[0.0; 5usize]; 30usize]; 3usize];
    pub static mut HeatCapacitySoilSolid: [f64; 40usize] = [0.0; 40usize];
    pub static mut HeatCondDrySoil: [f64; 40usize] = [0.0; 40usize];
    pub static mut HumusNitrogen: [[f64; 20usize]; 40usize] = [[0.0; 20usize]; 40usize];
    pub static mut HumusOrganicMatter: [[f64; 20usize]; 40usize] = [[0.0; 20usize]; 40usize];
    pub static mut LeafAge: [[[f64; 5usize]; 30usize]; 3usize] = [[[0.0; 5usize]; 30usize]; 3usize];
    pub static mut LeafAreaMainStem: [[f64; 30usize]; 3usize] = [[0.0; 30usize]; 3usize];
    pub static mut LeafArea: [f64; 20usize] = [0.0; 20usize];
    pub static mut LeafAreaIndexes: [f64; 20usize] = [0.0; 20usize];
    pub static mut LeafAreaNodes: [[[f64; 5usize]; 30usize]; 3usize] =
        [[[0.0; 5usize]; 30usize]; 3usize];
    pub static mut LeafAreaPreFru: [f64; 9usize] = [0.0; 9usize];
    pub static mut LeafWeightLayer: [f64; 20usize] = [0.0; 20usize];
    pub static mut LeafWeightMainStem: [[f64; 30usize]; 3usize] = [[0.0; 30usize]; 3usize];
    pub static mut LeafWeightNodes: [[[f64; 5usize]; 30usize]; 3usize] =
        [[[0.0; 5usize]; 30usize]; 3usize];
    pub static mut LeafWeightPreFru: [f64; 9usize] = [0.0; 9usize];
    pub static mut LightInterceptLayer: [f64; 20usize] = [0.0; 20usize];
    pub static mut LwpMinX: [f64; 3usize] = [0.0; 3usize];
    pub static mut LwpX: [f64; 3usize] = [0.0; 3usize];
    pub static mut MarginalWaterContent: [f64; 40usize] = [0.0; 40usize];
    pub static mut MaxWaterCapacity: [f64; 40usize] = [0.0; 40usize];
    pub static mut MulchTemp: [f64; 20usize] = [0.0; 20usize];
    pub static mut NO3FlowFraction: [f64; 40usize] = [0.0; 40usize];
    pub static mut PetioleWeightMainStem: [[f64; 30usize]; 3usize] = [[0.0; 30usize]; 3usize];
    pub static mut PetioleWeightNodes: [[[f64; 5usize]; 30usize]; 3usize] =
        [[[0.0; 5usize]; 30usize]; 3usize];
    pub static mut PetioleWeightPreFru: [f64; 9usize] = [0.0; 9usize];
    pub static mut PoreSpace: [f64; 40usize] = [0.0; 40usize];
    pub static mut PotGroBolls: [[[f64; 5usize]; 30usize]; 3usize] =
        [[[0.0; 5usize]; 30usize]; 3usize];
    pub static mut PotGroBurrs: [[[f64; 5usize]; 30usize]; 3usize] =
        [[[0.0; 5usize]; 30usize]; 3usize];
    pub static mut PotGroLeafAreaMainStem: [[f64; 30usize]; 3usize] = [[0.0; 30usize]; 3usize];
    pub static mut PotGroLeafAreaNodes: [[[f64; 5usize]; 30usize]; 3usize] =
        [[[0.0; 5usize]; 30usize]; 3usize];
    pub static mut PotGroLeafAreaPreFru: [f64; 9usize] = [0.0; 9usize];
    pub static mut PotGroLeafWeightMainStem: [[f64; 30usize]; 3usize] = [[0.0; 30usize]; 3usize];
    pub static mut PotGroLeafWeightNodes: [[[f64; 5usize]; 30usize]; 3usize] =
        [[[0.0; 5usize]; 30usize]; 3usize];
    pub static mut PotGroLeafWeightPreFru: [f64; 9usize] = [0.0; 9usize];
    pub static mut PotGroPetioleWeightMainStem: [[f64; 30usize]; 3usize] = [[0.0; 30usize]; 3usize];
    pub static mut PotGroPetioleWeightNodes: [[[f64; 5usize]; 30usize]; 3usize] =
        [[[0.0; 5usize]; 30usize]; 3usize];
    pub static mut PotGroPetioleWeightPreFru: [f64; 9usize] = [0.0; 9usize];
    pub static mut PotGroRoots: [[f64; 20usize]; 40usize] = [[0.0; 20usize]; 40usize];
    pub static mut PotGroSquares: [[[f64; 5usize]; 30usize]; 3usize] =
        [[[0.0; 5usize]; 30usize]; 3usize];
    pub static mut Radiation: [f64; 24usize] = [0.0; 24usize];
    pub static mut ReferenceETP: [f64; 24usize] = [0.0; 24usize];
    pub static mut RelativeHumidity: [f64; 24usize] = [0.0; 24usize];
    pub static mut rlat1: [f64; 40usize] = [0.0; 40usize];
    pub static mut rlat2: [f64; 40usize] = [0.0; 40usize];
    pub static mut RootAge: [[f64; 20usize]; 40usize] = [[0.0; 20usize]; 40usize];
    pub static mut RootGroFactor: [[f64; 20usize]; 40usize] = [[0.0; 20usize]; 40usize];
    pub static mut RootImpede: [[f64; 20usize]; 40usize] = [[0.0; 20usize]; 40usize];
    pub static mut RootWeight: [[[f64; 3usize]; 20usize]; 40usize] =
        [[[0.0; 3usize]; 20usize]; 40usize];
    pub static mut RootWtCapblUptake: [[f64; 20usize]; 40usize] = [[0.0; 20usize]; 40usize];
    pub static mut SandVolumeFraction: [f64; 40usize] = [0.0; 40usize];
    pub static mut SaturatedHydCond: [f64; 9usize] = [0.0; 9usize];
    pub static mut ShedByCarbonStress: [f64; 20usize] = [0.0; 20usize];
    pub static mut ShedByNitrogenStress: [f64; 20usize] = [0.0; 20usize];
    pub static mut ShedByWaterStress: [f64; 20usize] = [0.0; 20usize];
    pub static mut SoilPsi: [[f64; 20usize]; 40usize] = [[0.0; 20usize]; 40usize];
    pub static mut SoilTemp: [[f64; 20usize]; 40usize] = [[0.0; 20usize]; 40usize];
    pub static mut SoilTempDailyAvrg: [[f64; 20usize]; 40usize] = [[0.0; 20usize]; 40usize];
    pub static mut SquareWeight: [[[f64; 5usize]; 30usize]; 3usize] =
        [[[0.0; 5usize]; 30usize]; 3usize];
    pub static mut StemWeight: [f64; 365usize] = [0.0; 365usize];
    pub static mut thad: [f64; 40usize] = [0.0; 40usize];
    pub static mut thetar: [f64; 40usize] = [0.0; 40usize];
    pub static mut thetas: [f64; 9usize] = [0.0; 9usize];
    pub static mut thts: [f64; 40usize] = [0.0; 40usize];
    pub static mut VarPar: [f64; 61usize] = [0.0; 61usize];
    pub static mut VolNh4NContent: [[f64; 20usize]; 40usize] = [[0.0; 20usize]; 40usize];
    pub static mut VolNo3NContent: [[f64; 20usize]; 40usize] = [[0.0; 20usize]; 40usize];
    pub static mut VolUreaNContent: [[f64; 20usize]; 40usize] = [[0.0; 20usize]; 40usize];
    pub static mut VolWaterContent: [[f64; 20usize]; 40usize] = [[0.0; 20usize]; 40usize];
    pub static mut WindSpeed: [f64; 24usize] = [0.0; 24usize];
    pub static mut wk: [f64; 20usize] = [0.0; 20usize];
    pub fn total_leaf_weight() -> f64 {
        let legacy = crate::LegacyGlobalState::from_globals();
        let mut result = 0.0;
        if legacy.first_square <= 0 {
            result += 0.2;
        }

        for i in 0..legacy.num_pre_fru_nodes.max(0) as usize {
            result += legacy.leaf_weight_pre_fru[i];
        }

        for_each_fruiting_site(
            legacy.num_veg_branches,
            &legacy.num_fruit_branches,
            &legacy.num_nodes,
            |k, l, nodes| {
                result += legacy.leaf_weight_main_stem[[k, l]];
                result += legacy
                    .leaf_weight_nodes
                    .slice(ndarray::s![k, l, 0..nodes])
                    .sum();
            },
        );

        result
    }

    pub fn total_leaf_area() -> f64 {
        let legacy = crate::LegacyGlobalState::from_globals();
        let mut result = 0.0;
        if legacy.first_square <= 0 {
            result += 0.20 * 0.6;
        }

        for i in 0..legacy.num_pre_fru_nodes.max(0) as usize {
            result += legacy.leaf_area_pre_fru[i];
        }

        for_each_fruiting_site(
            legacy.num_veg_branches,
            &legacy.num_fruit_branches,
            &legacy.num_nodes,
            |k, l, nodes| {
                result += legacy.leaf_area_main_stem[[k, l]];
                result += legacy
                    .leaf_area_nodes
                    .slice(ndarray::s![k, l, 0..nodes])
                    .sum();
            },
        );

        result
    }
}

pub use global_rust_defs::*;

#[inline]
pub fn TotalLeafWeight() -> f64 {
    global_rust_defs::total_leaf_weight()
}

#[inline]
pub fn TotalLeafArea() -> f64 {
    global_rust_defs::total_leaf_area()
}
