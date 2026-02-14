use crate::atmosphere::{num_hours, Atmosphere};
use crate::general_functions::{wcond, GetFromClim};
use crate::model_state::{
    for_each_cell, for_each_cell_in_rows, for_each_fruiting_site, for_each_row_mut,
};
use crate::plant::growth::PlantGrowth;
use crate::plant::growth::{check_dry_matter_balance, defoliate, LeafResistance, PhysiologicalAge};
use crate::plant::phenology::cotton_phenology;
use crate::plant::Plant;
use crate::profile::{AgronomyOperation, FertilizationMethod, LightInterceptMethod, Profile};
use crate::soil::hydrology::{
    average_psi, capillary_flow, drain, soil_sum, ComputeIrrigation, WaterUptake,
};
use crate::soil::nitrogen::{soil_nitrogen, soil_nitrogen_average, soil_nitrogen_bal};
use crate::soil::Soil;
use crate::utils::{cell_distance, fmax, fmin, slab_horizontal_location, slab_vertical_location};
use crate::{
    addwtbl, bPollinSwitch, beta, dl, isw, maxl, nk, nl, noitr, thad, thts, wk,
    ActualTranspiration, AgeOfPreFruNode, AppliedWater, AverageLeafAge, AverageSoilPsi, Clim,
    Cotton2KError, CumFertilizerN, CumNitrogenUptake, CumTranspiration, CumWaterAdded,
    CumWaterDrained, DayEmerge, DayInc, DayOfSimulation, DayStart, DayStartPredIrrig,
    DayStopPredIrrig, Daynum, ElCondSatSoilToday, FirstSquare, Irrig, IrrigMethod, Kday, LeafAge,
    LeafArea, LeafAreaIndex, LeafAreaIndexes, LeafAreaMainStem, LeafAreaNodes, LeafAreaPreFru,
    LeafNitrogen, LightIntercept, LightInterceptLayer, LocationColumnDrip, LocationLayerDrip,
    LwpMax, LwpMin, MaxIrrigation, MaxWaterCapacity, ModelState, NO3FlowFraction, NodeLayer,
    NodeLayerPreFru, NumFruitBranches, NumIrrigations, NumLayersWithRoots, NumNodes,
    NumPreFruNodes, NumVegBranches, PerPlantArea, PlantHeight, PlantPopulation, PlantRowColumn,
    PoreSpace, ReferenceETP, RootColNumLeft, RootColNumRight, RootWeight, RootWtCapblUptake,
    RowSpace, SaturatedHydCond, SoilNitrogenLoss, SoilPsi, SupplyNH4N, SupplyNO3N, VolNh4NContent,
    VolNo3NContent, VolUreaNContent, VolWaterContent, WaterTableLayer, CLIMATE_METRIC_IRRD,
    CLIMATE_METRIC_RAIN,
};
use chrono::{Datelike, NaiveDate};
use ndarray::s;

fn column_shade_fraction(
    sw1: f64,
    plant_height: f64,
    light_intercept: f64,
    zint: f64,
    lmax: f64,
    leaf_area_index: f64,
) -> f64 {
    if sw1 >= plant_height {
        0.
    } else {
        let mut shade = 1. - (sw1 / plant_height).powi(2);
        if light_intercept < zint && leaf_area_index < lmax {
            shade *= light_intercept / zint;
        }
        shade
    }
}

fn update_column_radiation_fractions(
    rracol: &mut [f64],
    col_widths: &[f64],
    num_cols: usize,
    plant_col: usize,
    plant_height: f64,
    light_intercept: f64,
    leaf_area_index: f64,
    lmax: f64,
    zint: f64,
) {
    let mut width_accum = 0.; // sum of column widths
    let mut width_left = 0.; // sum of column widths up to location of plant row.
    for (k, &right_col_width) in col_widths.iter().take(num_cols).enumerate() {
        // distance from middle of a column to the plant row, cm.
        let (k0, sw1) = if k <= plant_col {
            // When the column is on the left of the plant row.
            let j = plant_col - k;
            let left_col_width = col_widths[j];
            width_accum += left_col_width;
            width_left = width_accum;
            (j, width_accum - left_col_width / 2.)
        } else {
            // When the column is on the right of the plant row.
            width_accum += right_col_width;
            (k, width_accum - width_left - right_col_width / 2.)
        };
        let shade = column_shade_fraction(
            sw1,
            plant_height,
            light_intercept,
            zint,
            lmax,
            leaf_area_index,
        );
        rracol[k0] = fmax(1. - shade, 0.05);
    }
}

fn total_leaf_weight_in_state(legacy: &crate::LegacyGlobalState) -> f64 {
    let mut total_leaf_weight = if legacy.first_square <= 0 { 0.2 } else { 0.0 };
    total_leaf_weight += legacy
        .leaf_weight_pre_fru
        .iter()
        .take(legacy.num_pre_fru_nodes.max(0) as usize)
        .sum::<f64>();
    for_each_fruiting_site(
        legacy.num_veg_branches,
        &legacy.num_fruit_branches,
        &legacy.num_nodes,
        |k, l, nodes| {
            total_leaf_weight += legacy.leaf_weight_main_stem[[k, l]];
            total_leaf_weight += legacy.leaf_weight_nodes.slice(s![k, l, 0..nodes]).sum();
        },
    );
    total_leaf_weight
}

fn leaf_resistance_stats(legacy: &crate::LegacyGlobalState) -> (usize, f64) {
    let mut num_leaves = legacy.num_pre_fru_nodes as usize;
    let mut resistance_sum = legacy
        .age_of_pre_fru_node
        .iter()
        .take(num_leaves)
        .map(|&leaf_age| LeafResistance(leaf_age))
        .sum::<f64>();
    for_each_fruiting_site(
        legacy.num_veg_branches,
        &legacy.num_fruit_branches,
        &legacy.num_nodes,
        |k, l, nodes| {
            num_leaves += nodes;
            resistance_sum += legacy
                .leaf_age
                .slice(s![k, l, 0..nodes])
                .iter()
                .map(|&leaf_age| LeafResistance(leaf_age))
                .sum::<f64>();
        },
    );
    (num_leaves, resistance_sum)
}

fn add_leaf_area_with_age(
    legacy: &mut crate::LegacyGlobalState,
    layer: usize,
    leaf_area: f64,
    leaf_age: f64,
) {
    legacy.leaf_area[layer] += leaf_area;
    legacy.average_leaf_age[layer] += leaf_area * leaf_age;
}

fn add_nutrients_to_soil_cell(
    legacy: &mut crate::LegacyGlobalState,
    layer: usize,
    column: usize,
    ammonium_mg: f64,
    nitrate_mg: f64,
    urea_mg: f64,
) {
    let cell_volume = legacy.dl[layer] * legacy.wk[column];
    legacy.vol_nh4_n_content[[layer, column]] += ammonium_mg / cell_volume;
    legacy.vol_no3_n_content[[layer, column]] += nitrate_mg / cell_volume;
    legacy.vol_urea_n_content[[layer, column]] += urea_mg / cell_volume;
}

fn add_uniform_nutrients_in_rows(
    legacy: &mut crate::LegacyGlobalState,
    row_start: usize,
    row_end: usize,
    num_cols: usize,
    nh4_delta: f64,
    no3_delta: f64,
    urea_delta: f64,
) {
    for_each_cell_in_rows(row_start, row_end, num_cols, |layer, column| {
        legacy.vol_nh4_n_content[[layer, column]] += nh4_delta;
        legacy.vol_no3_n_content[[layer, column]] += no3_delta;
        legacy.vol_urea_n_content[[layer, column]] += urea_delta;
    });
}

fn add_uniform_nutrients_in_row(
    legacy: &mut crate::LegacyGlobalState,
    row: usize,
    num_cols: usize,
    nh4_delta: f64,
    no3_delta: f64,
    urea_delta: f64,
) {
    add_uniform_nutrients_in_rows(
        legacy,
        row,
        row + 1,
        num_cols,
        nh4_delta,
        no3_delta,
        urea_delta,
    );
}

fn add_uniform_surface_water(legacy: &mut crate::LegacyGlobalState, water_delta: f64) {
    let num_cols = legacy.nk as usize;
    legacy
        .vol_water_content
        .slice_mut(s![0, 0..num_cols])
        .iter_mut()
        .for_each(|cell| *cell += water_delta);
}

fn find_water_table_layer(legacy: &crate::LegacyGlobalState, water_table_depth_cm: f64) -> i32 {
    let mut cumulative_depth = 0.;
    for (layer, &layer_depth) in legacy.dl.iter().take(legacy.nl as usize).enumerate() {
        cumulative_depth += layer_depth;
        if cumulative_depth > water_table_depth_cm {
            return layer as i32;
        }
    }
    1000
}

fn update_layer_water_table_saturation(
    water_row: &mut ndarray::ArrayViewMut1<'_, f64>,
    num_cols: usize,
    uplimit: f64,
    force_saturation: bool,
    layer_depth: f64,
    row_space: f64,
    col_widths: &[f64],
    added_water: &mut f64,
) {
    for (col, water_cell) in water_row.iter_mut().take(num_cols).enumerate() {
        let previous = *water_cell;
        if force_saturation || previous > uplimit {
            *water_cell = uplimit;
            *added_water +=
                10. * (*water_cell - previous) * layer_depth * col_widths[col] / row_space;
        }
    }
}

fn latered_gross_photosynthesis(
    legacy: &crate::LegacyGlobalState,
    pstand: f64,
    ptsred: f64,
    pnetcor: f64,
    ptnfac: f64,
) -> f64 {
    let mut pplant = 0.;
    let mut pstand_remain = pstand;
    for i in (0..20).rev() {
        if pstand_remain <= 0. {
            break;
        }
        if legacy.light_intercept_layer[i] <= 0. {
            continue;
        }
        let page = 1. - (legacy.average_leaf_age[i] / drop_leaf_age(legacy.leaf_area[i])).powi(2);
        let mut pplant_inc = 0.001
            * pstand_remain
            * legacy.light_intercept_layer[i]
            * legacy.per_plant_area
            * ptsred
            * pnetcor
            * ptnfac
            * page;
        if pplant_inc > pstand_remain {
            pplant_inc = pstand_remain;
        }
        pplant += pplant_inc;
        pstand_remain -= pplant_inc;
    }
    pplant
}

#[derive(Debug, Clone, Copy, Default)]
struct RootHydraulicSums {
    psinum: f64,
    rootvol: f64,
    rrlsum: f64,
    sumlv: f64,
    vh2sum: f64,
}

fn accumulate_root_zone_hydraulic_sums(
    legacy: &crate::LegacyGlobalState,
    vpsil: &[f64; 13],
    cmg: f64,
) -> RootHydraulicSums {
    let mut sums = RootHydraulicSums::default();
    for (l, k) in
        root_zone_cells_in_bounds(legacy, legacy.num_layers_with_roots as usize, |layer| {
            let left_col = legacy.root_col_num_left[layer];
            (left_col, left_col, false)
        })
    {
        let uptake_weight = legacy.root_wt_capbl_uptake[[l, k]];
        if uptake_weight >= vpsil[10] {
            let weighted_uptake = fmin(uptake_weight, vpsil[11]);
            sums.psinum += weighted_uptake;
            sums.sumlv += weighted_uptake * cmg;
            sums.rootvol += legacy.dl[l] * legacy.wk[k];
            let soil_psi = legacy.soil_psi[[l, k]];
            let rrl = if soil_psi <= vpsil[1] {
                vpsil[2] / cmg
            } else {
                (vpsil[3] - soil_psi * (vpsil[4] + vpsil[5] * soil_psi)) / cmg
            };
            sums.rrlsum += weighted_uptake / rrl;
            sums.vh2sum += legacy.vol_water_content[[l, k]] * weighted_uptake;
        }
    }
    sums
}

fn iterate_with_capillary_flow(
    legacy: &mut crate::LegacyGlobalState,
    iterations: usize,
    mut before_capillary: impl FnMut(&mut crate::LegacyGlobalState) -> Result<(), Cotton2KError>,
) -> Result<(), Cotton2KError> {
    for _ in 0..iterations {
        before_capillary(legacy)?;
        legacy.write_to_globals();
        capillary_flow();
        legacy.read_from_globals();
    }
    Ok(())
}

fn root_zone_cells_in_bounds(
    legacy: &crate::LegacyGlobalState,
    num_layers: usize,
    mut bounds: impl FnMut(usize) -> (i32, i32, bool),
) -> Vec<(usize, usize)> {
    let mut cells = Vec::new();
    if legacy.nk <= 0 {
        return cells;
    }
    let max_col = legacy.nk - 1;
    for layer in 0..num_layers {
        let (left_raw, right_raw, inclusive_right) = bounds(layer);
        let left = left_raw.max(0) as usize;
        let right = right_raw.min(max_col);
        if inclusive_right {
            if right < left as i32 {
                continue;
            }
            for col in left..=right as usize {
                cells.push((layer, col));
            }
        } else {
            let end = right.max(0) as usize;
            if end <= left {
                continue;
            }
            for col in left..end {
                cells.push((layer, col));
            }
        }
    }
    cells
}

fn sidedress_target_cells(
    legacy: &crate::LegacyGlobalState,
    layer: usize,
    column: usize,
) -> Vec<(usize, usize)> {
    let mut cells = vec![(layer, column)];
    if legacy.dl[layer] * legacy.wk[column] < 100. {
        if column < (legacy.nk - 1) as usize {
            cells.push((layer, column + 1));
        }
        if column > 0 {
            cells.push((layer, column - 1));
        }
        if layer < (legacy.nl - 1) as usize {
            cells.push((layer + 1, column));
        }
    }
    cells
}

fn broadcast_incorporation_layers(legacy: &crate::LegacyGlobalState, depth_cm: f64) -> usize {
    let mut layers = 0usize;
    let mut accumulated_depth = 0.;
    for (l, layer_depth) in legacy.dl.iter().take(legacy.nl as usize).enumerate() {
        accumulated_depth += layer_depth;
        if accumulated_depth >= depth_cm {
            layers = l + 1;
            break;
        }
    }
    layers
}

fn apply_broadcast_fertilization(
    legacy: &mut crate::LegacyGlobalState,
    ferc: f64,
    ammonium: f64,
    nitrate: f64,
    urea: f64,
) {
    let lplow = broadcast_incorporation_layers(legacy, 20.0);
    let fertdp: f64 = legacy.dl.iter().take(lplow).sum();
    let nh4_delta = ammonium * ferc / fertdp;
    let no3_delta = nitrate * ferc / fertdp;
    let urea_delta = urea * ferc / fertdp;
    add_uniform_nutrients_in_rows(
        legacy,
        0,
        lplow,
        legacy.nk as usize,
        nh4_delta,
        no3_delta,
        urea_delta,
    );
}

fn apply_foliar_fertilization(
    legacy: &mut crate::LegacyGlobalState,
    ferc: f64,
    ammonium: f64,
    nitrate: f64,
    urea: f64,
) {
    let interception = 0.70 * legacy.light_intercept;
    legacy.leaf_nitrogen += interception * (ammonium + urea) * 1000. / legacy.plant_population;
    let nh4_delta = ammonium * (1. - interception) * ferc / legacy.dl[0];
    let no3_delta = nitrate * ferc / legacy.dl[0];
    let urea_delta = urea * (1. - interception) * ferc / legacy.dl[0];
    add_uniform_nutrients_in_row(
        legacy,
        0,
        legacy.nk as usize,
        nh4_delta,
        no3_delta,
        urea_delta,
    );
}

fn apply_sidedress_fertilization(
    legacy: &mut crate::LegacyGlobalState,
    ferc: f64,
    ammonium: f64,
    nitrate: f64,
    urea: f64,
    drip_x: f64,
    drip_y: f64,
) -> Result<(), Cotton2KError> {
    let ksdr = slab_horizontal_location(drip_x, legacy.row_space)?;
    let lsdr = slab_vertical_location(drip_y)?;
    let target_cells = sidedress_target_cells(legacy, lsdr, ksdr);
    let cell_count = target_cells.len() as f64;
    let addamm = ammonium * ferc * legacy.row_space / cell_count;
    let addnit = nitrate * ferc * legacy.row_space / cell_count;
    let addnur = urea * ferc * legacy.row_space / cell_count;
    for (layer, column) in target_cells {
        add_nutrients_to_soil_cell(legacy, layer, column, addamm, addnit, addnur);
    }
    Ok(())
}

fn apply_drip_fertilization(
    legacy: &mut crate::LegacyGlobalState,
    ferc: f64,
    ammonium: f64,
    nitrate: f64,
    urea: f64,
) {
    add_nutrients_to_soil_cell(
        legacy,
        legacy.location_layer_drip as usize,
        legacy.location_column_drip as usize,
        ammonium * ferc * legacy.row_space,
        nitrate * ferc * legacy.row_space,
        urea * ferc * legacy.row_space,
    );
}

#[derive(Debug, Clone, Copy)]
pub struct State {
    pub date: NaiveDate,

    pub soil: Soil,
    pub plant: Plant,
    pub atmosphere: Atmosphere,
}

impl State {
    pub fn new(profile: &Profile, date: NaiveDate) -> Self {
        let atmosphere = Atmosphere::new(date, profile.longitude, profile.latitude);
        State {
            date,
            atmosphere,
            soil: Soil::new(profile),
            plant: Plant::new(),
        }
    }

    /// This function executes all the simulation computations in a day. It is called from [Profile::run()], and
    /// [Profile::adjust()].
    ///
    /// It calls the following functions:
    /// * [Profile::column_shading()]
    /// * [crate::soil::thermodynamics::SoilThermodynamics::simulate()]
    /// * [`State::soil_procedures()`]
    /// * [SoilNitrogen()]
    /// * [SoilSum()]
    /// * [PhysiologicalAge()]
    /// * [defoliate()]
    /// * [Profile::stress()]
    /// * [Profile::get_net_photosynthesis()]
    /// * [PlantGrowth::plant_growth()]
    /// * [CottonPhenology()]
    /// * [PlantNitrogen::plant_nitrogen()]
    /// * [check_dry_matter_balance()]
    /// * [PlantNitrogen::plant_nitrogen_balance()]
    /// * [SoilNitrogenBal()]
    /// * [SoilNitrogenAverage()]
    ///
    /// The following global variables are referenced here:
    /// * [DayEmerge]
    /// * [DayStart]
    /// * [Kday]
    /// * [LeafAreaIndex]
    /// The following global variables are set here:
    /// * [bEnd]
    /// * [DayInc]
    /// * [Daynum]
    /// * [DayOfSimulation]
    /// * [isw]
    /// * [Kday]
    pub fn simulate_this_day(
        &mut self,
        profile: &mut Profile,
        model_state: &mut ModelState,
    ) -> Result<(), Cotton2KError> {
        let root_tables = profile
            .root_impedance_tables
            .as_ref()
            .expect("root impedance tables must be initialized before simulation")
            .clone();
        {
            let legacy = &mut model_state.legacy;
            // Compute Daynum (day of year), Date, and DayOfSimulation (days from start of simulation).
            legacy.daynum += 1;
            legacy.day_of_simulation = legacy.daynum - legacy.day_start + 1;
            // Compute Kday (days from emergence).
            legacy.kday = if legacy.day_emerge <= 0 {
                0
            } else {
                legacy.daynum - legacy.day_emerge + 1
            };
            if legacy.kday < 0 {
                legacy.kday = 0;
            }
        }
        // The following functions are executed each day (also before emergence).
        self.column_shading(profile, model_state); // computes light interception and soil shading.
        self.atmosphere.meteorology(profile, model_state); // computes climate variables for today.
        model_state.legacy.write_to_globals();
        self.soil.thermodynamics.simulate(profile); // executes all modules of soil and canopy temperature.
        model_state.legacy.read_from_globals();
        self.soil_procedures(profile, model_state)?; // executes all other soil processes.
        model_state.legacy.write_to_globals();
        soil_nitrogen(); // computes nitrogen transformations in the soil.
        soil_sum(); // computes totals of water and N in the soil.
        model_state.legacy.read_from_globals();

        // The following is executed each day after plant emergence:
        if model_state.legacy.daynum >= model_state.legacy.day_emerge && model_state.legacy.isw > 0
        {
            // If this day is after emergence, assign to isw the value of 2.
            model_state.legacy.isw = 2;
            model_state.legacy.write_to_globals();
            model_state.legacy.day_inc = PhysiologicalAge(); // computes physiological age
            defoliate(); // effects of defoliants applied.
            model_state.legacy.read_from_globals();
            self.stress(profile, model_state); // computes water stress factors.
            self.get_net_photosynthesis(profile, model_state)?; // computes net photosynthesis.
            model_state.legacy.write_to_globals();
            self.plant
                .grow(self.atmosphere, &profile.agronomy_operations, &root_tables); // executes all modules of plant growth.
            cotton_phenology(); // executes all modules of plant phenology.
            self.plant.nitrogen.run(); // computes plant nitrogen allocation.
            check_dry_matter_balance(); // checks plant dry matter balance.
            model_state.legacy.read_from_globals();

            // If the relevant output flag is not zero, compute soil nitrogen balance and soil nitrogen averages by
            // layer, and write this information to files.
            if false {
                model_state.legacy.write_to_globals();
                self.plant.nitrogen.balance(); // checks plant nitrogen balance.
                soil_nitrogen_bal(); // checks soil nitrogen balance.
                soil_nitrogen_average(); // computes average soil nitrogen by layers.
                model_state.legacy.read_from_globals();
            }
        }
        // Check if the date to stop simulation has been reached, or if this is the last day with available weather
        // data. Simulation will also stop when no leaves remain on the plant.
        if self.date >= profile.last_day_weather_data {
            return Err(Cotton2KError {
                level: 0,
                message: String::from("No more weather data!"),
            });
        }
        if model_state.legacy.kday > 10 && model_state.legacy.leaf_area_index < 0.0002 {
            return Err(Cotton2KError {
                level: 0,
                message: String::from("Leaf area index is too small!"),
            });
        }
        model_state.legacy.write_to_globals();
        Ok(())
    }
    /// This function computes light interception by crop canopy and shading of soil columns by the plants. It is called from SimulateThisDay().
    ///
    /// The following global variables are referenced here:
    /// * [DayEmerge]
    /// * [Daynum]
    /// * [isw]
    /// * [LeafAreaIndex]
    /// * [PlantHeight]
    /// * [PlantRowColumn]
    /// * [nk]
    /// * [RowSpace]
    ///
    /// The following global variables are set here:
    /// * [LightIntercept]
    /// * [rracol]
    fn column_shading(&mut self, profile: &mut Profile, model_state: &mut ModelState) {
        let legacy = &mut model_state.legacy;
        // Before emergence: no light interception and no shading.
        // LightIntercept is assigned zero, and the rracol array is assigned 1.
        if legacy.daynum < legacy.day_emerge || legacy.isw <= 0 || legacy.day_emerge <= 0 {
            legacy.light_intercept = 0.;
            self.soil
                .thermodynamics
                .rracol
                .iter_mut()
                .take(legacy.nk as usize)
                .for_each(|rracol| *rracol = 1.0);
            legacy.write_to_globals();
            return;
        }
        // Compute the maximum leaf area index until this day (lmax).
        if legacy.leaf_area_index > profile.lmax {
            profile.lmax = legacy.leaf_area_index;
        }
        // Light interception is computed by two methods:
        //
        // 1. It is assumed to be proportional to the ratio of plant height to row spacing.
        // light interception computed from plant height.
        let zint = 1.0756 * legacy.plant_height / legacy.row_space;
        match profile.light_intercept_method {
            LightInterceptMethod::Latered => {
                let params = profile
                    .light_intercept_parameters
                    .as_ref()
                    .expect("light intercept parameters must be provided for Latered method");
                legacy.leaf_area.slice_mut(s![0..20]).fill(0.);
                legacy.average_leaf_age.slice_mut(s![0..20]).fill(0.);
                for m in 0..9 {
                    add_leaf_area_with_age(
                        legacy,
                        legacy.node_layer_pre_fru[m] as usize,
                        legacy.leaf_area_pre_fru[m],
                        legacy.age_of_pre_fru_node[m],
                    );
                }
                let mut branch_sites = Vec::new();
                for_each_fruiting_site(
                    legacy.num_veg_branches,
                    &legacy.num_fruit_branches,
                    &legacy.num_nodes,
                    |k, l, nodes| branch_sites.push((k, l, nodes)),
                );
                for (k, l, nodes) in branch_sites {
                    let layer = legacy.node_layer[[k, l]] as usize;
                    add_leaf_area_with_age(
                        legacy,
                        layer,
                        legacy.leaf_area_main_stem[[k, l]],
                        legacy.leaf_age[[k, l, 0]],
                    );
                    for m in 0..nodes {
                        add_leaf_area_with_age(
                            legacy,
                            layer,
                            legacy.leaf_area_nodes[[k, l, m]],
                            legacy.leaf_age[[k, l, m]],
                        );
                    }
                }
                if legacy.first_square <= 0 {
                    add_leaf_area_with_age(legacy, 0, 0.20 * 0.6, 0.0);
                }
                let mut light_through = 0.;
                for (
                    (((avg_leaf_age, leaf_area_index), &leaf_area), &param),
                    light_intercept_layer,
                ) in legacy
                    .average_leaf_age
                    .iter_mut()
                    .zip(legacy.leaf_area_indexes.iter_mut())
                    .zip(legacy.leaf_area.iter())
                    .zip(params.iter())
                    .zip(legacy.light_intercept_layer.iter_mut())
                {
                    *avg_leaf_age /= leaf_area;
                    *leaf_area_index = leaf_area / legacy.per_plant_area;
                    *light_intercept_layer = 1. - (param * *leaf_area_index).exp();
                    light_through += param * *leaf_area_index;
                }
                legacy.light_intercept = 1. - light_through.exp();
            }
            LightInterceptMethod::Fry1980 => {
                legacy.light_intercept = 0.39 * legacy.leaf_area_index.powf(0.68);
            }
            LightInterceptMethod::Original => {
                // 2. It is computed as a function of leaf area index. If LeafAreaIndex is not greater than 0.5 lfint is a linear function of it.

                // light interception computed from leaf area index.
                let lfint = if legacy.leaf_area_index <= 0.5 {
                    0.80 * legacy.leaf_area_index
                } else {
                    // If the leaf area index is greater than 0.5, lfint is computed as an exponential function of LeafAreaIndex.
                    1. - (0.07 - 1.16 * legacy.leaf_area_index).exp()
                };
                // If lfint is greater then zint, LightIntercept is their average value.
                // Otherwise, if the LeafAreaIndex is decreasing, it is lfint. Else it is zint.
                legacy.light_intercept = if lfint > zint {
                    0.5 * (zint + lfint)
                } else if legacy.leaf_area_index < profile.lmax {
                    lfint
                } else {
                    zint
                };
            }
        }
        // The value of LightIntercept is between zero and one.
        if legacy.light_intercept < 0. {
            legacy.light_intercept = 0.;
        }
        if legacy.light_intercept > 1. {
            legacy.light_intercept = 1.;
        }
        update_column_radiation_fractions(
            &mut self.soil.thermodynamics.rracol,
            legacy
                .wk
                .as_slice()
                .expect("wk should be contiguous in standard layout"),
            legacy.nk as usize,
            legacy.plant_row_column as usize,
            legacy.plant_height,
            legacy.light_intercept,
            legacy.leaf_area_index,
            profile.lmax,
            zint,
        );
        legacy.write_to_globals();
    }

    fn process_daily_rainfall_and_runoff(&mut self, legacy: &mut crate::LegacyGlobalState) -> f64 {
        let mut rain_today = GetFromClim(CLIMATE_METRIC_RAIN, self.date.ordinal() as i32);
        legacy.b_pollin_switch = rain_today < 2.5;
        let mut runoff_today = 0.;
        if rain_today >= 2. {
            runoff_today = self.soil.hydrology.runoff(rain_today);
            if runoff_today < rain_today {
                rain_today -= runoff_today;
            } else {
                rain_today = 0.;
            }
            let j = legacy.daynum - legacy.day_start;
            let mut clim = Clim.write().expect("Clim lock poisoned");
            clim[j as usize].Rain = rain_today;
        }
        self.soil.hydrology.runoff = runoff_today;
        GetFromClim(CLIMATE_METRIC_RAIN, legacy.daynum)
    }

    fn apply_predicted_irrigation(
        &mut self,
        legacy: &mut crate::LegacyGlobalState,
        water_to_apply: &mut f64,
        drip_water_amount: &mut f64,
    ) {
        if legacy.max_irrigation > 0.
            && legacy.daynum >= legacy.day_start_pred_irrig
            && legacy.daynum < legacy.day_stop_pred_irrig
        {
            legacy.write_to_globals();
            ComputeIrrigation();
            legacy.read_from_globals();
            if legacy.irrig_method == 2 {
                *drip_water_amount = legacy.applied_water;
            } else {
                *water_to_apply += legacy.applied_water;
            }
            legacy.applied_water = 0.;
        }
    }

    fn apply_scheduled_irrigation(
        &mut self,
        legacy: &mut crate::LegacyGlobalState,
        water_to_apply: &mut f64,
        drip_water_amount: &mut f64,
    ) {
        let irrig = Irrig.read().expect("Irrig lock poisoned");
        for i in 0..legacy.num_irrigations as usize {
            if legacy.daynum == irrig[i].day {
                if irrig[i].method == 2 {
                    *drip_water_amount += irrig[i].amount;
                    legacy.location_column_drip = irrig[i].LocationColumnDrip;
                    legacy.location_layer_drip = irrig[i].LocationLayerDrip;
                } else {
                    *water_to_apply += irrig[i].amount;
                }
                break;
            }
        }
    }
    /// This function manages all the soil related processes, and is executed once each day.
    ///
    /// It is called from [Profile::simulate_this_day()] and it calls the following functions:
    /// * [ApplyFertilizer()]
    /// * [AveragePsi()]
    /// * [CapillaryFlow()]
    /// * [ComputeIrrigation()]
    /// * [DripFlow()]
    /// * gravity-flow redistribution for rain/surface irrigation
    /// * [RootsCapableOfUptake()]
    /// * [WaterUptake()]
    /// * [WaterTable()]
    ///
    /// The following global variables are referenced here:
    /// * [ActualTranspiration]
    /// * [Clim]
    /// * [DayEmerge]
    /// * [Daynum]
    /// * [DayStartPredIrrig]
    /// * [DayStopPredIrrig]
    /// * [dl]
    /// * [Irrig]
    /// * [IrrigMethod]
    /// * [isw]
    /// * [Kday]
    /// * [MaxIrrigation]
    /// * [nk]
    /// * [nl]
    /// * [NumIrrigations]
    /// * [PerPlantArea]
    /// * [SoilPsi]
    /// * [RowSpace]
    /// * [SupplyNH4N]
    /// * [SupplyNO3N]
    /// * [VolWaterContent]
    /// * [wk]
    ///
    /// The following global variables are set here:
    /// * [AverageSoilPsi]
    /// * [CumNitrogenUptake]
    /// * [CumTranspiration]
    /// * [CumWaterAdded]
    /// * [LocationColumnDrip]
    /// * [LocationLayerDrip]
    /// * [noitr]
    fn soil_procedures(
        &mut self,
        profile: &mut Profile,
        model_state: &mut ModelState,
    ) -> Result<(), Cotton2KError> {
        // The following constant parameters are used:
        const cpardrip: f64 = 0.2;
        const cparelse: f64 = 0.4;
        // Call function ApplyFertilizer() for nitrogen fertilizer application.
        self.apply_fertilizer(profile, model_state)?;
        let legacy = &mut model_state.legacy;
        let mut drip_water_amount = 0f64; // amount of water applied by drip irrigation

        // amount of water applied by non-drip irrigation or rainfall
        let mut water_to_apply = self.process_daily_rainfall_and_runoff(legacy);

        // If irrigation is to be predicted for this day, call ComputeIrrigation() to compute the actual amount of irrigation.
        self.apply_predicted_irrigation(legacy, &mut water_to_apply, &mut drip_water_amount);

        // When water is added by an irrigation defined in the input: update the amount of applied water.
        self.apply_scheduled_irrigation(legacy, &mut water_to_apply, &mut drip_water_amount);
        legacy.cum_water_added += water_to_apply + drip_water_amount;

        // The following will be executed only after plant emergence
        if legacy.daynum >= legacy.day_emerge && legacy.isw > 0 {
            RootsCapableOfUptake(legacy); // function computes roots capable of uptake for each soil cell
            legacy.write_to_globals();
            legacy.average_soil_psi = average_psi(); // function computes the average matric soil water potential in the root zone, weighted by the roots-capable-of-uptake.
            WaterUptake(); // function  computes water and nitrogen uptake by plants.
            legacy.read_from_globals();

            // Update the cumulative sums of actual transpiration (CumTranspiration, mm) and total uptake of nitrogen (CumNitrogenUptake, mg N per slab, converted from total N supply, g per plant).
            legacy.cum_transpiration += legacy.actual_transpiration;
            legacy.cum_nitrogen_uptake +=
                (legacy.supply_no3_n + legacy.supply_nh4_n) * 10. * legacy.row_space
                    / legacy.per_plant_area;
        }

        // Call function WaterTable() for saturating soil below water table.
        if profile.num_watertable_data > 0 {
            self.watertable(profile, legacy)?;
        }

        if water_to_apply > 0. {
            // For rain or surface irrigation.
            // The number of iterations is computed from the thickness of the first soil layer.
            legacy.noitr = (cparelse * water_to_apply / (legacy.dl[0] + 2.) + 1.) as i32;
            // the amount of water applied, mm per iteration.
            let applywat = water_to_apply / legacy.noitr as f64;
            // The following redistribution steps are called noitr times per day:
            // surface/flood gravity flow update, followed by CapillaryFlow().
            iterate_with_capillary_flow(legacy, legacy.noitr as usize, |legacy| {
                let water_delta = 0.10 * applywat / legacy.dl[0];
                add_uniform_surface_water(legacy, water_delta);
                legacy.write_to_globals();
                let water_drained_out = drain();
                legacy.read_from_globals();
                if water_drained_out > 0. {
                    legacy.cum_water_drained += 10. * water_drained_out / legacy.row_space;
                }
                Ok(())
            })?;
        }

        if drip_water_amount > 0. {
            // For drip irrigation.
            // The number of iterations is computed from the volume of the soil cell in which the water is applied.
            legacy.noitr = (cpardrip * drip_water_amount
                / (legacy.dl[legacy.location_layer_drip as usize]
                    * legacy.wk[legacy.location_column_drip as usize])
                + 1.) as i32;
            // the amount of water applied, mm per iteration.
            let applywat = drip_water_amount / legacy.noitr as f64;
            // If water is applied, DripFlow() is called followed by CapillaryFlow().
            iterate_with_capillary_flow(legacy, legacy.noitr as usize, |legacy| {
                DripFlow(applywat, legacy)
            })?;
        }

        // When no water is added, there is only one iteration in this day.
        if water_to_apply + drip_water_amount <= 0. {
            legacy.noitr = 1;
            iterate_with_capillary_flow(legacy, 1, |_| Ok(()))?;
        }
        legacy.write_to_globals();
        Ok(())
    }

    /// This function simulates the application of nitrogen fertilizer on each date of application.
    /// It is called from [`State::soil_procedures()`].
    ///
    /// The following global variables are referenced here:
    /// * [Daynum]
    /// * [dl]
    /// * [LightIntercept]
    /// * [LocationColumnDrip]
    /// * [LocationLayerDrip]
    /// * [nk]
    /// * [nl]
    /// * [RowSpace]
    /// * [wk]
    ///
    /// The following global variables are set here:
    /// * [CumFertilizerN]
    /// * [LeafNitrogen]
    /// * [VolNh4NContent]
    /// * [VolNo3NContent]
    /// * [VolUreaNContent]
    fn apply_fertilizer(
        &mut self,
        profile: &Profile,
        model_state: &mut ModelState,
    ) -> Result<(), Cotton2KError> {
        let ferc = 0.01; // constant used to convert kgs per ha to mg cm-2
        let legacy = &mut model_state.legacy;

        for ao in &profile.agronomy_operations {
            match ao {
                AgronomyOperation::fertilization {
                    date,
                    urea,
                    nitrate,
                    ammonium,
                    method,
                    drip_x,
                    drip_y,
                } => {
                    if legacy.daynum as u32 == date.ordinal() {
                        legacy.cum_fertilizer_n +=
                            ferc * legacy.row_space * (ammonium + nitrate + urea);
                        match method {
                            FertilizationMethod::Broadcast => {
                                apply_broadcast_fertilization(
                                    legacy, ferc, *ammonium, *nitrate, *urea,
                                );
                            }
                            FertilizationMethod::Foliar => {
                                apply_foliar_fertilization(
                                    legacy, ferc, *ammonium, *nitrate, *urea,
                                );
                            }
                            FertilizationMethod::Sidedress => {
                                apply_sidedress_fertilization(
                                    legacy, ferc, *ammonium, *nitrate, *urea, *drip_x, *drip_y,
                                )?;
                            }
                            FertilizationMethod::Drip => {
                                apply_drip_fertilization(legacy, ferc, *ammonium, *nitrate, *urea);
                            }
                        }
                        legacy.write_to_globals();
                    }
                }
                _ => {}
            }
        }
        Ok(())
    }

    /// This function computes the water stress variables affecting the cotton plants.
    /// It is called by [Profile::simulate_this_day()] and calls [LeafWaterPotential()].
    ///
    /// The following global variables are referenced here:
    /// * [Kday]
    /// * [LwpMin]
    /// * [LwpMax]
    /// The following global variables are set here:
    /// * [AverageLwp]
    /// * [AverageLwpMin]
    /// * [LwpMinX]
    /// * [LwpX]
    /// * [Profile::ptsred]
    /// * [WaterStress]
    /// * [WaterStressStem]
    fn stress(&mut self, profile: &mut Profile, model_state: &mut ModelState) {
        // The following constant parameters are used:
        const vstrs: [f64; 9] = [-3.0, 3.229, 1.907, 0.321, -0.10, 1.230, 0.340, 0.30, 0.05];
        // Call LeafWaterPotential() to compute leaf water potentials.
        let legacy = &mut model_state.legacy;
        LeafWaterPotential(legacy);
        // The running averages, for the last three days, are computed: AverageLwpMin is the average of LwpMin, and
        // AverageLwp of LwpMin + LwpMax.
        legacy.average_lwp_min += (legacy.lwp_min - legacy.lwp_min_x[2]) / 3.;
        legacy.average_lwp += (legacy.lwp_min + legacy.lwp_max - legacy.lwp_x[2]) / 3.;
        for i in [2, 1] {
            legacy.lwp_min_x[i] = legacy.lwp_min_x[i - 1];
            legacy.lwp_x[i] = legacy.lwp_x[i - 1];
        }
        legacy.lwp_min_x[0] = legacy.lwp_min;
        legacy.lwp_x[0] = legacy.lwp_min + legacy.lwp_max;
        //     No stress effects before 5th day after emergence.
        if legacy.kday < 5 {
            profile.ptsred = 1.;
            legacy.water_stress = 1.;
            legacy.water_stress_stem = 1.;
            return;
        }
        // The computation of ptsred, the effect of moisture stress on the photosynthetic rate, is based on the
        // following work:
        // * Ephrath, J.E., Marani, A., Bravdo, B.A., 1990. Effects of moisture stress on stomatal resistance and
        //   photosynthetic rate in cotton (Gossypium hirsutum) 1. Controlled levels of stress. Field Crops Res.
        //   23:117-131.
        //
        // It is a function of AverageLwpMin (average LwpMin for the last three days).
        if legacy.average_lwp_min < vstrs[0] {
            legacy.average_lwp_min = vstrs[0];
        }
        profile.ptsred =
            vstrs[1] + legacy.average_lwp_min * (vstrs[2] + vstrs[3] * legacy.average_lwp_min);
        if profile.ptsred > 1. {
            profile.ptsred = 1.;
        }
        // The general moisture stress factor (WaterStress) is computed as an empirical function of AverageLwp. psilim,
        // the value of AverageLwp at the maximum value of the function, is used for truncating it.

        // The minimum value of WaterStress is 0.05, and the maximum is 1.
        let psilim = -0.5 * vstrs[5] / vstrs[6]; // limiting value of AverageLwp.
        if legacy.average_lwp > psilim {
            legacy.water_stress = 1.;
        } else {
            legacy.water_stress =
                vstrs[4] - legacy.average_lwp * (vstrs[5] + vstrs[6] * legacy.average_lwp);
            if legacy.water_stress > 1. {
                legacy.water_stress = 1.;
            }
            if legacy.water_stress < 0.05 {
                legacy.water_stress = 0.05;
            }
        }
        // Water stress affecting plant height and stem growth (WaterStressStem) is assumed to be more severe than
        // WaterStress, especially at low WaterStress values.
        legacy.water_stress_stem =
            legacy.water_stress * (1. + vstrs[7] * (2. - legacy.water_stress)) - vstrs[7];
        if legacy.water_stress_stem < vstrs[8] {
            legacy.water_stress_stem = vstrs[8];
        }

        legacy.write_to_globals();
    }
    /// This function simulates the net photosynthesis of cotton  plants. It is called daily by
    /// [Profile::simulate_this_day()]. This is essentially the routine of GOSSYM with minor changes.
    ///
    /// The following global and file scope variables are referenced here:
    /// * [BurrWeightOpenBolls]
    /// * [CottonWeightOpenBolls]
    /// * [Daynum]
    /// * [DayEmerge]
    /// * [DayTimeTemp]
    /// * [iyear]
    /// * [Kday]
    /// * [LeafNConc]
    /// * [LightIntercept]
    /// * [PerPlantArea]
    /// * [PlantWeight]
    /// * [Profile::ptsred]
    /// * [StemWeight]
    /// * [TotalLeafWeight()]
    ///
    /// The following global variables are set here:
    /// * [bEnd]
    /// * [CumNetPhotosynth]
    /// * [NetPhotosynthesis]
    ///
    /// References:
    ///
    /// * Baker et. al. (1972). Simulation of Growth and Yield in Cotton: I. Gross photosynthesis, respiration and
    ///   growth. Crop Sci. 12:431-435.
    /// * Harper et. al. (1973) Carbon dioxide and the photosynthesis of field crops. A metered carbon dioxide release
    ///   in cotton under field conditions.  Agron. J. 65:7-11.
    /// * Baker (1965) Effects of certain environmental factors on net assimilation in cotton. Crop Sci. 5:53-56 (Fig 5).
    fn get_net_photosynthesis(
        &mut self,
        profile: &Profile,
        model_state: &mut ModelState,
    ) -> Result<(), Cotton2KError> {
        //  constants:
        const gsubr: f64 = 0.375; // the growth resiration factor.
        const rsubo: f64 = 0.0032; // maintenance respiration factor.
        const vpnet: [f64; 4] = [1.30, 0.034, 0.010, 0.32];
        const co2parm: [f64; 45] =
            // parameters used to correct photosynthesis for ambient CO2 concentration.
            [
                1.0235, 1.0264, 1.0285, 1.0321, 1.0335, 1.0353, 1.0385, 1.0403, 1.0431, 1.0485,
                1.0538, 1.0595, 1.0627, 1.0663, 1.0716, 1.0752, 1.0784, 1.0823, 1.0880, 1.0923,
                1.0968, 1.1019, 1.1087, 1.1172, 1.1208, 1.1243, 1.1311, 1.1379, 1.1435, 1.1490,
                1.1545, 1.1601, 1.1656, 1.1712, 1.1767, 1.1823, 1.1878, 1.1934, 1.1990, 1.2045,
                1.2101, 1.2156, 1.2212, 1.2267, 1.2323,
            ];
        // Note: co2parm is for icrease in ambient CO2 concentration changes from 1959 (308 ppm).
        // The first 28 values (up to 1987) are from GOSSYM. The other values (up to 2004) are derived from data of the
        // Carbon Dioxide Information Analysis Center (CDIAC).
        //
        // Exit the function and end simulation if there are no leaves.

        let legacy = &mut model_state.legacy;
        let total_leaf_weight = total_leaf_weight_in_state(legacy);
        if total_leaf_weight <= 0. {
            return Err(Cotton2KError {
                level: 0,
                message: String::from("Leaf weight is less than 0!"),
            });
        }

        // Get the CO2 correction factor (pnetcor) for photosynthesis, using AmbientCO2Factor and a factor that may be
        // variety specific (vpnet[0]).

        // correction factor for gross photosynthesis.
        let pnetcor = profile.ambient_CO2_factor
            * vpnet[0]
            * profile.co2_enrichment.unwrap_or_default().factor;
        // Compute ptnfac, the effect of leaf N concentration on photosynthesis, using an empirical relationship.
        // correction factor for low nitrogen content in leaves.
        let mut ptnfac =
            vpnet[3] + (legacy.leaf_n_conc - vpnet[2]) * (1. - vpnet[3]) / (vpnet[1] - vpnet[2]);
        if ptnfac > 1. {
            ptnfac = 1.;
        }
        if ptnfac < vpnet[3] {
            ptnfac = vpnet[3];
        }
        // Convert the average daily short wave radiation from langley per day, to Watts per square meter (wattsm).
        // average daily global radiation, W m-2.
        let wattsm = GetFromClim(CLIMATE_METRIC_IRRD, legacy.daynum) * 697.45
            / (num_hours(self.atmosphere.daylength) * 60.);
        // Compute pstand as an empirical function of wattsm (based on Baker et al., 1972).
        // gross photosynthesis for a non-stressed full canopy.
        let pstand = 2.3908 + wattsm * (1.37379 - wattsm * 0.00054136);
        // Convert it to gross photosynthesis per plant (pplant), using PerPlantArea and corrections for light interception by canopy, ambient CO2 concentration, water stress and low N in the leaves.
        // actual gross photosynthetic rate, g per plant per day.
        let pplant = match profile.light_intercept_method {
            LightInterceptMethod::Latered => {
                latered_gross_photosynthesis(legacy, pstand, profile.ptsred, pnetcor, ptnfac)
            }
            _ => {
                0.001
                    * pstand
                    * legacy.light_intercept
                    * legacy.per_plant_area
                    * profile.ptsred
                    * pnetcor
                    * ptnfac
            }
        };
        // Compute the photorespiration factor (rsubl) as a linear function af average day time temperature.
        let rsubl = 0.0032125 + 0.0066875 * legacy.day_time_temp; // photorespiration factor.

        // Photorespiration (lytres) is computed as a proportion of gross photosynthetic rate.
        let lytres = rsubl * pplant; // rate of photorespiration, g per plant per day.

        // Old stems are those more than voldstm = 32 calendar days old.
        // Maintenance respiration is computed on the basis of plant dry weight, minus the old stems and the dry tissue of opened bolls.
        let voldstm = 32;
        let kkday = legacy.kday - voldstm; // day of least recent actively growing stems.

        // weight of old stems.
        let oldstmwt = if kkday < 1 {
            0.
        } else {
            legacy.stem_weight[kkday as usize]
        };
        // maintenance respiration, g per plant per day.
        let bmain = (legacy.plant_weight
            - legacy.cotton_weight_open_bolls
            - legacy.burr_weight_open_bolls
            - oldstmwt)
            * rsubo;
        // Net photosynthesis is computed by substracting photo-respiration and maintenance respiration from the gross rate of photosynthesis. To avoid computational problems, make sure that pts is positive and non-zero.
        // intermediate computation of NetPhotosynthesis.
        let mut pts = pplant - lytres - bmain;
        if pts < 0.00001 {
            pts = 0.00001;
        }
        // The growth respiration (gsubr) supplies energy for converting the supplied carbohydrates to plant tissue dry
        // matter. 0.68182 converts CO2 to CH2O. NetPhotosynthesis is the computed net photosynthesis, in g per plant
        // per day.

        legacy.net_photosynthesis = pts / (1. + gsubr) * 0.68182;
        // CumNetPhotosynth is the cumulative value of NetPhotosynthesis, from day of emergence.
        legacy.cum_net_photosynth += legacy.net_photosynthesis;
        legacy.write_to_globals();
        Ok(())
    }
    /// This function sets the water saturation of the soil layers below the water table, if it has been defined in the input.
    /// It is called from [Profile::soil_procedures()] if water table data have been input.
    ///
    /// The following global variables are referenced here:
    /// * [Daynum]
    /// * [dl]
    /// * [MaxWaterCapacity]
    /// * [nk]
    /// * [nl]
    /// * [PoreSpace]
    /// * [MaxWaterCapacity]
    /// * [RowSpace]
    /// * [wk]
    ///
    /// The following global variables are set here:
    /// * [addwtbl]
    /// * [ElCondSatSoilToday]
    /// * [WaterTableLayer]
    /// * [VolWaterContent]
    fn watertable(
        &mut self,
        profile: &Profile,
        legacy: &mut crate::LegacyGlobalState,
    ) -> Result<(), Cotton2KError> {
        if profile.num_watertable_data == 0 {
            return Ok(());
        }
        // Find the depth of water table for this day.
        let mut lwtable = 201f64; // level of water table on this day, cm
        let current_day = legacy.daynum as u32;
        for ao in &profile.agronomy_operations {
            if let AgronomyOperation::watertable(watertable) = ao {
                if current_day >= watertable.date.ordinal() {
                    lwtable = watertable.level;
                    legacy.el_cond_sat_soil_today = watertable.ecs;
                }
            }
        }

        // Find the number of the uppermost layer of water table.
        legacy.water_table_layer = find_water_table_layer(legacy, lwtable);

        // The total water entering the soil slab (addwtbl) is computed.
        // It is used to check the water balance in the soil.
        let row_space = legacy.row_space;
        let water_table_layer = legacy.water_table_layer;
        let num_cols = legacy.nk as usize;
        let nl_layers = legacy.nl as usize;
        let layer_depths = &legacy.dl;
        let col_widths = &legacy.wk;
        let pore_space = &legacy.pore_space;
        let max_water_capacity = &legacy.max_water_capacity;
        let added_water = &mut legacy.addwtbl;
        for_each_row_mut(
            &mut legacy.vol_water_content,
            nl_layers,
            |l, mut water_row| {
                let layer_depth = layer_depths[l];
                let uplimit = if l as i32 >= water_table_layer {
                    pore_space[l]
                } else {
                    max_water_capacity[l]
                };
                update_layer_water_table_saturation(
                    &mut water_row,
                    num_cols,
                    uplimit,
                    l as i32 >= water_table_layer,
                    layer_depth,
                    row_space,
                    col_widths
                        .as_slice()
                        .expect("wk should be contiguous in standard layout"),
                    added_water,
                );
            },
        );
        legacy.write_to_globals();
        Ok(())
    }
}

/// This function simulates the leaf water potential of cotton plants. It has been adapted from the model of Moshe Meron (The relation of cotton leaf water potential to soil water content in the irrigated management range. PhD dissertation, UC Davis, 1984).
///
///  It is called from [Profile::stress()]. It calls [wcond()] and [LeafResistance()].
///
/// The following global variables are referenced here:
/// * [AgeOfPreFruNode]
/// * [AverageSoilPsi]
/// * [beta]
/// * [dl]
/// * [Kday]
/// * [LeafAge]
/// * [NumFruitBranches]
/// * [NumLayersWithRoots]
/// * [NumNodes]
/// * [NumPreFruNodes]
/// * [NumVegBranches]
/// * [PlantHeight]
/// * [PoreSpace]
/// * [ReferenceETP]
/// * [RootColNumLeft]
/// * [RootColNumRight]
/// * [RootWtCapblUptake]
/// * [SaturatedHydCond]
/// * [SoilPsi]
/// * [thad]
/// * [thts]
/// * [VolWaterContent]
/// * [wk]
///
/// The following global variables are set here:
/// * [LwpMin]
/// * [LwpMax]
fn LeafWaterPotential(legacy: &mut crate::LegacyGlobalState) {
    // Constant parameters used:
    const cmg: f64 = 3200.; // length in cm per g dry weight of roots, based on an average root diameter of 0.06 cm, and a specific weight of 0.11 g  dw per cubic cm.
    const psild0: f64 = -1.32; // maximum values of LwpMin
    const psiln0: f64 = -0.40; // maximum values of LwpMax.
    const rtdiam: f64 = 0.06; // average root diameter in cm.
    const vpsil: [f64; 13] = [
        0.48, -5.0, 27000., 4000., 9200., 920., 0.000012, -0.15, -1.70, -3.5, 0.1e-9, 0.025, 0.80,
    ];
    // Leaf water potential is not computed during 10 days after emergence. Constant values are assumed for this period.
    if legacy.kday <= 10 {
        legacy.lwp_max = psiln0;
        legacy.lwp_min = psild0;
        return;
    }
    // Compute shoot resistance (rshoot) as a function of plant height.
    let rshoot: f64 = vpsil[0] * legacy.plant_height / 100.; // shoot resistance, Mpa hours per cm.

    let rroot; // root resistance, Mpa hours per cm.

    // Loop over all soil cells with roots. Check if RootWtCapblUptake is greater than vpsil[10].
    // All average values computed for the root zone, are weighted by RootWtCapblUptake (root weight capable of uptake), but the weight assigned will not be greater than vpsil[11].
    let sums = accumulate_root_zone_hydraulic_sums(legacy, &vpsil, cmg);
    // Compute average root resistance (rroot) and average soil water content (vh2).
    let dumyrs: f64; // intermediate variable for computing cond.
    let vh2: f64; // average of soil water content, for all soil soil cells with roots.
    if sums.psinum > 0. && sums.sumlv > 0. {
        rroot = sums.psinum / sums.rrlsum;
        vh2 = sums.vh2sum / sums.psinum;
        dumyrs = fmax(
            1.001,
            (1. / (std::f64::consts::PI * sums.sumlv / sums.rootvol)).sqrt() / rtdiam,
        );
    } else {
        rroot = 0.;
        vh2 = legacy.thad[0];
        dumyrs = 1.001;
    }
    // Compute hydraulic conductivity (cond), and soil resistance near the root surface (rsoil). soil hydraulic conductivity near the root surface.
    let cond = fmax(
        wcond(
            vh2,
            legacy.thad[0],
            legacy.thts[0],
            legacy.beta[0],
            legacy.saturated_hyd_cond[0],
            legacy.pore_space[0],
        ) / 24.
            * 2.
            * sums.sumlv
            / sums.rootvol
            / dumyrs.ln(),
        vpsil[6],
    );
    let rsoil = 0.0001 / (2. * std::f64::consts::PI * cond); // soil resistance, Mpa hours per cm.

    // Compute leaf resistance (LeafResistance) as the average of the resistances of all existing leaves.
    let (numl, sumrl) = leaf_resistance_stats(legacy);
    let rleaf = sumrl / numl as f64; // leaf resistance, Mpa hours per cm.

    // The total resistance to transpiration, MPa hours per cm, (rtotal) is computed.
    let rtotal = rsoil + rroot + rshoot + rleaf;
    // Compute maximum (early morning) leaf water potential, LwpMax, from soil water potential (AverageSoilPsi, converted from bars to MPa).
    // Check for minimum and maximum values.
    legacy.lwp_max = vpsil[7] + 0.1 * legacy.average_soil_psi;
    if legacy.lwp_max < vpsil[8] {
        legacy.lwp_max = vpsil[8];
    }
    if legacy.lwp_max > psiln0 {
        legacy.lwp_max = psiln0;
    }
    //     Compute minimum (at time of maximum transpiration rate) leaf water
    //     potential, LwpMin, from
    //  maximum transpiration rate (etmax) and total resistance to transpiration
    //  (rtotal).
    let etmax = legacy
        .reference_etp
        .iter()
        .take(24)
        .fold(0.0, |current_max, &hourly_etp| {
            fmax(current_max, hourly_etp)
        }); // the maximum hourly rate of evapotranspiration for this day.
    legacy.lwp_min = legacy.lwp_max - 0.1 * fmax(etmax, vpsil[12]) * rtotal;
    //     Check for minimum and maximum values.
    if legacy.lwp_min < vpsil[9] {
        legacy.lwp_min = vpsil[9];
    }
    if legacy.lwp_min > psild0 {
        legacy.lwp_min = psild0;
    }
}
/// This function computes the water redistribution in the soil after irrigation by a drip system.
/// It also computes the resulting redistribution of nitrate and urea N.
/// It is called by [`State::soil_procedures()`] `noitr` times per day.
/// It calls function CellDistrib().
///
/// The following argument is used:
/// Drip - amount of irrigation applied by the drip method, mm.
///
/// The following global variables are referenced:
/// * [dl]
/// * [LocationColumnDrip]
/// * [LocationLayerDrip]
/// * [MaxWaterCapacity]
/// * [nk]
/// * [nl]
/// * [NO3FlowFraction]
/// * [PoreSpace]
/// * [RowSpace]
/// * [wk]
///
/// The following global variables are set:
/// * [CumWaterDrained]
/// * [SoilNitrogenLoss]
/// * [VolWaterContent]
/// * [VolNo3NContent]
/// * [VolUreaNContent]
fn DripFlow(Drip: f64, legacy: &mut crate::LegacyGlobalState) -> Result<(), Cotton2KError> {
    let dims = GridDims::from_legacy(legacy);
    let mut ring_flow = DripRingFlowBuffer::new(dims.max_rings); // amount of water and nitrogen transferred between rings.

    // Incoming flow of water (Drip, in mm) is converted to dripw(0), in cm^3 per slab.
    ring_flow[0].water = Drip * legacy.row_space * 0.10;
    let l0 = legacy.location_layer_drip as usize; //  layer where the drip emitter is situated
    let k0 = legacy.location_column_drip as usize; //  column where the drip emitter is situated

    if initialize_drip_source_cell(legacy, l0, k0, &mut ring_flow) {
        return Ok(());
    }

    // the difference between the maximum water capacity (at a water content of uplimit) of this ring of soil cell, and the actual water content, cm3.
    let mut h2odef;
    // Loop of concentric rings of cells, starting from ring 1.
    // Assign zero to the sums sv, st, sn, sn1, su and su1.
    for kr in 1..dims.max_rings {
        let radius = (6 * kr) as f64; // radius (cm) of the wetting ring
        let (ring_cells, ring_totals) = collect_drip_ring_cells_and_totals(legacy, l0, k0, radius)?;
        // Compute the amount of water needed to saturate all the cells in this ring (h2odef).
        h2odef = ring_totals.st - ring_totals.sv;
        // Test if the amount of incoming flow, dripw(kr), is greater than h2odef.
        if ring_flow[kr].water <= h2odef {
            // In this case, this will be the last wetted ring.
            // Update VolWaterContent in this ring, by wetting each cell in proportion to its defcit.
            // Update VolNo3NContent and VolUreaNContent of the cells in this ring by the same proportion.
            // This is executed for all the cells in the ring.
            apply_drip_ring_partial_flow(
                legacy,
                &ring_cells,
                ring_flow[kr].water,
                ring_flow[kr].no3,
                ring_flow[kr].urea,
                h2odef,
            );
            return Ok(());
        }
        let ring_result = compute_drip_ring_result(ring_flow[kr], ring_totals, h2odef);
        // For all the cells in the ring, as in the 1st cell, saturate VolWaterContent to uplimit, and update VolNo3NContent and VolUreaNContent.
        saturate_drip_ring_cells(
            legacy,
            &ring_cells,
            ring_result.outflow.no3_loss,
            ring_totals.sn,
            ring_result.concentration.no3,
            ring_result.outflow.urea_loss,
            ring_totals.su,
            ring_result.concentration.urea,
        );
        if finalize_or_advance_drip_ring(
            legacy,
            &mut ring_flow,
            kr,
            l0,
            dims,
            ring_result.outflow.flow,
        ) {
            return Ok(());
        }
        // Repeat all these procedures for the next ring.
    }
    Ok(())
}

#[derive(Debug, Clone, Copy)]
struct RingCell {
    layer: usize,
    column: usize,
    uplimit: f64,
    deficit: f64,
}

#[derive(Debug, Clone, Copy, Default)]
struct DripRingFlow {
    water: f64,
    no3: f64,
    urea: f64,
}

#[derive(Debug, Clone, Copy)]
struct DripRingComputation {
    concentration: DripRingConcentration,
    outflow: DripRingOutflow,
}

#[derive(Debug, Clone, Copy)]
struct DripRingConcentration {
    no3: f64,
    urea: f64,
}

#[derive(Debug, Clone, Copy)]
struct DripRingOutflow {
    flow: DripRingFlow,
    no3_loss: f64,
    urea_loss: f64,
}

#[derive(Debug, Clone, Copy)]
struct GridDims {
    nl: usize,
    nk: usize,
    max_rings: usize,
}

impl GridDims {
    fn from_legacy(legacy: &crate::LegacyGlobalState) -> Self {
        Self {
            nl: legacy.nl as usize,
            nk: legacy.nk as usize,
            max_rings: maxl as usize,
        }
    }
}

#[derive(Debug, Clone)]
struct DripRingFlowBuffer {
    by_ring: Vec<DripRingFlow>,
}

impl DripRingFlowBuffer {
    fn new(size: usize) -> Self {
        Self {
            by_ring: vec![DripRingFlow::default(); size],
        }
    }
}

impl std::ops::Index<usize> for DripRingFlowBuffer {
    type Output = DripRingFlow;

    fn index(&self, index: usize) -> &Self::Output {
        &self.by_ring[index]
    }
}

impl std::ops::IndexMut<usize> for DripRingFlowBuffer {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.by_ring[index]
    }
}

fn initialize_drip_source_cell(
    legacy: &mut crate::LegacyGlobalState,
    layer: usize,
    column: usize,
    ring_flow: &mut DripRingFlowBuffer,
) -> bool {
    let h2odef = (legacy.max_water_capacity[layer] - legacy.vol_water_content[[layer, column]])
        * legacy.dl[layer]
        * legacy.wk[column];
    if ring_flow[0].water <= h2odef {
        legacy.vol_water_content[[layer, column]] +=
            ring_flow[0].water / (legacy.dl[layer] * legacy.wk[column]);
        return true;
    }
    ring_flow[1].water = ring_flow[0].water - h2odef;
    (ring_flow[1].no3, legacy.vol_no3_n_content[[layer, column]]) =
        compute_drip_cell_nutrient_outflow(
            legacy,
            layer,
            column,
            ring_flow[0].water,
            ring_flow[1].water,
            legacy.vol_no3_n_content[[layer, column]],
        );
    (
        ring_flow[1].urea,
        legacy.vol_urea_n_content[[layer, column]],
    ) = compute_drip_cell_nutrient_outflow(
        legacy,
        layer,
        column,
        ring_flow[0].water,
        ring_flow[1].water,
        legacy.vol_urea_n_content[[layer, column]],
    );
    legacy.vol_water_content[[layer, column]] = legacy.max_water_capacity[layer];
    false
}

fn compute_drip_ring_result(
    incoming: DripRingFlow,
    totals: DripRingTotals,
    h2odef: f64,
) -> DripRingComputation {
    let concentration = DripRingConcentration {
        no3: (totals.sn + incoming.no3) / (totals.sv + incoming.water),
        urea: (totals.su + incoming.urea) / (totals.sv + incoming.water),
    };
    let water_out = incoming.water - h2odef;
    let (no3_out, no3_loss) =
        limit_ring_nutrient_outflow(water_out * concentration.no3, incoming.no3, totals.sn1);
    let (urea_out, xuloss) =
        limit_ring_nutrient_outflow(water_out * concentration.urea, incoming.urea, totals.su1);
    DripRingComputation {
        concentration,
        outflow: DripRingOutflow {
            flow: DripRingFlow {
                water: water_out,
                no3: no3_out,
                urea: urea_out,
            },
            no3_loss,
            urea_loss: xuloss,
        },
    }
}

fn finalize_or_advance_drip_ring(
    legacy: &mut crate::LegacyGlobalState,
    ring_flow: &mut DripRingFlowBuffer,
    ring_idx: usize,
    source_layer: usize,
    dims: GridDims,
    outflow: DripRingFlow,
) -> bool {
    if ring_idx < (dims.nl - source_layer - 1) && ring_idx < dims.max_rings - 1 {
        ring_flow[ring_idx + 1] = outflow;
        return false;
    }
    legacy.cum_water_drained += 10. * outflow.water / legacy.row_space;
    legacy.soil_nitrogen_loss += outflow.no3 + outflow.urea;
    true
}

#[derive(Debug, Clone, Copy, Default)]
struct DripRingTotals {
    // sum of actual water content in a ring of cells, cm3
    sv: f64,
    // sum of total water capacity in a ring of cells, cm3
    st: f64,
    // sum of nitrate N content in a ring of cells, mg
    sn: f64,
    // sum of movable nitrate N content in a ring of cells, mg
    sn1: f64,
    // sum of urea N content in a ring of cells, mg
    su: f64,
    // sum of movable urea N content in a ring of cells, mg
    su1: f64,
}

fn compute_drip_cell_nutrient_outflow(
    legacy: &crate::LegacyGlobalState,
    layer: usize,
    column: usize,
    incoming_water: f64,
    outgoing_water: f64,
    nutrient_content: f64,
) -> (f64, f64) {
    if nutrient_content <= 1.0e-30 {
        return (0., nutrient_content);
    }
    let cell_volume = legacy.dl[layer] * legacy.wk[column];
    let total_water_after_inflow =
        legacy.vol_water_content[[layer, column]] + incoming_water / cell_volume;
    let concentration = nutrient_content / total_water_after_inflow;
    if concentration * legacy.max_water_capacity[layer]
        < legacy.no3_flow_fraction[layer] * nutrient_content
    {
        let outflow = legacy.no3_flow_fraction[layer] * nutrient_content * cell_volume;
        let updated_content = (1. - legacy.no3_flow_fraction[layer]) * nutrient_content;
        (outflow, updated_content)
    } else {
        let outflow = outgoing_water * concentration;
        let updated_content = legacy.max_water_capacity[layer] * concentration;
        (outflow, updated_content)
    }
}

fn limit_ring_nutrient_outflow(raw_outflow: f64, inflow: f64, movable_cap: f64) -> (f64, f64) {
    if raw_outflow <= inflow {
        return (raw_outflow, 0.);
    }
    let mut loss = raw_outflow - inflow;
    if loss > movable_cap {
        loss = movable_cap;
    }
    (inflow + loss, loss)
}

fn apply_drip_ring_partial_flow(
    legacy: &mut crate::LegacyGlobalState,
    ring_cells: &[RingCell],
    dripw: f64,
    dripn: f64,
    dripu: f64,
    h2odef: f64,
) {
    for_each_drip_ring_cell_mut(legacy, ring_cells, |legacy, cell| {
        legacy.vol_water_content[[cell.layer, cell.column]] += dripw * cell.deficit / h2odef;
        legacy.vol_no3_n_content[[cell.layer, cell.column]] += dripn * cell.deficit / h2odef;
        legacy.vol_urea_n_content[[cell.layer, cell.column]] += dripu * cell.deficit / h2odef;
    });
}

fn saturate_drip_ring_cells(
    legacy: &mut crate::LegacyGlobalState,
    ring_cells: &[RingCell],
    xnloss: f64,
    sn: f64,
    cnw: f64,
    xuloss: f64,
    su: f64,
    cuw: f64,
) {
    for_each_drip_ring_cell_mut(legacy, ring_cells, |legacy, cell| {
        legacy.vol_water_content[[cell.layer, cell.column]] = cell.uplimit;
        legacy.vol_no3_n_content[[cell.layer, cell.column]] = if xnloss <= 0. {
            cell.uplimit * cnw
        } else {
            legacy.vol_no3_n_content[[cell.layer, cell.column]] * (1. - xnloss / sn)
        };
        legacy.vol_urea_n_content[[cell.layer, cell.column]] = if xuloss <= 0. {
            cell.uplimit * cuw
        } else {
            legacy.vol_urea_n_content[[cell.layer, cell.column]] * (1. - xuloss / su)
        };
    });
}

fn for_each_drip_ring_cell_mut(
    legacy: &mut crate::LegacyGlobalState,
    ring_cells: &[RingCell],
    mut f: impl FnMut(&mut crate::LegacyGlobalState, RingCell),
) {
    for &cell in ring_cells {
        f(legacy, cell);
    }
}

fn collect_drip_ring_cells_and_totals(
    legacy: &crate::LegacyGlobalState,
    center_layer: usize,
    center_col: usize,
    radius: f64,
) -> Result<(Vec<RingCell>, DripRingTotals), Cotton2KError> {
    let dims = GridDims::from_legacy(legacy);
    let mut cells = Vec::new();
    let mut totals = DripRingTotals::default();
    let mut distance_error = None;
    for_each_cell_in_rows(1, dims.nl, dims.nk, |l, k| {
        if distance_error.is_some() {
            return;
        }
        let uplimit = if l >= legacy.water_table_layer as usize {
            legacy.pore_space[l]
        } else {
            legacy.max_water_capacity[l]
        };
        match cell_distance(l, k, center_layer, center_col, legacy.row_space) {
            Ok(dist) => {
                if dist <= radius && dist > (radius - 6.) {
                    let deficit = uplimit - legacy.vol_water_content[[l, k]];
                    cells.push(RingCell {
                        layer: l,
                        column: k,
                        uplimit,
                        deficit,
                    });
                    let area = legacy.dl[l] * legacy.wk[k];
                    totals.sv += legacy.vol_water_content[[l, k]] * area;
                    totals.st += uplimit * area;
                    totals.sn += legacy.vol_no3_n_content[[l, k]] * area;
                    totals.sn1 +=
                        legacy.vol_no3_n_content[[l, k]] * area * legacy.no3_flow_fraction[l];
                    totals.su += legacy.vol_urea_n_content[[l, k]] * area;
                    totals.su1 +=
                        legacy.vol_urea_n_content[[l, k]] * area * legacy.no3_flow_fraction[l];
                }
            }
            Err(err) => {
                distance_error = Some(err);
            }
        }
    });
    if let Some(err) = distance_error {
        return Err(err);
    }
    Ok((cells, totals))
}

fn drop_leaf_age(lai: f64) -> f64 {
    140. - 1. * lai
}
//     This function computes the weight of roots capable of uptake for all soil
//     cells.
//  It is called from [`State::soil_procedures()`].
//
//     The following global variables are referenced here:
//       nk, nl, NumLayersWithRoots, RootColNumLeft, RootColNumRight,
//       RootWeight.
//     The following global variable is set here:     RootWtCapblUptake
fn RootsCapableOfUptake(legacy: &mut crate::LegacyGlobalState) {
    // the indices for the relative capability of uptake (between 0 and 1) of water and nutrients by root age classes.
    const cuind: [f64; 3] = [1., 0.5, 0.];
    for_each_cell(legacy.nl as usize, legacy.nk as usize, |l, k| {
        legacy.root_wt_capbl_uptake[[l, k]] = 0.;
    });
    //     Loop for all soil soil cells with roots. compute for each soil cell
    //     root-weight capable
    //  of uptake (RootWtCapblUptake) as the sum of products of root weight and
    //  capability of uptake index (cuind) for each root class in it.
    for (l, k) in
        root_zone_cells_in_bounds(legacy, legacy.num_layers_with_roots as usize, |layer| {
            (
                legacy.root_col_num_left[layer],
                legacy.root_col_num_right[layer],
                true,
            )
        })
    {
        legacy.root_wt_capbl_uptake[[l, k]] = legacy
            .root_weight
            .slice(s![l, k, 0..cuind.len()])
            .iter()
            .zip(cuind.iter())
            .filter(|(&weight, _)| weight > 1.0e-15)
            .map(|(&weight, &coeff)| weight * coeff)
            .sum();
    }
}
