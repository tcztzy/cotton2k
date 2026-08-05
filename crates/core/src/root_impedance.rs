//! Lookup tables for soil mechanical impedance during root growth.
//!
//! [`RootImpedanceTables`] stores monotonic water-content and bulk-density axes
//! together with the corresponding impedance curves. The plant crate owns the
//! interpolation algorithm; this core type owns only caller-provided data and
//! reports emptiness through [`RootImpedanceTables::is_empty`].

#[derive(Debug, Clone, Default, PartialEq)]
/// Axes and curves for interpolating soil mechanical impedance.
pub struct RootImpedanceTables {
    /// Soil volumetric-water-content axis.
    pub water_content: Vec<f64>,
    /// Soil bulk-density axis.
    pub bulk_density: Vec<f64>,
    /// Impedance curves indexed by bulk-density axis.
    pub impedance: Vec<Vec<f64>>,
}

impl RootImpedanceTables {
    /// Builds a lookup table from caller-owned axes and curves.
    pub fn new(water_content: Vec<f64>, bulk_density: Vec<f64>, impedance: Vec<Vec<f64>>) -> Self {
        Self {
            water_content,
            bulk_density,
            impedance,
        }
    }

    /// Reports whether any required lookup axis or curve is absent.
    pub fn is_empty(&self) -> bool {
        self.water_content.is_empty() || self.bulk_density.is_empty() || self.impedance.is_empty()
    }
}
