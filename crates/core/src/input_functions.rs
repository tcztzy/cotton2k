//! Small scalar parameterization helpers used while initializing model inputs.
//!
//! The functions here encode the compact legacy equations that turn cultivar or
//! site coefficients into initialized values. They are pure calculations with
//! no I/O, shared-state access, or allocation.

/// Evaluates the legacy three-coefficient initialization equation.
pub fn form(c0: f64, d0: f64, g0: f64) -> f64 {
    (2.0 / (1.0 + (c0 / d0 - 1.0) * g0) + 1.0 / (1.0 + (c0 / d0 - 1.0) * (1.0 - 2.0 * g0))) / 3.0
}
