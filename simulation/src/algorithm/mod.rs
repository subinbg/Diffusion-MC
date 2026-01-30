//! DMC algorithm implementations.
//!
//! This module provides both Pure DMC and Importance-Sampled DMC algorithms.
//!
//! # Algorithm Comparison
//!
//! | Feature | Pure DMC | Importance Sampled |
//! |---------|----------|-------------------|
//! | Branching | W = exp(-δτ(V-E_T)) | W = exp(-δτ(E_L-E_T)) |
//! | Stability | Unstable at Coulomb singularities | Stable with cusp condition |
//! | Drift | None (isotropic diffusion) | v_D = (ℏ/m)∇ln Ψ_T |
//! | Detailed Balance | Automatic | Requires Metropolis correction |
//!
//! See [README: Summary Table](https://github.com/subinbg/Diffusion-MC#summary-pure-dmc-vs-importance-sampled-dmc)

mod pure_dmc;
mod importance;

pub use pure_dmc::PureDmc;
pub use importance::ImportanceSampledDmc;

use crate::config::SimulationConfig;
use crate::walker::Population;

/// Result of a single DMC step.
#[derive(Clone, Debug)]
pub struct StepResult {
    /// Energy estimate from this step.
    pub energy_estimate: f64,
    /// Current population size.
    pub population_size: usize,
    /// Acceptance ratio (for importance sampling only).
    pub acceptance_ratio: Option<f64>,
    /// Local energies of all walkers (for statistics).
    pub local_energies: Vec<f64>,
}

/// Trait for DMC algorithm variants.
pub trait DmcAlgorithm: Send + Sync {
    /// Perform a single imaginary time step.
    ///
    /// # Equation (Trotter-Suzuki Decomposition)
    ///
    /// ```text
    /// |Ψ(τ+δτ)⟩ = exp(-δτ(V̂-E_T)/ℏ) exp(-δτT̂/ℏ) |Ψ(τ)⟩ + O(δτ²)
    /// ```
    ///
    /// This splits the evolution into:
    /// 1. Kinetic step (diffusion)
    /// 2. Potential step (branching)
    ///
    /// See [README: Trotter-Suzuki Decomposition](https://github.com/subinbg/Diffusion-MC#trotter-suzuki-decomposition)
    fn step(
        &mut self,
        population: &mut Population,
        trial_energy: f64,
        config: &SimulationConfig,
    ) -> StepResult;

    /// Name of the algorithm (for output labeling).
    fn name(&self) -> &'static str;
}
