//! Pure Diffusion Monte Carlo algorithm.
//!
//! # Algorithm
//!
//! Pure DMC directly simulates the imaginary-time Schrödinger equation:
//!
//! 1. **Diffusion**: x' = x + √(δτ) × ξ, where ξ ~ N(0,1)
//! 2. **Branching**: M = ⌊W + u⌋, where W = exp(-δτ(V(x')-E_T))
//!
//! # Warning
//!
//! Pure DMC is numerically unstable for systems with Coulomb singularities
//! (V → -∞). The branching weight W → ∞ near nuclei, causing population explosion.
//!
//! Use Importance-Sampled DMC for atomic/molecular systems.
//!
//! See [README: Pure Diffusion Monte Carlo](https://github.com/subinbg/Diffusion-MC#pure-diffusion-monte-carlo)

use crate::algorithm::{DmcAlgorithm, StepResult};
use crate::config::{SimulationConfig, SystemType};
use crate::physics::branching_weight_pure;
use crate::system::{H2Ion, H2Molecule, Hydrogen, QuantumSystem};
use crate::walker::Population;

/// Pure DMC algorithm implementation.
///
/// # Algorithm Steps
///
/// For each step:
/// 1. Diffuse all walkers (Gaussian random walk)
/// 2. Calculate potential V(x) for each walker
/// 3. Compute branching weights W = exp(-δτ(V-E_T))
/// 4. Apply stochastic branching (birth/death)
pub struct PureDmc {
    /// Maximum offspring per walker (prevents explosion).
    pub max_offspring: usize,
    /// The quantum system being simulated.
    system: Box<dyn QuantumSystem>,
}

impl PureDmc {
    /// Create a new Pure DMC algorithm for the given system.
    pub fn new(config: &SimulationConfig) -> Self {
        let max_offspring = match &config.algorithm {
            crate::config::AlgorithmType::Pure { max_offspring } => *max_offspring,
            _ => 3,
        };

        let system: Box<dyn QuantumSystem> = match &config.system {
            SystemType::Hydrogen => Box::new(Hydrogen::new()),
            SystemType::H2Ion { bond_length } => Box::new(H2Ion::new(*bond_length)),
            SystemType::H2Molecule { bond_length } => Box::new(H2Molecule::new(*bond_length)),
        };

        Self {
            max_offspring,
            system,
        }
    }

    /// Get reference to the quantum system.
    pub fn system(&self) -> &dyn QuantumSystem {
        self.system.as_ref()
    }
}

impl DmcAlgorithm for PureDmc {
    fn step(
        &mut self,
        population: &mut Population,
        trial_energy: f64,
        config: &SimulationConfig,
    ) -> StepResult {
        let time_step = config.time_step;

        // Step 1: Diffusion
        // x' = x + √(δτ) × ξ
        // See README: Step 2: Diffusion (Kinetic Update)
        population.diffuse(time_step);

        // Step 2: Calculate potentials and local energies
        let mut potentials = Vec::with_capacity(population.size());
        for walker in population.iter() {
            let v = self.system.potential(&walker.positions);
            potentials.push(v);
        }

        // Step 3: Calculate branching weights
        // W = exp(-δτ(V(x) - E_T))
        // See README: Step 3: Weighting and Branching
        let weights: Vec<f64> = potentials
            .iter()
            .map(|&v| branching_weight_pure(v, trial_energy, time_step))
            .collect();

        // Step 4: Branching
        // M = ⌊W + u⌋
        population.branch(&weights, self.max_offspring);
        population.increment_ages();

        // Energy estimate is average potential
        let avg_potential = potentials.iter().sum::<f64>() / potentials.len() as f64;

        StepResult {
            energy_estimate: avg_potential,
            population_size: population.size(),
            acceptance_ratio: None, // Pure DMC has no Metropolis step
            local_energies: potentials,
        }
    }

    fn name(&self) -> &'static str {
        "Pure DMC"
    }
}
