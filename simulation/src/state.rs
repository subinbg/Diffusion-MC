//! Simulation state and iterator for generator-based DMC.
//!
//! This module provides a streaming/generator API for DMC simulations,
//! allowing step-by-step control and real-time monitoring.
//!
//! # Example
//!
//! ```rust,no_run
//! use simulation::prelude::*;
//!
//! let config = SimulationConfig::builder()
//!     .num_walkers(1000)
//!     .time_step(0.01)
//!     .total_steps(10_000)
//!     .system(SystemType::Hydrogen)
//!     .algorithm(AlgorithmType::ImportanceSampled { max_offspring: 3 })
//!     .trial_wavefunction(TrialWfParams::Hydrogen { alpha: 1.0 })
//!     .build()
//!     .unwrap();
//!
//! let mut sim = Simulation::new(config).unwrap();
//!
//! for step in &mut sim {
//!     println!("Step {}: E = {:.6} Ha", step.step, step.energy_estimate);
//! }
//! ```

use crate::algorithm::{DmcAlgorithm, ImportanceSampledDmc, PureDmc};
use crate::config::{AlgorithmType, SimulationConfig, SystemType};
use crate::physics::update_trial_energy;
use crate::walker::Population;
use crate::Vec3;
use serde::{Deserialize, Serialize};

/// Data yielded for each simulation step.
///
/// Contains all information about a single imaginary-time step,
/// suitable for real-time visualization and statistics collection.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SimulationStep {
    /// Current step number (0-indexed).
    pub step: usize,

    /// Imaginary time τ = step × δτ.
    pub imaginary_time: f64,

    /// Energy estimate from this step (average of local energies).
    pub energy_estimate: f64,

    /// Current trial energy E_T (used for population control).
    pub trial_energy: f64,

    /// Current population size (number of walkers).
    pub population_size: usize,

    /// Acceptance ratio for this step (importance sampling only).
    pub acceptance_ratio: Option<f64>,

    /// Local energies of all walkers (for variance calculation).
    pub local_energies: Vec<f64>,

    /// Walker positions (for density visualization).
    /// Each inner Vec contains positions for one walker's electrons.
    pub walker_positions: Vec<Vec<Vec3>>,

    /// Whether this step is past the equilibration phase.
    pub is_equilibrated: bool,
}

/// Factory for creating simulation generators.
pub struct Simulation;

impl Simulation {
    /// Create a new simulation generator from configuration.
    ///
    /// Returns a `SimulationState` that implements `Iterator<Item = SimulationStep>`.
    pub fn new(config: SimulationConfig) -> Result<SimulationState, SimulationError> {
        SimulationState::new(config)
    }
}

/// Error type for simulation initialization.
#[derive(Debug, thiserror::Error)]
pub enum SimulationError {
    #[error("Configuration error: {0}")]
    Config(String),

    #[error("Initialization error: {0}")]
    Init(String),
}

/// Mutable simulation state that yields `SimulationStep` on each iteration.
///
/// Implements `Iterator` for convenient streaming access.
pub struct SimulationState {
    /// Current step number.
    step: usize,

    /// Current trial energy E_T.
    trial_energy: f64,

    /// Walker population.
    population: Population,

    /// DMC algorithm instance.
    algorithm: Box<dyn DmcAlgorithm>,

    /// Simulation configuration.
    config: SimulationConfig,

    /// Number of electrons in the system.
    num_electrons: usize,

    /// Exact energy for the system (if known).
    exact_energy: Option<f64>,

    /// System name for display.
    system_name: String,
}

impl SimulationState {
    /// Create a new simulation state.
    fn new(config: SimulationConfig) -> Result<Self, SimulationError> {
        // Create algorithm
        let algorithm: Box<dyn DmcAlgorithm> = match &config.algorithm {
            AlgorithmType::Pure { .. } => Box::new(PureDmc::new(&config)),
            AlgorithmType::ImportanceSampled { .. } => Box::new(ImportanceSampledDmc::new(&config)),
        };

        // Get system info
        let num_electrons = match &config.system {
            SystemType::Hydrogen => 1,
            SystemType::H2Ion { .. } => 1,
            SystemType::H2Molecule { .. } => 2,
        };

        let exact_energy = match &config.system {
            SystemType::Hydrogen => Some(-0.5),
            SystemType::H2Ion { .. } => Some(-0.6026),
            SystemType::H2Molecule { .. } => Some(-1.1745),
        };

        let system_name = match &config.system {
            SystemType::Hydrogen => "Hydrogen (H)".to_string(),
            SystemType::H2Ion { bond_length } => format!("H₂⁺ (R={:.2})", bond_length),
            SystemType::H2Molecule { bond_length } => format!("H₂ (R={:.2})", bond_length),
        };

        // Initialize population
        let mut population = Population::new(config.num_walkers, config.seed);
        population.initialize_gaussian(num_electrons, 1.0);

        // Initialize trial energy
        let trial_energy = exact_energy.unwrap_or(-0.5);

        Ok(Self {
            step: 0,
            trial_energy,
            population,
            algorithm,
            config,
            num_electrons,
            exact_energy,
            system_name,
        })
    }

    /// Get the current step number.
    pub fn current_step(&self) -> usize {
        self.step
    }

    /// Get the current trial energy.
    pub fn current_trial_energy(&self) -> f64 {
        self.trial_energy
    }

    /// Get the current population size.
    pub fn population_size(&self) -> usize {
        self.population.size()
    }

    /// Check if the simulation is past equilibration.
    pub fn is_equilibrated(&self) -> bool {
        self.step >= self.config.equilibration_steps
    }

    /// Check if the simulation is finished.
    pub fn is_finished(&self) -> bool {
        self.step >= self.config.total_steps
    }

    /// Get the total number of steps.
    pub fn total_steps(&self) -> usize {
        self.config.total_steps
    }

    /// Get the equilibration steps count.
    pub fn equilibration_steps(&self) -> usize {
        self.config.equilibration_steps
    }

    /// Get the exact energy for this system (if known).
    pub fn exact_energy(&self) -> Option<f64> {
        self.exact_energy
    }

    /// Get the system name.
    pub fn system_name(&self) -> &str {
        &self.system_name
    }

    /// Get the algorithm name.
    pub fn algorithm_name(&self) -> &str {
        self.algorithm.name()
    }

    /// Get the number of electrons.
    pub fn num_electrons(&self) -> usize {
        self.num_electrons
    }

    /// Get the time step.
    pub fn time_step(&self) -> f64 {
        self.config.time_step
    }

    /// Get the target population size.
    pub fn target_population(&self) -> usize {
        self.config.num_walkers
    }

    /// Get a reference to the configuration.
    pub fn config(&self) -> &SimulationConfig {
        &self.config
    }

    /// Get a reference to the population (for advanced access).
    pub fn population(&self) -> &Population {
        &self.population
    }

    /// Reset the simulation to initial state.
    pub fn reset(&mut self) {
        self.step = 0;
        self.trial_energy = self.exact_energy.unwrap_or(-0.5);
        self.population = Population::new(self.config.num_walkers, self.config.seed);
        self.population.initialize_gaussian(self.num_electrons, 1.0);
    }

    /// Run a single step and return the result.
    ///
    /// This is called internally by the iterator but can also be used
    /// for manual step-by-step control.
    pub fn advance(&mut self) -> Option<SimulationStep> {
        if self.step >= self.config.total_steps {
            return None;
        }

        // Perform one DMC step
        let result = self
            .algorithm
            .step(&mut self.population, self.trial_energy, &self.config);

        // Collect walker positions before updating trial energy
        let walker_positions: Vec<Vec<Vec3>> = self
            .population
            .iter()
            .map(|w| w.positions.clone())
            .collect();

        // Update trial energy using population control
        self.trial_energy = update_trial_energy(
            result.energy_estimate,
            result.population_size,
            self.config.num_walkers,
            self.config.feedback_alpha,
        );

        let is_equilibrated = self.step >= self.config.equilibration_steps;

        let step_data = SimulationStep {
            step: self.step,
            imaginary_time: self.step as f64 * self.config.time_step,
            energy_estimate: result.energy_estimate,
            trial_energy: self.trial_energy,
            population_size: result.population_size,
            acceptance_ratio: result.acceptance_ratio,
            local_energies: result.local_energies,
            walker_positions,
            is_equilibrated,
        };

        self.step += 1;

        Some(step_data)
    }
}

impl Iterator for SimulationState {
    type Item = SimulationStep;

    fn next(&mut self) -> Option<Self::Item> {
        self.advance()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::{SimulationConfigBuilder, TrialWfParams};

    #[test]
    fn simulation_iterator_yields_correct_steps() {
        let config = SimulationConfigBuilder::default()
            .num_walkers(100)
            .time_step(0.01)
            .total_steps(100)
            .equilibration_steps(10)
            .seed(42)
            .system(SystemType::Hydrogen)
            .algorithm(AlgorithmType::ImportanceSampled { max_offspring: 3 })
            .trial_wavefunction(TrialWfParams::Hydrogen { alpha: 1.0 })
            .build()
            .unwrap();

        let mut sim = Simulation::new(config).unwrap();
        let mut count = 0;
        let mut last_step = 0;

        for step in &mut sim {
            assert_eq!(step.step, count);
            assert!(!step.walker_positions.is_empty());
            last_step = step.step;
            count += 1;
        }

        assert_eq!(count, 100);
        assert_eq!(last_step, 99);
    }

    #[test]
    fn simulation_reset_works() {
        let config = SimulationConfigBuilder::default()
            .num_walkers(100)
            .time_step(0.01)
            .total_steps(50)
            .equilibration_steps(5)
            .seed(42)
            .system(SystemType::Hydrogen)
            .algorithm(AlgorithmType::ImportanceSampled { max_offspring: 3 })
            .trial_wavefunction(TrialWfParams::Hydrogen { alpha: 1.0 })
            .build()
            .unwrap();

        let mut sim = Simulation::new(config).unwrap();

        // Run some steps
        for _ in 0..25 {
            sim.next();
        }
        assert_eq!(sim.current_step(), 25);

        // Reset
        sim.reset();
        assert_eq!(sim.current_step(), 0);
        assert!(!sim.is_finished());

        // Can run again
        let step = sim.next().unwrap();
        assert_eq!(step.step, 0);
    }
}
