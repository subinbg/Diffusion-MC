//! Configuration for DMC simulations.
//!
//! Supports both TOML file configuration and programmatic builder API.
//!
//! # TOML Example
//!
//! ```toml
//! [simulation]
//! num_walkers = 10000
//! time_step = 0.01
//! total_steps = 100000
//! equilibration_steps = 1000
//! seed = 42
//!
//! [algorithm]
//! type = "importance_sampled"
//! max_offspring = 3
//! feedback_alpha = 1.0
//!
//! [system]
//! type = "hydrogen"
//! ```
//!
//! # Builder Example
//!
//! ```rust,no_run
//! use simulation::config::{SimulationConfig, SystemType, AlgorithmType, TrialWfParams};
//!
//! let config = SimulationConfig::builder()
//!     .num_walkers(10_000)
//!     .time_step(0.01)
//!     .total_steps(100_000)
//!     .system(SystemType::Hydrogen)
//!     .algorithm(AlgorithmType::ImportanceSampled { max_offspring: 3 })
//!     .trial_wavefunction(TrialWfParams::Hydrogen { alpha: 1.0 })
//!     .build()
//!     .unwrap();
//! ```

use serde::{Deserialize, Serialize};
use std::path::Path;
use thiserror::Error;

/// Configuration errors.
#[derive(Error, Debug)]
pub enum ConfigError {
    #[error("Missing required field: {0}")]
    MissingField(&'static str),

    #[error("Invalid value: {0}")]
    InvalidValue(&'static str),

    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),

    #[error("TOML parse error: {0}")]
    TomlParse(#[from] toml::de::Error),
}

/// Complete simulation configuration.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SimulationConfig {
    /// Number of walkers (N₀ in population control).
    pub num_walkers: usize,

    /// Imaginary time step δτ in atomic units.
    ///
    /// # Equation
    /// Diffusion variance: σ² = ℏδτ/m = δτ (in atomic units)
    pub time_step: f64,

    /// Total number of simulation steps.
    pub total_steps: usize,

    /// Steps before collecting statistics (equilibration).
    pub equilibration_steps: usize,

    /// Random seed for reproducibility.
    pub seed: u64,

    /// Algorithm type (Pure DMC or Importance Sampled).
    pub algorithm: AlgorithmType,

    /// Quantum system to simulate.
    pub system: SystemType,

    /// Trial wavefunction parameters (for importance sampling).
    pub trial_wavefunction: Option<TrialWfParams>,

    /// Population control feedback parameter α.
    ///
    /// # Equation
    /// E_T = ⟨E_L⟩ - α ln(N/N₀)
    pub feedback_alpha: f64,
}

impl SimulationConfig {
    /// Create a new builder for SimulationConfig.
    pub fn builder() -> SimulationConfigBuilder {
        SimulationConfigBuilder::default()
    }

    /// Load configuration from a TOML file.
    pub fn from_toml<P: AsRef<Path>>(path: P) -> Result<Self, ConfigError> {
        let contents = std::fs::read_to_string(path)?;
        let toml_config: TomlConfig = toml::from_str(&contents)?;
        toml_config.into_simulation_config()
    }
}

/// Algorithm type selection.
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(tag = "type", rename_all = "snake_case")]
pub enum AlgorithmType {
    /// Pure DMC: Gaussian diffusion + V-based branching.
    ///
    /// # Equation
    /// W = exp(-δτ(V(x) - E_T))
    Pure {
        /// Maximum offspring per walker.
        max_offspring: usize,
    },

    /// Importance-sampled DMC: drift-diffusion + E_L-based branching.
    ///
    /// # Equation
    /// W = exp(-δτ(E_L(x) - E_T))
    ImportanceSampled {
        /// Maximum offspring per walker.
        max_offspring: usize,
    },
}

impl Default for AlgorithmType {
    fn default() -> Self {
        Self::Pure { max_offspring: 3 }
    }
}

/// Quantum system type selection.
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(tag = "type", rename_all = "snake_case")]
pub enum SystemType {
    /// Hydrogen atom (1 electron, 1 nucleus).
    Hydrogen,

    /// Hydrogen molecule ion H₂⁺ (1 electron, 2 nuclei).
    H2Ion {
        /// Internuclear distance in Bohr radii.
        bond_length: f64,
    },

    /// Hydrogen molecule H₂ (2 electrons, 2 nuclei).
    H2Molecule {
        /// Internuclear distance in Bohr radii.
        bond_length: f64,
    },
}

impl Default for SystemType {
    fn default() -> Self {
        Self::Hydrogen
    }
}

/// Trial wavefunction parameters.
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(tag = "type", rename_all = "snake_case")]
pub enum TrialWfParams {
    /// Hydrogen: Ψ_T = exp(-αr)
    Hydrogen {
        /// Variational parameter (cusp: α = 1).
        alpha: f64,
    },

    /// H₂⁺: Ψ_T = exp(-αr₁) + exp(-αr₂)
    H2Ion {
        /// Variational parameter.
        alpha: f64,
    },

    /// H₂: Heitler-London × Jastrow
    H2Molecule {
        /// Variational parameter for atomic orbitals.
        alpha: f64,
        /// Jastrow correlation parameter.
        jastrow_b: f64,
    },
}

/// Builder for SimulationConfig.
#[derive(Default)]
pub struct SimulationConfigBuilder {
    num_walkers: Option<usize>,
    time_step: Option<f64>,
    total_steps: Option<usize>,
    equilibration_steps: Option<usize>,
    seed: Option<u64>,
    algorithm: Option<AlgorithmType>,
    system: Option<SystemType>,
    trial_wavefunction: Option<TrialWfParams>,
    feedback_alpha: Option<f64>,
}

impl SimulationConfigBuilder {
    /// Set number of walkers (N₀).
    ///
    /// # Population Control
    ///
    /// This is the target population size used in the feedback law:
    /// ```text
    /// E_T(τ) = ⟨E_L⟩ - α ln(N/N₀)
    /// ```
    pub fn num_walkers(mut self, n: usize) -> Self {
        self.num_walkers = Some(n);
        self
    }

    /// Set imaginary time step (δτ).
    ///
    /// # Diffusion Variance
    ///
    /// The time step determines the Gaussian diffusion width:
    /// ```text
    /// σ² = ℏδτ/m = δτ  (in atomic units)
    /// ```
    ///
    /// Typical values: 0.001 - 0.1 atomic units.
    pub fn time_step(mut self, dt: f64) -> Self {
        self.time_step = Some(dt);
        self
    }

    /// Set total number of simulation steps.
    pub fn total_steps(mut self, n: usize) -> Self {
        self.total_steps = Some(n);
        self
    }

    /// Set equilibration steps (before collecting statistics).
    pub fn equilibration_steps(mut self, n: usize) -> Self {
        self.equilibration_steps = Some(n);
        self
    }

    /// Set random seed for reproducibility.
    pub fn seed(mut self, s: u64) -> Self {
        self.seed = Some(s);
        self
    }

    /// Set algorithm type.
    pub fn algorithm(mut self, alg: AlgorithmType) -> Self {
        self.algorithm = Some(alg);
        self
    }

    /// Set quantum system.
    pub fn system(mut self, sys: SystemType) -> Self {
        self.system = Some(sys);
        self
    }

    /// Set trial wavefunction parameters.
    pub fn trial_wavefunction(mut self, wf: TrialWfParams) -> Self {
        self.trial_wavefunction = Some(wf);
        self
    }

    /// Set population control feedback parameter α.
    ///
    /// # Stability
    ///
    /// The feedback law E_T = ⟨E_L⟩ - α ln(N/N₀)
    /// has local stability with decay rate δ ∝ exp(-ατ/ℏ).
    pub fn feedback_alpha(mut self, alpha: f64) -> Self {
        self.feedback_alpha = Some(alpha);
        self
    }

    /// Build the final configuration, validating all parameters.
    pub fn build(self) -> Result<SimulationConfig, ConfigError> {
        let num_walkers = self
            .num_walkers
            .ok_or(ConfigError::MissingField("num_walkers"))?;
        let time_step = self
            .time_step
            .ok_or(ConfigError::MissingField("time_step"))?;

        if time_step <= 0.0 {
            return Err(ConfigError::InvalidValue("time_step must be positive"));
        }

        if num_walkers == 0 {
            return Err(ConfigError::InvalidValue("num_walkers must be > 0"));
        }

        let algorithm = self.algorithm.unwrap_or_default();

        // Importance sampling requires a trial wavefunction
        if matches!(algorithm, AlgorithmType::ImportanceSampled { .. })
            && self.trial_wavefunction.is_none()
        {
            return Err(ConfigError::InvalidValue(
                "Importance-sampled DMC requires a trial wavefunction",
            ));
        }

        Ok(SimulationConfig {
            num_walkers,
            time_step,
            total_steps: self.total_steps.unwrap_or(10_000),
            equilibration_steps: self.equilibration_steps.unwrap_or(1_000),
            seed: self.seed.unwrap_or_else(rand::random),
            algorithm,
            system: self.system.unwrap_or_default(),
            trial_wavefunction: self.trial_wavefunction,
            feedback_alpha: self.feedback_alpha.unwrap_or(1.0),
        })
    }
}

/// Internal TOML configuration structure.
#[derive(Deserialize)]
struct TomlConfig {
    simulation: TomlSimulation,
    algorithm: TomlAlgorithm,
    system: TomlSystem,
}

#[derive(Deserialize)]
struct TomlSimulation {
    num_walkers: usize,
    time_step: f64,
    #[serde(default = "default_total_steps")]
    total_steps: usize,
    #[serde(default = "default_equilibration_steps")]
    equilibration_steps: usize,
    #[serde(default)]
    seed: Option<u64>,
}

fn default_total_steps() -> usize {
    10_000
}
fn default_equilibration_steps() -> usize {
    1_000
}

#[derive(Deserialize)]
struct TomlAlgorithm {
    #[serde(rename = "type")]
    alg_type: String,
    #[serde(default = "default_max_offspring")]
    max_offspring: usize,
    #[serde(default = "default_feedback_alpha")]
    feedback_alpha: f64,
}

fn default_max_offspring() -> usize {
    3
}
fn default_feedback_alpha() -> f64 {
    1.0
}

#[derive(Deserialize)]
struct TomlSystem {
    #[serde(rename = "type")]
    sys_type: String,
    #[serde(default)]
    bond_length: Option<f64>,
    #[serde(default)]
    trial_wavefunction: Option<TomlTrialWf>,
}

#[derive(Deserialize)]
struct TomlTrialWf {
    #[serde(default = "default_alpha")]
    alpha: f64,
    #[serde(default)]
    jastrow_b: Option<f64>,
}

fn default_alpha() -> f64 {
    1.0
}

impl TomlConfig {
    fn into_simulation_config(self) -> Result<SimulationConfig, ConfigError> {
        let algorithm = match self.algorithm.alg_type.as_str() {
            "pure" => AlgorithmType::Pure {
                max_offspring: self.algorithm.max_offspring,
            },
            "importance_sampled" => AlgorithmType::ImportanceSampled {
                max_offspring: self.algorithm.max_offspring,
            },
            _ => {
                return Err(ConfigError::InvalidValue(
                    "algorithm type must be 'pure' or 'importance_sampled'",
                ))
            }
        };

        let system = match self.system.sys_type.as_str() {
            "hydrogen" => SystemType::Hydrogen,
            "h2_ion" => SystemType::H2Ion {
                bond_length: self.system.bond_length.unwrap_or(2.0),
            },
            "h2_molecule" => SystemType::H2Molecule {
                bond_length: self.system.bond_length.unwrap_or(1.4),
            },
            _ => {
                return Err(ConfigError::InvalidValue(
                    "system type must be 'hydrogen', 'h2_ion', or 'h2_molecule'",
                ))
            }
        };

        let trial_wavefunction = self.system.trial_wavefunction.map(|wf| match &system {
            SystemType::Hydrogen => TrialWfParams::Hydrogen { alpha: wf.alpha },
            SystemType::H2Ion { .. } => TrialWfParams::H2Ion { alpha: wf.alpha },
            SystemType::H2Molecule { .. } => TrialWfParams::H2Molecule {
                alpha: wf.alpha,
                jastrow_b: wf.jastrow_b.unwrap_or(0.5),
            },
        });

        Ok(SimulationConfig {
            num_walkers: self.simulation.num_walkers,
            time_step: self.simulation.time_step,
            total_steps: self.simulation.total_steps,
            equilibration_steps: self.simulation.equilibration_steps,
            seed: self.simulation.seed.unwrap_or_else(rand::random),
            algorithm,
            system,
            trial_wavefunction,
            feedback_alpha: self.algorithm.feedback_alpha,
        })
    }
}
