//! # dmc-core: Diffusion Monte Carlo Library
//!
//! A Rust library for quantum ground state calculations using Diffusion Monte Carlo.
//!
//! ## Overview
//!
//! This library implements both **Pure DMC** and **Importance-Sampled DMC** algorithms
//! for calculating ground state energies of simple quantum systems:
//!
//! | System | Electrons | Nuclei | Exact Energy (Ha) |
//! |--------|-----------|--------|-------------------|
//! | Hydrogen (H) | 1 | 1 | -0.5 |
//! | H₂⁺ ion | 1 | 2 | ≈ -0.6026 |
//! | H₂ molecule | 2 | 2 | ≈ -1.1745 |
//!
//! ## Generator API (Recommended)
//!
//! The generator API provides step-by-step control over the simulation,
//! enabling real-time visualization and custom output handling.
//!
//! ```rust,no_run
//! use dmc_core::prelude::*;
//!
//! let config = SimulationConfig::builder()
//!     .num_walkers(10_000)
//!     .time_step(0.01)
//!     .total_steps(50_000)
//!     .system(SystemType::Hydrogen)
//!     .algorithm(AlgorithmType::ImportanceSampled { max_offspring: 3 })
//!     .trial_wavefunction(TrialWfParams::Hydrogen { alpha: 1.0 })
//!     .build()
//!     .unwrap();
//!
//! // Create simulation generator
//! let mut sim = Simulation::new(config).unwrap();
//!
//! // Iterate step by step
//! for step in &mut sim {
//!     if step.is_equilibrated {
//!         println!("Step {}: E = {:.6} Ha, N = {}",
//!                  step.step, step.energy_estimate, step.population_size);
//!     }
//! }
//! ```
//!
//! ## From TOML Configuration
//!
//! ```rust,no_run
//! use dmc_core::prelude::*;
//!
//! let config = SimulationConfig::from_toml("dmc_config.toml").unwrap();
//! let mut sim = Simulation::new(config).unwrap();
//!
//! for step in &mut sim {
//!     // Process each step...
//! }
//! ```
//!
//! ## Theoretical Background
//!
//! See the project [README](https://github.com/subinbg/Diffusion-MC) for comprehensive
//! theoretical derivations of all equations implemented in this library.
//!
//! ### Key Equations
//!
//! - **Imaginary-time Schrödinger equation**:
//!   ∂Ψ/∂τ = -(1/ℏ)(Ĥ - E_T)Ψ
//!
//! - **Gaussian Green's function** (diffusion):
//!   G(x←y) = (m/2πℏδτ)^(3/2) exp(-m|x-y|²/2ℏδτ)
//!
//! - **Local energy** (importance sampling):
//!   E_L(x) = ĤΨ_T(x)/Ψ_T(x) = -(ℏ²/2m)(∇²Ψ_T/Ψ_T) + V(x)
//!
//! - **Drift velocity**:
//!   v_D = (ℏ/m)∇ln|Ψ_T|

pub mod algorithm;
pub mod config;
pub mod physics;
pub mod state;
pub mod system;
pub mod walker;

/// Convenient re-exports for common usage.
pub mod prelude {
    pub use crate::config::{
        AlgorithmType, ConfigError, SimulationConfig, SimulationConfigBuilder, SystemType,
        TrialWfParams,
    };
    pub use crate::state::{Simulation, SimulationError, SimulationState, SimulationStep};
    pub use crate::system::{QuantumSystem, TrialWavefunction};
    pub use crate::walker::{Population, Walker};
}

/// 3D position vector type (atomic units: Bohr radii).
pub type Vec3 = nalgebra::Vector3<f64>;
