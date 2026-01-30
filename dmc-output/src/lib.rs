//! # dmc-output: Output Writers and Statistics for DMC
//!
//! This crate provides output handling for Diffusion Monte Carlo simulations.
//!
//! ## Features
//!
//! - Multiple output formats: CSV, JSON, TOML
//! - Console output with progress bar
//! - Statistical analysis (block averaging, error estimation)
//! - Aggregator for combining simulation steps into results
//!
//! ## Quick Start
//!
//! ```rust,no_run
//! use dmc_core::prelude::*;
//! use dmc_output::prelude::*;
//! use std::time::Instant;
//!
//! // Configure simulation
//! let config = SimulationConfig::builder()
//!     .num_walkers(5_000)
//!     .time_step(0.01)
//!     .total_steps(10_000)
//!     .equilibration_steps(1_000)
//!     .system(SystemType::Hydrogen)
//!     .algorithm(AlgorithmType::ImportanceSampled { max_offspring: 3 })
//!     .trial_wavefunction(TrialWfParams::Hydrogen { alpha: 1.0 })
//!     .build()
//!     .unwrap();
//!
//! // Configure output
//! let output_config = OutputConfig {
//!     format: OutputFormat::All,
//!     directory: "./output".to_string(),
//!     prefix: "hydrogen".to_string(),
//!     energy_interval: 10,
//! };
//!
//! // Create simulation and aggregator
//! let mut sim = Simulation::new(config).unwrap();
//! let mut aggregator = OutputAggregator::new(
//!     &output_config,
//!     sim.system_name(),
//!     sim.total_steps(),
//!     sim.exact_energy().unwrap_or(-0.5),
//!     sim.target_population(),
//! ).unwrap();
//!
//! let start = Instant::now();
//!
//! // Run simulation
//! for step in &mut sim {
//!     aggregator.process_step(&step);
//! }
//!
//! // Compute and write results
//! let result = aggregator.finalize(
//!     sim.algorithm_name(),
//!     sim.system_name(),
//!     sim.exact_energy(),
//!     start.elapsed().as_secs_f64(),
//! );
//!
//! print_result_summary(&result);
//! print_output_info(&output_config);
//! ```

pub mod statistics;
pub mod types;
pub mod writers;

mod aggregator;

pub use aggregator::OutputAggregator;

/// Convenient re-exports for common usage.
pub mod prelude {
    pub use crate::aggregator::OutputAggregator;
    pub use crate::statistics::{block_average, final_statistics};
    pub use crate::types::{BlockData, ConvergenceData, OutputConfig, OutputFormat, SimulationResult};
    pub use crate::writers::{
        print_output_info, print_result_summary, print_simulation_header, ConsoleWriter, CsvWriter,
        FileWriter, JsonWriter, OutputError, OutputWriter, TomlWriter,
    };
}
