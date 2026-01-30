//! # examples: Example Simulations for DMC
//!
//! This crate provides example simulations and output handling for
//! Diffusion Monte Carlo simulations.
//!
//! ## Available Binaries
//!
//! - `hydrogen` - Hydrogen atom simulation
//! - `h2_ion` - H2+ ion simulation
//! - `h2_molecule` - H2 molecule simulation
//!
//! ## Features
//!
//! - Multiple output formats: CSV, JSON, TOML
//! - Console output with progress bar
//! - Statistical analysis (block averaging, error estimation)
//! - Aggregator for combining simulation steps into results
//!
//! ## Usage
//!
//! ```bash
//! # Run hydrogen simulation
//! cargo run --bin hydrogen
//!
//! # Run H2+ ion simulation
//! cargo run --bin h2_ion
//!
//! # Run H2 molecule simulation
//! cargo run --bin h2_molecule
//!
//! # Run H2 molecule with custom config
//! cargo run --bin h2_molecule -- config.toml
//! ```

pub mod statistics;
pub mod types;
pub mod writers;

mod aggregator;

pub use aggregator::OutputAggregator;

// Re-export simulation types for convenience
pub use simulation::prelude::*;

/// Convenient re-exports for common usage.
pub mod prelude {
    pub use crate::aggregator::OutputAggregator;
    pub use crate::statistics::{block_average, final_statistics};
    pub use crate::types::{BlockData, ConvergenceData, OutputConfig, OutputFormat, SimulationResult};
    pub use crate::writers::{
        print_output_info, print_result_summary, print_simulation_header, ConsoleWriter, CsvWriter,
        FileWriter, JsonWriter, OutputError, OutputWriter, TomlWriter,
    };
    // Also re-export simulation prelude
    pub use simulation::prelude::*;
}
