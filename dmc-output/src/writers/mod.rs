//! Output writers for DMC simulations.
//!
//! This module provides various writers for outputting simulation results
//! in different formats (CSV, JSON, TOML) and to different destinations
//! (files, console).

mod console;
mod csv;
mod file;
mod json;
mod toml_writer;

pub use console::{print_output_info, print_result_summary, print_simulation_header, ConsoleWriter};
pub use csv::CsvWriter;
pub use file::FileWriter;
pub use json::JsonWriter;
pub use toml_writer::TomlWriter;

use crate::types::SimulationResult;
use dmc_core::state::SimulationStep;

/// Trait for output writers.
///
/// Writers can process simulation steps and finalize output at the end.
pub trait OutputWriter {
    /// Process a single simulation step.
    fn write_step(&mut self, step: &SimulationStep) -> Result<(), OutputError>;

    /// Finalize the output with the complete simulation result.
    fn finalize(&mut self, result: &SimulationResult) -> Result<(), OutputError>;

    /// Flush any buffered output.
    fn flush(&mut self) -> Result<(), OutputError>;
}

/// Output error type.
#[derive(Debug, thiserror::Error)]
pub enum OutputError {
    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),

    #[error("Serialization error: {0}")]
    Serialization(String),
}
