//! JSON output writer for DMC simulations.

use super::{OutputError, OutputWriter};
use crate::types::SimulationResult;
use simulation::state::SimulationStep;
use std::fs::{self, File};
use std::path::Path;

/// JSON output writer.
///
/// Writes complete simulation results to JSON format.
pub struct JsonWriter {
    directory: String,
    prefix: String,
}

impl JsonWriter {
    /// Create a new JSON writer.
    ///
    /// # Arguments
    ///
    /// * `directory` - Output directory
    /// * `prefix` - File name prefix
    pub fn new(directory: &str, prefix: &str) -> Result<Self, OutputError> {
        fs::create_dir_all(directory)?;

        Ok(Self {
            directory: directory.to_string(),
            prefix: prefix.to_string(),
        })
    }
}

impl OutputWriter for JsonWriter {
    fn write_step(&mut self, _step: &SimulationStep) -> Result<(), OutputError> {
        // JSON writer only outputs at finalization
        Ok(())
    }

    fn finalize(&mut self, result: &SimulationResult) -> Result<(), OutputError> {
        let path = Path::new(&self.directory).join(format!("{}_result.json", self.prefix));
        let file = File::create(path)?;
        serde_json::to_writer_pretty(file, result)
            .map_err(|e| OutputError::Serialization(e.to_string()))?;
        Ok(())
    }

    fn flush(&mut self) -> Result<(), OutputError> {
        Ok(())
    }
}
