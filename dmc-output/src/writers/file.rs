//! Combined file writer that manages multiple output formats.

use super::{CsvWriter, JsonWriter, OutputError, OutputWriter, TomlWriter};
use crate::types::{OutputConfig, OutputFormat, SimulationResult};
use dmc_core::state::SimulationStep;

/// Combined file writer that manages multiple output formats.
///
/// Automatically creates appropriate writers based on the output configuration.
pub struct FileWriter {
    writers: Vec<Box<dyn OutputWriter>>,
}

impl FileWriter {
    /// Create a new file writer from output configuration.
    pub fn new(config: &OutputConfig) -> Result<Self, OutputError> {
        let mut writers: Vec<Box<dyn OutputWriter>> = Vec::new();

        if config.format.includes_csv() {
            writers.push(Box::new(CsvWriter::new(
                &config.directory,
                &config.prefix,
                config.energy_interval,
            )?));
        }

        if config.format.includes_json() {
            writers.push(Box::new(JsonWriter::new(&config.directory, &config.prefix)?));
        }

        if config.format.includes_toml() {
            writers.push(Box::new(TomlWriter::new(&config.directory, &config.prefix)?));
        }

        Ok(Self { writers })
    }

    /// Create a file writer with specific formats.
    pub fn with_formats(
        directory: &str,
        prefix: &str,
        format: OutputFormat,
        energy_interval: usize,
    ) -> Result<Self, OutputError> {
        let config = OutputConfig {
            format,
            directory: directory.to_string(),
            prefix: prefix.to_string(),
            energy_interval,
        };
        Self::new(&config)
    }
}

impl OutputWriter for FileWriter {
    fn write_step(&mut self, step: &SimulationStep) -> Result<(), OutputError> {
        for writer in &mut self.writers {
            writer.write_step(step)?;
        }
        Ok(())
    }

    fn finalize(&mut self, result: &SimulationResult) -> Result<(), OutputError> {
        for writer in &mut self.writers {
            writer.finalize(result)?;
        }
        Ok(())
    }

    fn flush(&mut self) -> Result<(), OutputError> {
        for writer in &mut self.writers {
            writer.flush()?;
        }
        Ok(())
    }
}
