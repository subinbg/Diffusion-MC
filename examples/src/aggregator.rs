//! Output aggregator for collecting simulation steps into results.

use crate::statistics::{block_average, final_statistics};
use crate::types::{ConvergenceData, OutputConfig, SimulationResult};
use crate::writers::{ConsoleWriter, FileWriter, OutputError, OutputWriter};
use simulation::state::SimulationStep;

/// Aggregator for collecting simulation steps and computing statistics.
///
/// Combines console output (progress bar) with file output, and collects
/// energy/population traces for statistical analysis.
pub struct OutputAggregator {
    console: ConsoleWriter,
    file_writer: Option<FileWriter>,
    energy_trace: Vec<f64>,
    population_trace: Vec<usize>,
    equilibration_steps: usize,
    total_steps: usize,
}

impl OutputAggregator {
    /// Create a new output aggregator.
    ///
    /// # Arguments
    ///
    /// * `config` - Output configuration
    /// * `system_name` - Name of the quantum system
    /// * `total_steps` - Total number of simulation steps
    /// * `initial_energy` - Initial trial energy
    /// * `initial_population` - Initial population size
    pub fn new(
        config: &OutputConfig,
        system_name: &str,
        total_steps: usize,
        initial_energy: f64,
        initial_population: usize,
    ) -> Result<Self, OutputError> {
        let console = ConsoleWriter::new(
            system_name,
            total_steps,
            initial_energy,
            initial_population,
            50, // Update interval
        );

        let file_writer = Some(FileWriter::new(config)?);

        Ok(Self {
            console,
            file_writer,
            energy_trace: Vec::with_capacity(total_steps),
            population_trace: Vec::with_capacity(total_steps),
            equilibration_steps: 0,
            total_steps,
        })
    }

    /// Create a silent aggregator (no console output, no file output).
    pub fn silent(total_steps: usize) -> Self {
        Self {
            console: ConsoleWriter::silent(),
            file_writer: None,
            energy_trace: Vec::with_capacity(total_steps),
            population_trace: Vec::with_capacity(total_steps),
            equilibration_steps: 0,
            total_steps,
        }
    }

    /// Create an aggregator with only console output (no files).
    pub fn console_only(
        system_name: &str,
        total_steps: usize,
        initial_energy: f64,
        initial_population: usize,
    ) -> Self {
        Self {
            console: ConsoleWriter::new(
                system_name,
                total_steps,
                initial_energy,
                initial_population,
                50,
            ),
            file_writer: None,
            energy_trace: Vec::with_capacity(total_steps),
            population_trace: Vec::with_capacity(total_steps),
            equilibration_steps: 0,
            total_steps,
        }
    }

    /// Process a single simulation step.
    pub fn process_step(&mut self, step: &SimulationStep) {
        // Collect statistics after equilibration
        if step.is_equilibrated {
            self.energy_trace.push(step.energy_estimate);
            self.population_trace.push(step.population_size);

            if self.equilibration_steps == 0 {
                self.equilibration_steps = step.step;
            }
        }

        // Write to console and file writers
        let _ = self.console.write_step(step);
        if let Some(ref mut writer) = self.file_writer {
            let _ = writer.write_step(step);
        }
    }

    /// Finalize and compute statistics.
    ///
    /// Returns the complete simulation result.
    pub fn finalize(
        mut self,
        algorithm: &str,
        system: &str,
        exact_energy: Option<f64>,
        wall_time_seconds: f64,
    ) -> SimulationResult {
        // Calculate block statistics
        let block_size = 100.min(self.energy_trace.len() / 10).max(1);
        let blocks = if self.energy_trace.len() >= block_size * 2 {
            block_average(&self.energy_trace, &self.population_trace, block_size)
        } else {
            Vec::new()
        };

        // Calculate final energy and error
        let (energy, energy_error) = if !blocks.is_empty() {
            final_statistics(&blocks)
        } else if !self.energy_trace.is_empty() {
            let mean = self.energy_trace.iter().sum::<f64>() / self.energy_trace.len() as f64;
            let var: f64 = self
                .energy_trace
                .iter()
                .map(|e| (e - mean).powi(2))
                .sum::<f64>()
                / (self.energy_trace.len() - 1).max(1) as f64;
            (mean, (var / self.energy_trace.len() as f64).sqrt())
        } else {
            (0.0, 0.0)
        };

        let result = SimulationResult {
            algorithm: algorithm.to_string(),
            system: system.to_string(),
            total_steps: self.total_steps as u64,
            equilibration_steps: self.equilibration_steps as u64,
            energy,
            energy_error,
            exact_energy,
            energy_trace: self.energy_trace,
            blocks,
            convergence: ConvergenceData::default(),
            wall_time_seconds,
        };

        // Finalize writers
        let _ = self.console.finalize(&result);
        if let Some(ref mut writer) = self.file_writer {
            let _ = writer.finalize(&result);
        }

        result
    }

    /// Get the collected energy trace.
    pub fn energy_trace(&self) -> &[f64] {
        &self.energy_trace
    }

    /// Get the collected population trace.
    pub fn population_trace(&self) -> &[usize] {
        &self.population_trace
    }
}
