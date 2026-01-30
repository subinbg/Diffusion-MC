//! CSV output writer for DMC simulations.

use super::{OutputError, OutputWriter};
use crate::types::SimulationResult;
use dmc_core::state::SimulationStep;
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::Path;

/// CSV output writer.
///
/// Writes step-by-step energy data to CSV files.
pub struct CsvWriter {
    energy_file: Option<BufWriter<File>>,
    directory: String,
    prefix: String,
    interval: usize,
    step_count: usize,
}

impl CsvWriter {
    /// Create a new CSV writer.
    ///
    /// # Arguments
    ///
    /// * `directory` - Output directory
    /// * `prefix` - File name prefix
    /// * `interval` - Write interval (1 = every step)
    pub fn new(directory: &str, prefix: &str, interval: usize) -> Result<Self, OutputError> {
        fs::create_dir_all(directory)?;

        let path = Path::new(directory).join(format!("{}_energy.csv", prefix));
        let file = File::create(path)?;
        let mut writer = BufWriter::new(file);
        writeln!(
            writer,
            "step,time,energy,trial_energy,population,acceptance_ratio"
        )?;

        Ok(Self {
            energy_file: Some(writer),
            directory: directory.to_string(),
            prefix: prefix.to_string(),
            interval: interval.max(1),
            step_count: 0,
        })
    }

    /// Write block data to a separate CSV file.
    fn write_blocks(&self, result: &SimulationResult) -> Result<(), OutputError> {
        if result.blocks.is_empty() {
            return Ok(());
        }

        let path = Path::new(&self.directory).join(format!("{}_blocks.csv", self.prefix));
        let mut file = File::create(path)?;
        writeln!(
            file,
            "block,start_step,end_step,mean_energy,std_error,mean_population"
        )?;

        for block in &result.blocks {
            writeln!(
                file,
                "{},{},{},{:.8},{:.8},{:.2}",
                block.block_index,
                block.start_step,
                block.end_step,
                block.mean_energy,
                block.std_error,
                block.mean_population
            )?;
        }

        Ok(())
    }

    /// Write summary to a separate CSV file.
    fn write_summary(&self, result: &SimulationResult) -> Result<(), OutputError> {
        let path = Path::new(&self.directory).join(format!("{}_summary.csv", self.prefix));
        let mut file = File::create(path)?;
        writeln!(file, "parameter,value")?;
        writeln!(file, "algorithm,{}", result.algorithm)?;
        writeln!(file, "system,{}", result.system)?;
        writeln!(file, "total_steps,{}", result.total_steps)?;
        writeln!(file, "equilibration_steps,{}", result.equilibration_steps)?;
        writeln!(file, "energy,{:.8}", result.energy)?;
        writeln!(file, "energy_error,{:.8}", result.energy_error)?;
        if let Some(exact) = result.exact_energy {
            writeln!(file, "exact_energy,{:.8}", exact)?;
        }
        writeln!(file, "wall_time_seconds,{:.2}", result.wall_time_seconds)?;

        Ok(())
    }
}

impl OutputWriter for CsvWriter {
    fn write_step(&mut self, step: &SimulationStep) -> Result<(), OutputError> {
        if !step.is_equilibrated {
            return Ok(());
        }

        self.step_count += 1;
        if self.step_count % self.interval != 0 {
            return Ok(());
        }

        if let Some(ref mut writer) = self.energy_file {
            let accept_str = step
                .acceptance_ratio
                .map(|a| format!("{:.4}", a))
                .unwrap_or_default();
            writeln!(
                writer,
                "{},{:.6},{:.8},{:.8},{},{}",
                step.step,
                step.imaginary_time,
                step.energy_estimate,
                step.trial_energy,
                step.population_size,
                accept_str
            )?;
        }

        Ok(())
    }

    fn finalize(&mut self, result: &SimulationResult) -> Result<(), OutputError> {
        self.flush()?;
        self.write_blocks(result)?;
        self.write_summary(result)?;
        Ok(())
    }

    fn flush(&mut self) -> Result<(), OutputError> {
        if let Some(ref mut writer) = self.energy_file {
            writer.flush()?;
        }
        Ok(())
    }
}
