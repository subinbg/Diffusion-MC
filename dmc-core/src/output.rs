//! Output data structures and writers for DMC simulations.
//!
//! Supports multiple output formats:
//! - **CSV**: For spreadsheets and data analysis tools
//! - **JSON**: For programmatic access
//! - **TOML**: Human-readable, structured summary
//!
//! # Output Files
//!
//! - `{prefix}_energy.csv`: Step-by-step energy trace
//! - `{prefix}_blocks.csv`: Block-averaged statistics
//! - `{prefix}_result.json`: Complete results in JSON
//! - `{prefix}_result.toml`: Human-readable summary in TOML

use crate::config::OutputConfig;
use console::style;
use serde::{Deserialize, Serialize};
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::Path;

/// Single step output data.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct StepData {
    /// Step number.
    pub step: u64,
    /// Imaginary time τ = step × δτ.
    pub imaginary_time: f64,
    /// Energy estimate from this step.
    pub energy_estimate: f64,
    /// Current trial energy E_T.
    pub trial_energy: f64,
    /// Current population size.
    pub population_size: usize,
    /// Acceptance ratio (importance sampling only).
    pub acceptance_ratio: Option<f64>,
}

/// Block-averaged statistics.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct BlockData {
    /// Block index.
    pub block_index: u64,
    /// First step in block.
    pub start_step: u64,
    /// Last step in block.
    pub end_step: u64,
    /// Mean energy over block.
    pub mean_energy: f64,
    /// Standard error of the mean.
    pub std_error: f64,
    /// Mean population over block.
    pub mean_population: f64,
}

/// Convergence diagnostics.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct ConvergenceData {
    /// Autocorrelation time (in steps).
    pub autocorrelation_time: f64,
    /// Effective number of independent samples.
    pub effective_samples: f64,
    /// Is the simulation converged?
    pub is_converged: bool,
}

/// Final simulation results.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SimulationResult {
    /// Algorithm name.
    pub algorithm: String,
    /// System name.
    pub system: String,

    /// Total simulation steps.
    pub total_steps: u64,
    /// Equilibration steps.
    pub equilibration_steps: u64,

    /// Final energy estimate (Hartree).
    pub energy: f64,
    /// Statistical error (standard error of the mean).
    pub energy_error: f64,
    /// Exact energy for comparison (if known).
    pub exact_energy: Option<f64>,

    /// Time series of energy estimates.
    pub energy_trace: Vec<f64>,
    /// Block-averaged data.
    pub blocks: Vec<BlockData>,

    /// Convergence diagnostics.
    pub convergence: ConvergenceData,

    /// Wall-clock time in seconds.
    pub wall_time_seconds: f64,
}

impl SimulationResult {
    /// Calculate relative error compared to exact energy.
    pub fn relative_error(&self) -> Option<f64> {
        self.exact_energy.map(|exact| {
            ((self.energy - exact) / exact).abs()
        })
    }
}

/// Output writer for DMC simulations.
pub struct OutputWriter {
    config: OutputConfig,
    energy_file: Option<BufWriter<File>>,
    step_count: u64,
}

impl OutputWriter {
    /// Create a new output writer.
    pub fn new(config: &OutputConfig) -> std::io::Result<Self> {
        // Create output directory
        fs::create_dir_all(&config.directory)?;

        let energy_file = match config.format {
            crate::config::OutputFormat::Csv
            | crate::config::OutputFormat::Both
            | crate::config::OutputFormat::All => {
                let path = Path::new(&config.directory)
                    .join(format!("{}_energy.csv", config.prefix));
                let file = File::create(path)?;
                let mut writer = BufWriter::new(file);
                writeln!(
                    writer,
                    "step,time,energy,trial_energy,population,acceptance_ratio"
                )?;
                Some(writer)
            }
            crate::config::OutputFormat::Json | crate::config::OutputFormat::Toml => None,
        };

        Ok(Self {
            config: config.clone(),
            energy_file,
            step_count: 0,
        })
    }

    /// Write a single step's data.
    pub fn write_step(&mut self, data: &StepData) -> std::io::Result<()> {
        self.step_count += 1;

        if self.config.energy_interval > 0
            && self.step_count % self.config.energy_interval as u64 != 0
        {
            return Ok(());
        }

        if let Some(ref mut writer) = self.energy_file {
            let accept_str = data
                .acceptance_ratio
                .map(|a| format!("{:.4}", a))
                .unwrap_or_default();
            writeln!(
                writer,
                "{},{:.6},{:.8},{:.8},{},{}",
                data.step,
                data.imaginary_time,
                data.energy_estimate,
                data.trial_energy,
                data.population_size,
                accept_str
            )?;
        }

        Ok(())
    }

    /// Flush all buffers.
    pub fn flush(&mut self) -> std::io::Result<()> {
        if let Some(ref mut writer) = self.energy_file {
            writer.flush()?;
        }
        Ok(())
    }

    /// Write final simulation results.
    pub fn write_result(&self, result: &SimulationResult) -> std::io::Result<()> {
        use crate::config::OutputFormat;

        // Write JSON if requested
        if matches!(
            self.config.format,
            OutputFormat::Json | OutputFormat::Both | OutputFormat::All
        ) {
            let path = Path::new(&self.config.directory)
                .join(format!("{}_result.json", self.config.prefix));
            let file = File::create(path)?;
            serde_json::to_writer_pretty(file, result)?;
        }

        // Write TOML summary if requested
        if matches!(
            self.config.format,
            OutputFormat::Toml | OutputFormat::All
        ) {
            self.write_toml_result(result)?;
        }

        // Write summary CSV if only CSV format
        if matches!(self.config.format, OutputFormat::Csv) {
            let path = Path::new(&self.config.directory)
                .join(format!("{}_summary.csv", self.config.prefix));
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
        }

        // Write blocks CSV if using CSV format
        if matches!(
            self.config.format,
            OutputFormat::Csv | OutputFormat::Both | OutputFormat::All
        ) && !result.blocks.is_empty()
        {
            let path = Path::new(&self.config.directory)
                .join(format!("{}_blocks.csv", self.config.prefix));
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
        }

        Ok(())
    }

    /// Write human-readable TOML summary.
    fn write_toml_result(&self, result: &SimulationResult) -> std::io::Result<()> {
        let path = Path::new(&self.config.directory)
            .join(format!("{}_result.toml", self.config.prefix));
        let mut file = File::create(path)?;

        writeln!(file, "# DMC Simulation Result")?;
        writeln!(file, "# Generated by dmc-core")?;
        writeln!(file)?;

        writeln!(file, "[simulation]")?;
        writeln!(file, "algorithm = \"{}\"", result.algorithm)?;
        writeln!(file, "system = \"{}\"", result.system)?;
        writeln!(file, "total_steps = {}", result.total_steps)?;
        writeln!(file, "equilibration_steps = {}", result.equilibration_steps)?;
        writeln!(file, "wall_time_seconds = {:.2}", result.wall_time_seconds)?;
        writeln!(file)?;

        writeln!(file, "[energy]")?;
        writeln!(file, "# Final energy estimate in Hartree")?;
        writeln!(file, "value = {:.8}", result.energy)?;
        writeln!(file, "error = {:.8}", result.energy_error)?;
        if let Some(exact) = result.exact_energy {
            writeln!(file, "reference = {:.8}", exact)?;
            let deviation = result.energy - exact;
            let percent = 100.0 * deviation.abs() / exact.abs();
            writeln!(file, "deviation = {:.8}", deviation)?;
            writeln!(file, "deviation_percent = {:.4}", percent)?;
        }
        writeln!(file)?;

        if !result.blocks.is_empty() {
            writeln!(file, "[statistics]")?;
            writeln!(file, "num_blocks = {}", result.blocks.len())?;
            let block_energies: Vec<f64> = result.blocks.iter().map(|b| b.mean_energy).collect();
            let min_e = block_energies.iter().cloned().fold(f64::INFINITY, f64::min);
            let max_e = block_energies.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
            writeln!(file, "energy_min = {:.8}", min_e)?;
            writeln!(file, "energy_max = {:.8}", max_e)?;
            writeln!(file, "energy_range = {:.8}", max_e - min_e)?;
        }

        Ok(())
    }
}

/// Calculate block-averaged statistics from energy trace.
///
/// # Arguments
///
/// * `energies` - Vector of energy samples
/// * `populations` - Vector of population sizes
/// * `block_size` - Number of samples per block
///
/// # Returns
///
/// Vector of block statistics.
pub fn block_average(
    energies: &[f64],
    populations: &[usize],
    block_size: usize,
) -> Vec<BlockData> {
    let num_blocks = energies.len() / block_size;
    let mut blocks = Vec::with_capacity(num_blocks);

    for i in 0..num_blocks {
        let start = i * block_size;
        let end = start + block_size;

        let block_energies = &energies[start..end];
        let block_pops = &populations[start..end];

        let mean_energy = block_energies.iter().sum::<f64>() / block_size as f64;
        let mean_pop = block_pops.iter().sum::<usize>() as f64 / block_size as f64;

        // Standard error within block
        let variance: f64 = block_energies
            .iter()
            .map(|e| (e - mean_energy).powi(2))
            .sum::<f64>()
            / (block_size - 1) as f64;
        let std_error = (variance / block_size as f64).sqrt();

        blocks.push(BlockData {
            block_index: i as u64,
            start_step: start as u64,
            end_step: end as u64,
            mean_energy,
            std_error,
            mean_population: mean_pop,
        });
    }

    blocks
}

/// Calculate final energy estimate and error from block data.
///
/// # Returns
///
/// (mean_energy, standard_error)
pub fn final_statistics(blocks: &[BlockData]) -> (f64, f64) {
    if blocks.is_empty() {
        return (0.0, 0.0);
    }

    let n = blocks.len() as f64;
    let mean: f64 = blocks.iter().map(|b| b.mean_energy).sum::<f64>() / n;

    // Standard error of block means
    let variance: f64 = blocks
        .iter()
        .map(|b| (b.mean_energy - mean).powi(2))
        .sum::<f64>()
        / (n - 1.0);
    let std_error = (variance / n).sqrt();

    (mean, std_error)
}

// ============================================================================
// CLI Display Utilities
// ============================================================================

/// Print a professional, formatted result summary to the terminal.
///
/// Use this for consistent output across all examples.
pub fn print_result_summary(result: &SimulationResult) {
    println!();
    println!(
        "{}",
        style("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")
            .dim()
    );
    println!(
        "  {}",
        style("SIMULATION RESULTS").bold().cyan()
    );
    println!(
        "{}",
        style("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")
            .dim()
    );
    println!();

    // System info
    println!(
        "  {} {}",
        style("System:").dim(),
        style(&result.system).white().bold()
    );
    println!(
        "  {} {}",
        style("Algorithm:").dim(),
        result.algorithm
    );
    println!();

    // Energy result - the main output
    println!(
        "  {} {} {} {} {}",
        style("Energy:").bold(),
        style(format!("{:>12.6}", result.energy)).cyan().bold(),
        style("±").dim(),
        style(format!("{:.6}", result.energy_error)).cyan(),
        style("Ha").dim()
    );

    // Comparison with exact value if available
    if let Some(exact) = result.exact_energy {
        let deviation = result.energy - exact;
        let percent = 100.0 * deviation.abs() / exact.abs();

        println!(
            "  {} {:>12.6} {}",
            style("Reference:").dim(),
            exact,
            style("Ha").dim()
        );

        let deviation_style = if percent < 0.5 {
            style(format!("{:+.6}", deviation)).green()
        } else if percent < 2.0 {
            style(format!("{:+.6}", deviation)).yellow()
        } else {
            style(format!("{:+.6}", deviation)).red()
        };

        println!(
            "  {} {} {} ({:.2}%)",
            style("Deviation:").dim(),
            deviation_style,
            style("Ha").dim(),
            percent
        );
    }

    println!();

    // Statistics
    println!(
        "  {} {} steps ({} equilibration)",
        style("Steps:").dim(),
        result.total_steps,
        result.equilibration_steps
    );

    if !result.blocks.is_empty() {
        println!(
            "  {} {} blocks analyzed",
            style("Blocks:").dim(),
            result.blocks.len()
        );
    }

    println!(
        "  {} {:.2} seconds",
        style("Time:").dim(),
        result.wall_time_seconds
    );

    println!();
    println!(
        "{}",
        style("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")
            .dim()
    );
}

/// Print a formatted header for the simulation.
pub fn print_simulation_header(title: &str, config: &crate::config::SimulationConfig) {
    println!();
    println!(
        "{}",
        style("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")
            .dim()
    );
    println!("  {}", style(title).bold().cyan());
    println!(
        "{}",
        style("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")
            .dim()
    );
    println!();
    println!(
        "  {} {}",
        style("Walkers:").dim(),
        config.num_walkers
    );
    println!(
        "  {} {} a.u.",
        style("Time step:").dim(),
        config.time_step
    );
    println!(
        "  {} {} ({} equilibration)",
        style("Steps:").dim(),
        config.total_steps,
        config.equilibration_steps
    );
    println!();
}

/// Print output file locations.
pub fn print_output_info(directory: &str, prefix: &str, format: &crate::config::OutputFormat) {
    use crate::config::OutputFormat;

    println!();
    println!("  {}", style("Output files:").dim());

    match format {
        OutputFormat::Csv => {
            println!("    {} {}/{}_energy.csv", style("•").dim(), directory, prefix);
            println!("    {} {}/{}_summary.csv", style("•").dim(), directory, prefix);
        }
        OutputFormat::Json => {
            println!("    {} {}/{}_result.json", style("•").dim(), directory, prefix);
        }
        OutputFormat::Toml => {
            println!("    {} {}/{}_result.toml", style("•").dim(), directory, prefix);
        }
        OutputFormat::Both => {
            println!("    {} {}/{}_energy.csv", style("•").dim(), directory, prefix);
            println!("    {} {}/{}_result.json", style("•").dim(), directory, prefix);
        }
        OutputFormat::All => {
            println!("    {} {}/{}_energy.csv", style("•").dim(), directory, prefix);
            println!("    {} {}/{}_result.json", style("•").dim(), directory, prefix);
            println!("    {} {}/{}_result.toml", style("•").dim(), directory, prefix);
        }
    }
    println!();
}
