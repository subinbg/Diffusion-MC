//! Console output with progress bar for DMC simulations.

use super::{OutputError, OutputWriter};
use crate::types::{OutputConfig, SimulationResult};
use console::style;
use simulation::config::SimulationConfig;
use simulation::state::SimulationStep;
use indicatif::{ProgressBar, ProgressStyle};

/// Console writer with progress bar.
pub struct ConsoleWriter {
    progress_bar: Option<ProgressBar>,
    total_steps: usize,
    update_interval: usize,
    running_energy_sum: f64,
    running_count: usize,
    last_acceptance: f64,
}

impl ConsoleWriter {
    /// Create a new console writer with progress bar.
    ///
    /// # Arguments
    ///
    /// * `system_name` - Name of the quantum system
    /// * `total_steps` - Total number of simulation steps
    /// * `initial_energy` - Initial trial energy
    /// * `initial_population` - Initial population size
    /// * `update_interval` - Progress bar update interval
    pub fn new(
        system_name: &str,
        total_steps: usize,
        initial_energy: f64,
        initial_population: usize,
        update_interval: usize,
    ) -> Self {
        let pb = ProgressBar::new(total_steps as u64);
        pb.set_style(
            ProgressStyle::default_bar()
                .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({percent}%) {msg}")
                .unwrap()
                .progress_chars("==-"),
        );
        pb.set_message(format!(
            "{}  E: {:.6} Ha  N: {}",
            style(system_name).bold(),
            initial_energy,
            initial_population
        ));

        Self {
            progress_bar: Some(pb),
            total_steps,
            update_interval: update_interval.max(1),
            running_energy_sum: 0.0,
            running_count: 0,
            last_acceptance: 0.0,
        }
    }

    /// Create a silent console writer (no progress bar).
    pub fn silent() -> Self {
        Self {
            progress_bar: None,
            total_steps: 0,
            update_interval: 1,
            running_energy_sum: 0.0,
            running_count: 0,
            last_acceptance: 0.0,
        }
    }
}

impl OutputWriter for ConsoleWriter {
    fn write_step(&mut self, step: &SimulationStep) -> Result<(), OutputError> {
        // Update running statistics
        if step.is_equilibrated {
            self.running_energy_sum += step.energy_estimate;
            self.running_count += 1;
        }

        if let Some(ar) = step.acceptance_ratio {
            self.last_acceptance = ar;
        }

        // Update progress bar
        if let Some(ref pb) = self.progress_bar {
            if step.step % self.update_interval == 0 || step.step == self.total_steps - 1 {
                let phase = if step.is_equilibrated {
                    style("Sampling").green()
                } else {
                    style("Equilibrating").yellow()
                };

                let avg_energy = if self.running_count > 0 {
                    self.running_energy_sum / self.running_count as f64
                } else {
                    step.energy_estimate
                };

                let acceptance_str = if self.last_acceptance > 0.0 {
                    format!("  Accept: {:.1}%", self.last_acceptance * 100.0)
                } else {
                    String::new()
                };

                pb.set_position(step.step as u64);
                pb.set_message(format!(
                    "{}  E: {} Ha  N: {}{}",
                    phase,
                    style(format!("{:.6}", avg_energy)).cyan(),
                    style(step.population_size).magenta(),
                    acceptance_str,
                ));
            }
        }

        Ok(())
    }

    fn finalize(&mut self, result: &SimulationResult) -> Result<(), OutputError> {
        if let Some(pb) = self.progress_bar.take() {
            pb.finish_with_message(format!(
                "{}  E = {} +/- {} Ha  ({:.1}s)",
                style("Complete").green().bold(),
                style(format!("{:.6}", result.energy)).cyan().bold(),
                style(format!("{:.6}", result.energy_error)).cyan(),
                result.wall_time_seconds,
            ));
        }
        Ok(())
    }

    fn flush(&mut self) -> Result<(), OutputError> {
        Ok(())
    }
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
        style("------------------------------------------------------------").dim()
    );
    println!("  {}", style("SIMULATION RESULTS").bold().cyan());
    println!(
        "{}",
        style("------------------------------------------------------------").dim()
    );
    println!();

    // System info
    println!(
        "  {} {}",
        style("System:").dim(),
        style(&result.system).white().bold()
    );
    println!("  {} {}", style("Algorithm:").dim(), result.algorithm);
    println!();

    // Energy result - the main output
    println!(
        "  {} {} {} {} {}",
        style("Energy:").bold(),
        style(format!("{:>12.6}", result.energy)).cyan().bold(),
        style("+/-").dim(),
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
        style("------------------------------------------------------------").dim()
    );
}

/// Print a formatted header for the simulation.
pub fn print_simulation_header(title: &str, config: &SimulationConfig) {
    println!();
    println!(
        "{}",
        style("------------------------------------------------------------").dim()
    );
    println!("  {}", style(title).bold().cyan());
    println!(
        "{}",
        style("------------------------------------------------------------").dim()
    );
    println!();
    println!("  {} {}", style("Walkers:").dim(), config.num_walkers);
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
pub fn print_output_info(config: &OutputConfig) {
    println!();
    println!("  {}", style("Output files:").dim());

    let format = &config.format;

    if format.includes_csv() {
        println!(
            "    {} {}/{}_energy.csv",
            style("*").dim(),
            config.directory,
            config.prefix
        );
        println!(
            "    {} {}/{}_summary.csv",
            style("*").dim(),
            config.directory,
            config.prefix
        );
    }

    if format.includes_json() {
        println!(
            "    {} {}/{}_result.json",
            style("*").dim(),
            config.directory,
            config.prefix
        );
    }

    if format.includes_toml() {
        println!(
            "    {} {}/{}_result.toml",
            style("*").dim(),
            config.directory,
            config.prefix
        );
    }

    println!();
}
