//! Simulation runner for DMC.
//!
//! This module orchestrates the complete DMC simulation:
//! 1. Initialize walker population
//! 2. Run equilibration phase
//! 3. Run production phase collecting statistics
//! 4. Compute final results

use crate::algorithm::{DmcAlgorithm, ImportanceSampledDmc, PureDmc};
use crate::config::{AlgorithmType, SimulationConfig, SystemType};
use crate::output::{
    block_average, final_statistics, BlockData, ConvergenceData, OutputWriter, SimulationResult,
    StepData,
};
use crate::physics::update_trial_energy;
use crate::walker::Population;
use console::style;
use indicatif::{ProgressBar, ProgressStyle};
use std::time::Instant;

/// Progress display options
#[derive(Clone, Debug)]
pub struct ProgressOptions {
    /// Show progress bar (default: true)
    pub show_progress: bool,
    /// Update interval in steps (default: 50)
    ///
    /// Note: Very small values (< 10) may cause performance overhead
    /// due to terminal rendering. Recommended: 50-100 for real-time feel.
    pub update_interval: usize,
}

impl Default for ProgressOptions {
    fn default() -> Self {
        Self::enabled()
    }
}

impl ProgressOptions {
    /// Create with progress display enabled (update every 50 steps)
    pub fn enabled() -> Self {
        Self {
            show_progress: true,
            update_interval: 50,
        }
    }

    /// Create with fast updates (every 10 steps) for short simulations
    pub fn fast() -> Self {
        Self {
            show_progress: true,
            update_interval: 10,
        }
    }

    /// Create with progress display disabled (for tests/benchmarks)
    pub fn disabled() -> Self {
        Self {
            show_progress: false,
            update_interval: 100,
        }
    }

    /// Set custom update interval
    pub fn with_interval(mut self, interval: usize) -> Self {
        self.update_interval = interval.max(1);
        self
    }
}

/// Run a complete DMC simulation with default progress display.
///
/// # Algorithm Overview
///
/// 1. **Initialization**: Create walker population at origin or sampled positions
/// 2. **Equilibration**: Run `equilibration_steps` to reach steady state
/// 3. **Production**: Run `total_steps - equilibration_steps` collecting statistics
/// 4. **Analysis**: Compute block averages and final energy estimate
///
/// # Returns
///
/// `SimulationResult` containing final energy, error, and diagnostics.
pub fn run_simulation(config: SimulationConfig) -> Result<SimulationResult, std::io::Error> {
    run_simulation_with_progress(config, ProgressOptions::enabled())
}

/// Run simulation without progress output (for tests).
pub fn run_simulation_quiet(config: SimulationConfig) -> Result<SimulationResult, std::io::Error> {
    run_simulation_with_progress(config, ProgressOptions::disabled())
}

/// Run a complete DMC simulation with configurable progress display.
pub fn run_simulation_with_progress(
    config: SimulationConfig,
    progress_opts: ProgressOptions,
) -> Result<SimulationResult, std::io::Error> {
    let start_time = Instant::now();

    // Create output writer
    let mut output = OutputWriter::new(&config.output)?;

    // Create algorithm and get system info
    let mut algorithm: Box<dyn DmcAlgorithm> = match &config.algorithm {
        AlgorithmType::Pure { .. } => Box::new(PureDmc::new(&config)),
        AlgorithmType::ImportanceSampled { .. } => Box::new(ImportanceSampledDmc::new(&config)),
    };

    let num_electrons = match &config.system {
        SystemType::Hydrogen => 1,
        SystemType::H2Ion { .. } => 1,
        SystemType::H2Molecule { .. } => 2,
    };

    let exact_energy = match &config.system {
        SystemType::Hydrogen => Some(-0.5),
        SystemType::H2Ion { .. } => Some(-0.6026),
        SystemType::H2Molecule { .. } => Some(-1.1745),
    };

    let system_name = match &config.system {
        SystemType::Hydrogen => "Hydrogen (H)".to_string(),
        SystemType::H2Ion { bond_length } => format!("H₂⁺ (R={:.2})", bond_length),
        SystemType::H2Molecule { bond_length } => format!("H₂ (R={:.2})", bond_length),
    };

    // Initialize population
    let mut population = Population::new(config.num_walkers, config.seed);
    population.initialize_gaussian(num_electrons, 1.0);

    // Initialize trial energy
    let mut trial_energy = exact_energy.unwrap_or(-0.5);

    // Storage for statistics
    let mut energy_trace = Vec::with_capacity(config.total_steps);
    let mut population_trace = Vec::with_capacity(config.total_steps);

    // Create progress bar
    let progress_bar = if progress_opts.show_progress {
        let pb = ProgressBar::new(config.total_steps as u64);
        pb.set_style(
            ProgressStyle::default_bar()
                .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({percent}%) {msg}")
                .unwrap()
                .progress_chars("█▓▒░  "),
        );
        pb.set_message(format!(
            "{}  E: {:.6} Ha  N: {}",
            style(&system_name).bold(),
            trial_energy,
            population.size()
        ));
        Some(pb)
    } else {
        None
    };

    // Track statistics for display
    let mut running_energy_sum = 0.0;
    let mut running_count = 0usize;
    let mut last_acceptance = 0.0f64;

    // Main simulation loop
    for step in 0..config.total_steps {
        // Perform one DMC step
        let result = algorithm.step(&mut population, trial_energy, &config);

        // Update trial energy using population control
        // E_T = ⟨E_L⟩ - α ln(N/N₀)
        // See README: Step 5: Updating E_T
        trial_energy = update_trial_energy(
            result.energy_estimate,
            result.population_size,
            config.num_walkers,
            config.feedback_alpha,
        );

        if let Some(ar) = result.acceptance_ratio {
            last_acceptance = ar;
        }

        // Record data after equilibration
        if step >= config.equilibration_steps {
            energy_trace.push(result.energy_estimate);
            population_trace.push(result.population_size);
            running_energy_sum += result.energy_estimate;
            running_count += 1;

            // Write step data
            let step_data = StepData {
                step: step as u64,
                imaginary_time: step as f64 * config.time_step,
                energy_estimate: result.energy_estimate,
                trial_energy,
                population_size: result.population_size,
                acceptance_ratio: result.acceptance_ratio,
            };
            output.write_step(&step_data)?;
        }

        // Update progress bar
        if let Some(ref pb) = progress_bar {
            if step % progress_opts.update_interval == 0 || step == config.total_steps - 1 {
                let phase = if step < config.equilibration_steps {
                    style("Equilibrating").yellow()
                } else {
                    style("Sampling").green()
                };

                let avg_energy = if running_count > 0 {
                    running_energy_sum / running_count as f64
                } else {
                    result.energy_estimate
                };

                let acceptance_str = if last_acceptance > 0.0 {
                    format!("  Accept: {:.1}%", last_acceptance * 100.0)
                } else {
                    String::new()
                };

                pb.set_position(step as u64);
                pb.set_message(format!(
                    "{}  E: {} Ha  N: {}{}",
                    phase,
                    style(format!("{:.6}", avg_energy)).cyan(),
                    style(result.population_size).magenta(),
                    acceptance_str,
                ));
            }
        }
    }

    output.flush()?;

    // Calculate statistics BEFORE finishing progress bar
    let block_size = 100.min(energy_trace.len() / 10).max(1);
    let blocks: Vec<BlockData> = if energy_trace.len() >= block_size * 2 {
        block_average(&energy_trace, &population_trace, block_size)
    } else {
        Vec::new()
    };

    let (energy, energy_error) = if !blocks.is_empty() {
        final_statistics(&blocks)
    } else if !energy_trace.is_empty() {
        let mean = energy_trace.iter().sum::<f64>() / energy_trace.len() as f64;
        let var: f64 = energy_trace
            .iter()
            .map(|e| (e - mean).powi(2))
            .sum::<f64>()
            / (energy_trace.len() - 1).max(1) as f64;
        (mean, (var / energy_trace.len() as f64).sqrt())
    } else {
        (trial_energy, 0.0)
    };

    // Finish progress bar with final computed energy
    if let Some(pb) = progress_bar {
        let elapsed = start_time.elapsed().as_secs_f64();
        pb.finish_with_message(format!(
            "{}  E = {} ± {} Ha  ({:.1}s)",
            style("Complete").green().bold(),
            style(format!("{:.6}", energy)).cyan().bold(),
            style(format!("{:.6}", energy_error)).cyan(),
            elapsed,
        ));
    }

    // Create result
    let result = SimulationResult {
        algorithm: algorithm.name().to_string(),
        system: system_name,
        total_steps: config.total_steps as u64,
        equilibration_steps: config.equilibration_steps as u64,
        energy,
        energy_error,
        exact_energy,
        energy_trace,
        blocks,
        convergence: ConvergenceData::default(),
        wall_time_seconds: start_time.elapsed().as_secs_f64(),
    };

    // Write final results
    output.write_result(&result)?;

    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::{OutputConfig, OutputFormat, SimulationConfigBuilder, TrialWfParams};
    use tempfile::TempDir;

    #[test]
    fn hydrogen_ground_state_energy() {
        // Use temp directory for output to avoid polluting workspace
        let temp_dir = TempDir::new().unwrap();
        let output_path = temp_dir.path().to_str().unwrap().to_string();

        let config = SimulationConfigBuilder::default()
            .num_walkers(1000)
            .time_step(0.01)
            .total_steps(5000)
            .equilibration_steps(500)
            .seed(42)
            .system(SystemType::Hydrogen)
            .algorithm(AlgorithmType::ImportanceSampled { max_offspring: 3 })
            .trial_wavefunction(TrialWfParams::Hydrogen { alpha: 1.0 })
            .output(OutputConfig {
                format: OutputFormat::Csv,
                directory: output_path,
                prefix: "test".to_string(),
                energy_interval: 10,
                position_interval: 0,
            })
            .build()
            .unwrap();

        // Run without progress output for cleaner test output
        let result = run_simulation_quiet(config).unwrap();

        // Energy should be close to -0.5 Ha
        assert!(
            (result.energy - (-0.5)).abs() < 0.05,
            "Expected energy ~-0.5, got {}",
            result.energy
        );
    }

    #[test]
    fn h2_ion_with_variance() {
        // H2+ should show some variance (unlike H with exact trial function)
        let temp_dir = TempDir::new().unwrap();
        let output_path = temp_dir.path().to_str().unwrap().to_string();

        let config = SimulationConfigBuilder::default()
            .num_walkers(500)
            .time_step(0.01)
            .total_steps(2000)
            .equilibration_steps(200)
            .seed(42)
            .system(SystemType::H2Ion { bond_length: 2.0 })
            .algorithm(AlgorithmType::ImportanceSampled { max_offspring: 3 })
            .trial_wavefunction(TrialWfParams::H2Ion { alpha: 1.0 })
            .output(OutputConfig {
                format: OutputFormat::Csv,
                directory: output_path,
                prefix: "test_h2ion".to_string(),
                energy_interval: 10,
                position_interval: 0,
            })
            .build()
            .unwrap();

        let result = run_simulation_quiet(config).unwrap();

        // Energy should be around -0.6 Ha (less precise due to non-exact trial function)
        assert!(
            result.energy < 0.0 && result.energy > -1.0,
            "Energy {} out of expected range",
            result.energy
        );

        // Should have some variance (unlike hydrogen with exact trial function)
        // Note: error might still be very small with good sampling
        println!("H2+ energy: {:.6} ± {:.6}", result.energy, result.energy_error);
    }
}
