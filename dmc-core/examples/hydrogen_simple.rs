//! Simple hydrogen atom DMC simulation.
//!
//! This example demonstrates importance-sampled DMC for the hydrogen atom.
//!
//! # Expected Result
//!
//! Ground state energy: -0.5 Hartree (exact)
//!
//! Note: With the exact trial wavefunction (α = Z = 1), the local energy
//! E_L is constant everywhere, so there is zero variance in the energy
//! estimate. This is known as the "zero-variance property".

use dmc_core::prelude::*;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Configure simulation
    let config = SimulationConfig::builder()
        .num_walkers(5_000)
        .time_step(0.01)
        .total_steps(10_000)
        .equilibration_steps(1_000)
        .seed(42)
        .system(SystemType::Hydrogen)
        .algorithm(AlgorithmType::ImportanceSampled { max_offspring: 3 })
        .trial_wavefunction(TrialWfParams::Hydrogen { alpha: 1.0 })
        .feedback_alpha(1.0)
        .output(OutputConfig {
            format: OutputFormat::Toml,
            directory: "./output".to_string(),
            prefix: "hydrogen".to_string(),
            energy_interval: 10,
            position_interval: 0,
        })
        .build()?;

    // Print header
    print_simulation_header("Hydrogen Atom DMC", &config);

    // Run simulation
    let result = run_simulation(config.clone())?;

    // Print standardized results
    print_result_summary(&result);

    // Print output file info
    print_output_info(
        &config.output.directory,
        &config.output.prefix,
        &config.output.format,
    );

    Ok(())
}
