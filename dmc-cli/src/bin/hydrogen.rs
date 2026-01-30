//! Simple hydrogen atom DMC simulation.
//!
//! This example demonstrates importance-sampled DMC for the hydrogen atom.
//!
//! # Expected Result
//!
//! Ground state energy: -0.5 Hartree (exact)
//!
//! Note: With the exact trial wavefunction (alpha = Z = 1), the local energy
//! E_L is constant everywhere, so there is zero variance in the energy
//! estimate. This is known as the "zero-variance property".

use dmc_core::prelude::*;
use dmc_output::prelude::*;
use std::time::Instant;

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
        .build()?;

    // Configure output
    let output_config = OutputConfig {
        format: OutputFormat::Toml,
        directory: "./output".to_string(),
        prefix: "hydrogen".to_string(),
        energy_interval: 10,
    };

    // Print header
    print_simulation_header("Hydrogen Atom DMC", &config);

    // Create simulation and aggregator
    let mut sim = Simulation::new(config)?;
    let mut aggregator = OutputAggregator::new(
        &output_config,
        sim.system_name(),
        sim.total_steps(),
        sim.exact_energy().unwrap_or(-0.5),
        sim.target_population(),
    )?;

    let start = Instant::now();

    // Run simulation
    for step in &mut sim {
        aggregator.process_step(&step);
    }

    // Finalize and print results
    let result = aggregator.finalize(
        sim.algorithm_name(),
        sim.system_name(),
        sim.exact_energy(),
        start.elapsed().as_secs_f64(),
    );

    print_result_summary(&result);
    print_output_info(&output_config);

    Ok(())
}
