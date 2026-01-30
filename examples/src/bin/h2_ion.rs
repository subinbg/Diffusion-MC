//! H2+ ion DMC simulation.
//!
//! This example demonstrates importance-sampled DMC for the hydrogen
//! molecule ion (H2+).
//!
//! # Expected Result
//!
//! Ground state energy: approximately -0.6026 Hartree
//!
//! Note: Unlike hydrogen with the exact trial wavefunction, H2+ uses
//! an LCAO approximation that is not exact. This results in statistical
//! variance in the energy estimate.

use examples::prelude::*;
use std::time::Instant;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Configure simulation
    // Using equilibrium bond length of ~2.0 Bohr
    let config = SimulationConfig::builder()
        .num_walkers(5_000)
        .time_step(0.005) // Smaller time step for molecular system
        .total_steps(20_000)
        .equilibration_steps(2_000)
        .seed(42)
        .system(SystemType::H2Ion { bond_length: 2.0 })
        .algorithm(AlgorithmType::ImportanceSampled { max_offspring: 3 })
        .trial_wavefunction(TrialWfParams::H2Ion { alpha: 1.0 })
        .feedback_alpha(1.0)
        .build()?;

    // Configure output
    let output_config = OutputConfig {
        format: OutputFormat::Toml,
        directory: "./output".to_string(),
        prefix: "h2_ion".to_string(),
        energy_interval: 10,
    };

    // Print header
    print_simulation_header("H2+ Ion DMC Simulation", &config);

    // Create simulation and aggregator
    let mut sim = Simulation::new(config)?;
    let mut aggregator = OutputAggregator::new(
        &output_config,
        sim.system_name(),
        sim.total_steps(),
        sim.exact_energy().unwrap_or(-0.6026),
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
