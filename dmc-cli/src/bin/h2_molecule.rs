//! H2 molecule DMC simulation.
//!
//! This example demonstrates importance-sampled DMC for the hydrogen
//! molecule (H2), loaded optionally from TOML configuration.
//!
//! # Usage
//!
//! ```bash
//! cargo run --bin h2_molecule
//! cargo run --bin h2_molecule -- config.toml
//! ```
//!
//! # Expected Result
//!
//! Ground state energy: approximately -1.1745 Hartree (experimental)
//!
//! Note: The Heitler-London + Jastrow trial wavefunction is approximate,
//! so expect ~1-2% deviation from the exact value due to:
//! - Fixed-node approximation
//! - Time step bias
//! - Non-optimal variational parameters

use dmc_core::prelude::*;
use dmc_output::prelude::*;
use std::env;
use std::time::Instant;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = env::args().collect();

    let config = if args.len() > 1 {
        // Load from TOML file
        SimulationConfig::from_toml(&args[1])?
    } else {
        // Use optimized H2 configuration for better accuracy
        // Smaller time step and more walkers reduce systematic errors
        SimulationConfig::builder()
            .num_walkers(10_000)
            .time_step(0.002) // Smaller time step for 2-electron system
            .total_steps(50_000)
            .equilibration_steps(5_000)
            .seed(123)
            .system(SystemType::H2Molecule { bond_length: 1.4 })
            .algorithm(AlgorithmType::ImportanceSampled { max_offspring: 3 })
            .trial_wavefunction(TrialWfParams::H2Molecule {
                alpha: 1.2,     // Slightly optimized
                jastrow_b: 0.5, // Correlation parameter
            })
            .feedback_alpha(1.0)
            .build()?
    };

    // Configure output
    let output_config = OutputConfig {
        format: OutputFormat::All, // Generate all output formats
        directory: "./output".to_string(),
        prefix: "h2_molecule".to_string(),
        energy_interval: 10,
    };

    // Print header
    print_simulation_header("H2 Molecule DMC Simulation", &config);

    // Create simulation and aggregator
    let mut sim = Simulation::new(config)?;
    let mut aggregator = OutputAggregator::new(
        &output_config,
        sim.system_name(),
        sim.total_steps(),
        sim.exact_energy().unwrap_or(-1.1745),
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
