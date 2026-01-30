//! H₂ molecule simulation loaded from TOML configuration.
//!
//! This example demonstrates loading DMC configuration from a TOML file.
//!
//! # Usage
//!
//! ```bash
//! cargo run --example h2_from_toml -- config.toml
//! ```
//!
//! # Expected Result
//!
//! Ground state energy: ≈ -1.1745 Hartree (experimental)
//!
//! Note: The Heitler-London + Jastrow trial wavefunction is approximate,
//! so expect ~1-2% deviation from the exact value due to:
//! - Fixed-node approximation
//! - Time step bias
//! - Non-optimal variational parameters

use dmc_core::prelude::*;
use std::env;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = env::args().collect();

    let config = if args.len() > 1 {
        // Load from TOML file
        SimulationConfig::from_toml(&args[1])?
    } else {
        // Use optimized H₂ configuration for better accuracy
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
                alpha: 1.2,    // Slightly optimized
                jastrow_b: 0.5,
            })
            .feedback_alpha(1.0)
            .output(OutputConfig {
                format: OutputFormat::All, // Generate all output formats
                directory: "./output".to_string(),
                prefix: "h2_molecule".to_string(),
                energy_interval: 10,
                position_interval: 0,
            })
            .build()?
    };

    // Print header
    print_simulation_header("H₂ Molecule DMC Simulation", &config);

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
