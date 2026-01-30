//! H₂⁺ ion simulation demonstrating energy fluctuations.
//!
//! Unlike hydrogen with an exact trial wavefunction (α = Z = 1),
//! the H₂⁺ LCAO trial function is not exact, so you'll see
//! actual statistical fluctuations in the energy.
//!
//! # Expected Result
//!
//! Ground state energy at R=2.0 Bohr: ≈ -0.6026 Hartree

use dmc_core::prelude::*;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Configure H2+ simulation
    let config = SimulationConfig::builder()
        .num_walkers(5_000)
        .time_step(0.005) // Smaller timestep for molecule
        .total_steps(20_000)
        .equilibration_steps(2_000)
        .seed(12345)
        .system(SystemType::H2Ion { bond_length: 2.0 }) // R = 2.0 Bohr
        .algorithm(AlgorithmType::ImportanceSampled { max_offspring: 3 })
        .trial_wavefunction(TrialWfParams::H2Ion { alpha: 1.0 })
        .feedback_alpha(1.0)
        .output(OutputConfig {
            format: OutputFormat::Toml,
            directory: "./output".to_string(),
            prefix: "h2_ion".to_string(),
            energy_interval: 10,
            position_interval: 0,
        })
        .build()?;

    // Print header
    print_simulation_header("H₂⁺ Molecule Ion DMC", &config);

    // Run simulation
    let result = run_simulation(config.clone())?;

    // Print standardized results
    print_result_summary(&result);

    // Show energy trace statistics (demonstrates fluctuations)
    if !result.energy_trace.is_empty() {
        let min_e = result
            .energy_trace
            .iter()
            .cloned()
            .fold(f64::INFINITY, f64::min);
        let max_e = result
            .energy_trace
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max);

        println!();
        println!("  Energy trace (shows fluctuations unlike H with exact Ψ_T):");
        println!("    Min: {:.6} Ha", min_e);
        println!("    Max: {:.6} Ha", max_e);
        println!("    Range: {:.6} Ha", max_e - min_e);
    }

    // Print output file info
    print_output_info(
        &config.output.directory,
        &config.output.prefix,
        &config.output.format,
    );

    Ok(())
}
