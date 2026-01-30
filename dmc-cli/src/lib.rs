//! # dmc-cli: Command-Line Interface for DMC
//!
//! This crate provides example simulations for Diffusion Monte Carlo.
//!
//! ## Available Binaries
//!
//! - `hydrogen` - Hydrogen atom simulation
//! - `h2_ion` - H2+ ion simulation
//! - `h2_molecule` - H2 molecule simulation
//!
//! ## Usage
//!
//! ```bash
//! # Run hydrogen simulation
//! cargo run --bin hydrogen
//!
//! # Run H2+ ion simulation
//! cargo run --bin h2_ion
//!
//! # Run H2 molecule simulation
//! cargo run --bin h2_molecule
//!
//! # Run H2 molecule with custom config
//! cargo run --bin h2_molecule -- config.toml
//! ```

// Re-export commonly used items
pub use dmc_core::prelude::*;
pub use dmc_output::prelude::*;
