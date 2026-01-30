//! CLI entry point for web site builder.
//!
//! # Usage
//!
//! ```bash
//! # Build WASM module
//! cargo run -p web -- wasm
//!
//! # Build the complete site
//! cargo run -p web -- build
//!
//! # Build and serve with dev server
//! cargo run -p web -- serve
//!
//! # Serve on custom port
//! cargo run -p web -- serve --port 3000
//! ```

#[cfg(not(target_arch = "wasm32"))]
mod cli;
#[cfg(not(target_arch = "wasm32"))]
mod site;

#[cfg(not(target_arch = "wasm32"))]
fn main() -> Result<(), Box<dyn std::error::Error>> {
    cli::run()
}

#[cfg(target_arch = "wasm32")]
fn main() {
    // This binary is not meant to run in WASM
    unreachable!("This binary should not be compiled for WASM");
}
