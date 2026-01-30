//! CLI entry point for web site builder.
//!
//! # Usage
//!
//! ```bash
//! # Build WASM module
//! cargo run -p web --features cli -- wasm
//!
//! # Build the complete site
//! cargo run -p web --features cli -- build
//!
//! # Build and serve with dev server
//! cargo run -p web --features cli -- serve
//!
//! # Serve on custom port
//! cargo run -p web --features cli -- serve --port 3000
//! ```

#[cfg(feature = "cli")]
mod cli;
#[cfg(feature = "cli")]
mod site;

#[cfg(feature = "cli")]
fn main() -> Result<(), Box<dyn std::error::Error>> {
    cli::run()
}

#[cfg(not(feature = "cli"))]
fn main() {
    eprintln!("CLI not available. Build with --features cli");
    eprintln!("Example: cargo run -p web --features cli -- build");
    std::process::exit(1);
}
