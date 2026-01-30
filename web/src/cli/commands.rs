//! CLI commands for building and serving the site.

use clap::{Parser, Subcommand};
use std::fs;
use std::path::Path;
use std::process::Command;

use crate::site::{assets, generate_documentation_page, generate_simulator_page};
use super::server::serve_site;

#[derive(Parser)]
#[command(name = "dmc-site")]
#[command(about = "Build and serve the Diffusion MC website")]
struct Cli {
    #[command(subcommand)]
    command: Option<Commands>,
}

#[derive(Subcommand)]
enum Commands {
    /// Build WASM module using wasm-pack
    Wasm,
    /// Build the complete site (documentation + simulator)
    Build,
    /// Build and serve the site with a development server
    Serve {
        /// Port to serve on
        #[arg(short, long, default_value = "8080")]
        port: u16,
    },
}

/// Run the CLI.
pub fn run() -> Result<(), Box<dyn std::error::Error>> {
    let cli = Cli::parse();

    match cli.command {
        Some(Commands::Wasm) => {
            build_wasm()?;
        }
        Some(Commands::Build) => {
            build_site()?;
        }
        Some(Commands::Serve { port }) => {
            build_site()?;
            serve_site(port)?;
        }
        None => {
            // Default: just build the site
            build_site()?;
        }
    }

    Ok(())
}

/// Build WASM module using wasm-pack.
fn build_wasm() -> Result<(), Box<dyn std::error::Error>> {
    println!("Building WASM module...");

    // Check if wasm-pack is installed
    let wasm_pack_check = Command::new("wasm-pack").arg("--version").output();

    if wasm_pack_check.is_err() {
        println!("wasm-pack not found. Installing...");
        let status = Command::new("cargo")
            .args(["install", "wasm-pack"])
            .status()?;

        if !status.success() {
            return Err("Failed to install wasm-pack".into());
        }
    }

    // Build the WASM module - output to dist/simulator/pkg/
    let status = Command::new("wasm-pack")
        .args(["build", "--target", "web", "--out-dir", "../dist/simulator/pkg", "web"])
        .status()?;

    if !status.success() {
        return Err("wasm-pack build failed".into());
    }

    println!("WASM build complete! Output in dist/simulator/pkg/");
    Ok(())
}

/// Build the complete site.
fn build_site() -> Result<(), Box<dyn std::error::Error>> {
    // Read the README.md file
    let readme_path = Path::new("README.md");
    let readme_content = fs::read_to_string(readme_path)?;

    // Generate the documentation page
    let doc_html = generate_documentation_page(&readme_content)?;

    // Create dist directory and write output
    let dist_path = Path::new("dist");
    fs::create_dir_all(dist_path)?;
    fs::write(dist_path.join("index.html"), doc_html)?;

    // Build simulator page
    let simulator_dest = dist_path.join("simulator");
    fs::create_dir_all(&simulator_dest)?;

    // Check if WASM has been built
    let pkg_path = simulator_dest.join("pkg");
    if !pkg_path.exists() || !pkg_path.join("web.js").exists() {
        println!("WASM not built yet. Building now...");
        build_wasm()?;
    }

    // Generate simulator HTML from maud template
    let simulator_html = generate_simulator_page();
    fs::write(simulator_dest.join("index.html"), simulator_html)?;

    // Write the embedded JavaScript file
    fs::write(simulator_dest.join("main.js"), assets::SIMULATOR_JS)?;

    println!("Site built successfully! Output in dist/");
    Ok(())
}
