//! Data structures for DMC simulation output.
//!
//! These types represent the results and intermediate data from DMC simulations.

use serde::{Deserialize, Serialize};

/// Block-averaged statistics.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct BlockData {
    /// Block index.
    pub block_index: u64,
    /// First step in block.
    pub start_step: u64,
    /// Last step in block.
    pub end_step: u64,
    /// Mean energy over block.
    pub mean_energy: f64,
    /// Standard error of the mean.
    pub std_error: f64,
    /// Mean population over block.
    pub mean_population: f64,
}

/// Convergence diagnostics.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct ConvergenceData {
    /// Autocorrelation time (in steps).
    pub autocorrelation_time: f64,
    /// Effective number of independent samples.
    pub effective_samples: f64,
    /// Is the simulation converged?
    pub is_converged: bool,
}

/// Final simulation results.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SimulationResult {
    /// Algorithm name.
    pub algorithm: String,
    /// System name.
    pub system: String,

    /// Total simulation steps.
    pub total_steps: u64,
    /// Equilibration steps.
    pub equilibration_steps: u64,

    /// Final energy estimate (Hartree).
    pub energy: f64,
    /// Statistical error (standard error of the mean).
    pub energy_error: f64,
    /// Exact energy for comparison (if known).
    pub exact_energy: Option<f64>,

    /// Time series of energy estimates.
    pub energy_trace: Vec<f64>,
    /// Block-averaged data.
    pub blocks: Vec<BlockData>,

    /// Convergence diagnostics.
    pub convergence: ConvergenceData,

    /// Wall-clock time in seconds.
    pub wall_time_seconds: f64,
}

impl SimulationResult {
    /// Calculate relative error compared to exact energy.
    pub fn relative_error(&self) -> Option<f64> {
        self.exact_energy
            .map(|exact| ((self.energy - exact) / exact).abs())
    }
}

/// Output configuration.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct OutputConfig {
    /// Output format.
    #[serde(default)]
    pub format: OutputFormat,

    /// Output directory.
    #[serde(default = "default_directory")]
    pub directory: String,

    /// File prefix.
    #[serde(default)]
    pub prefix: String,

    /// Interval for writing energy data (0 = every step).
    #[serde(default = "default_energy_interval")]
    pub energy_interval: usize,
}

fn default_directory() -> String {
    "./output".to_string()
}

fn default_energy_interval() -> usize {
    1
}

/// Output format selection.
#[derive(Clone, Debug, Default, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "lowercase")]
pub enum OutputFormat {
    /// CSV format for spreadsheets and data analysis tools
    #[default]
    Csv,
    /// JSON format for programmatic access
    Json,
    /// TOML format for human-readable, structured output
    Toml,
    /// Both CSV and JSON
    Both,
    /// All formats (CSV, JSON, TOML)
    All,
}

impl OutputFormat {
    /// Check if CSV output is enabled.
    pub fn includes_csv(&self) -> bool {
        matches!(
            self,
            OutputFormat::Csv | OutputFormat::Both | OutputFormat::All
        )
    }

    /// Check if JSON output is enabled.
    pub fn includes_json(&self) -> bool {
        matches!(
            self,
            OutputFormat::Json | OutputFormat::Both | OutputFormat::All
        )
    }

    /// Check if TOML output is enabled.
    pub fn includes_toml(&self) -> bool {
        matches!(self, OutputFormat::Toml | OutputFormat::All)
    }
}
