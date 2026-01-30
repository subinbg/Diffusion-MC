//! # dmc-web: WebAssembly Interface for DMC
//!
//! This crate provides a WebAssembly interface for running DMC simulations
//! in the browser with real-time visualization.

mod density;

use dmc_core::config::{AlgorithmType, SimulationConfig, SystemType, TrialWfParams};
use dmc_core::state::{Simulation, SimulationState, SimulationStep};
use serde::{Deserialize, Serialize};
use wasm_bindgen::prelude::*;

pub use density::{compute_density_2d, Plane};

/// Initialize panic hook for better error messages in browser console.
#[wasm_bindgen(start)]
pub fn init() {
    console_error_panic_hook::set_once();
}

/// Web-compatible simulation configuration.
#[derive(Serialize, Deserialize)]
pub struct WebConfig {
    pub system: String,           // "hydrogen", "h2_ion", "h2_molecule"
    pub algorithm: String,        // "pure", "importance_sampled"
    pub num_walkers: usize,
    pub time_step: f64,
    pub total_steps: usize,
    pub equilibration_steps: usize,
    pub seed: Option<u64>,
    pub bond_length: Option<f64>, // For H2+ and H2
    pub alpha: Option<f64>,       // Trial wavefunction parameter
    pub jastrow_b: Option<f64>,   // For H2 molecule
}

impl Default for WebConfig {
    fn default() -> Self {
        Self {
            system: "hydrogen".to_string(),
            algorithm: "importance_sampled".to_string(),
            num_walkers: 1000,
            time_step: 0.01,
            total_steps: 5000,
            equilibration_steps: 500,
            seed: Some(42),
            bond_length: None,
            alpha: Some(1.0),
            jastrow_b: None,
        }
    }
}

/// Web-compatible step result.
#[derive(Serialize, Deserialize)]
pub struct WebStepResult {
    pub step: usize,
    pub imaginary_time: f64,
    pub energy_estimate: f64,
    pub trial_energy: f64,
    pub population_size: usize,
    pub acceptance_ratio: Option<f64>,
    pub is_equilibrated: bool,
}

impl From<&SimulationStep> for WebStepResult {
    fn from(step: &SimulationStep) -> Self {
        Self {
            step: step.step,
            imaginary_time: step.imaginary_time,
            energy_estimate: step.energy_estimate,
            trial_energy: step.trial_energy,
            population_size: step.population_size,
            acceptance_ratio: step.acceptance_ratio,
            is_equilibrated: step.is_equilibrated,
        }
    }
}

/// Web-compatible statistics.
#[derive(Serialize, Deserialize, Default)]
pub struct WebStatistics {
    pub mean_energy: f64,
    pub energy_error: f64,
    pub sample_count: usize,
    pub exact_energy: Option<f64>,
}

/// WebAssembly-compatible simulation wrapper.
#[wasm_bindgen]
pub struct WebSimulation {
    state: Option<SimulationState>,
    config: WebConfig,
    update_interval: usize,
    energy_trace: Vec<f64>,
    population_trace: Vec<usize>,
}

#[wasm_bindgen]
impl WebSimulation {
    /// Create a new simulation from JSON configuration.
    #[wasm_bindgen(constructor)]
    pub fn new(config_json: &str) -> Result<WebSimulation, JsValue> {
        let config: WebConfig = serde_json::from_str(config_json)
            .map_err(|e| JsValue::from_str(&format!("Config parse error: {}", e)))?;

        let sim_config = build_config(&config)?;
        let state = Simulation::new(sim_config)
            .map_err(|e| JsValue::from_str(&format!("Simulation error: {}", e)))?;

        Ok(Self {
            state: Some(state),
            config,
            update_interval: 1,
            energy_trace: Vec::new(),
            population_trace: Vec::new(),
        })
    }

    /// Create a simulation with default configuration.
    #[wasm_bindgen]
    pub fn new_default() -> Result<WebSimulation, JsValue> {
        let config = WebConfig::default();
        let config_json = serde_json::to_string(&config).unwrap();
        Self::new(&config_json)
    }

    /// Run a single step and return the result as JSON.
    #[wasm_bindgen]
    pub fn step(&mut self) -> Result<JsValue, JsValue> {
        let state = self
            .state
            .as_mut()
            .ok_or_else(|| JsValue::from_str("Simulation not initialized"))?;

        match state.advance() {
            Some(step) => {
                if step.is_equilibrated {
                    self.energy_trace.push(step.energy_estimate);
                    self.population_trace.push(step.population_size);
                }
                let result = WebStepResult::from(&step);
                serde_wasm_bindgen::to_value(&result)
                    .map_err(|e| JsValue::from_str(&format!("Serialization error: {}", e)))
            }
            None => Ok(JsValue::NULL),
        }
    }

    /// Run multiple steps and return the last result.
    #[wasm_bindgen]
    pub fn run_steps(&mut self, n: usize) -> Result<JsValue, JsValue> {
        let mut last_result = None;

        for _ in 0..n {
            let state = self
                .state
                .as_mut()
                .ok_or_else(|| JsValue::from_str("Simulation not initialized"))?;

            match state.advance() {
                Some(step) => {
                    if step.is_equilibrated {
                        self.energy_trace.push(step.energy_estimate);
                        self.population_trace.push(step.population_size);
                    }
                    last_result = Some(WebStepResult::from(&step));
                }
                None => break,
            }
        }

        match last_result {
            Some(result) => serde_wasm_bindgen::to_value(&result)
                .map_err(|e| JsValue::from_str(&format!("Serialization error: {}", e))),
            None => Ok(JsValue::NULL),
        }
    }

    /// Get 2D electron density on a plane.
    ///
    /// Returns a flat array of density values for a grid.
    #[wasm_bindgen]
    pub fn get_density_2d(
        &self,
        plane: &str,
        slice_pos: f64,
        grid_size: usize,
        bounds: f64,
        sigma: f64,
    ) -> Result<Vec<f64>, JsValue> {
        let state = self
            .state
            .as_ref()
            .ok_or_else(|| JsValue::from_str("Simulation not initialized"))?;

        let plane = match plane {
            "xy" | "XY" => Plane::XY,
            "xz" | "XZ" => Plane::XZ,
            "yz" | "YZ" => Plane::YZ,
            _ => return Err(JsValue::from_str("Invalid plane: use xy, xz, or yz")),
        };

        let density = compute_density_2d(
            state.population(),
            plane,
            slice_pos,
            (-bounds, bounds),
            grid_size,
            sigma,
        );

        Ok(density)
    }

    /// Get current statistics as JSON.
    #[wasm_bindgen]
    pub fn get_statistics(&self) -> Result<JsValue, JsValue> {
        let state = self
            .state
            .as_ref()
            .ok_or_else(|| JsValue::from_str("Simulation not initialized"))?;

        let stats = if !self.energy_trace.is_empty() {
            let n = self.energy_trace.len() as f64;
            let mean = self.energy_trace.iter().sum::<f64>() / n;
            let variance: f64 = self
                .energy_trace
                .iter()
                .map(|e| (e - mean).powi(2))
                .sum::<f64>()
                / (n - 1.0).max(1.0);
            let error = (variance / n).sqrt();

            WebStatistics {
                mean_energy: mean,
                energy_error: error,
                sample_count: self.energy_trace.len(),
                exact_energy: state.exact_energy(),
            }
        } else {
            WebStatistics {
                exact_energy: state.exact_energy(),
                ..Default::default()
            }
        };

        serde_wasm_bindgen::to_value(&stats)
            .map_err(|e| JsValue::from_str(&format!("Serialization error: {}", e)))
    }

    /// Set the update interval (for UI throttling).
    #[wasm_bindgen]
    pub fn set_update_interval(&mut self, interval: usize) {
        self.update_interval = interval.max(1);
    }

    /// Get current step number.
    #[wasm_bindgen]
    pub fn current_step(&self) -> usize {
        self.state
            .as_ref()
            .map(|s| s.current_step())
            .unwrap_or(0)
    }

    /// Get current energy estimate.
    #[wasm_bindgen]
    pub fn current_energy(&self) -> f64 {
        self.state
            .as_ref()
            .map(|s| s.current_trial_energy())
            .unwrap_or(0.0)
    }

    /// Get current population size.
    #[wasm_bindgen]
    pub fn population_size(&self) -> usize {
        self.state
            .as_ref()
            .map(|s| s.population_size())
            .unwrap_or(0)
    }

    /// Check if simulation is past equilibration.
    #[wasm_bindgen]
    pub fn is_equilibrated(&self) -> bool {
        self.state
            .as_ref()
            .map(|s| s.is_equilibrated())
            .unwrap_or(false)
    }

    /// Check if simulation is finished.
    #[wasm_bindgen]
    pub fn is_finished(&self) -> bool {
        self.state
            .as_ref()
            .map(|s| s.is_finished())
            .unwrap_or(true)
    }

    /// Get total steps.
    #[wasm_bindgen]
    pub fn total_steps(&self) -> usize {
        self.state
            .as_ref()
            .map(|s| s.total_steps())
            .unwrap_or(0)
    }

    /// Get equilibration steps.
    #[wasm_bindgen]
    pub fn equilibration_steps(&self) -> usize {
        self.state
            .as_ref()
            .map(|s| s.equilibration_steps())
            .unwrap_or(0)
    }

    /// Get system name.
    #[wasm_bindgen]
    pub fn system_name(&self) -> String {
        self.state
            .as_ref()
            .map(|s| s.system_name().to_string())
            .unwrap_or_default()
    }

    /// Get exact energy (if known).
    #[wasm_bindgen]
    pub fn exact_energy(&self) -> Option<f64> {
        self.state.as_ref().and_then(|s| s.exact_energy())
    }

    /// Reset the simulation.
    #[wasm_bindgen]
    pub fn reset(&mut self) -> Result<(), JsValue> {
        let sim_config = build_config(&self.config)?;
        let state = Simulation::new(sim_config)
            .map_err(|e| JsValue::from_str(&format!("Simulation error: {}", e)))?;
        self.state = Some(state);
        self.energy_trace.clear();
        self.population_trace.clear();
        Ok(())
    }

    /// Get the energy trace as JSON array.
    #[wasm_bindgen]
    pub fn get_energy_trace(&self) -> Result<JsValue, JsValue> {
        serde_wasm_bindgen::to_value(&self.energy_trace)
            .map_err(|e| JsValue::from_str(&format!("Serialization error: {}", e)))
    }

    /// Get the population trace as JSON array.
    #[wasm_bindgen]
    pub fn get_population_trace(&self) -> Result<JsValue, JsValue> {
        serde_wasm_bindgen::to_value(&self.population_trace)
            .map_err(|e| JsValue::from_str(&format!("Serialization error: {}", e)))
    }
}

/// Build SimulationConfig from WebConfig.
fn build_config(config: &WebConfig) -> Result<SimulationConfig, JsValue> {
    let system = match config.system.as_str() {
        "hydrogen" | "h" => SystemType::Hydrogen,
        "h2_ion" | "h2+" => SystemType::H2Ion {
            bond_length: config.bond_length.unwrap_or(2.0),
        },
        "h2_molecule" | "h2" => SystemType::H2Molecule {
            bond_length: config.bond_length.unwrap_or(1.4),
        },
        _ => return Err(JsValue::from_str("Invalid system type")),
    };

    let algorithm = match config.algorithm.as_str() {
        "pure" => AlgorithmType::Pure { max_offspring: 3 },
        "importance_sampled" | "importance" => {
            AlgorithmType::ImportanceSampled { max_offspring: 3 }
        }
        _ => return Err(JsValue::from_str("Invalid algorithm type")),
    };

    let trial_wf = match &system {
        SystemType::Hydrogen => Some(TrialWfParams::Hydrogen {
            alpha: config.alpha.unwrap_or(1.0),
        }),
        SystemType::H2Ion { .. } => Some(TrialWfParams::H2Ion {
            alpha: config.alpha.unwrap_or(1.0),
        }),
        SystemType::H2Molecule { .. } => Some(TrialWfParams::H2Molecule {
            alpha: config.alpha.unwrap_or(1.2),
            jastrow_b: config.jastrow_b.unwrap_or(0.5),
        }),
    };

    let mut builder = SimulationConfig::builder()
        .num_walkers(config.num_walkers)
        .time_step(config.time_step)
        .total_steps(config.total_steps)
        .equilibration_steps(config.equilibration_steps)
        .system(system)
        .algorithm(algorithm)
        .feedback_alpha(1.0);

    if let Some(seed) = config.seed {
        builder = builder.seed(seed);
    }

    if let Some(wf) = trial_wf {
        builder = builder.trial_wavefunction(wf);
    }

    builder
        .build()
        .map_err(|e| JsValue::from_str(&format!("Config error: {}", e)))
}

// Set up console error panic hook
mod console_error_panic_hook {
    use std::panic;

    pub fn set_once() {
        static SET: std::sync::Once = std::sync::Once::new();
        SET.call_once(|| {
            panic::set_hook(Box::new(|info| {
                let msg = info.to_string();
                web_sys::console::error_1(&wasm_bindgen::JsValue::from_str(&msg));
            }));
        });
    }
}
