//! Importance-Sampled Diffusion Monte Carlo algorithm.
//!
//! # Algorithm
//!
//! Importance sampling introduces a trial wavefunction Ψ_T to:
//! - Guide walkers toward physically relevant regions
//! - Replace divergent V(x) with bounded E_L(x)
//!
//! Steps:
//! 1. **Drift-Diffusion**: x' = x + v_D(x)δτ + √(δτ) × ξ
//! 2. **Metropolis**: Accept with probability A (restores detailed balance)
//! 3. **Branching**: M = ⌊W + u⌋, where W = exp(-δτ(E_L-E_T))
//!
//! See [README: Importance Sampled DMC](https://github.com/subinbg/Diffusion-MC#importance-sampled-diffusion-monte-carlo)

use crate::algorithm::{DmcAlgorithm, StepResult};
use crate::config::{SimulationConfig, SystemType, TrialWfParams};
use crate::physics::{branching_weight_importance, drift_velocity, local_energy, metropolis_acceptance};
use crate::system::{H2Ion, H2Molecule, Hydrogen, QuantumSystem, TrialWavefunction};
use crate::walker::Population;
use rand::Rng;
use rand_distr::StandardNormal;

/// Importance-Sampled DMC algorithm implementation.
///
/// # Stability
///
/// Unlike pure DMC, importance sampling is stable at Coulomb singularities
/// because the local energy E_L(x) remains bounded if the trial wavefunction
/// satisfies the cusp condition.
pub struct ImportanceSampledDmc {
    /// Maximum offspring per walker.
    pub max_offspring: usize,
    /// The quantum system being simulated.
    system: Box<dyn QuantumSystem>,
    /// The trial wavefunction.
    trial_wf: Box<dyn TrialWavefunction>,
}

impl ImportanceSampledDmc {
    /// Create a new Importance-Sampled DMC algorithm.
    pub fn new(config: &SimulationConfig) -> Self {
        let max_offspring = match &config.algorithm {
            crate::config::AlgorithmType::ImportanceSampled { max_offspring } => *max_offspring,
            _ => 3,
        };

        let (system, trial_wf): (Box<dyn QuantumSystem>, Box<dyn TrialWavefunction>) =
            match (&config.system, &config.trial_wavefunction) {
                (SystemType::Hydrogen, Some(TrialWfParams::Hydrogen { alpha })) => {
                    let h = Hydrogen::new().with_trial_wavefunction(*alpha);
                    let wf = crate::system::HydrogenTrialWf::new(*alpha);
                    (Box::new(h), Box::new(wf))
                }
                (SystemType::H2Ion { bond_length }, Some(TrialWfParams::H2Ion { alpha })) => {
                    let h2 = H2Ion::new(*bond_length).with_trial_wavefunction(*alpha);
                    let wf = crate::system::H2IonTrialWf::new(
                        *alpha,
                        h2.nuclear_config().positions.clone(),
                    );
                    (Box::new(h2), Box::new(wf))
                }
                (
                    SystemType::H2Molecule { bond_length },
                    Some(TrialWfParams::H2Molecule { alpha, jastrow_b }),
                ) => {
                    let h2 = H2Molecule::new(*bond_length).with_trial_wavefunction(*alpha, *jastrow_b);
                    let nuc = h2.nuclear_config().positions.clone();
                    let wf = crate::system::H2TrialWf::new(*alpha, *jastrow_b, [nuc[0], nuc[1]]);
                    (Box::new(h2), Box::new(wf))
                }
                // Default: Hydrogen with α = 1
                (SystemType::Hydrogen, None) => {
                    let h = Hydrogen::new().with_trial_wavefunction(1.0);
                    let wf = crate::system::HydrogenTrialWf::new(1.0);
                    (Box::new(h), Box::new(wf))
                }
                _ => panic!("Mismatched system and trial wavefunction types"),
            };

        Self {
            max_offspring,
            system,
            trial_wf,
        }
    }

    /// Get reference to the quantum system.
    pub fn system(&self) -> &dyn QuantumSystem {
        self.system.as_ref()
    }

    /// Get reference to the trial wavefunction.
    pub fn trial_wf(&self) -> &dyn TrialWavefunction {
        self.trial_wf.as_ref()
    }
}

impl DmcAlgorithm for ImportanceSampledDmc {
    fn step(
        &mut self,
        population: &mut Population,
        trial_energy: f64,
        config: &SimulationConfig,
    ) -> StepResult {
        let time_step = config.time_step;
        let sigma = time_step.sqrt();
        let num_walkers = population.size();

        // First, collect electron counts for each walker
        let electron_counts: Vec<usize> = population.iter().map(|w| w.positions.len()).collect();

        // Pre-generate all random numbers we need
        let rng = population.rng();
        let mut gauss_randoms: Vec<Vec<[f64; 3]>> = Vec::with_capacity(num_walkers);
        let mut uniform_randoms: Vec<f64> = Vec::with_capacity(num_walkers);

        for &num_electrons in &electron_counts {
            let mut walker_randoms = Vec::with_capacity(num_electrons);
            for _ in 0..num_electrons {
                walker_randoms.push([
                    rng.sample::<f64, _>(StandardNormal),
                    rng.sample::<f64, _>(StandardNormal),
                    rng.sample::<f64, _>(StandardNormal),
                ]);
            }
            gauss_randoms.push(walker_randoms);
            uniform_randoms.push(rng.gen());
        }

        let mut accepted = 0usize;
        let mut local_energies_new = Vec::with_capacity(num_walkers);
        let mut weights = Vec::with_capacity(num_walkers);
        let mut new_positions_list = Vec::with_capacity(num_walkers);

        // Process each walker (read-only iteration to collect updates)
        for (i, walker) in population.iter().enumerate() {
            // Calculate old state properties
            let psi_old = self.trial_wf.value(&walker.positions);
            let drift_old = drift_velocity(&*self.trial_wf, &walker.positions);
            let e_l_old = local_energy(
                &*self.system,
                &*self.trial_wf,
                &walker.positions,
            );

            // Step 1: Propose new position via drift-diffusion
            // x' = x + v_D(x)δτ + √(δτ) × ξ
            // See README: Step 2: Drift-Diffusion Update
            let mut new_positions = walker.positions.clone();
            for (j, (pos, drift)) in new_positions.iter_mut().zip(drift_old.iter()).enumerate() {
                let randoms = &gauss_randoms[i][j];
                pos.x += drift.x * time_step + randoms[0] * sigma;
                pos.y += drift.y * time_step + randoms[1] * sigma;
                pos.z += drift.z * time_step + randoms[2] * sigma;
            }

            // Calculate new state properties
            let psi_new = self.trial_wf.value(&new_positions);
            let drift_new = drift_velocity(&*self.trial_wf, &new_positions);
            let e_l_new = local_energy(
                &*self.system,
                &*self.trial_wf,
                &new_positions,
            );

            // Step 2: Metropolis acceptance/rejection
            // A = min(1, |Ψ_T(x')|² G(x←x') / |Ψ_T(x)|² G(x'←x))
            // See README: Step 3: Metropolis Acceptance/Rejection
            let accept_prob = metropolis_acceptance(
                psi_old,
                psi_new,
                &walker.positions,
                &new_positions,
                &drift_old,
                &drift_new,
                time_step,
            );

            let (final_positions, final_e_l) = if uniform_randoms[i] < accept_prob {
                accepted += 1;
                (new_positions, e_l_new)
            } else {
                (walker.positions.clone(), e_l_old)
            };

            local_energies_new.push(final_e_l);

            // Step 3: Calculate branching weight
            // W = exp(-δτ((E_L(x') + E_L(x))/2 - E_T))
            // See README: Step 4: Branching
            let weight =
                branching_weight_importance(e_l_old, final_e_l, trial_energy, time_step);
            weights.push(weight);

            new_positions_list.push(final_positions);
        }

        // Apply position updates
        for (walker, new_pos) in population.iter_mut().zip(new_positions_list.into_iter()) {
            walker.positions = new_pos;
        }

        // Step 4: Apply branching
        population.branch(&weights, self.max_offspring);
        population.increment_ages();

        // Calculate average local energy
        let avg_energy = if !local_energies_new.is_empty() {
            local_energies_new.iter().sum::<f64>() / local_energies_new.len() as f64
        } else {
            trial_energy
        };

        StepResult {
            energy_estimate: avg_energy,
            population_size: population.size(),
            acceptance_ratio: Some(accepted as f64 / num_walkers as f64),
            local_energies: local_energies_new,
        }
    }

    fn name(&self) -> &'static str {
        "Importance-Sampled DMC"
    }
}
