//! Physics calculations for DMC.
//!
//! This module contains the core physics functions:
//! - Local energy calculation
//! - Green's functions
//! - Branching weights

use crate::system::{QuantumSystem, TrialWavefunction};
use crate::Vec3;

/// Calculate the local energy at given electron positions.
///
/// # Mathematical Definition
///
/// The local energy is defined as the action of the Hamiltonian on the
/// trial wavefunction, divided by the trial wavefunction:
///
/// ```text
/// E_L(x) = ĤΨ_T(x)/Ψ_T(x) = -(ℏ²/2m)(∇²Ψ_T/Ψ_T) + V(x)
/// ```
///
/// In atomic units (ℏ = m = 1):
/// ```text
/// E_L(x) = -½(∇²Ψ_T/Ψ_T) + V(x)
/// ```
///
/// # Cusp Condition
///
/// If the trial wavefunction satisfies the Kato cusp condition (α = Z),
/// the singularities in kinetic and potential terms cancel, making E_L bounded.
///
/// For hydrogen with α = Z = 1:
/// ```text
/// E_L(r) = -α²/2 + (α - Z)/r → E_L = -½ (constant)
/// ```
///
/// # References
///
/// - [README: The Local Energy](https://github.com/subinbg/Diffusion-MC#the-local-energy)
/// - [README: The Cusp Condition](https://github.com/subinbg/Diffusion-MC#the-cusp-condition-stability)
///
/// # Arguments
///
/// * `system` - The quantum system providing the potential V(x)
/// * `trial_wf` - Trial wavefunction providing ∇²Ψ_T/Ψ_T
/// * `electrons` - Current electron positions
///
/// # Returns
///
/// Local energy in Hartree atomic units.
pub fn local_energy<S, T>(system: &S, trial_wf: &T, electrons: &[Vec3]) -> f64
where
    S: QuantumSystem + ?Sized,
    T: TrialWavefunction + ?Sized,
{
    // Kinetic energy: -½ ∇²Ψ_T/Ψ_T (in atomic units)
    let kinetic = -0.5 * trial_wf.laplacian_ratio(electrons);

    // Potential energy: V(x)
    let potential = system.potential(electrons);

    kinetic + potential
}

/// Calculate drift velocity from trial wavefunction gradient.
///
/// # Equation
///
/// ```text
/// v_D = (ℏ/m) ∇ln|Ψ_T|
/// ```
///
/// In atomic units (ℏ = m = 1):
/// ```text
/// v_D = ∇ln|Ψ_T|
/// ```
///
/// # References
///
/// - [README: The Drift-Diffusion Green's Function](https://github.com/subinbg/Diffusion-MC#the-drift-diffusion-greens-function)
///
/// # Arguments
///
/// * `trial_wf` - Trial wavefunction
/// * `electrons` - Current electron positions
///
/// # Returns
///
/// Drift velocity for each electron (Vec3 per electron).
pub fn drift_velocity<T>(trial_wf: &T, electrons: &[Vec3]) -> Vec<Vec3>
where
    T: TrialWavefunction + ?Sized,
{
    // In atomic units, v_D = ∇ln|Ψ_T|
    trial_wf.gradient_log(electrons)
}

/// Calculate branching weight for pure DMC.
///
/// # Equation
///
/// ```text
/// W = exp(-δτ(V(x) - E_T)/ℏ)
/// ```
///
/// In atomic units (ℏ = 1):
/// ```text
/// W = exp(-δτ(V(x) - E_T))
/// ```
///
/// # References
///
/// - [README: Step 3: Weighting and Branching](https://github.com/subinbg/Diffusion-MC#step-3-weighting-and-branching-potential-update)
///
/// # Arguments
///
/// * `potential` - V(x) at current position
/// * `trial_energy` - Current trial energy E_T
/// * `time_step` - Imaginary time step δτ
///
/// # Returns
///
/// Branching weight W.
pub fn branching_weight_pure(potential: f64, trial_energy: f64, time_step: f64) -> f64 {
    (-(potential - trial_energy) * time_step).exp()
}

/// Calculate branching weight for importance-sampled DMC.
///
/// # Equation
///
/// Uses the average local energy at old and new positions:
///
/// ```text
/// W = exp(-δτ((E_L(x') + E_L(x))/2 - E_T))
/// ```
///
/// This averaging reduces time-step error to O(δτ²).
///
/// # References
///
/// - [README: Step 4: Branching](https://github.com/subinbg/Diffusion-MC#step-4-branching-the-potential-step)
///
/// # Arguments
///
/// * `local_energy_old` - E_L at old position
/// * `local_energy_new` - E_L at new position
/// * `trial_energy` - Current trial energy E_T
/// * `time_step` - Imaginary time step δτ
///
/// # Returns
///
/// Branching weight W.
pub fn branching_weight_importance(
    local_energy_old: f64,
    local_energy_new: f64,
    trial_energy: f64,
    time_step: f64,
) -> f64 {
    let e_l_avg = 0.5 * (local_energy_old + local_energy_new);
    (-(e_l_avg - trial_energy) * time_step).exp()
}

/// Calculate Metropolis acceptance probability for importance-sampled DMC.
///
/// # Equation
///
/// The acceptance probability corrects for the asymmetry of the Green's function:
///
/// ```text
/// A = min(1, |Ψ_T(x')|² G(x←x') / |Ψ_T(x)|² G(x'←x))
/// ```
///
/// Explicitly:
/// ```text
/// A = min(1, |Ψ_T(x')|²/|Ψ_T(x)|² × exp[-(|x-x'-v_D(x')δτ|² - |x'-x-v_D(x)δτ|²)/(2σ²)])
/// ```
///
/// where σ² = δτ (in atomic units).
///
/// # References
///
/// - [README: Step 3: Metropolis Acceptance/Rejection](https://github.com/subinbg/Diffusion-MC#step-3-metropolis-acceptancerejection)
/// - [README: Asymmetry of the Green's Function](https://github.com/subinbg/Diffusion-MC#asymmetry-of-the-greens-function)
///
/// # Arguments
///
/// * `psi_old` - Ψ_T(x) at old position
/// * `psi_new` - Ψ_T(x') at proposed new position
/// * `old_positions` - Old electron positions
/// * `new_positions` - Proposed new electron positions
/// * `drift_old` - Drift velocity at old position
/// * `drift_new` - Drift velocity at new position
/// * `time_step` - Imaginary time step δτ
///
/// # Returns
///
/// Acceptance probability A ∈ [0, 1].
pub fn metropolis_acceptance(
    psi_old: f64,
    psi_new: f64,
    old_positions: &[Vec3],
    new_positions: &[Vec3],
    drift_old: &[Vec3],
    drift_new: &[Vec3],
    time_step: f64,
) -> f64 {
    // |Ψ_T(x')|² / |Ψ_T(x)|²
    let psi_ratio_squared = (psi_new / psi_old).powi(2);

    // Green's function ratio: G(x←x') / G(x'←x)
    // G(x←x') ∝ exp(-|x - x' - v_D(x')δτ|² / (2σ²))
    // G(x'←x) ∝ exp(-|x' - x - v_D(x)δτ|² / (2σ²))
    let sigma_sq = time_step; // σ² = δτ in atomic units

    let mut exponent = 0.0;
    for i in 0..old_positions.len() {
        // Forward move: x → x'
        let forward_diff = new_positions[i] - old_positions[i] - drift_old[i] * time_step;
        let forward_dist_sq = forward_diff.norm_squared();

        // Backward move: x' → x
        let backward_diff = old_positions[i] - new_positions[i] - drift_new[i] * time_step;
        let backward_dist_sq = backward_diff.norm_squared();

        exponent += (forward_dist_sq - backward_dist_sq) / (2.0 * sigma_sq);
    }

    let green_ratio = exponent.exp();

    (psi_ratio_squared * green_ratio).min(1.0)
}

/// Update trial energy using population control feedback.
///
/// # Equation
///
/// ```text
/// E_T(τ) = ⟨E_L⟩ - α ln(N/N₀)
/// ```
///
/// where:
/// - ⟨E_L⟩ is the average local energy (or average potential for pure DMC)
/// - α is the feedback parameter
/// - N is current population, N₀ is target population
///
/// # Stability
///
/// This feedback law ensures exponential stability:
/// ```text
/// δN/δτ ≈ -α δN / ℏ
/// ```
///
/// Leading to δ(τ) ∝ exp(-ατ/ℏ).
///
/// # References
///
/// - [README: Population Control Bias and Feedback Law](https://github.com/subinbg/Diffusion-MC#population-control-bias-and-feedback-law)
/// - [README: Step 5: Updating E_T](https://github.com/subinbg/Diffusion-MC#step-5-updating-e_t-population-control)
///
/// # Arguments
///
/// * `avg_energy` - Average local energy (or potential) of current walkers
/// * `current_population` - Current number of walkers N
/// * `target_population` - Target number of walkers N₀
/// * `feedback_alpha` - Feedback parameter α
///
/// # Returns
///
/// Updated trial energy E_T.
pub fn update_trial_energy(
    avg_energy: f64,
    current_population: usize,
    target_population: usize,
    feedback_alpha: f64,
) -> f64 {
    let ratio = current_population as f64 / target_population as f64;
    avg_energy - feedback_alpha * ratio.ln()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::system::{Hydrogen, HydrogenTrialWf};
    use approx::assert_relative_eq;

    #[test]
    fn local_energy_hydrogen_cusp() {
        // With α = 1 (cusp condition), E_L should be -0.5 everywhere
        let h = Hydrogen::new();
        let wf = HydrogenTrialWf::new(1.0);

        for r in [0.1, 0.5, 1.0, 2.0, 5.0] {
            let pos = vec![Vec3::new(r, 0.0, 0.0)];
            let e_l = local_energy(&h, &wf, &pos);
            assert_relative_eq!(e_l, -0.5, epsilon = 1e-6);
        }
    }

    #[test]
    fn branching_weight_unity_at_equilibrium() {
        // If V = E_T, weight should be 1.0
        let w = branching_weight_pure(-0.5, -0.5, 0.01);
        assert_relative_eq!(w, 1.0, epsilon = 1e-10);
    }

    #[test]
    fn branching_weight_less_than_one_high_energy() {
        // If V > E_T, walker should die (W < 1)
        let w = branching_weight_pure(0.0, -0.5, 0.01);
        assert!(w < 1.0);
    }

    #[test]
    fn branching_weight_greater_than_one_low_energy() {
        // If V < E_T, walker should multiply (W > 1)
        let w = branching_weight_pure(-1.0, -0.5, 0.01);
        assert!(w > 1.0);
    }

    #[test]
    fn trial_energy_feedback() {
        // If N > N₀, E_T should decrease
        let e_t_1 = update_trial_energy(-0.5, 1100, 1000, 1.0);
        let e_t_2 = update_trial_energy(-0.5, 1000, 1000, 1.0);
        assert!(e_t_1 < e_t_2);

        // If N < N₀, E_T should increase
        let e_t_3 = update_trial_energy(-0.5, 900, 1000, 1.0);
        assert!(e_t_3 > e_t_2);
    }
}
