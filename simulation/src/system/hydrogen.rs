//! Hydrogen atom (H): 1 electron, 1 nucleus.
//!
//! # Exact Solution
//!
//! - Ground state energy: E₀ = -0.5 Hartree = -13.6 eV
//! - Ground state wavefunction: Ψ₀(r) = (1/√π) exp(-r) (in atomic units)
//!
//! # Trial Wavefunction
//!
//! ```text
//! Ψ_T(r) = exp(-αr)
//! ```
//!
//! The Kato cusp condition requires α = Z = 1 for bounded local energy.

use crate::system::{NuclearConfig, QuantumSystem, TrialWavefunction};
use crate::Vec3;

/// Hydrogen atom system.
///
/// Single electron in the Coulomb field of a proton at the origin.
#[derive(Clone, Debug)]
pub struct Hydrogen {
    nuclear_config: NuclearConfig,
    trial_wf: Option<HydrogenTrialWf>,
}

impl Default for Hydrogen {
    fn default() -> Self {
        Self::new()
    }
}

impl Hydrogen {
    /// Create a new hydrogen atom with nucleus at the origin.
    pub fn new() -> Self {
        Self {
            nuclear_config: NuclearConfig {
                positions: vec![Vec3::zeros()],
                charges: vec![1.0],
            },
            trial_wf: None,
        }
    }

    /// Add a trial wavefunction for importance sampling.
    ///
    /// # Arguments
    ///
    /// * `alpha` - Variational parameter. Use α = 1.0 for cusp condition.
    pub fn with_trial_wavefunction(mut self, alpha: f64) -> Self {
        self.trial_wf = Some(HydrogenTrialWf::new(alpha));
        self
    }
}

impl QuantumSystem for Hydrogen {
    fn num_electrons(&self) -> usize {
        1
    }

    fn nuclear_config(&self) -> &NuclearConfig {
        &self.nuclear_config
    }

    /// Coulomb potential V(r) = -Z/r = -1/r for hydrogen.
    ///
    /// # Equation
    ///
    /// ```text
    /// V(r) = -1/r
    /// ```
    ///
    /// See [README: Step 3](https://github.com/subinbg/Diffusion-MC#step-3-weighting-and-branching-potential-update)
    fn potential(&self, electrons: &[Vec3]) -> f64 {
        let r = electrons[0].norm();
        if r < 1e-10 {
            return -1e10; // Regularization for numerical stability
        }
        -1.0 / r
    }

    fn trial_wavefunction(&self) -> Option<&dyn TrialWavefunction> {
        self.trial_wf.as_ref().map(|wf| wf as &dyn TrialWavefunction)
    }

    fn exact_energy(&self) -> Option<f64> {
        Some(-0.5)
    }
}

/// Trial wavefunction for hydrogen: Ψ_T(r) = exp(-αr).
///
/// # Cusp Condition
///
/// For the local energy E_L to be bounded, α must equal Z = 1.
///
/// With α = Z, the kinetic singularity (+α/r) exactly cancels the
/// potential singularity (-Z/r), giving:
///
/// ```text
/// E_L(r) = -α²/2 + (α - Z)/r  →  E_L = -Z²/2 = -0.5 Ha (constant)
/// ```
///
/// See [README: The Cusp Condition](https://github.com/subinbg/Diffusion-MC#the-cusp-condition-stability)
#[derive(Clone, Debug)]
pub struct HydrogenTrialWf {
    /// Variational parameter (cusp condition: α = Z = 1).
    pub alpha: f64,
}

impl HydrogenTrialWf {
    /// Create a new hydrogen trial wavefunction.
    ///
    /// # Arguments
    ///
    /// * `alpha` - Variational parameter. Use 1.0 for optimal (cusp condition).
    pub fn new(alpha: f64) -> Self {
        Self { alpha }
    }
}

impl TrialWavefunction for HydrogenTrialWf {
    /// Ψ_T(r) = exp(-αr)
    fn value(&self, electrons: &[Vec3]) -> f64 {
        let r = electrons[0].norm();
        (-self.alpha * r).exp()
    }

    fn log_value(&self, electrons: &[Vec3]) -> f64 {
        let r = electrons[0].norm();
        -self.alpha * r
    }

    /// ∇ln Ψ_T = ∇(-αr) = -α r̂
    ///
    /// # Derivation
    ///
    /// ```text
    /// ln Ψ_T = -αr
    /// ∇(ln Ψ_T) = -α ∇r = -α (r/|r|) = -α r̂
    /// ```
    fn gradient_log(&self, electrons: &[Vec3]) -> Vec<Vec3> {
        let r_vec = electrons[0];
        let r_norm = r_vec.norm();
        if r_norm < 1e-10 {
            return vec![Vec3::zeros()];
        }
        // -α * r̂ = -α * r/|r|
        vec![-self.alpha * r_vec / r_norm]
    }

    /// ∇²Ψ_T/Ψ_T = α² - 2α/r
    ///
    /// # Derivation
    ///
    /// For Ψ_T = exp(-αr) in 3D spherical coordinates:
    ///
    /// ```text
    /// ∇² = d²/dr² + (2/r)(d/dr)
    ///
    /// dΨ_T/dr = -α Ψ_T
    /// d²Ψ_T/dr² = α² Ψ_T
    ///
    /// ∇²Ψ_T = α² Ψ_T - (2α/r) Ψ_T
    /// ∇²Ψ_T/Ψ_T = α² - 2α/r
    /// ```
    ///
    /// See [README: Calculating the Kinetic Term](https://github.com/subinbg/Diffusion-MC#calculating-the-kinetic-term)
    fn laplacian_ratio(&self, electrons: &[Vec3]) -> f64 {
        let r = electrons[0].norm();
        if r < 1e-10 {
            // At origin, return just the constant term
            return self.alpha * self.alpha;
        }
        self.alpha * self.alpha - 2.0 * self.alpha / r
    }

    fn parameters(&self) -> Vec<f64> {
        vec![self.alpha]
    }

    fn set_parameters(&mut self, params: &[f64]) {
        if !params.is_empty() {
            self.alpha = params[0];
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn hydrogen_potential_at_bohr_radius() {
        let h = Hydrogen::new();
        let r = vec![Vec3::new(1.0, 0.0, 0.0)];
        assert_relative_eq!(h.potential(&r), -1.0, epsilon = 1e-10);
    }

    #[test]
    fn hydrogen_potential_at_half_bohr() {
        let h = Hydrogen::new();
        let r = vec![Vec3::new(0.5, 0.0, 0.0)];
        assert_relative_eq!(h.potential(&r), -2.0, epsilon = 1e-10);
    }

    #[test]
    fn trial_wf_value_at_origin() {
        let wf = HydrogenTrialWf::new(1.0);
        let r = vec![Vec3::zeros()];
        assert_relative_eq!(wf.value(&r), 1.0, epsilon = 1e-10);
    }

    #[test]
    fn trial_wf_value_at_bohr_radius() {
        let wf = HydrogenTrialWf::new(1.0);
        let r = vec![Vec3::new(1.0, 0.0, 0.0)];
        assert_relative_eq!(wf.value(&r), (-1.0_f64).exp(), epsilon = 1e-10);
    }

    #[test]
    fn gradient_log_direction() {
        let wf = HydrogenTrialWf::new(1.0);
        let r = vec![Vec3::new(1.0, 0.0, 0.0)];
        let grad = wf.gradient_log(&r);
        // Should point toward origin (negative x direction)
        assert!(grad[0].x < 0.0);
        assert_relative_eq!(grad[0].y, 0.0, epsilon = 1e-10);
        assert_relative_eq!(grad[0].z, 0.0, epsilon = 1e-10);
    }

    #[test]
    fn cusp_condition_local_energy_constant() {
        // With α = 1, E_L should be -0.5 everywhere
        let h = Hydrogen::new().with_trial_wavefunction(1.0);
        let wf = h.trial_wavefunction().unwrap();

        for r in [0.1, 0.5, 1.0, 2.0, 5.0] {
            let pos = vec![Vec3::new(r, 0.0, 0.0)];
            // E_L = -0.5 * laplacian_ratio + V
            let kinetic = -0.5 * wf.laplacian_ratio(&pos);
            let potential = h.potential(&pos);
            let e_l = kinetic + potential;
            assert_relative_eq!(e_l, -0.5, epsilon = 1e-6);
        }
    }
}
