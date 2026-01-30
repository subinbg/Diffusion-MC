//! Hydrogen molecule ion (H₂⁺): 1 electron, 2 nuclei.
//!
//! # System Description
//!
//! One electron shared between two protons separated by distance R.
//!
//! # Exact Energy
//!
//! At equilibrium bond length R ≈ 2.0 Bohr: E₀ ≈ -0.6026 Hartree
//!
//! # Trial Wavefunction (LCAO)
//!
//! Linear Combination of Atomic Orbitals:
//! ```text
//! Ψ_T = exp(-α r₁) + exp(-α r₂)
//! ```
//!
//! where r₁, r₂ are distances to the two nuclei.

use crate::system::{NuclearConfig, QuantumSystem, TrialWavefunction};
use crate::Vec3;

/// Hydrogen molecule ion (H₂⁺).
///
/// # Geometry
///
/// Two protons placed symmetrically along the x-axis:
/// - Nucleus 1 at (-R/2, 0, 0)
/// - Nucleus 2 at (+R/2, 0, 0)
#[derive(Clone, Debug)]
pub struct H2Ion {
    nuclear_config: NuclearConfig,
    bond_length: f64,
    trial_wf: Option<H2IonTrialWf>,
}

impl H2Ion {
    /// Create a new H₂⁺ ion with given internuclear distance.
    ///
    /// # Arguments
    ///
    /// * `bond_length` - Internuclear distance R in Bohr radii.
    ///   Equilibrium value is approximately 2.0 Bohr.
    pub fn new(bond_length: f64) -> Self {
        let half_r = bond_length / 2.0;
        Self {
            nuclear_config: NuclearConfig {
                positions: vec![
                    Vec3::new(-half_r, 0.0, 0.0),
                    Vec3::new(half_r, 0.0, 0.0),
                ],
                charges: vec![1.0, 1.0],
            },
            bond_length,
            trial_wf: None,
        }
    }

    /// Add a trial wavefunction for importance sampling.
    ///
    /// # Arguments
    ///
    /// * `alpha` - Variational parameter for the atomic orbitals.
    pub fn with_trial_wavefunction(mut self, alpha: f64) -> Self {
        self.trial_wf = Some(H2IonTrialWf::new(
            alpha,
            self.nuclear_config.positions.clone(),
        ));
        self
    }

    /// Get the bond length.
    pub fn bond_length(&self) -> f64 {
        self.bond_length
    }
}

impl QuantumSystem for H2Ion {
    fn num_electrons(&self) -> usize {
        1
    }

    fn nuclear_config(&self) -> &NuclearConfig {
        &self.nuclear_config
    }

    /// Potential energy for H₂⁺.
    ///
    /// # Equation
    ///
    /// ```text
    /// V(r) = -1/r₁ - 1/r₂ + 1/R
    /// ```
    ///
    /// where:
    /// - r₁ = |r - R₁| (electron-nucleus 1 distance)
    /// - r₂ = |r - R₂| (electron-nucleus 2 distance)
    /// - R = |R₁ - R₂| (internuclear distance, constant)
    fn potential(&self, electrons: &[Vec3]) -> f64 {
        let r_vec = electrons[0];
        let nuc = &self.nuclear_config.positions;

        let r1 = (r_vec - nuc[0]).norm().max(1e-10);
        let r2 = (r_vec - nuc[1]).norm().max(1e-10);

        // Electron-nuclear attraction
        let v_en = -1.0 / r1 - 1.0 / r2;

        // Nuclear-nuclear repulsion (constant for fixed nuclei)
        let v_nn = 1.0 / self.bond_length;

        v_en + v_nn
    }

    fn trial_wavefunction(&self) -> Option<&dyn TrialWavefunction> {
        self.trial_wf.as_ref().map(|wf| wf as &dyn TrialWavefunction)
    }

    fn exact_energy(&self) -> Option<f64> {
        // Depends on bond length; this is approximate for R = 2.0 Bohr
        Some(-0.6026)
    }
}

/// LCAO trial wavefunction for H₂⁺.
///
/// # Form
///
/// ```text
/// Ψ_T = exp(-α r₁) + exp(-α r₂)
/// ```
///
/// where r₁, r₂ are distances to the two nuclei.
#[derive(Clone, Debug)]
pub struct H2IonTrialWf {
    /// Variational parameter.
    pub alpha: f64,
    /// Nuclear positions (stored for gradient/laplacian calculations).
    nuclear_positions: Vec<Vec3>,
}

impl H2IonTrialWf {
    /// Create a new H₂⁺ trial wavefunction.
    pub fn new(alpha: f64, nuclear_positions: Vec<Vec3>) -> Self {
        Self {
            alpha,
            nuclear_positions,
        }
    }
}

impl TrialWavefunction for H2IonTrialWf {
    /// Ψ_T = exp(-α r₁) + exp(-α r₂)
    fn value(&self, electrons: &[Vec3]) -> f64 {
        let r_vec = electrons[0];
        let r1 = (r_vec - self.nuclear_positions[0]).norm();
        let r2 = (r_vec - self.nuclear_positions[1]).norm();

        (-self.alpha * r1).exp() + (-self.alpha * r2).exp()
    }

    /// ∇ln Ψ_T for H₂⁺ LCAO.
    ///
    /// # Derivation
    ///
    /// Let φ₁ = exp(-αr₁), φ₂ = exp(-αr₂), Ψ_T = φ₁ + φ₂
    ///
    /// ```text
    /// ∇ln Ψ_T = ∇Ψ_T / Ψ_T = (∇φ₁ + ∇φ₂) / (φ₁ + φ₂)
    ///
    /// ∇φ₁ = -α (r - R₁)/r₁ × φ₁
    /// ∇φ₂ = -α (r - R₂)/r₂ × φ₂
    /// ```
    fn gradient_log(&self, electrons: &[Vec3]) -> Vec<Vec3> {
        let r_vec = electrons[0];

        let d1 = r_vec - self.nuclear_positions[0];
        let d2 = r_vec - self.nuclear_positions[1];
        let r1 = d1.norm().max(1e-10);
        let r2 = d2.norm().max(1e-10);

        let phi1 = (-self.alpha * r1).exp();
        let phi2 = (-self.alpha * r2).exp();
        let psi = phi1 + phi2;

        if psi.abs() < 1e-10 {
            return vec![Vec3::zeros()];
        }

        // ∇φ_i = -α r̂_i φ_i
        let grad_phi1 = -self.alpha * (d1 / r1) * phi1;
        let grad_phi2 = -self.alpha * (d2 / r2) * phi2;

        vec![(grad_phi1 + grad_phi2) / psi]
    }

    /// ∇²Ψ_T/Ψ_T for H₂⁺ LCAO.
    ///
    /// # Derivation
    ///
    /// For each atomic orbital φ_i = exp(-αr_i):
    /// ```text
    /// ∇²φ_i/φ_i = α² - 2α/r_i
    /// ```
    ///
    /// For the sum:
    /// ```text
    /// ∇²Ψ_T/Ψ_T = [∇²φ₁ + ∇²φ₂] / [φ₁ + φ₂]
    /// ```
    fn laplacian_ratio(&self, electrons: &[Vec3]) -> f64 {
        let r_vec = electrons[0];

        let d1 = r_vec - self.nuclear_positions[0];
        let d2 = r_vec - self.nuclear_positions[1];
        let r1 = d1.norm().max(1e-10);
        let r2 = d2.norm().max(1e-10);

        let phi1 = (-self.alpha * r1).exp();
        let phi2 = (-self.alpha * r2).exp();
        let psi = phi1 + phi2;

        if psi.abs() < 1e-10 {
            return self.alpha * self.alpha;
        }

        // ∇²φ_i = (α² - 2α/r_i) φ_i
        let lap_phi1 = (self.alpha * self.alpha - 2.0 * self.alpha / r1) * phi1;
        let lap_phi2 = (self.alpha * self.alpha - 2.0 * self.alpha / r2) * phi2;

        (lap_phi1 + lap_phi2) / psi
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
    fn h2_ion_potential_at_midpoint() {
        let h2 = H2Ion::new(2.0); // R = 2.0 Bohr
        // Electron at midpoint between nuclei
        let r = vec![Vec3::zeros()];
        // V = -1/1 - 1/1 + 1/2 = -1.5
        assert_relative_eq!(h2.potential(&r), -1.5, epsilon = 1e-10);
    }

    #[test]
    fn h2_ion_potential_symmetry() {
        let h2 = H2Ion::new(2.0);
        let r1 = vec![Vec3::new(0.5, 0.3, 0.0)];
        let r2 = vec![Vec3::new(-0.5, 0.3, 0.0)];
        assert_relative_eq!(h2.potential(&r1), h2.potential(&r2), epsilon = 1e-10);
    }

    #[test]
    fn trial_wf_symmetry() {
        let h2 = H2Ion::new(2.0).with_trial_wavefunction(1.0);
        let wf = h2.trial_wavefunction().unwrap();

        let r1 = vec![Vec3::new(0.5, 0.3, 0.0)];
        let r2 = vec![Vec3::new(-0.5, 0.3, 0.0)];
        assert_relative_eq!(wf.value(&r1), wf.value(&r2), epsilon = 1e-10);
    }
}
