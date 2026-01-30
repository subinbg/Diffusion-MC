//! Hydrogen molecule (H₂): 2 electrons, 2 nuclei.
//!
//! # System Description
//!
//! Two electrons shared between two protons, forming a covalent bond.
//!
//! # Experimental Ground State
//!
//! - Energy: E₀ ≈ -1.1745 Hartree
//! - Equilibrium bond length: R_eq ≈ 1.4 Bohr (0.74 Å)
//!
//! # Trial Wavefunction
//!
//! Heitler-London form with Jastrow correlation factor:
//!
//! ```text
//! Ψ_T = [exp(-α(r₁ₐ + r₂ᵦ)) + exp(-α(r₁ᵦ + r₂ₐ))] × J(r₁₂)
//! ```
//!
//! where the Jastrow factor handles the electron-electron cusp:
//!
//! ```text
//! J(r₁₂) = exp(r₁₂ / (2(1 + b·r₁₂)))
//! ```

use crate::system::{NuclearConfig, QuantumSystem, TrialWavefunction};
use crate::Vec3;

/// Hydrogen molecule (H₂).
///
/// # Geometry
///
/// Two protons placed symmetrically along the x-axis:
/// - Nucleus A at (-R/2, 0, 0)
/// - Nucleus B at (+R/2, 0, 0)
#[derive(Clone, Debug)]
pub struct H2Molecule {
    nuclear_config: NuclearConfig,
    bond_length: f64,
    trial_wf: Option<H2TrialWf>,
}

impl H2Molecule {
    /// Create a new H₂ molecule with given internuclear distance.
    ///
    /// # Arguments
    ///
    /// * `bond_length` - Internuclear distance R in Bohr radii.
    ///   Equilibrium value is approximately 1.4 Bohr.
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
    /// * `alpha` - Variational parameter for atomic orbitals.
    /// * `jastrow_b` - Jastrow correlation parameter.
    pub fn with_trial_wavefunction(mut self, alpha: f64, jastrow_b: f64) -> Self {
        self.trial_wf = Some(H2TrialWf::new(
            alpha,
            jastrow_b,
            [
                self.nuclear_config.positions[0],
                self.nuclear_config.positions[1],
            ],
        ));
        self
    }

    /// Get the bond length.
    pub fn bond_length(&self) -> f64 {
        self.bond_length
    }
}

impl QuantumSystem for H2Molecule {
    fn num_electrons(&self) -> usize {
        2
    }

    fn nuclear_config(&self) -> &NuclearConfig {
        &self.nuclear_config
    }

    /// Potential energy for H₂.
    ///
    /// # Equation
    ///
    /// ```text
    /// V(r₁, r₂) = V_en + V_ee + V_nn
    ///
    /// V_en = -1/r₁ₐ - 1/r₁ᵦ - 1/r₂ₐ - 1/r₂ᵦ  (electron-nuclear)
    /// V_ee = +1/r₁₂                            (electron-electron)
    /// V_nn = +1/R                              (nuclear-nuclear)
    /// ```
    ///
    /// This matches the potential() function in dqmc.F90 (lines 81-131).
    fn potential(&self, electrons: &[Vec3]) -> f64 {
        let nuc = &self.nuclear_config.positions;

        // Electron-nuclear attraction
        let mut v_en = 0.0;
        for e in electrons {
            let r_a = (e - nuc[0]).norm().max(1e-10);
            let r_b = (e - nuc[1]).norm().max(1e-10);
            v_en -= 1.0 / r_a + 1.0 / r_b;
        }

        // Electron-electron repulsion
        let r12 = (electrons[0] - electrons[1]).norm().max(1e-10);
        let v_ee = 1.0 / r12;

        // Nuclear-nuclear repulsion
        let v_nn = 1.0 / self.bond_length;

        v_en + v_ee + v_nn
    }

    fn trial_wavefunction(&self) -> Option<&dyn TrialWavefunction> {
        self.trial_wf.as_ref().map(|wf| wf as &dyn TrialWavefunction)
    }

    fn exact_energy(&self) -> Option<f64> {
        // Experimental value at equilibrium (R = 1.4 Bohr)
        Some(-1.1745)
    }
}

/// Heitler-London + Jastrow trial wavefunction for H₂.
///
/// # Form
///
/// ```text
/// Ψ_T = Ψ_HL × J
///
/// Ψ_HL = exp(-α(r₁ₐ + r₂ᵦ)) + exp(-α(r₁ᵦ + r₂ₐ))
/// J(r₁₂) = exp(r₁₂ / (2(1 + b·r₁₂)))
/// ```
///
/// # Physical Interpretation
///
/// - **Heitler-London part**: Symmetric spatial wavefunction representing
///   covalent bonding (electron 1 on nucleus A, electron 2 on nucleus B,
///   plus exchange term).
///
/// - **Jastrow factor**: Handles electron-electron cusp condition.
///   As r₁₂ → 0, J → exp(r₁₂/2) ensures the correct cusp behavior.
#[derive(Clone, Debug)]
pub struct H2TrialWf {
    /// Variational parameter for atomic orbitals.
    pub alpha: f64,
    /// Jastrow correlation parameter.
    pub jastrow_b: f64,
    /// Nuclear positions [A, B].
    nuclear_positions: [Vec3; 2],
}

impl H2TrialWf {
    /// Create a new H₂ trial wavefunction.
    pub fn new(alpha: f64, jastrow_b: f64, nuclear_positions: [Vec3; 2]) -> Self {
        Self {
            alpha,
            jastrow_b,
            nuclear_positions,
        }
    }

    /// Heitler-London part: Ψ_HL = exp(-α(r₁ₐ + r₂ᵦ)) + exp(-α(r₁ᵦ + r₂ₐ))
    fn heitler_london(&self, electrons: &[Vec3]) -> f64 {
        let r1a = (electrons[0] - self.nuclear_positions[0]).norm();
        let r1b = (electrons[0] - self.nuclear_positions[1]).norm();
        let r2a = (electrons[1] - self.nuclear_positions[0]).norm();
        let r2b = (electrons[1] - self.nuclear_positions[1]).norm();

        (-self.alpha * (r1a + r2b)).exp() + (-self.alpha * (r1b + r2a)).exp()
    }

    /// Jastrow factor: J(r₁₂) = exp(r₁₂ / (2(1 + b·r₁₂)))
    fn jastrow(&self, r12: f64) -> f64 {
        (r12 / (2.0 * (1.0 + self.jastrow_b * r12))).exp()
    }

    /// du/dr₁₂ = 1 / (2(1 + b·r₁₂)²)
    fn jastrow_du(&self, r12: f64) -> f64 {
        let denom = 1.0 + self.jastrow_b * r12;
        1.0 / (2.0 * denom * denom)
    }

    /// d²u/dr₁₂² = -2b / (2(1 + b·r₁₂)³)
    fn jastrow_d2u(&self, r12: f64) -> f64 {
        let denom = 1.0 + self.jastrow_b * r12;
        -self.jastrow_b / (denom * denom * denom)
    }
}

impl TrialWavefunction for H2TrialWf {
    /// Ψ_T = Ψ_HL × J
    fn value(&self, electrons: &[Vec3]) -> f64 {
        let r12 = (electrons[0] - electrons[1]).norm();
        self.heitler_london(electrons) * self.jastrow(r12)
    }

    /// ∇ln Ψ_T = ∇ln Ψ_HL + ∇ln J
    ///
    /// This is computed for each electron separately.
    fn gradient_log(&self, electrons: &[Vec3]) -> Vec<Vec3> {
        let r1a_vec = electrons[0] - self.nuclear_positions[0];
        let r1b_vec = electrons[0] - self.nuclear_positions[1];
        let r2a_vec = electrons[1] - self.nuclear_positions[0];
        let r2b_vec = electrons[1] - self.nuclear_positions[1];
        let r12_vec = electrons[0] - electrons[1];

        let r1a = r1a_vec.norm().max(1e-10);
        let r1b = r1b_vec.norm().max(1e-10);
        let r2a = r2a_vec.norm().max(1e-10);
        let r2b = r2b_vec.norm().max(1e-10);
        let r12 = r12_vec.norm().max(1e-10);

        // Unit vectors
        let r1a_hat = r1a_vec / r1a;
        let r1b_hat = r1b_vec / r1b;
        let r2a_hat = r2a_vec / r2a;
        let r2b_hat = r2b_vec / r2b;
        let r12_hat = r12_vec / r12;

        // Heitler-London contributions
        let phi1 = (-self.alpha * (r1a + r2b)).exp(); // electron 1 on A, 2 on B
        let phi2 = (-self.alpha * (r1b + r2a)).exp(); // electron 1 on B, 2 on A
        let psi_hl = phi1 + phi2;

        if psi_hl.abs() < 1e-10 {
            return vec![Vec3::zeros(), Vec3::zeros()];
        }

        // ∇₁ln Ψ_HL = (∇₁Ψ_HL) / Ψ_HL
        // ∇₁φ₁ = -α r̂₁ₐ φ₁, ∇₁φ₂ = -α r̂₁ᵦ φ₂
        let grad1_hl = (-self.alpha * r1a_hat * phi1 - self.alpha * r1b_hat * phi2) / psi_hl;

        // ∇₂ln Ψ_HL
        // ∇₂φ₁ = -α r̂₂ᵦ φ₁, ∇₂φ₂ = -α r̂₂ₐ φ₂
        let grad2_hl = (-self.alpha * r2b_hat * phi1 - self.alpha * r2a_hat * phi2) / psi_hl;

        // Jastrow contribution: ∇ᵢln J = ∇ᵢu
        // u = r₁₂ / (2(1 + b·r₁₂))
        // ∇₁u = (du/dr₁₂) × ∇₁r₁₂ = (du/dr₁₂) × r̂₁₂
        // ∇₂u = (du/dr₁₂) × ∇₂r₁₂ = (du/dr₁₂) × (-r̂₁₂)
        let du = self.jastrow_du(r12);
        let grad1_j = du * r12_hat;
        let grad2_j = -du * r12_hat;

        vec![grad1_hl + grad1_j, grad2_hl + grad2_j]
    }

    /// ∇²Ψ_T/Ψ_T for H₂.
    ///
    /// This is complex due to the product form Ψ_T = Ψ_HL × J.
    ///
    /// Using ∇²(fg)/fg = ∇²f/f + ∇²g/g + 2(∇f/f)·(∇g/g)
    fn laplacian_ratio(&self, electrons: &[Vec3]) -> f64 {
        let r1a_vec = electrons[0] - self.nuclear_positions[0];
        let r1b_vec = electrons[0] - self.nuclear_positions[1];
        let r2a_vec = electrons[1] - self.nuclear_positions[0];
        let r2b_vec = electrons[1] - self.nuclear_positions[1];
        let r12_vec = electrons[0] - electrons[1];

        let r1a = r1a_vec.norm().max(1e-10);
        let r1b = r1b_vec.norm().max(1e-10);
        let r2a = r2a_vec.norm().max(1e-10);
        let r2b = r2b_vec.norm().max(1e-10);
        let r12 = r12_vec.norm().max(1e-10);

        let r1a_hat = r1a_vec / r1a;
        let r1b_hat = r1b_vec / r1b;
        let r2a_hat = r2a_vec / r2a;
        let r2b_hat = r2b_vec / r2b;
        let r12_hat = r12_vec / r12;

        // Heitler-London terms
        let phi1 = (-self.alpha * (r1a + r2b)).exp();
        let phi2 = (-self.alpha * (r1b + r2a)).exp();
        let psi_hl = phi1 + phi2;

        if psi_hl.abs() < 1e-10 {
            return self.alpha * self.alpha * 4.0; // Approximate
        }

        // ∇²φ₁/φ₁ for electron 1: α² - 2α/r₁ₐ
        // ∇²φ₁/φ₁ for electron 2: α² - 2α/r₂ᵦ
        let lap1_phi1_over_phi1 = self.alpha * self.alpha - 2.0 * self.alpha / r1a;
        let lap2_phi1_over_phi1 = self.alpha * self.alpha - 2.0 * self.alpha / r2b;

        // ∇²φ₂/φ₂ for electron 1: α² - 2α/r₁ᵦ
        // ∇²φ₂/φ₂ for electron 2: α² - 2α/r₂ₐ
        let lap1_phi2_over_phi2 = self.alpha * self.alpha - 2.0 * self.alpha / r1b;
        let lap2_phi2_over_phi2 = self.alpha * self.alpha - 2.0 * self.alpha / r2a;

        // ∇²Ψ_HL/Ψ_HL = (∇²φ₁ + ∇²φ₂) / (φ₁ + φ₂)
        let lap1_hl = (lap1_phi1_over_phi1 * phi1 + lap1_phi2_over_phi2 * phi2) / psi_hl;
        let lap2_hl = (lap2_phi1_over_phi1 * phi1 + lap2_phi2_over_phi2 * phi2) / psi_hl;

        // Jastrow Laplacian contribution
        // ∇²ln J = ∇²u = d²u/dr₁₂² + (2/r₁₂)(du/dr₁₂)
        // For electron 1: adds this term
        // For electron 2: adds this term (same contribution due to |r₁₂|)
        let du = self.jastrow_du(r12);
        let d2u = self.jastrow_d2u(r12);
        let lap_j_per_electron = d2u + 2.0 * du / r12;

        // Cross term: 2(∇HL/HL)·(∇J/J)
        // For electron 1: 2 * grad1_hl · (du * r̂₁₂)
        let grad1_hl = (-self.alpha * r1a_hat * phi1 - self.alpha * r1b_hat * phi2) / psi_hl;
        let grad2_hl = (-self.alpha * r2b_hat * phi1 - self.alpha * r2a_hat * phi2) / psi_hl;

        let cross1 = 2.0 * grad1_hl.dot(&(du * r12_hat));
        let cross2 = 2.0 * grad2_hl.dot(&(-du * r12_hat));

        // Total Laplacian ratio
        lap1_hl + lap2_hl + 2.0 * lap_j_per_electron + cross1 + cross2
    }

    fn parameters(&self) -> Vec<f64> {
        vec![self.alpha, self.jastrow_b]
    }

    fn set_parameters(&mut self, params: &[f64]) {
        if params.len() >= 2 {
            self.alpha = params[0];
            self.jastrow_b = params[1];
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn h2_potential_at_equilibrium() {
        let h2 = H2Molecule::new(1.4);
        // Two electrons at nuclei positions
        let electrons = vec![
            Vec3::new(-0.7, 0.0, 0.0), // On nucleus A
            Vec3::new(0.7, 0.0, 0.0),  // On nucleus B
        ];
        // This should give a finite (negative) potential
        let v = h2.potential(&electrons);
        assert!(v < 0.0);
    }

    #[test]
    fn h2_potential_electron_symmetry() {
        let h2 = H2Molecule::new(1.4);
        let e1 = vec![
            Vec3::new(0.0, 0.5, 0.0),
            Vec3::new(0.0, -0.5, 0.0),
        ];
        let e2 = vec![
            Vec3::new(0.0, -0.5, 0.0),
            Vec3::new(0.0, 0.5, 0.0),
        ];
        // Swapping electrons should give same potential
        assert_relative_eq!(h2.potential(&e1), h2.potential(&e2), epsilon = 1e-10);
    }

    #[test]
    fn h2_trial_wf_symmetry() {
        let h2 = H2Molecule::new(1.4).with_trial_wavefunction(1.0, 0.5);
        let wf = h2.trial_wavefunction().unwrap();

        let e1 = vec![
            Vec3::new(0.0, 0.5, 0.0),
            Vec3::new(0.0, -0.5, 0.0),
        ];
        let e2 = vec![
            Vec3::new(0.0, -0.5, 0.0),
            Vec3::new(0.0, 0.5, 0.0),
        ];
        // Symmetric wavefunction
        assert_relative_eq!(wf.value(&e1), wf.value(&e2), epsilon = 1e-10);
    }
}
