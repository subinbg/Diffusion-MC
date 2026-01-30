//! Quantum system abstractions for DMC simulations.
//!
//! This module defines the core traits [`QuantumSystem`] and [`TrialWavefunction`]
//! that abstract over different quantum systems (H, H₂⁺, H₂).
//!
//! # Mathematical Background
//!
//! The DMC method solves the imaginary-time Schrödinger equation:
//!
//! ```text
//! ∂Ψ/∂τ = -(1/ℏ)(Ĥ - E_T)Ψ
//! ```
//!
//! See [README: Theoretical Background](https://github.com/subinbg/Diffusion-MC#theoretical-background)

mod hydrogen;
mod h2_ion;
mod h2_molecule;

pub use hydrogen::{Hydrogen, HydrogenTrialWf};
pub use h2_ion::{H2Ion, H2IonTrialWf};
pub use h2_molecule::{H2Molecule, H2TrialWf};

use crate::Vec3;

/// Nuclear configuration: fixed positions and charges.
///
/// In the Born-Oppenheimer approximation, nuclei are treated as fixed point charges.
#[derive(Clone, Debug)]
pub struct NuclearConfig {
    /// Nuclear positions in atomic units (Bohr radii).
    pub positions: Vec<Vec3>,
    /// Nuclear charges in units of e (proton charge = 1.0).
    pub charges: Vec<f64>,
}

/// Trait for quantum systems that can be simulated with DMC.
///
/// # Equation Reference
///
/// The potential energy in position representation:
///
/// ```text
/// V(x) = -Σ_{i,I} Z_I/|r_i - R_I| + Σ_{i<j} 1/|r_i - r_j| + Σ_{I<J} Z_I Z_J/|R_I - R_J|
/// ```
///
/// where:
/// - First term: electron-nuclear attraction
/// - Second term: electron-electron repulsion
/// - Third term: nuclear-nuclear repulsion (constant for fixed nuclei)
///
/// See [README: Step 3](https://github.com/subinbg/Diffusion-MC#step-3-weighting-and-branching-potential-update)
pub trait QuantumSystem: Send + Sync {
    /// Number of electrons in the system.
    fn num_electrons(&self) -> usize;

    /// Nuclear configuration (positions and charges).
    fn nuclear_config(&self) -> &NuclearConfig;

    /// Calculate potential energy V(x) at given electron positions.
    ///
    /// # Arguments
    ///
    /// * `electrons` - Slice of electron position vectors
    ///
    /// # Returns
    ///
    /// Potential energy in Hartree atomic units.
    fn potential(&self, electrons: &[Vec3]) -> f64;

    /// Get the associated trial wavefunction (for importance sampling).
    ///
    /// Returns `None` for pure DMC without importance sampling.
    fn trial_wavefunction(&self) -> Option<&dyn TrialWavefunction>;

    /// Exact ground state energy, if known analytically.
    ///
    /// Used for validation and error reporting.
    fn exact_energy(&self) -> Option<f64>;
}

/// Trial wavefunction for importance-sampled DMC.
///
/// # Theory
///
/// The trial wavefunction Ψ_T(x) guides walkers and defines:
///
/// - **Drift velocity**: `v_D = (ℏ/m) ∇ln|Ψ_T|`
/// - **Local energy**: `E_L = ĤΨ_T/Ψ_T`
///
/// # Cusp Condition
///
/// For numerical stability at Coulomb singularities, Ψ_T must satisfy the
/// Kato cusp condition. For a hydrogen-like atom with nuclear charge Z:
///
/// ```text
/// Ψ_T(r) = exp(-αr)  with α = Z
/// ```
///
/// This ensures the singularity in the kinetic term cancels the potential singularity,
/// making E_L(x) bounded everywhere.
///
/// See [README: The Cusp Condition](https://github.com/subinbg/Diffusion-MC#the-cusp-condition-stability)
pub trait TrialWavefunction: Send + Sync {
    /// Evaluate Ψ_T(x) at given electron positions.
    fn value(&self, electrons: &[Vec3]) -> f64;

    /// Evaluate ln|Ψ_T(x)|.
    ///
    /// Often more numerically stable than computing value and taking log.
    fn log_value(&self, electrons: &[Vec3]) -> f64 {
        self.value(electrons).abs().ln()
    }

    /// Gradient of ln|Ψ_T| for each electron.
    ///
    /// # Equation
    ///
    /// Returns ∇_i ln|Ψ_T(x)| for each electron i.
    ///
    /// The drift velocity is then:
    /// ```text
    /// v_{D,i} = (ℏ/m) ∇_i ln|Ψ_T(x)|
    /// ```
    ///
    /// In atomic units (ℏ = m = 1), v_D = ∇ln|Ψ_T|.
    ///
    /// See [README: The Drift-Diffusion Green's Function](https://github.com/subinbg/Diffusion-MC#the-drift-diffusion-greens-function)
    fn gradient_log(&self, electrons: &[Vec3]) -> Vec<Vec3>;

    /// Laplacian of Ψ_T divided by Ψ_T.
    ///
    /// # Equation
    ///
    /// ```text
    /// ∇²Ψ_T/Ψ_T = Σ_i [(∇_i ln Ψ_T)² + ∇_i² ln Ψ_T]
    /// ```
    ///
    /// This is used in the local energy calculation:
    /// ```text
    /// E_L(x) = -(ℏ²/2m)(∇²Ψ_T/Ψ_T) + V(x)
    /// ```
    ///
    /// See [README: The Local Energy](https://github.com/subinbg/Diffusion-MC#the-local-energy)
    fn laplacian_ratio(&self, electrons: &[Vec3]) -> f64;

    /// Variational parameters of the trial wavefunction.
    fn parameters(&self) -> Vec<f64>;

    /// Update variational parameters.
    fn set_parameters(&mut self, params: &[f64]);
}
