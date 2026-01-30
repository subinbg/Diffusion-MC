//! Walker population management for DMC.
//!
//! In DMC, the wavefunction is represented by a population of "walkers"
//! (random walkers in configuration space) that undergo diffusion and branching.
//!
//! # Mathematical Interpretation
//!
//! The walker density ρ(x, τ) approximates the wavefunction Ψ(x, τ):
//!
//! ```text
//! Ψ(x, τ) ≈ (1/N) Σᵢ δ(x - xᵢ(τ))
//! ```
//!
//! As τ → ∞, this converges to the ground state Φ₀(x).

use crate::Vec3;
use rand::prelude::*;
use rand_distr::StandardNormal;

/// State of a single walker (replica/pseudoparticle).
#[derive(Clone, Debug)]
pub struct Walker {
    /// Electron positions (one Vec3 per electron).
    pub positions: Vec<Vec3>,
    /// Weight (for weighted averages; usually 1.0 with integer branching).
    pub weight: f64,
    /// Age in simulation steps (for diagnostics).
    pub age: u64,
}

impl Walker {
    /// Create a new walker with given electron positions.
    pub fn new(positions: Vec<Vec3>) -> Self {
        Self {
            positions,
            weight: 1.0,
            age: 0,
        }
    }

    /// Create a walker with all electrons at the origin.
    pub fn at_origin(num_electrons: usize) -> Self {
        Self::new(vec![Vec3::zeros(); num_electrons])
    }

    /// Create a copy (clone) for branching.
    pub fn spawn(&self) -> Self {
        Self {
            positions: self.positions.clone(),
            weight: 1.0,
            age: 0,
        }
    }
}

/// Population of walkers with birth/death dynamics.
///
/// # Branching
///
/// After each step, walkers are replicated or killed based on their weight W:
///
/// ```text
/// M = ⌊W + u⌋
/// ```
///
/// where u ~ Uniform(0,1). This preserves ⟨M⟩ = W on average.
///
/// See [README: Stochastic Branching](https://github.com/subinbg/Diffusion-MC#step-3-weighting-and-branching-potential-update)
pub struct Population {
    walkers: Vec<Walker>,
    target_size: usize,
    rng: StdRng,
}

impl Population {
    /// Create a new population with given target size and random seed.
    pub fn new(target_size: usize, seed: u64) -> Self {
        Self {
            walkers: Vec::with_capacity(target_size * 2),
            target_size,
            rng: StdRng::seed_from_u64(seed),
        }
    }

    /// Initialize walkers using a sampler function.
    ///
    /// The sampler is called `target_size` times to generate initial positions.
    pub fn initialize<F>(&mut self, mut sampler: F)
    where
        F: FnMut(&mut StdRng) -> Vec<Vec3>,
    {
        self.walkers.clear();
        for _ in 0..self.target_size {
            let positions = sampler(&mut self.rng);
            self.walkers.push(Walker::new(positions));
        }
    }

    /// Initialize all walkers at the origin.
    pub fn initialize_at_origin(&mut self, num_electrons: usize) {
        self.walkers.clear();
        for _ in 0..self.target_size {
            self.walkers.push(Walker::at_origin(num_electrons));
        }
    }

    /// Initialize walkers with Gaussian-distributed positions.
    ///
    /// Each electron position component is sampled from N(0, σ²).
    pub fn initialize_gaussian(&mut self, num_electrons: usize, sigma: f64) {
        self.initialize(|rng| {
            (0..num_electrons)
                .map(|_| {
                    Vec3::new(
                        rng.sample::<f64, _>(StandardNormal) * sigma,
                        rng.sample::<f64, _>(StandardNormal) * sigma,
                        rng.sample::<f64, _>(StandardNormal) * sigma,
                    )
                })
                .collect()
        });
    }

    /// Current number of walkers.
    pub fn size(&self) -> usize {
        self.walkers.len()
    }

    /// Target population size (N₀).
    pub fn target_size(&self) -> usize {
        self.target_size
    }

    /// Iterate over walkers.
    pub fn iter(&self) -> impl Iterator<Item = &Walker> {
        self.walkers.iter()
    }

    /// Mutable iteration over walkers.
    pub fn iter_mut(&mut self) -> impl Iterator<Item = &mut Walker> {
        self.walkers.iter_mut()
    }

    /// Get mutable reference to the RNG.
    pub fn rng(&mut self) -> &mut StdRng {
        &mut self.rng
    }

    /// Apply branching based on weights.
    ///
    /// # Algorithm
    ///
    /// For each walker with weight W:
    /// 1. Compute M = ⌊W + u⌋ where u ~ Uniform(0,1)
    /// 2. If M = 0: walker dies
    /// 3. If M = 1: walker survives
    /// 4. If M > 1: walker is cloned (M-1) times
    ///
    /// # Equation
    ///
    /// ```text
    /// M_i = min(⌊W_i + u⌋, max_offspring)
    /// ```
    ///
    /// See [README: Stochastic Branching](https://github.com/subinbg/Diffusion-MC#step-3-weighting-and-branching-potential-update)
    pub fn branch(&mut self, weights: &[f64], max_offspring: usize) {
        let mut new_walkers = Vec::with_capacity(self.walkers.len() * 2);

        for (walker, &weight) in self.walkers.iter().zip(weights.iter()) {
            let u: f64 = self.rng.gen();
            let offspring = ((weight + u) as usize).min(max_offspring);

            for _ in 0..offspring {
                let mut new_walker = walker.spawn();
                new_walker.age = walker.age;
                new_walkers.push(new_walker);
            }
        }

        self.walkers = new_walkers;
    }

    /// Increment age of all walkers.
    pub fn increment_ages(&mut self) {
        for walker in &mut self.walkers {
            walker.age += 1;
        }
    }

    /// Apply Gaussian diffusion to all walkers.
    ///
    /// # Equation
    ///
    /// For each electron position:
    /// ```text
    /// x' = x + √(δτ) × ξ
    /// ```
    ///
    /// where ξ ~ N(0, 1) for each component (in atomic units with ℏ = m = 1).
    ///
    /// See [README: Step 2: Diffusion](https://github.com/subinbg/Diffusion-MC#step-2-diffusion-kinetic-update)
    pub fn diffuse(&mut self, time_step: f64) {
        let sigma = time_step.sqrt();

        for walker in &mut self.walkers {
            for pos in &mut walker.positions {
                pos.x += self.rng.sample::<f64, _>(StandardNormal) * sigma;
                pos.y += self.rng.sample::<f64, _>(StandardNormal) * sigma;
                pos.z += self.rng.sample::<f64, _>(StandardNormal) * sigma;
            }
        }
    }

    /// Apply drift-diffusion to all walkers with given drift velocities.
    ///
    /// # Equation
    ///
    /// For each electron position with drift velocity v_D:
    /// ```text
    /// x' = x + v_D × δτ + √(δτ) × ξ
    /// ```
    ///
    /// See [README: Step 2: Drift-Diffusion Update](https://github.com/subinbg/Diffusion-MC#step-2-drift-diffusion-update-the-kinetic-step)
    pub fn drift_diffuse(&mut self, time_step: f64, drift_velocities: &[Vec<Vec3>]) {
        let sigma = time_step.sqrt();

        for (walker, drifts) in self.walkers.iter_mut().zip(drift_velocities.iter()) {
            for (pos, drift) in walker.positions.iter_mut().zip(drifts.iter()) {
                pos.x += drift.x * time_step + self.rng.sample::<f64, _>(StandardNormal) * sigma;
                pos.y += drift.y * time_step + self.rng.sample::<f64, _>(StandardNormal) * sigma;
                pos.z += drift.z * time_step + self.rng.sample::<f64, _>(StandardNormal) * sigma;
            }
        }
    }

    /// Get positions of all walkers (for output).
    pub fn all_positions(&self) -> Vec<&[Vec3]> {
        self.walkers
            .iter()
            .map(|w| w.positions.as_slice())
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn population_initialization() {
        let mut pop = Population::new(100, 42);
        pop.initialize_at_origin(1);
        assert_eq!(pop.size(), 100);
    }

    #[test]
    fn branching_preserves_average() {
        let mut pop = Population::new(1000, 42);
        pop.initialize_at_origin(1);

        // Weights averaging to 1.0 should preserve population size approximately
        let weights: Vec<f64> = (0..1000).map(|i| 0.5 + (i as f64) / 1000.0).collect();
        let total_weight: f64 = weights.iter().sum();

        pop.branch(&weights, 3);

        // Population should be close to total weight
        let size = pop.size() as f64;
        assert!(
            (size - total_weight).abs() < 100.0,
            "Expected ~{}, got {}",
            total_weight,
            size
        );
    }

    #[test]
    fn diffusion_changes_positions() {
        let mut pop = Population::new(10, 42);
        pop.initialize_at_origin(1);

        let initial: Vec<Vec3> = pop.iter().map(|w| w.positions[0]).collect();

        pop.diffuse(0.01);

        let final_: Vec<Vec3> = pop.iter().map(|w| w.positions[0]).collect();

        // Positions should have changed
        let any_changed = initial
            .iter()
            .zip(final_.iter())
            .any(|(a, b)| (a - b).norm() > 1e-10);
        assert!(any_changed, "Diffusion should change positions");
    }
}
