//! 2D electron density computation for visualization.

use dmc_core::walker::Population;

/// Plane for 2D density projection.
#[derive(Clone, Copy, Debug)]
pub enum Plane {
    /// XY plane (z = slice_position)
    XY,
    /// XZ plane (y = slice_position)
    XZ,
    /// YZ plane (x = slice_position)
    YZ,
}

/// Compute 2D electron density on a plane using Gaussian smearing.
///
/// # Arguments
///
/// * `population` - Walker population
/// * `plane` - Projection plane (XY, XZ, YZ)
/// * `slice_position` - Position of the slice along the perpendicular axis
/// * `grid_bounds` - (min, max) bounds for the grid
/// * `grid_size` - Number of grid points in each dimension
/// * `sigma` - Gaussian smearing width
///
/// # Returns
///
/// Flattened array of density values (row-major order).
pub fn compute_density_2d(
    population: &Population,
    plane: Plane,
    slice_position: f64,
    grid_bounds: (f64, f64),
    grid_size: usize,
    sigma: f64,
) -> Vec<f64> {
    let (min_bound, max_bound) = grid_bounds;
    let grid_spacing = (max_bound - min_bound) / grid_size as f64;
    let sigma2 = sigma * sigma;
    let norm = 1.0 / (2.0 * std::f64::consts::PI * sigma2);

    // Initialize density grid
    let mut density = vec![0.0; grid_size * grid_size];

    // For each walker
    for walker in population.iter() {
        // For each electron in the walker
        for pos in &walker.positions {
            // Extract coordinates based on plane
            let (u, v, w) = match plane {
                Plane::XY => (pos.x, pos.y, pos.z),
                Plane::XZ => (pos.x, pos.z, pos.y),
                Plane::YZ => (pos.y, pos.z, pos.x),
            };

            // Skip if electron is far from the slice
            let dw = (w - slice_position).abs();
            if dw > 3.0 * sigma {
                continue;
            }

            // Gaussian weight based on distance from slice
            let slice_weight = (-0.5 * dw * dw / sigma2).exp();

            // Find grid cells that could be affected (within 3 sigma)
            let i_min = ((u - 3.0 * sigma - min_bound) / grid_spacing)
                .floor()
                .max(0.0) as usize;
            let i_max = ((u + 3.0 * sigma - min_bound) / grid_spacing)
                .ceil()
                .min(grid_size as f64) as usize;
            let j_min = ((v - 3.0 * sigma - min_bound) / grid_spacing)
                .floor()
                .max(0.0) as usize;
            let j_max = ((v + 3.0 * sigma - min_bound) / grid_spacing)
                .ceil()
                .min(grid_size as f64) as usize;

            // Add Gaussian contribution to nearby grid points
            for i in i_min..i_max {
                let grid_u = min_bound + (i as f64 + 0.5) * grid_spacing;
                let du = grid_u - u;

                for j in j_min..j_max {
                    let grid_v = min_bound + (j as f64 + 0.5) * grid_spacing;
                    let dv = grid_v - v;

                    let r2 = du * du + dv * dv;
                    let weight = slice_weight * norm * (-0.5 * r2 / sigma2).exp();

                    density[j * grid_size + i] += weight;
                }
            }
        }
    }

    // Normalize by number of walkers
    let n_walkers = population.size() as f64;
    if n_walkers > 0.0 {
        for d in &mut density {
            *d /= n_walkers;
        }
    }

    density
}

#[cfg(test)]
mod tests {
    use super::*;
    use dmc_core::walker::Population;

    #[test]
    fn density_non_negative() {
        let mut pop = Population::new(100, 42);
        pop.initialize_gaussian(1, 1.0);

        let density = compute_density_2d(&pop, Plane::XY, 0.0, (-5.0, 5.0), 50, 0.5);

        assert_eq!(density.len(), 50 * 50);
        for d in &density {
            assert!(*d >= 0.0);
        }
    }

    #[test]
    fn density_has_maximum_near_center() {
        let mut pop = Population::new(1000, 42);
        pop.initialize_gaussian(1, 0.5); // Small spread

        let density = compute_density_2d(&pop, Plane::XY, 0.0, (-5.0, 5.0), 50, 0.5);

        // Find maximum
        let _center_idx = 25 * 50 + 25; // Center of grid
        let center_region: f64 = density[24 * 50 + 24..27 * 50 + 27]
            .iter()
            .sum::<f64>();
        let edge_region: f64 = density[0..3 * 50].iter().sum::<f64>();

        assert!(center_region > edge_region);
    }
}
