//! Statistical analysis functions for DMC simulations.
//!
//! Provides block averaging and error estimation.

use crate::types::BlockData;

/// Calculate block-averaged statistics from energy trace.
///
/// # Arguments
///
/// * `energies` - Vector of energy samples
/// * `populations` - Vector of population sizes
/// * `block_size` - Number of samples per block
///
/// # Returns
///
/// Vector of block statistics.
///
/// # Example
///
/// ```rust
/// use dmc_output::statistics::block_average;
///
/// let energies = vec![-0.5; 1000];
/// let populations = vec![100; 1000];
/// let blocks = block_average(&energies, &populations, 100);
///
/// assert_eq!(blocks.len(), 10);
/// ```
pub fn block_average(
    energies: &[f64],
    populations: &[usize],
    block_size: usize,
) -> Vec<BlockData> {
    if block_size == 0 {
        return Vec::new();
    }

    let num_blocks = energies.len() / block_size;
    let mut blocks = Vec::with_capacity(num_blocks);

    for i in 0..num_blocks {
        let start = i * block_size;
        let end = start + block_size;

        let block_energies = &energies[start..end];
        let block_pops = &populations[start..end];

        let mean_energy = block_energies.iter().sum::<f64>() / block_size as f64;
        let mean_pop = block_pops.iter().sum::<usize>() as f64 / block_size as f64;

        // Standard error within block
        let variance: f64 = block_energies
            .iter()
            .map(|e| (e - mean_energy).powi(2))
            .sum::<f64>()
            / (block_size - 1).max(1) as f64;
        let std_error = (variance / block_size as f64).sqrt();

        blocks.push(BlockData {
            block_index: i as u64,
            start_step: start as u64,
            end_step: end as u64,
            mean_energy,
            std_error,
            mean_population: mean_pop,
        });
    }

    blocks
}

/// Calculate final energy estimate and error from block data.
///
/// # Returns
///
/// (mean_energy, standard_error)
///
/// # Example
///
/// ```rust
/// use dmc_output::statistics::{block_average, final_statistics};
///
/// let energies = vec![-0.5; 1000];
/// let populations = vec![100; 1000];
/// let blocks = block_average(&energies, &populations, 100);
/// let (mean, error) = final_statistics(&blocks);
///
/// assert!((mean - (-0.5)).abs() < 0.01);
/// ```
pub fn final_statistics(blocks: &[BlockData]) -> (f64, f64) {
    if blocks.is_empty() {
        return (0.0, 0.0);
    }

    let n = blocks.len() as f64;
    let mean: f64 = blocks.iter().map(|b| b.mean_energy).sum::<f64>() / n;

    if blocks.len() == 1 {
        return (mean, blocks[0].std_error);
    }

    // Standard error of block means
    let variance: f64 = blocks
        .iter()
        .map(|b| (b.mean_energy - mean).powi(2))
        .sum::<f64>()
        / (n - 1.0);
    let std_error = (variance / n).sqrt();

    (mean, std_error)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn block_average_simple() {
        let energies = vec![-0.5; 100];
        let populations = vec![100; 100];
        let blocks = block_average(&energies, &populations, 10);

        assert_eq!(blocks.len(), 10);
        for block in &blocks {
            assert!((block.mean_energy - (-0.5)).abs() < 1e-10);
            assert!(block.std_error < 1e-10);
        }
    }

    #[test]
    fn final_statistics_simple() {
        let energies = vec![-0.5; 100];
        let populations = vec![100; 100];
        let blocks = block_average(&energies, &populations, 10);
        let (mean, error) = final_statistics(&blocks);

        assert!((mean - (-0.5)).abs() < 1e-10);
        assert!(error < 1e-10);
    }

    #[test]
    fn empty_blocks() {
        let (mean, error) = final_statistics(&[]);
        assert_eq!(mean, 0.0);
        assert_eq!(error, 0.0);
    }
}
