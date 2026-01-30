//! Site generation module.
//!
//! This module handles generating the documentation and simulator HTML pages.

pub mod assets;
pub mod markdown;
pub mod templates;

pub use templates::{generate_documentation_page, generate_simulator_page};
