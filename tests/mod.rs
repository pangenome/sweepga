#![allow(clippy::uninlined_format_args)]

/// Test module declarations
mod fastga_integration;
// mod synthetic_genomes; // Loaded in fastga_integration.rs
mod test_utils;
mod unit_tests;

// Re-export for convenience
// pub use synthetic_genomes::*; // Not available at crate root
pub use test_utils::*;
