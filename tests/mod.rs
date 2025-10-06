#![allow(clippy::uninlined_format_args)]

/// Test module declarations
mod fastga_integration;
mod synthetic_genomes;
mod test_utils;
mod unit_tests;

// Re-export for convenience
pub use synthetic_genomes::*;
pub use test_utils::*;
