/// Build script for compiling ONElib as a static library.
use std::env;
use std::path::PathBuf;

fn main() {
    println!("cargo:rerun-if-changed=build.rs");
    println!("cargo:rerun-if-changed=deps/FASTGA");
    println!("cargo:rerun-if-changed=rust_onelib.c");

    let manifest_dir = PathBuf::from(env::var("CARGO_MANIFEST_DIR").unwrap());
    let fastga_dir = manifest_dir.join("deps").join("FASTGA");

    // Build ONElib and our wrapper functions
    cc::Build::new()
        .file("rust_onelib.c")
        .file(fastga_dir.join("ONElib.c"))
        .file(fastga_dir.join("GDB.c"))
        .file(fastga_dir.join("gene_core.c"))
        .file(fastga_dir.join("alncode.c"))
        .include(&fastga_dir)
        .include(&manifest_dir) // For rust_onelib.c
        .flag("-O3")
        .flag("-fno-strict-aliasing")
        .warnings(false)
        .define("_GNU_SOURCE", None)
        .compile("onelib");

    // Link required system libraries
    println!("cargo:rustc-link-lib=pthread");
    println!("cargo:rustc-link-lib=m");
    println!("cargo:rustc-link-lib=z");
}
