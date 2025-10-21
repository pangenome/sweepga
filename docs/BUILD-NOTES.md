# Build Notes for Debian 10 + Guix Mixed Environment

## Issue

This system has both Debian 10 (glibc 2.28) and Guix (glibc 2.35) installed. When building Rust projects that use C dependencies (via bindgen/cc), the mixed glibc versions cause linker errors:

```
ld: /usr/lib/x86_64-linux-gnu/librt.so: undefined reference to `__pthread_barrier_wait@GLIBC_PRIVATE'
ld: /usr/lib/x86_64-linux-gnu/librt.so: undefined reference to `__libc_dlsym@GLIBC_PRIVATE'
...
```

This happens because:
1. The system's `librt.so.1` (from Debian glibc 2.28) depends on GLIBC_PRIVATE symbols
2. Guix's libraries in PATH/LD_LIBRARY_PATH provide glibc 2.35 which has different GLIBC_PRIVATE symbols
3. The symbols are incompatible → linker failure

## Solution

Build with a pristine environment that excludes Guix paths:

### Method 1: Use the build script (recommended)

```bash
# Build only
./scripts/build-clean.sh

# Build and install to ~/.cargo/bin
./scripts/build-clean.sh --install
```

### Method 2: Manual build command

```bash
env -i \
  PATH="/home/erikg/.cargo/bin:/usr/local/bin:/usr/bin:/bin" \
  HOME="/home/erikg" \
  USER="$USER" \
  RUSTUP_HOME="/export/local/home/erikg/.rustup" \
  CARGO_HOME="/home/erikg/.cargo" \
  CC="/usr/bin/gcc" \
  CXX="/usr/bin/g++" \
  BINDGEN_EXTRA_CLANG_ARGS="-I/usr/lib/gcc/x86_64-linux-gnu/8/include -I/usr/local/include -I/usr/lib/gcc/x86_64-linux-gnu/8/include-fixed -I/usr/include/x86_64-linux-gnu -I/usr/include" \
  cargo build --release
```

## Key Points

1. **Remove Guix from environment**: `env -i` clears the environment, then we selectively add only system paths
2. **Force system GCC**: `CC=/usr/bin/gcc` ensures we use Debian's GCC 8, not Guix's GCC 11
3. **Provide system includes**: `BINDGEN_EXTRA_CLANG_ARGS` tells bindgen where to find system headers
4. **Keep cargo/rustup paths**: Still need access to Rust toolchain

## Why Not Fix Guix/Debian Conflict?

The conflict is inherent when mixing package managers with different glibc versions. Options:
- Use only Debian packages (remove Guix)
- Use only Guix packages (use Guix's rust/cargo)
- Build in isolated environment (this solution)

We chose the third option as it's non-invasive and works reliably.

## Related Issues

- onecode-rs and fastga-rs both use bindgen → both affected
- Any Rust crate with C dependencies may hit this issue
- Standard Rust-only crates build fine without this workaround

## Changes Made

1. **sweepga src/main.rs**: Replaced unstable `is_multiple_of()` with stable `% 2 == 0`
2. **sweepga src/bin/alnstats.rs**: Same fix for is_multiple_of()
3. **Cargo.toml**: Changed onecode dependency to local path (temporary for testing)
   - Revert to git dependency before committing

## Testing

```bash
./build-clean.sh
./target/release/sweepga --version  # Should print: sweepga 0.1.0
```
