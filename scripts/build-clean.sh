#!/bin/bash
# Build sweepga in a pristine environment without Guix contamination
#
# This script removes Guix paths from the environment to avoid glibc version
# conflicts (Debian 10 system glibc 2.28 vs Guix glibc 2.35) which cause
# "undefined reference to __pthread_barrier_wait@GLIBC_PRIVATE" linker errors.

set -e

echo "Building sweepga with clean system environment..."

# Build with minimal clean environment
env -i \
  PATH="/home/erikg/.cargo/bin:/usr/local/bin:/usr/bin:/bin" \
  HOME="/home/erikg" \
  USER="${USER:-erikg}" \
  RUSTUP_HOME="${RUSTUP_HOME:-/export/local/home/erikg/.rustup}" \
  CARGO_HOME="${CARGO_HOME:-/home/erikg/.cargo}" \
  CC="/usr/bin/gcc" \
  CXX="/usr/bin/g++" \
  BINDGEN_EXTRA_CLANG_ARGS="-I/usr/lib/gcc/x86_64-linux-gnu/8/include -I/usr/local/include -I/usr/lib/gcc/x86_64-linux-gnu/8/include-fixed -I/usr/include/x86_64-linux-gnu -I/usr/include" \
  bash -c "cargo build --release"

echo ""
echo "Build successful!"
echo "Binary location: ./target/release/sweepga"
