#!/bin/bash
# Build sweepga in a pristine environment without Guix contamination
#
# This script removes Guix paths from the environment to avoid glibc version
# conflicts (Debian 10 system glibc 2.28 vs Guix glibc 2.35) which cause
# "undefined reference to __pthread_barrier_wait@GLIBC_PRIVATE" linker errors.
#
# Usage: ./scripts/build-clean.sh [--install]
#   --install: Also install to ~/.cargo/bin after building

set -e

# Parse arguments
INSTALL=false
if [[ "$1" == "--install" ]]; then
    INSTALL=true
elif [[ "$1" == "--help" || "$1" == "-h" ]]; then
    echo "Usage: $0 [--install]"
    echo ""
    echo "Build sweepga in a clean environment (without Guix paths)"
    echo ""
    echo "Options:"
    echo "  --install    Also install to ~/.cargo/bin after building"
    echo "  --help, -h   Show this help message"
    exit 0
fi

# Detect paths (make more portable)
HOME_DIR="${HOME:-/home/erikg}"
CARGO_BIN="${HOME_DIR}/.cargo/bin"
RUSTUP_HOME="${RUSTUP_HOME:-${HOME_DIR}/.rustup}"
CARGO_HOME="${CARGO_HOME:-${HOME_DIR}/.cargo}"

# Detect system GCC include paths
GCC_VERSION=$(ls /usr/lib/gcc/x86_64-linux-gnu/ 2>/dev/null | sort -V | tail -1)
GCC_INCLUDE="/usr/lib/gcc/x86_64-linux-gnu/${GCC_VERSION}/include"

echo "Building sweepga with clean system environment..."

# Build command
CARGO_CMD="cargo build --release"
if [[ "$INSTALL" == true ]]; then
    CARGO_CMD="cargo install --path . --force"
fi

# Build with minimal clean environment
env -i \
  PATH="${CARGO_BIN}:/usr/local/bin:/usr/bin:/bin" \
  HOME="${HOME_DIR}" \
  USER="${USER}" \
  RUSTUP_HOME="${RUSTUP_HOME}" \
  CARGO_HOME="${CARGO_HOME}" \
  CC="/usr/bin/gcc" \
  CXX="/usr/bin/g++" \
  BINDGEN_EXTRA_CLANG_ARGS="-I${GCC_INCLUDE} -I/usr/local/include -I${GCC_INCLUDE}-fixed -I/usr/include/x86_64-linux-gnu -I/usr/include" \
  bash -c "${CARGO_CMD}"

echo ""
if [[ "$INSTALL" == true ]]; then
    echo "Installation successful!"
    echo "Installed: ${CARGO_BIN}/sweepga"
    echo "           ${CARGO_BIN}/alnstats"
else
    echo "Build successful!"
    echo "Binary location: ./target/release/sweepga"
    echo ""
    echo "To install, run: $0 --install"
fi
