#!/usr/bin/env bash
set -e

# Install script for sweepga with embedded FastGA binaries
# This ensures FastGA helper tools are installed alongside sweepga

echo "=== Installing sweepga with embedded FastGA binaries ==="

# Determine CARGO_HOME
CARGO_HOME="${CARGO_HOME:-$HOME/.cargo}"
echo "CARGO_HOME: $CARGO_HOME"

# Build the project first to ensure FastGA binaries are compiled
echo ""
echo "Building sweepga..."
cargo build --release

# Find the FastGA binaries
echo ""
echo "Locating FastGA binaries..."
FASTGA_BIN=$(find target/release/build -type f -name "FastGA" | head -1)
ALNTOPAF_BIN=$(find target/release/build -type f -name "ALNtoPAF" | head -1)
PAFTOALN_BIN=$(find target/release/build -type f -name "PAFtoALN" | head -1)

if [ -z "$FASTGA_BIN" ] || [ -z "$ALNTOPAF_BIN" ] || [ -z "$PAFTOALN_BIN" ]; then
    echo "ERROR: Could not find all FastGA binaries after build"
    echo "  FastGA: ${FASTGA_BIN:-NOT FOUND}"
    echo "  ALNtoPAF: ${ALNTOPAF_BIN:-NOT FOUND}"
    echo "  PAFtoALN: ${PAFTOALN_BIN:-NOT FOUND}"
    exit 1
fi

echo "  FastGA: $FASTGA_BIN"
echo "  ALNtoPAF: $ALNTOPAF_BIN"
echo "  PAFtoALN: $PAFTOALN_BIN"

# Install main binaries using cargo install
echo ""
echo "Installing sweepga binaries to $CARGO_HOME/bin/..."
cargo install --path . --force

# Create lib directory for FastGA tools
LIB_DIR="$CARGO_HOME/lib/sweepga"
echo ""
echo "Creating library directory: $LIB_DIR"
mkdir -p "$LIB_DIR"

# Copy FastGA binaries
echo ""
echo "Installing FastGA helper binaries..."
cp "$FASTGA_BIN" "$LIB_DIR/FastGA"
cp "$ALNTOPAF_BIN" "$LIB_DIR/ALNtoPAF"
cp "$PAFTOALN_BIN" "$LIB_DIR/PAFtoALN"

# Make them executable
chmod +x "$LIB_DIR/FastGA"
chmod +x "$LIB_DIR/ALNtoPAF"
chmod +x "$LIB_DIR/PAFtoALN"

echo ""
echo "=== Installation complete! ==="
echo ""
echo "Installed binaries:"
echo "  sweepga:  $CARGO_HOME/bin/sweepga"
echo "  alnstats: $CARGO_HOME/bin/alnstats"
echo ""
echo "FastGA helper binaries (not in PATH):"
echo "  FastGA:   $LIB_DIR/FastGA"
echo "  ALNtoPAF: $LIB_DIR/ALNtoPAF"
echo "  PAFtoALN: $LIB_DIR/PAFtoALN"
echo ""
echo "Test the installation:"
echo "  sweepga --help"
echo "  sweepga data/scerevisiae8.fa.gz -b 5000 | head"
