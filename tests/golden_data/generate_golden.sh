#!/bin/bash
set -euo pipefail

# Script to generate golden reference files for regression testing
# WARNING: Only run this when intentionally changing output behavior!

echo "=== Generating Golden Reference Files ==="
echo "WARNING: This will overwrite existing golden files!"
read -p "Continue? (y/N) " -n 1 -r
echo
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo "Aborted."
    exit 1
fi

cd "$(dirname "$0")"
GOLDEN_DIR="$(pwd)"
PROJECT_ROOT="$(cd ../.. && pwd)"

echo "Project root: $PROJECT_ROOT"
echo "Golden dir: $GOLDEN_DIR"

# Build release binary
echo "Building release binary..."
cd "$PROJECT_ROOT"
cargo build --release --quiet

# Create test input (deterministic - use existing test data)
echo "Creating test input..."
# Use fixed temp directory to ensure deterministic output
# (mktemp creates random paths which can affect .1aln binary format)
TEMP_DIR="/tmp/sweepga_golden_gen"
rm -rf "$TEMP_DIR"
mkdir -p "$TEMP_DIR"
trap "rm -rf $TEMP_DIR" EXIT

# Use existing test data if available
if [ -f "$PROJECT_ROOT/data/scerevisiae8.fa.gz" ]; then
    echo "Using scerevisiae8.fa.gz for golden data..."
    cp "$PROJECT_ROOT/data/scerevisiae8.fa.gz" "$TEMP_DIR/test_input.fa.gz"
    gunzip "$TEMP_DIR/test_input.fa.gz"
else
    # Fallback: generate synthetic data (deterministic)
    echo "Generating synthetic test data..."
    python3 << 'PYTHON'
import random
random.seed(42)  # Deterministic

def generate_seq(length, seed_offset=0):
    random.seed(42 + seed_offset)
    bases = ['A', 'C', 'G', 'T']
    return ''.join(random.choice(bases) for _ in range(length))

# Generate 3 related sequences (3kb each)
base = generate_seq(3000, 0)
seq1 = base
# Mutate 1%
seq2 = list(base)
for i in random.sample(range(len(seq2)), 30):
    seq2[i] = random.choice(['A', 'C', 'G', 'T'])
seq2 = ''.join(seq2)

with open('/tmp/test_input_golden.fa', 'w') as f:
    f.write(f'>seq1\n{seq1}\n')
    f.write(f'>seq2\n{seq2}\n')
PYTHON
    cp /tmp/test_input_golden.fa "$TEMP_DIR/test_input.fa"
fi

# Generate .1aln output (use binary directly to match test environment)
echo "Generating .1aln golden file..."
"$PROJECT_ROOT/target/release/sweepga" "$TEMP_DIR/test_input.fa" > "$GOLDEN_DIR/golden_output.1aln" 2>/dev/null || true

# Generate filtered .1aln (1:1 mode)
echo "Generating filtered .1aln golden file..."
if [ -f "$GOLDEN_DIR/golden_output.1aln" ]; then
    "$PROJECT_ROOT/target/release/sweepga" "$GOLDEN_DIR/golden_output.1aln" -n 1:1 > "$GOLDEN_DIR/golden_filtered_1to1.1aln" 2>/dev/null || true
fi

# Generate PAF output
echo "Generating PAF golden file..."
"$PROJECT_ROOT/target/release/sweepga" "$TEMP_DIR/test_input.fa" --paf > "$GOLDEN_DIR/golden_output.paf" 2>/dev/null || true

# Generate checksums
echo "Generating checksums..."
cd "$GOLDEN_DIR"
rm -f checksums.txt

# For .1aln files: use normalized checksums (ONEview + filter provenance + sort)
# This removes non-deterministic provenance records ('!' lines with timestamps/paths)
# and sorts output to handle non-deterministic record ordering from FastGA multithreading
for file in golden_*.1aln; do
    if [ -f "$file" ]; then
        echo "Computing normalized checksum for $file..."
        CHECKSUM=$(ONEview "$file" | grep -v '^!' | sort | sha256sum | awk '{print $1}')
        echo "$CHECKSUM  $file" >> checksums.txt
    fi
done

# For PAF files: use regular checksums (deterministic)
for file in golden_*.paf; do
    if [ -f "$file" ]; then
        sha256sum "$file" >> checksums.txt
    fi
done

echo ""
echo "=== Generated Files ==="
ls -lh golden_*
echo ""
echo "=== Checksums ==="
cat checksums.txt
echo ""
echo "✓ Golden files generated successfully!"
echo "Review the changes, then commit with: git add tests/golden_data/"
