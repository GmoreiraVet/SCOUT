#!/bin/bash
#
# Example workflows for tick species identification analysis
# Shows both CONFIG-BASED and COMMAND-LINE approaches
#

echo "╔═══════════════════════════════════════════════════════════════════════╗"
echo "║         Tick Species Identification - Example Workflows               ║"
echo "║              (Both CONFIG and COMMAND-LINE approaches)               ║"
echo "╚═══════════════════════════════════════════════════════════════════════╝"
echo ""

show_workflow() {
    local num=$1
    local name=$2
    
    echo ""
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "WORKFLOW $num: $name"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo ""
}

# ============================================================================
# WORKFLOW 1: Basic Analysis - Standard Thresholds
# ============================================================================

show_workflow 1 "Basic Analysis (Standard Thresholds)"

echo "✓ APPROACH 1: CONFIG-BASED (Recommended for most users)"
echo ""
echo "  1. Open tick_species_metagenome_identifier.py in a text editor"
echo "  2. Find the CONFIG section (around line 30-40)"
echo "  3. Edit these settings:"
echo ""
cat << 'EOF'
CONFIG = {
    'input_fastq': 'your_sample.fastq',           # ← Change this
    'reference_fasta': 'Tick_Specifier_Ticks.fasta',
    'output_dir': './tick_analysis',
    'threads': 8,
    'min_identity': 0.85,
    'min_coverage': 0.80,
}
EOF
echo ""
echo "  4. Save and run:"
echo "     python tick_species_metagenome_identifier.py"
echo ""

echo "✓ APPROACH 2: COMMAND-LINE (For scripting/automation)"
echo ""
echo "  python tick_species_metagenome_identifier.py \\"
echo "    -r Tick_Specifier_Ticks.fasta \\"
echo "    -i your_sample.fastq \\"
echo "    -o ./tick_analysis \\"
echo "    -t 8"
echo ""

echo "✓ APPROACH 3: MIXED (Config + command-line override)"
echo ""
echo "  # Edit CONFIG in script, then override with command-line:"
echo "  python tick_species_metagenome_identifier.py \\"
echo "    --threads 16 \\"
echo "    --min-identity 0.90"
echo ""

# ============================================================================
# WORKFLOW 2: High-Confidence Mode (Strict QC)
# ============================================================================

show_workflow 2 "High-Confidence Analysis (Strict QC)"

echo "✓ APPROACH 1: CONFIG-BASED"
echo ""
echo "  In the CONFIG section, uncomment or modify:"
echo ""
cat << 'EOF'
CONFIG = {
    'input_fastq': 'your_sample.fastq',
    'reference_fasta': 'Tick_Specifier_Ticks.fasta',
    'min_identity': 0.95,          # ← Stricter threshold
    'min_coverage': 0.90,          # ← Stricter threshold
    'threads': 8,
}
EOF
echo ""
echo "  Then run: python tick_species_metagenome_identifier.py"
echo ""

echo "✓ APPROACH 2: COMMAND-LINE"
echo ""
echo "  python tick_species_metagenome_identifier.py \\"
echo "    -r Tick_Specifier_Ticks.fasta \\"
echo "    -i your_sample.fastq \\"
echo "    --min-identity 0.95 \\"
echo "    --min-coverage 0.90"
echo ""

# ============================================================================
# WORKFLOW 3: Sensitive Mode (Detect Novel Species)
# ============================================================================

show_workflow 3 "Sensitive Analysis (Novel Species Detection)"

echo "✓ APPROACH 1: CONFIG-BASED"
echo ""
echo "  Use relaxed thresholds:"
echo ""
cat << 'EOF'
CONFIG = {
    'input_fastq': 'your_sample.fastq',
    'min_identity': 0.80,          # ← More permissive
    'min_coverage': 0.75,          # ← More permissive
    'blast_word_size': 11,         # ← More sensitive BLAST
}
EOF
echo ""

echo "✓ APPROACH 2: COMMAND-LINE"
echo ""
echo "  python tick_species_metagenome_identifier.py \\"
echo "    -r Tick_Specifier_Ticks.fasta \\"
echo "    -i your_sample.fastq \\"
echo "    --min-identity 0.80 \\"
echo "    --min-coverage 0.75"
echo ""

# ============================================================================
# WORKFLOW 4: Multi-Sample Analysis
# ============================================================================

show_workflow 4 "Multi-Sample Analysis with Comparison"

echo "✓ APPROACH 1: CONFIG-BASED (Recommended)"
echo ""
echo "  Create a Python script: batch_analysis.py"
echo ""
cat << 'EOF'
import subprocess
from pathlib import Path
import shutil

# Process each sample
for fastq in Path('.').glob('*.fastq'):
    sample_name = fastq.stem
    
    # Copy template script
    shutil.copy('tick_species_metagenome_identifier.py', f'run_{sample_name}.py')
    
    # Read and modify config
    with open(f'run_{sample_name}.py', 'r') as f:
        content = f.read()
    
    # Replace config values
    content = content.replace(
        "'input_fastq': 'your_metagenome_reads.fastq'",
        f"'input_fastq': '{fastq}'"
    )
    content = content.replace(
        "'output_dir': './tick_analysis'",
        f"'output_dir': './results_{sample_name}'"
    )
    
    # Write modified script
    with open(f'run_{sample_name}.py', 'w') as f:
        f.write(content)
    
    # Run it
    print(f"Processing {fastq}...")
    subprocess.run(['python', f'run_{sample_name}.py'])
    
    # Clean up
    Path(f'run_{sample_name}.py').unlink()

print("✓ All samples processed")

# Compare results
print("\nNow you can compare results:")
print("  python tick_analysis_utils.py")
EOF
echo ""
echo "  Run it: python batch_analysis.py"
echo ""

echo "✓ APPROACH 2: COMMAND-LINE LOOP"
echo ""
echo "  Create: batch_analysis.sh"
echo ""
cat << 'EOF'
#!/bin/bash

for fastq in *.fastq; do
    sample="${fastq%.fastq}"
    echo "Processing: $fastq"
    
    python tick_species_metagenome_identifier.py \
        -r Tick_Specifier_Ticks.fasta \
        -i "$fastq" \
        -o "results_${sample}" \
        -t 8
done

echo "✓ All samples processed"
EOF
echo ""
echo "  Run it: bash batch_analysis.sh"
echo ""

# ============================================================================
# WORKFLOW 5: Filter Results
# ============================================================================

show_workflow 5 "Filter Results by Confidence"

echo "✓ Using tick_analysis_utils.py"
echo ""
echo "  1. Open tick_analysis_utils.py"
echo "  2. Find OPERATION_FILTER_RESULTS (near line 20)"
echo "  3. Edit the config:"
echo ""
cat << 'EOF'
OPERATION_FILTER_RESULTS = """
config = {
    'input_json': 'tick_analysis/species_identification_results.json',
    'output_file': 'tick_analysis/filtered_results.json',
    'min_mapq': 40,              # Minimum mapping quality
    'min_identity': 0.90,         # Minimum identity
    'require_blast_validation': True,  # Only BLAST-validated
}
"""
EOF
echo ""
echo "  4. Uncomment the operation in main():"
echo "     exec(OPERATION_FILTER_RESULTS)"
echo "     TickAnalysisUtils.filter_by_confidence(...)"
echo ""
echo "  5. Run: python tick_analysis_utils.py"
echo ""

# ============================================================================
# WORKFLOW 6: Generate CSV
# ============================================================================

show_workflow 6 "Export Results to CSV (for Excel)"

echo "✓ Using tick_analysis_utils.py"
echo ""
echo "  1. Open tick_analysis_utils.py"
echo "  2. Find OPERATION_GENERATE_CSV"
echo "  3. Edit:"
echo ""
cat << 'EOF'
OPERATION_GENERATE_CSV = """
config = {
    'input_json': 'tick_analysis/species_identification_results.json',
    'output_csv': 'results.csv',
}
"""
EOF
echo ""
echo "  4. Uncomment in main(): exec(OPERATION_GENERATE_CSV)"
echo "  5. Run: python tick_analysis_utils.py"
echo "  6. Open results.csv in Excel"
echo ""

# ============================================================================
# WORKFLOW 7: Compare Multiple Samples
# ============================================================================

show_workflow 7 "Compare Multiple Samples"

echo "✓ Using tick_analysis_utils.py"
echo ""
echo "  1. Open tick_analysis_utils.py"
echo "  2. Find OPERATION_COMPARE_SAMPLES"
echo "  3. Edit:"
echo ""
cat << 'EOF'
OPERATION_COMPARE_SAMPLES = """
config = {
    'json_files': [
        'sample1_results/species_identification_results.json',
        'sample2_results/species_identification_results.json',
        'sample3_results/species_identification_results.json',
    ],
    'output_file': 'comparison.json',
}
"""
EOF
echo ""
echo "  4. Uncomment in main()"
echo "  5. Run: python tick_analysis_utils.py"
echo ""

# ============================================================================
# Summary
# ============================================================================

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "SUMMARY"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "You have TWO options for configuring and running:"
echo ""
echo "1. CONFIG-BASED (Edit CONFIG dict in script) - RECOMMENDED"
echo "   ✓ Easiest to use"
echo "   ✓ No command-line knowledge needed"
echo "   ✓ Easy to repeat same analysis"
echo "   ✓ Great for reproducibility"
echo ""
echo "2. COMMAND-LINE ARGUMENTS"
echo "   ✓ Good for scripting/automation"
echo "   ✓ Can mix with config file"
echo "   ✓ Override specific settings"
echo ""
echo "3. MIXED APPROACH"
echo "   ✓ Use config file as base"
echo "   ✓ Override specific args from command-line"
echo ""
echo "For detailed help:"
echo "  python tick_species_metagenome_identifier.py --help"
echo "  cat README.md"
echo "  cat QUICK_REFERENCE.md"
echo ""
