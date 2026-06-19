#!/bin/bash
#
# Example workflows for tick species identification analysis
# Shows common use cases and analysis patterns
#

set -e

# Configuration
REFERENCE="Tick_Specifier_Ticks.fasta"
THREADS=8
OUTPUT_DIR="./analysis_results"

echo "╔════════════════════════════════════════════════════════════╗"
echo "║    Tick Species Identification - Example Workflows          ║"
echo "╚════════════════════════════════════════════════════════════╝"
echo ""

# Function to run a workflow
run_workflow() {
    local workflow_name=$1
    local fastq_file=$2
    
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "WORKFLOW: $workflow_name"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo ""
}

# ============================================================================
# WORKFLOW 1: Basic Analysis - Standard Thresholds
# ============================================================================

workflow_basic() {
    local fastq=$1
    local output="${OUTPUT_DIR}/basic_analysis"
    
    run_workflow "Basic Analysis (Standard Thresholds)" "$fastq"
    
    echo "Running standard analysis..."
    echo "Command:"
    echo "  python tick_species_metagenome_identifier.py \\"
    echo "    -r $REFERENCE \\"
    echo "    -i $fastq \\"
    echo "    -o $output \\"
    echo "    -t $THREADS"
    echo ""
    
    # Uncomment to run:
    # python tick_species_metagenome_identifier.py \
    #     -r "$REFERENCE" \
    #     -i "$fastq" \
    #     -o "$output" \
    #     -t "$THREADS"
    
    echo "Output files:"
    echo "  • $output/species_identification_report.txt"
    echo "  • $output/species_identification_results.json"
    echo ""
}

# ============================================================================
# WORKFLOW 2: High-Confidence Analysis
# ============================================================================

workflow_high_confidence() {
    local fastq=$1
    local output="${OUTPUT_DIR}/high_confidence_analysis"
    
    run_workflow "High-Confidence Analysis (Strict Thresholds)" "$fastq"
    
    echo "Using strict thresholds for only high-confidence species assignments:"
    echo ""
    echo "Command:"
    echo "  python tick_species_metagenome_identifier.py \\"
    echo "    -r $REFERENCE \\"
    echo "    -i $fastq \\"
    echo "    -o $output \\"
    echo "    --min-identity 0.95 \\"
    echo "    --min-coverage 0.90 \\"
    echo "    -t $THREADS"
    echo ""
    
    echo "Then filter by BLAST validation:"
    echo "  python tick_analysis_utils.py filter $output/species_identification_results.json \\"
    echo "    --min-mapq 40 \\"
    echo "    --min-identity 0.95 \\"
    echo "    --require-blast \\"
    echo "    -o $output/high_confidence_results.json"
    echo ""
}

# ============================================================================
# WORKFLOW 3: Sensitive Analysis (Novel Species Detection)
# ============================================================================

workflow_sensitive() {
    local fastq=$1
    local output="${OUTPUT_DIR}/sensitive_analysis"
    
    run_workflow "Sensitive Analysis (Novel Species Detection)" "$fastq"
    
    echo "Using relaxed thresholds to detect novel or divergent species:"
    echo ""
    echo "Command:"
    echo "  python tick_species_metagenome_identifier.py \\"
    echo "    -r $REFERENCE \\"
    echo "    -i $fastq \\"
    echo "    -o $output \\"
    echo "    --min-identity 0.80 \\"
    echo "    --min-coverage 0.75 \\"
    echo "    -t $THREADS"
    echo ""
    
    echo "Identify potentially novel reads:"
    echo "  python tick_analysis_utils.py mixed $output/species_identification_results.json \\"
    echo "    -o $output/novel_candidates.json"
    echo ""
}

# ============================================================================
# WORKFLOW 4: Multi-Sample Comparison
# ============================================================================

workflow_multi_sample() {
    local output="${OUTPUT_DIR}/multi_sample_comparison"
    
    run_workflow "Multi-Sample Comparison" ""
    
    echo "Analyzing multiple samples and comparing composition:"
    echo ""
    echo "Step 1: Analyze each sample individually"
    echo "  for fastq in sample_*.fastq; do"
    echo "    python tick_species_metagenome_identifier.py \\"
    echo "      -r $REFERENCE \\"
    echo "      -i \"\$fastq\" \\"
    echo "      -o \"${OUTPUT_DIR}/\${fastq%.fastq}\" \\"
    echo "      -t $THREADS"
    echo "  done"
    echo ""
    
    echo "Step 2: Compare results across samples"
    echo "  python tick_analysis_utils.py compare \\"
    echo "    ${OUTPUT_DIR}/sample_1/species_identification_results.json \\"
    echo "    ${OUTPUT_DIR}/sample_2/species_identification_results.json \\"
    echo "    ${OUTPUT_DIR}/sample_3/species_identification_results.json \\"
    echo "    -o $output/comparison_results.json"
    echo ""
    
    echo "Step 3: Export to CSV for Excel"
    echo "  for json in ${OUTPUT_DIR}/*/species_identification_results.json; do"
    echo "    python tick_analysis_utils.py csv \"\$json\" \\"
    echo "      -o \"\${json%/*}/results.csv\""
    echo "  done"
    echo ""
}

# ============================================================================
# WORKFLOW 5: Quality Control and Diagnostics
# ============================================================================

workflow_qc() {
    local fastq=$1
    local output="${OUTPUT_DIR}/qc_analysis"
    
    run_workflow "Quality Control and Diagnostics" "$fastq"
    
    echo "Comprehensive QC analysis:"
    echo ""
    
    echo "Step 1: Run standard analysis"
    echo "  python tick_species_metagenome_identifier.py \\"
    echo "    -r $REFERENCE \\"
    echo "    -i $fastq \\"
    echo "    -o $output \\"
    echo "    -t $THREADS"
    echo ""
    
    echo "Step 2: Extract species statistics"
    echo "  python tick_analysis_utils.py species $output/species_identification_results.json \\"
    echo "    -o $output/species_stats.json"
    echo ""
    
    echo "Step 3: Identify mixed/chimeric reads"
    echo "  python tick_analysis_utils.py mixed $output/species_identification_results.json \\"
    echo "    -o $output/mixed_reads_report.json"
    echo ""
    
    echo "Step 4: Export detailed CSV"
    echo "  python tick_analysis_utils.py csv $output/species_identification_results.json \\"
    echo "    -o $output/detailed_results.csv"
    echo ""
    
    echo "Output files:"
    echo "  • $output/species_identification_report.txt (summary)"
    echo "  • $output/species_stats.json (per-species statistics)"
    echo "  • $output/mixed_reads_report.json (potential chimeras)"
    echo "  • $output/detailed_results.csv (spreadsheet format)"
    echo ""
}

# ============================================================================
# WORKFLOW 6: Batch Processing
# ============================================================================

workflow_batch() {
    local output="${OUTPUT_DIR}/batch_processing"
    
    run_workflow "Batch Processing Multiple Files" ""
    
    echo "Process all FASTQ files in a directory:"
    echo ""
    echo "Script:"
    cat << 'EOF'
#!/bin/bash

REFERENCE="Tick_Specifier_Ticks.fasta"
OUTPUT_BASE="analysis_results"
THREADS=8

# Create output directory
mkdir -p "$OUTPUT_BASE"

# Process each FASTQ file
for fastq_file in *.fastq *.fastq.gz; do
    [ -e "$fastq_file" ] || continue  # Skip if no matches
    
    sample_name="${fastq_file%.*}"
    output_dir="$OUTPUT_BASE/$sample_name"
    
    echo "Processing: $fastq_file"
    
    python tick_species_metagenome_identifier.py \
        -r "$REFERENCE" \
        -i "$fastq_file" \
        -o "$output_dir" \
        -t "$THREADS"
    
    # Generate CSV
    python tick_analysis_utils.py csv \
        "$output_dir/species_identification_results.json" \
        -o "$output_dir/results.csv"
    
    echo "✓ Completed: $sample_name"
done

# Compare all samples
echo "Comparing all samples..."
python tick_analysis_utils.py compare \
    $OUTPUT_BASE/*/species_identification_results.json \
    -o "$OUTPUT_BASE/cross_sample_comparison.json"

echo "✓ All samples processed"
EOF
    
    echo ""
    echo "Save above script and run:"
    echo "  bash batch_analysis.sh"
    echo ""
}

# ============================================================================
# WORKFLOW 7: Advanced Filtering and Analysis
# ============================================================================

workflow_advanced_filtering() {
    local fastq=$1
    local output="${OUTPUT_DIR}/advanced_filtering"
    
    run_workflow "Advanced Filtering and Analysis" "$fastq"
    
    echo "Step-by-step filtering pipeline:"
    echo ""
    
    echo "Step 1: Initial analysis"
    echo "  python tick_species_metagenome_identifier.py \\"
    echo "    -r $REFERENCE -i $fastq -o $output -t $THREADS"
    echo ""
    
    echo "Step 2: Filter by BLAST validation only"
    echo "  python tick_analysis_utils.py filter $output/species_identification_results.json \\"
    echo "    --require-blast -o $output/blast_validated.json"
    echo ""
    
    echo "Step 3: Filter by high MAPQ only"
    echo "  python tick_analysis_utils.py filter $output/species_identification_results.json \\"
    echo "    --min-mapq 50 -o $output/mapq_filtered.json"
    echo ""
    
    echo "Step 4: Filter by high identity only"
    echo "  python tick_analysis_utils.py filter $output/species_identification_results.json \\"
    echo "    --min-identity 0.98 -o $output/high_identity.json"
    echo ""
    
    echo "Step 5: Combine all filters (most stringent)"
    echo "  python tick_analysis_utils.py filter $output/species_identification_results.json \\"
    echo "    --min-mapq 50 --min-identity 0.95 --require-blast \\"
    echo "    -o $output/stringent_filtered.json"
    echo ""
}

# ============================================================================
# Display menu
# ============================================================================

if [ $# -eq 0 ]; then
    echo "Usage: bash example_workflows.sh <workflow_number> [fastq_file]"
    echo ""
    echo "Available workflows:"
    echo "  1. Basic Analysis (Standard Thresholds)"
    echo "  2. High-Confidence Analysis (Strict Thresholds)"
    echo "  3. Sensitive Analysis (Novel Species Detection)"
    echo "  4. Multi-Sample Comparison"
    echo "  5. Quality Control and Diagnostics"
    echo "  6. Batch Processing"
    echo "  7. Advanced Filtering and Analysis"
    echo ""
    echo "Example:"
    echo "  bash example_workflows.sh 1 sample.fastq"
    echo ""
    exit 0
fi

WORKFLOW=$1
FASTQ=${2:-""}

case $WORKFLOW in
    1)
        [ -z "$FASTQ" ] && echo "Error: FASTQ file required" && exit 1
        workflow_basic "$FASTQ"
        ;;
    2)
        [ -z "$FASTQ" ] && echo "Error: FASTQ file required" && exit 1
        workflow_high_confidence "$FASTQ"
        ;;
    3)
        [ -z "$FASTQ" ] && echo "Error: FASTQ file required" && exit 1
        workflow_sensitive "$FASTQ"
        ;;
    4)
        workflow_multi_sample
        ;;
    5)
        [ -z "$FASTQ" ] && echo "Error: FASTQ file required" && exit 1
        workflow_qc "$FASTQ"
        ;;
    6)
        workflow_batch
        ;;
    7)
        [ -z "$FASTQ" ] && echo "Error: FASTQ file required" && exit 1
        workflow_advanced_filtering "$FASTQ"
        ;;
    *)
        echo "Unknown workflow: $WORKFLOW"
        exit 1
        ;;
esac

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
