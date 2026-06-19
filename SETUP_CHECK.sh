#!/usr/bin/env bash
# Quick Verification Checklist
# Verify your setup and understand available usage modes

echo "=========================================="
echo "TICK SPECIES IDENTIFIER - SETUP CHECK"
echo "=========================================="
echo ""

# Check Python
echo "✓ Checking Python version..."
python3 --version
echo ""

# Check files exist
echo "✓ Checking required files..."
files=(
    "tick_species_metagenome_identifier.py"
    "tick_analysis_utils.py"
    "setup.sh"
    "Tick_Specifier_Ticks.fasta"
    "README.md"
)

all_exist=true
for file in "${files[@]}"; do
    if [ -f "$file" ]; then
        echo "  ✓ $file"
    else
        echo "  ✗ $file (MISSING)"
        all_exist=false
    fi
done
echo ""

if [ "$all_exist" = true ]; then
    echo "✓ All files present!"
else
    echo "✗ Some files are missing. Run setup.sh first."
    exit 1
fi

echo ""
echo "=========================================="
echo "USAGE OPTIONS AVAILABLE"
echo "=========================================="
echo ""

echo "✓ OPTION 1: CONFIG-BASED (Edit script)"
echo "  Edit CONFIG dict at top of script"
echo "  $ nano tick_species_metagenome_identifier.py"
echo "  $ python tick_species_metagenome_identifier.py"
echo ""

echo "✓ OPTION 2: COMMAND-LINE (Use arguments)"
echo "  Use command-line flags"
echo "  $ python tick_species_metagenome_identifier.py \\"
echo "      -r Tick_Specifier_Ticks.fasta \\"
echo "      -i your_sample.fastq \\"
echo "      -t 8"
echo ""

echo "✓ OPTION 3: MIXED (CONFIG + override)"
echo "  Edit CONFIG, override with CLI args"
echo "  $ python tick_species_metagenome_identifier.py --threads 16"
echo ""

echo "=========================================="
echo "DOCUMENTATION FILES"
echo "=========================================="
echo ""
echo "  README.md           - Complete documentation"
echo "  QUICK_REFERENCE.md  - Quick lookup for common tasks"
echo "  USAGE_MODES.md      - Detailed usage modes guide"
echo "  example_workflows_dual.sh  - Real-world examples"
echo ""

echo "=========================================="
echo "NEXT STEPS"
echo "=========================================="
echo ""
echo "1. Run setup to install dependencies:"
echo "   $ bash setup.sh"
echo ""
echo "2. Choose your preferred usage approach:"
echo "   - CONFIG-BASED: Easiest for beginners"
echo "   - COMMAND-LINE: Best for automation"
echo "   - MIXED: Most flexible"
echo ""
echo "3. Read the appropriate documentation:"
echo "   $ cat README.md              # Full guide"
echo "   $ cat QUICK_REFERENCE.md     # Quick start"
echo "   $ cat USAGE_MODES.md         # All approaches explained"
echo ""
echo "4. Run the pipeline with your sample data"
echo ""
echo "=========================================="
echo ""
echo "✓ Setup check complete! You're ready to go."
echo ""
