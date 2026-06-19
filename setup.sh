#!/bin/bash
#
# Setup script for tick species identification pipeline
# Installs dependencies and validates installation
#

set -e  # Exit on error

COLOR_GREEN='\033[0;32m'
COLOR_YELLOW='\033[1;33m'
COLOR_RED='\033[0;31m'
COLOR_NC='\033[0m' # No Color

echo -e "${COLOR_GREEN}"
echo "╔════════════════════════════════════════════════════════════╗"
echo "║  Tick Species Identification Pipeline - Setup Script       ║"
echo "╚════════════════════════════════════════════════════════════╝"
echo -e "${COLOR_NC}\n"

# Check Python version
echo -e "${COLOR_YELLOW}Checking Python version...${COLOR_NC}"
PYTHON_VERSION=$(python3 --version 2>&1 | awk '{print $2}')
PYTHON_MAJOR=$(echo $PYTHON_VERSION | cut -d. -f1)
PYTHON_MINOR=$(echo $PYTHON_VERSION | cut -d. -f2)

if [ "$PYTHON_MAJOR" -lt 3 ] || [ "$PYTHON_MAJOR" -eq 3 -a "$PYTHON_MINOR" -lt 8 ]; then
    echo -e "${COLOR_RED}✗ Python 3.8+ required (found $PYTHON_VERSION)${COLOR_NC}"
    exit 1
fi
echo -e "${COLOR_GREEN}✓ Python $PYTHON_VERSION${COLOR_NC}"

# Detect system and package manager
echo -e "\n${COLOR_YELLOW}Detecting system...${COLOR_NC}"
if command -v conda &> /dev/null; then
    echo -e "${COLOR_GREEN}✓ Found conda${COLOR_NC}"
    PACKAGE_MANAGER="conda"
elif command -v apt-get &> /dev/null; then
    echo -e "${COLOR_GREEN}✓ Found apt-get (Ubuntu/Debian)${COLOR_NC}"
    PACKAGE_MANAGER="apt"
elif command -v brew &> /dev/null; then
    echo -e "${COLOR_GREEN}✓ Found Homebrew (macOS)${COLOR_NC}"
    PACKAGE_MANAGER="brew"
else
    echo -e "${COLOR_RED}✗ No recognized package manager found${COLOR_NC}"
    echo "  Please install conda, apt, or Homebrew"
    exit 1
fi

# Install dependencies
echo -e "\n${COLOR_YELLOW}Installing dependencies...${COLOR_NC}"

if [ "$PACKAGE_MANAGER" = "conda" ]; then
    echo "Installing with conda..."
    
    # Create environment if it doesn't exist
    if ! conda env list | grep -q "tick_analysis"; then
        echo "Creating conda environment 'tick_analysis'..."
        conda create -n tick_analysis -c bioconda minimap2 blast -y
        echo "Activate with: conda activate tick_analysis"
    fi
    
    conda run -n tick_analysis echo "Environment ready"
    INSTALL_PREFIX="conda run -n tick_analysis"
    
elif [ "$PACKAGE_MANAGER" = "apt" ]; then
    echo "Installing with apt-get..."
    sudo apt-get update
    sudo apt-get install -y minimap2 ncbi-blast+
    INSTALL_PREFIX=""
    
elif [ "$PACKAGE_MANAGER" = "brew" ]; then
    echo "Installing with Homebrew..."
    brew install minimap2 blast
    INSTALL_PREFIX=""
fi

# Verify installations
echo -e "\n${COLOR_YELLOW}Verifying installations...${COLOR_NC}"

for cmd in minimap2 blastn makeblastdb; do
    if $INSTALL_PREFIX command -v $cmd &> /dev/null; then
        VERSION=$($INSTALL_PREFIX $cmd -h 2>&1 | head -1)
        echo -e "${COLOR_GREEN}✓ $cmd${COLOR_NC}"
    else
        echo -e "${COLOR_RED}✗ $cmd not found${COLOR_NC}"
        exit 1
    fi
done

# Make scripts executable
echo -e "\n${COLOR_YELLOW}Setting up scripts...${COLOR_NC}"
chmod +x tick_species_metagenome_identifier.py 2>/dev/null || true
chmod +x tick_analysis_utils.py 2>/dev/null || true
echo -e "${COLOR_GREEN}✓ Scripts ready${COLOR_NC}"

# Create test data if reference exists
echo -e "\n${COLOR_YELLOW}Checking reference database...${COLOR_NC}"
if [ -f "Tick_Specifier_Ticks.fasta" ]; then
    FASTA_LINES=$(wc -l < Tick_Specifier_Ticks.fasta)
    echo -e "${COLOR_GREEN}✓ Reference database found ($FASTA_LINES lines)${COLOR_NC}"
else
    echo -e "${COLOR_YELLOW}⚠ Reference database not found${COLOR_NC}"
    echo "  Expected: Tick_Specifier_Ticks.fasta"
fi

# Test small example if provided
if [ $# -ge 1 ] && [ -f "$1" ]; then
    echo -e "\n${COLOR_YELLOW}Running test analysis...${COLOR_NC}"
    TEST_DIR="test_analysis_$(date +%s)"
    mkdir -p "$TEST_DIR"
    
    echo "Input file: $1"
    echo "Output directory: $TEST_DIR"
    
    $INSTALL_PREFIX python3 tick_species_metagenome_identifier.py \
        -r Tick_Specifier_Ticks.fasta \
        -i "$1" \
        -o "$TEST_DIR" \
        -t 4 2>&1 | head -50
    
    echo -e "\n${COLOR_GREEN}✓ Test run completed${COLOR_NC}"
    echo "  Results in: $TEST_DIR/"
    echo "  View report: cat $TEST_DIR/species_identification_report.txt"
fi

echo -e "\n${COLOR_GREEN}"
echo "╔════════════════════════════════════════════════════════════╗"
echo "║  Setup Complete! Ready to analyze tick metagenomes.        ║"
echo "╚════════════════════════════════════════════════════════════╝"
echo -e "${COLOR_NC}"

echo -e "\n${COLOR_YELLOW}Quick start:${COLOR_NC}"
echo "  OPTION 1: Edit CONFIG in script (Recommended)"
echo "    1. nano tick_species_metagenome_identifier.py"
echo "    2. Edit CONFIG section (lines 30-50)"
echo "    3. python tick_species_metagenome_identifier.py"
echo ""
echo "  OPTION 2: Command-line arguments"
echo "    python tick_species_metagenome_identifier.py -r Tick_Specifier_Ticks.fasta -i your_reads.fastq"
echo ""
echo "  OPTION 3: Mix both (config + override with args)"
echo "    python tick_species_metagenome_identifier.py --min-identity 0.95 --threads 16"
echo ""
echo "  View results in: tick_analysis/species_identification_report.txt"

echo -e "\n${COLOR_YELLOW}Documentation:${COLOR_NC}"
echo "  See README.md for detailed usage information"

echo -e "\n${COLOR_YELLOW}For help:${COLOR_NC}"
echo "  python tick_species_metagenome_identifier.py --help"
echo "  python tick_analysis_utils.py --help"
