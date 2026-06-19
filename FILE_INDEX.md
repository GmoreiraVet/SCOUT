# Tick Species Identification System - File Index

**Last Updated**: System Complete ✅  
**Status**: Production Ready ✅  
**Compatibility**: Linux, macOS, WSL ✅

---

## 📋 Quick Navigation

### Getting Started (Read These First)
1. **[SYSTEM_COMPLETE.md](SYSTEM_COMPLETE.md)** ← **START HERE**
   - Overview of what's available
   - Three usage approaches explained
   - Quick start guide (3 minutes)

2. **[SETUP_CHECK.sh](SETUP_CHECK.sh)**
   - Verify your installation
   - Check which approaches are available
   - List next steps

### Documentation (Choose One)
- **[README.md](README.md)** - Complete reference documentation
- **[QUICK_REFERENCE.md](QUICK_REFERENCE.md)** - Quick lookup for common tasks
- **[USAGE_MODES.md](USAGE_MODES.md)** - Detailed comparison of approaches

### Examples
- **[example_workflows_dual.sh](example_workflows_dual.sh)** - 7 real-world workflows

---

## 🔧 Core Scripts

### Main Pipeline
**[tick_species_metagenome_identifier.py](tick_species_metagenome_identifier.py)**
- Identifies tick species from ONT metagenomes
- Uses minimap2 + BLAST validation
- **Usage**: CONFIG-based, command-line, or mixed
- **CONFIG**: Lines 30-50 (edit here for simple use)
- **CLI**: Full argparse support with all parameters as flags
- **Output**: JSON results, text report, optional CSV

### Analysis Utilities
**[tick_analysis_utils.py](tick_analysis_utils.py)**
- Post-analysis operations (filtering, CSV export, comparison)
- 5 operations: filter, species, csv, compare, mixed
- **Usage**: CONFIG-based or command-line
- **CLI**: Subcommands for each operation (filter, csv, compare, etc.)
- **Output**: JSON, CSV, or printed statistics

### Installation
**[setup.sh](setup.sh)**
- Automated dependency installation
- Detects system (conda/apt/brew)
- Checks Python 3.8+
- Installs minimap2, BLAST+
- Shows all usage options
- **Run once**: `bash setup.sh`

---

## 📚 Reference Data

**[Tick_Specifier_Ticks.fasta](Tick_Specifier_Ticks.fasta)**
- Reference COI gene database for tick species
- FASTA format with headers: `>ID_Genus_species_additional_info`
- Used by default in pipeline
- Can be replaced with your own reference

---

## 📖 Documentation Files

### For Beginners
**[README.md](README.md)**
- Complete guide with 3 distinct usage sections
- Parameter descriptions
- Predefined configuration templates
- Troubleshooting guide
- Performance notes

### For Quick Lookup
**[QUICK_REFERENCE.md](QUICK_REFERENCE.md)**
- Common tasks at a glance
- Quick start examples for all 3 approaches
- Preset configurations
- Batch processing templates

### For Understanding Approaches
**[USAGE_MODES.md](USAGE_MODES.md)**
- Detailed explanation of CONFIG, CLI, and MIXED approaches
- Comparison table
- Examples for each approach
- Help choosing the best approach for your needs

### System Overview
**[SYSTEM_COMPLETE.md](SYSTEM_COMPLETE.md)**
- Summary of what's been built
- Feature checklist
- Getting started guide
- Design decisions explained

---

## 🚀 Three Usage Approaches (All Available!)

### ✓ Approach 1: CONFIG-BASED
Edit the script config, run without arguments
- **Best for**: Beginners, non-technical users
- **File**: `tick_species_metagenome_identifier.py` (lines 30-50)
- **Usage**: Edit CONFIG dict, save, run script

### ✓ Approach 2: COMMAND-LINE
Use command-line arguments
- **Best for**: Scripting, automation, Unix users
- **File**: Both scripts have full argparse support
- **Usage**: Use flags like `-r`, `-i`, `-t`, etc.

### ✓ Approach 3: MIXED
CONFIG file + command-line overrides
- **Best for**: Testing, flexibility, gradual transitions
- **Usage**: Edit CONFIG as base, override specific flags

---

## 📊 File Organization

```
TICK_DB/
├── tick_species_metagenome_identifier.py    [Main pipeline]
├── tick_analysis_utils.py                   [Post-analysis tools]
├── setup.sh                                  [Installation]
├── SETUP_CHECK.sh                            [Verification]
│
├── Tick_Specifier_Ticks.fasta               [Reference database]
├── COI_TICKS_DWC/                           [Output directory]
│
├── SYSTEM_COMPLETE.md                       [This is the status summary]
├── README.md                                [Complete documentation]
├── QUICK_REFERENCE.md                       [Quick lookup]
├── USAGE_MODES.md                           [Approach comparison]
├── FILE_INDEX.md                            [This file]
│
├── example_workflows_dual.sh                [Real-world examples]
├── example_workflows.sh                     [Original examples]
└── ParseTick_DWC.py                         [Original utility]
```

---

## 🎯 Getting Started

### 1️⃣ Quick Start (3 minutes)
```bash
# Run setup
bash setup.sh

# Edit the script
nano tick_species_metagenome_identifier.py

# Change CONFIG section with your file paths
# Save and run
python tick_species_metagenome_identifier.py
```

### 2️⃣ Or Use Command-Line
```bash
python tick_species_metagenome_identifier.py \
  -r Tick_Specifier_Ticks.fasta \
  -i your_sample.fastq \
  -t 8
```

### 3️⃣ Check Results
```bash
cat tick_analysis/species_identification_report.txt
cat tick_analysis/species_identification_results.json
```

---

## ❓ Finding Answers

| Question | Where to Look |
|----------|---------------|
| How do I get started? | [SYSTEM_COMPLETE.md](SYSTEM_COMPLETE.md) |
| How do I install? | [setup.sh](setup.sh) or [README.md](README.md) |
| What parameters exist? | [README.md](README.md) or [QUICK_REFERENCE.md](QUICK_REFERENCE.md) |
| Which approach should I use? | [USAGE_MODES.md](USAGE_MODES.md) |
| Can you show me examples? | [example_workflows_dual.sh](example_workflows_dual.sh) or [QUICK_REFERENCE.md](QUICK_REFERENCE.md) |
| How do I troubleshoot? | [README.md](README.md) → Troubleshooting section |
| What's available in commands? | `python script.py --help` |
| What file format is used? | [README.md](README.md) → File Formats section |

---

## ✅ What's Working

- ✅ Tick species identification from ONT metagenomes
- ✅ minimap2 alignment (fast, ONT-optimized)
- ✅ BLAST+ validation (high-confidence species)
- ✅ CONFIG-based approach (edit script)
- ✅ Command-line approach (use arguments)
- ✅ Mixed approach (CONFIG + overrides)
- ✅ Quality filtering
- ✅ Multi-sample analysis
- ✅ CSV export
- ✅ Comprehensive documentation
- ✅ Installation automation
- ✅ Backward compatibility

---

## 🔄 Latest Changes

**This session**: 
- ✅ Added command-line support to tick_analysis_utils.py
- ✅ Updated README.md with three usage sections
- ✅ Updated QUICK_REFERENCE.md with all approaches
- ✅ Created USAGE_MODES.md with detailed comparison
- ✅ Created SYSTEM_COMPLETE.md overview
- ✅ Created SETUP_CHECK.sh for verification
- ✅ Created FILE_INDEX.md (this file)

**All changes are backward compatible** - old scripts still work!

---

## 📌 Key Features Summary

| Feature | Method | Status |
|---------|--------|--------|
| Tick species identification | minimap2 + BLAST | ✅ |
| Configuration | CONFIG dict | ✅ |
| Configuration | Command-line args | ✅ |
| Configuration | Mixed (both) | ✅ |
| Quality filtering | Confidence thresholds | ✅ |
| Output formats | JSON, text, CSV | ✅ |
| Multi-sample | Batch processing | ✅ |
| Mixed read detection | Chimera detection | ✅ |
| Installation | Automated script | ✅ |
| Documentation | 4 guides | ✅ |

---

## 🎓 Recommended Reading Order

1. **[SYSTEM_COMPLETE.md](SYSTEM_COMPLETE.md)** (5 min) - Overview
2. **[SETUP_CHECK.sh](SETUP_CHECK.sh)** (3 min) - Verify setup
3. **Choose your approach**:
   - Beginner → [README.md](README.md) - Full guide
   - Quick start → [QUICK_REFERENCE.md](QUICK_REFERENCE.md) - Short version
   - Choosing → [USAGE_MODES.md](USAGE_MODES.md) - Comparison
4. **[example_workflows_dual.sh](example_workflows_dual.sh)** (10 min) - Real examples
5. **Start analyzing!**

---

## 🚀 You're Ready!

Everything is set up and ready to use. Pick your preferred approach and get started:

- **Prefer editing files?** → Use CONFIG-BASED approach
- **Prefer command-line?** → Use COMMAND-LINE approach
- **Want both?** → Use MIXED approach

**All three work perfectly!**

---

**Questions?** Check README.md, QUICK_REFERENCE.md, or USAGE_MODES.md
