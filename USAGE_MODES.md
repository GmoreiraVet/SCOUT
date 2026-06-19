# Tick Species Identification - Usage Modes

## ✓ All Three Approaches Now Supported

Your scripts now support **three different usage modes** - pick whichever you prefer!

---

## Mode 1: CONFIG-BASED (Recommended for Most Users)

Edit the CONFIG dictionary in the script itself, then run.

### Advantages:
- ✓ No command-line knowledge needed
- ✓ Easy to see all parameters at once
- ✓ Simple to repeat the same analysis
- ✓ Great for reproducibility
- ✓ Perfect for beginners

### Usage:

```bash
# 1. Edit the CONFIG section (top of script)
nano tick_species_metagenome_identifier.py

# 2. Change your settings
CONFIG = {
    'input_fastq': 'my_sample.fastq',          # ← Your file
    'reference_fasta': 'Tick_Specifier_Ticks.fasta',
    'output_dir': './my_results',              # ← Output location
    'threads': 8,                               # ← CPU cores
    'min_identity': 0.85,                      # ← Thresholds
    'min_coverage': 0.80,
}

# 3. Save and run
python tick_species_metagenome_identifier.py
```

---

## Mode 2: COMMAND-LINE (For Scripting/Automation)

Use command-line arguments instead of editing files.

### Advantages:
- ✓ Great for scripting/automation
- ✓ Can chain multiple analyses
- ✓ Easy to pass arguments in loops
- ✓ Familiar to Unix/Linux users
- ✓ Perfect for pipelines

### Usage:

```bash
# Basic example
python tick_species_metagenome_identifier.py \
  -r Tick_Specifier_Ticks.fasta \
  -i sample.fastq \
  -o my_results \
  -t 8

# With thresholds
python tick_species_metagenome_identifier.py \
  -r Tick_Specifier_Ticks.fasta \
  -i sample.fastq \
  --min-identity 0.95 \
  --min-coverage 0.90

# In a loop
for fastq in *.fastq; do
  python tick_species_metagenome_identifier.py \
    -r Tick_Specifier_Ticks.fasta \
    -i "$fastq" \
    -o "results_${fastq%.fastq}"
done
```

### Available Arguments:

```
-r, --reference          Path to reference FASTA
-i, --input             Path to input FASTQ
-o, --output            Output directory
-t, --threads           Number of CPU threads
--min-query-len         Minimum read length
--min-identity          Minimum identity (0-1)
--min-coverage          Minimum coverage (0-1)
--blast-evalue          BLAST E-value threshold
--blast-word-size       BLAST word size
--force-reindex         Force minimap2 re-indexing
```

---

## Mode 3: MIXED (CONFIG + Command-Line Override)

Best of both worlds! Use CONFIG as base, override specific values.

### Advantages:
- ✓ CONFIG file as permanent defaults
- ✓ Override specific settings when needed
- ✓ Quick parameter testing
- ✓ Less editing of script files

### Usage:

```bash
# CONFIG in script has defaults:
CONFIG = {
    'input_fastq': 'default_sample.fastq',
    'threads': 4,
    'min_identity': 0.85,
}

# Override from command-line:
python tick_species_metagenome_identifier.py \
  --input my_sample.fastq \
  --threads 16 \
  --min-identity 0.95

# Only some parameters change from command-line
python tick_species_metagenome_identifier.py \
  --threads 8 \
  --force-reindex
```

---

## Utilities - All Approaches Also Supported

### MODE 1: CONFIG-BASED

```bash
# Edit the script
nano tick_analysis_utils.py

# Configure operation near top
OPERATION_FILTER_RESULTS = """
config = {
    'input_json': 'results.json',
    'output_file': 'filtered.json',
    'min_identity': 0.90,
}
"""

# Uncomment in main() and run
python tick_analysis_utils.py
```

### MODE 2: COMMAND-LINE

```bash
# Filter results
python tick_analysis_utils.py filter -i results.json --min-mapq 40

# Generate CSV
python tick_analysis_utils.py csv -i results.json -o output.csv

# Compare samples
python tick_analysis_utils.py compare sample1.json sample2.json -o comparison.json

# Extract species stats
python tick_analysis_utils.py species -i results.json

# Find mixed reads
python tick_analysis_utils.py mixed -i results.json -o mixed.json
```

### MODE 3: MIXED

```bash
# CONFIG has defaults, command-line overrides
python tick_analysis_utils.py filter --min-identity 0.95
```

---

## Quick Comparison Table

| Feature | CONFIG-BASED | COMMAND-LINE | MIXED |
|---------|--------------|--------------|-------|
| Ease of use | ★★★★★ | ★★★☆☆ | ★★★★☆ |
| For scripting | ★★☆☆☆ | ★★★★★ | ★★★★☆ |
| Repeatability | ★★★★★ | ★★★☆☆ | ★★★★★ |
| Learning curve | Easy | Medium | Easy |
| Best for | Beginners | Automation | Everyone |

---

## Setup Script Instructions

```bash
bash setup.sh
```

This will:
1. Check Python version ✓
2. Detect your system (conda/apt/brew)
3. Install dependencies (minimap2, BLAST)
4. Verify installation
5. Show all three usage options

---

## Example: Choosing Your Approach

### Example 1: New User, Single Sample
→ **Use CONFIG-BASED Mode**
```bash
nano tick_species_metagenome_identifier.py  # Edit once
python tick_species_metagenome_identifier.py  # Run
```

### Example 2: Experienced User, Automation
→ **Use COMMAND-LINE Mode**
```bash
for f in *.fastq; do
  python tick_species_metagenome_identifier.py -r ref.fasta -i "$f"
done
```

### Example 3: Testing Different Parameters
→ **Use MIXED Mode**
```bash
# CONFIG has your defaults
python tick_species_metagenome_identifier.py --min-identity 0.95  # Override
python tick_species_metagenome_identifier.py --min-identity 0.80  # Test another
```

---

## Getting Help

```bash
# See available command-line options
python tick_species_metagenome_identifier.py --help
python tick_analysis_utils.py --help
python tick_analysis_utils.py filter --help
python tick_analysis_utils.py csv --help

# Read documentation
cat README.md
cat QUICK_REFERENCE.md
cat USAGE_MODES.md  # This file
```

---

## Backward Compatibility

**Don't worry!** If you're used to command-line tools:
- Old command-line approach **still works** ✓
- All previous scripts are compatible ✓
- No need to learn anything new ✓

If you prefer editing config files:
- New CONFIG-based approach **is available** ✓
- No command-line knowledge needed ✓
- Simpler for most users ✓

---

## Summary

✅ **CONFIG-BASED**: Edit script, run (Easiest)  
✅ **COMMAND-LINE**: Use arguments (Best for automation)  
✅ **MIXED**: CONFIG + overrides (Most flexible)  

**Pick whichever method you prefer - they all work!**
