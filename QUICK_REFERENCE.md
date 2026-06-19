# Tick Species Identification - Quick Reference

## Two Usage Approaches Available

### ✓ Approach 1: CONFIG-BASED (Recommended)
Edit the script configuration, no command-line knowledge needed

### ✓ Approach 2: COMMAND-LINE 
Use command-line arguments (great for scripting/automation)

### ✓ Approach 3: MIXED
Use config file as base + command-line overrides

---

| File | Purpose |
|------|---------|
| `tick_species_metagenome_identifier.py` | Main pipeline - alignment + BLAST analysis |
| `tick_analysis_utils.py` | Utility tools for filtering and analysis |
| `setup.sh` | Install dependencies and verify setup |
| `example_workflows.sh` | Common analysis workflows |
| `README.md` | Full documentation |
| `Tick_Specifier_Ticks.fasta` | Reference tick COI database |

## Installation (One-Time Setup)

```bash
# Automatic setup (recommended)
bash setup.sh

# Manual installation
conda create -n tick_analysis -c bioconda minimap2 blast -y
conda activate tick_analysis
```

## Quick Start

### Method 1: Edit CONFIG (Recommended)
```bash
# Edit the config at top of script
nano tick_species_metagenome_identifier.py

# Change CONFIG section:
CONFIG = {
    'input_fastq': 'your_sample.fastq',
    'reference_fasta': 'Tick_Specifier_Ticks.fasta',
    'threads': 8,
}

# Save and run
python tick_species_metagenome_identifier.py
```

### Method 2: Command-Line Only
```bash
python tick_species_metagenome_identifier.py \
  -r Tick_Specifier_Ticks.fasta \
  -i your_sample.fastq \
  -t 8
```

### Method 3: Mixed (CONFIG + Override)
```bash
# Edit config in script first, then override with args:
python tick_species_metagenome_identifier.py \
  --threads 16 \
  --min-identity 0.95
```

## Common Configurations

### 1. Standard Analysis (Default)
```python
CONFIG = {
    'input_fastq': 'reads.fastq',
    'reference_fasta': 'Tick_Specifier_Ticks.fasta',
    'threads': 8,
    'min_identity': 0.85,
    'min_coverage': 0.80,
}
```

### 2. High-Confidence Mode (Strict QC)
```python
CONFIG = {
    'input_fastq': 'reads.fastq',
    'reference_fasta': 'Tick_Specifier_Ticks.fasta',
    'threads': 8,
    'min_identity': 0.95,        # Stricter
    'min_coverage': 0.90,         # Stricter
}
```

### 3. Sensitive Mode (Novel Species)
```python
CONFIG = {
    'input_fastq': 'reads.fastq',
    'reference_fasta': 'Tick_Specifier_Ticks.fasta',
    'threads': 8,
    'min_identity': 0.80,         # More permissive
    'min_coverage': 0.75,         # More permissive
}
```

### 4. Fast Mode (Low Memory)
```python
CONFIG = {
    'input_fastq': 'reads.fastq',
    'reference_fasta': 'Tick_Specifier_Ticks.fasta',
    'threads': 2,                 # Fewer threads
    'blast_word_size': 15,        # Faster BLAST
}
```

## Utility Operations

### Filter by Confidence
Edit `tick_analysis_utils.py`:
```python
# Near top, configure:
OPERATION_FILTER_RESULTS = """
config = {
    'input_json': 'tick_analysis/species_identification_results.json',
    'output_file': 'filtered.json',
    'min_mapq': 40,
    'min_identity': 0.90,
    'require_blast_validation': True,
}
"""

# In main(), uncomment:
# exec(OPERATION_FILTER_RESULTS)
# TickAnalysisUtils.filter_by_confidence(...)
```

### Generate CSV
```python
# Configure:
OPERATION_GENERATE_CSV = """
config = {
    'input_json': 'tick_analysis/species_identification_results.json',
    'output_csv': 'results.csv',
}
"""

# Uncomment in main():
# exec(OPERATION_GENERATE_CSV)
# TickAnalysisUtils.generate_csv_report(...)
```

### Compare Samples
```python
# Configure:
OPERATION_COMPARE_SAMPLES = """
config = {
    'json_files': [
        'sample1_analysis/species_identification_results.json',
        'sample2_analysis/species_identification_results.json',
    ],
    'output_file': 'comparison.json',
}
"""

# Uncomment in main():
# exec(OPERATION_COMPARE_SAMPLES)
# comparison = TickAnalysisUtils.compare_samples(...)
```

### Extract Species Stats
```python
# Configure:
OPERATION_EXTRACT_SPECIES = """
config = {
    'input_json': 'tick_analysis/species_identification_results.json',
    'output_file': 'stats.json',
}
"""

# Uncomment in main():
# exec(OPERATION_EXTRACT_SPECIES)
# stats = TickAnalysisUtils.extract_species_list(...)
```

### Find Mixed Reads
```python
# Configure:
OPERATION_IDENTIFY_MIXED = """
config = {
    'input_json': 'tick_analysis/species_identification_results.json',
    'output_file': 'mixed_reads.json',
}
"""

# Uncomment in main():
# exec(OPERATION_IDENTIFY_MIXED)
# mixed = TickAnalysisUtils.identify_mixed_reads(...)
```

## Understanding Output

### Report.txt - Key Metrics

| Metric | What It Means |
|--------|--------------|
| Species | Identified tick species |
| Reads | Number of reads assigned |
| % | Relative abundance |
| Validated | Reads passing BLAST validation |
| Val. % | Validation rate (confidence) |
| MAPQ | Mapping quality (30=good, 50+=excellent) |
| Identity | % sequence match to reference |

### Species Abundance Interpretation

**High Validation Rate (>90%)**
- ✓ Confident species assignment
- ✓ Good sequence quality
- ✓ Clear match to reference

**Low Validation Rate (<70%)**
- ⚠ Divergent from reference
- ⚠ Possible novel variant
- ⚠ Check with more sensitive analysis

## Performance Tuning

| Issue | Solution |
|-------|----------|
| Slow processing | Increase threads: `-t 16` |
| Missing species | Lower thresholds: `--min-identity 0.80` |
| Too many low-quality hits | Increase thresholds: `--min-identity 0.95` |
| Out of memory | Reduce threads, process smaller batches |

## Troubleshooting

```bash
# Check if minimap2 is installed
which minimap2

# Check if BLAST is installed
which blastn

# Verify reference database
head Tick_Specifier_Ticks.fasta

# Check FASTQ format
head -4 your_metagenome.fastq

# Count reads
grep -c "^@" your_metagenome.fastq
```

## Batch Processing Template

Create a file called `batch_analysis.sh`:

```bash
#!/bin/bash

REFERENCE="Tick_Specifier_Ticks.fasta"
THREADS=8

for fastq in *.fastq; do
    sample="${fastq%.fastq}"
    echo "Processing: $fastq"
    
    # Create a copy of the script for this sample
    cp tick_species_metagenome_identifier.py "run_${sample}.py"
    
    # Edit the config in the copy
    sed -i "s/'input_fastq': '[^']*'/'input_fastq': '${fastq}'/" "run_${sample}.py"
    sed -i "s/'output_dir': '[^']*'/'output_dir': '.\/results_${sample}'/" "run_${sample}.py"
    
    # Run it
    python "run_${sample}.py"
    
    # Clean up
    rm "run_${sample}.py"
done

echo "✓ All samples processed"
```

Run it:
```bash
bash batch_analysis.sh
```

## Parameter Reference

| Parameter | Default | Range | Notes |
|-----------|---------|-------|-------|
| `--min-identity` | 0.85 | 0.00-1.00 | Lower = more sensitive |
| `--min-coverage` | 0.80 | 0.00-1.00 | Coverage of query sequence |
| `--min-query-len` | 500 | 100-5000 | Minimum read length (bp) |
| `-t/--threads` | 4 | 1-32 | CPU cores to use |

## Tips for Best Results

1. **Pre-process reads**: Use `fastqc` to assess read quality
2. **Use appropriate thresholds**: 
   - 0.95+ for known species
   - 0.85-0.90 for species detection
   - 0.80 for novel/divergent detection
3. **Check validation rates**: >90% = high confidence
4. **Run high-confidence filter**: `--require-blast` for critical results
5. **Compare MAPQ and identity**: Both should be high for good assignments
6. **Batch similar analyses**: Combine samples for better statistics
7. **Export to CSV**: Easier to review in Excel with sorting/filtering

## Output Files Explained

```
tick_analysis/
├── reference.mmi                              # Minimap2 index
├── reference_blast.nhr/nin/nsq               # BLAST database
├── *.sam                                      # Alignments (SAM format)
├── *_blast.tsv                                # BLAST results
├── species_identification_report.txt          # Summary report ← READ THIS
└── species_identification_results.json        # Detailed results
```

## Advanced Workflows

### Multi-Step Analysis
```bash
# 1. Run analysis
python tick_species_metagenome_identifier.py \
  -r Tick_Specifier_Ticks.fasta \
  -i reads.fastq -o analysis

# 2. Filter by BLAST validation
python tick_analysis_utils.py filter \
  analysis/species_identification_results.json \
  --require-blast -o analysis/blast_only.json

# 3. Export to CSV
python tick_analysis_utils.py csv \
  analysis/blast_only.json \
  -o analysis/final_results.csv

# 4. Identify outliers
python tick_analysis_utils.py mixed \
  analysis/species_identification_results.json \
  -o analysis/outliers.json
```

## System Requirements

- **Python**: 3.8+
- **RAM**: 2 GB minimum (8 GB recommended)
- **Storage**: ~100 MB for reference + 2x input size for temp files
- **CPU**: 4+ cores recommended

## Common Issues & Solutions

**Error: "minimap2 not found"**
```bash
conda install -c bioconda minimap2
```

**Error: "BLAST database creation failed"**
```bash
chmod -R u+w output_directory/
```

**Low species detection**
```bash
# Try more sensitive thresholds
python tick_species_metagenome_identifier.py \
  -r Tick_Specifier_Ticks.fasta \
  -i reads.fastq \
  --min-identity 0.80 \
  --min-coverage 0.75
```

**Out of memory**
```bash
# Reduce threads
python tick_species_metagenome_identifier.py \
  -r Tick_Specifier_Ticks.fasta \
  -i reads.fastq \
  -t 2
```

## Examples

### Example 1: Simple Analysis

**CONFIG-BASED:**
```bash
# Edit CONFIG, then:
python tick_species_metagenome_identifier.py
```

**COMMAND-LINE:**
```bash
python tick_species_metagenome_identifier.py \
  -r Tick_Specifier_Ticks.fasta \
  -i sample.fastq
```

### Example 2: Strict Quality Control

**CONFIG-BASED:**
```python
CONFIG = {
    'input_fastq': 'sample.fastq',
    'reference_fasta': 'Tick_Specifier_Ticks.fasta',
    'min_identity': 0.95,
    'min_coverage': 0.90,
    'threads': 8,
}
```

**COMMAND-LINE:**
```bash
python tick_species_metagenome_identifier.py \
  -r Tick_Specifier_Ticks.fasta \
  -i sample.fastq \
  --min-identity 0.95 \
  --min-coverage 0.90
```

**MIXED (CONFIG + OVERRIDE):**
```bash
# Edit CONFIG in script with defaults, then:
python tick_species_metagenome_identifier.py \
  --min-identity 0.95 \
  --min-coverage 0.90
```

### Example 3: Multi-Sample Analysis

**CONFIG-BASED BATCH SCRIPT:**

Create `batch_analysis.py`:
```python
import subprocess
from pathlib import Path
import shutil

for fastq in Path('.').glob('*.fastq'):
    sample_name = fastq.stem
    print(f"\nProcessing {fastq}...")
    
    # Copy and modify script for this sample
    shutil.copy('tick_species_metagenome_identifier.py', f'run_{sample_name}.py')
    
    with open(f'run_{sample_name}.py', 'r') as f:
        content = f.read()
    
    content = content.replace(
        "'input_fastq': 'your_metagenome_reads.fastq'",
        f"'input_fastq': '{fastq}'"
    )
    
    with open(f'run_{sample_name}.py', 'w') as f:
        f.write(content)
    
    subprocess.run(['python', f'run_{sample_name}.py'])
    Path(f'run_{sample_name}.py').unlink()

print("\n✓ All done!")
```

Run it:
```bash
python batch_analysis.py
```

**COMMAND-LINE BATCH SCRIPT:**

Create `batch_analysis.sh`:
```bash
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
```

Run it:
```bash
bash batch_analysis.sh
```

## Getting Help

- Full documentation: See `README.md`
- Workflow examples: `bash example_workflows.sh`
- Script help: `python tick_species_metagenome_identifier.py --help`
- Utility help: `python tick_analysis_utils.py --help`

## Citation

If using this toolkit in research, please cite:
- minimap2: Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences
- BLAST+: Camacho, C., et al. (2009). BLAST+: architecture and applications

---

**Last Updated**: 2024-06-19  
**Version**: 1.0
