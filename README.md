# SCOUT - Species Classification Using COI Tags

## Overview

This toolkit identifies tick species from Oxford Nanopore Technology (ONT) metagenome samples using a dual-approach pipeline:

1. **Minimap2 Alignment**: Fast, sensitive mapping of ONT reads to a COI reference database
2. **BLAST Validation**: High-confidence species assignment through sequence similarity search

The pipeline produces detailed species composition reports, abundance estimates, and validation metrics.

---

## System Requirements

### Required Dependencies
- **Python 3.8+**
- **minimap2** - ONT-optimized long-read aligner
- **BLAST+** - NCBI BLAST suite for sequence validation
- **Python Libraries**: json, re, subprocess, pathlib, collections (standard library)

### Installation

#### Option 1: Conda (Recommended)
```bash
# Create a new environment
conda create -n tick_analysis -c bioconda minimap2 blast python=3.10

# Activate environment
conda activate tick_analysis
```

#### Option 2: Manual Installation
```bash
# Ubuntu/Debian
sudo apt-get install minimap2 ncbi-blast+

# macOS with Homebrew
brew install minimap2 blast
```

### Verify Installation
```bash
minimap2 -h          # Should show minimap2 version
blastn -h            # Should show BLAST version
python3 --version    # Should show Python 3.8+
```

---

## Input Data

### Required Files

1. **Reference Database** (`Tick_Specifier_Ticks.fasta`)
   - FASTA file with tick COI gene sequences
   - Headers format: `>ID_Genus_species_additional_info`
   - Example:
     ```
     >ARAK096-13_Rhipicephalus_maculatus_2011-01-29_KWS-ICIPE_101
     AACAATATATTTAATTTTCGGCGTATGATCTGGGATATTAGGATTAAGAATAAGAATATTGATTCGATTAGAATTAGGGCAACCAGG...
     ```

2. **ONT Metagenome Reads** (your input FASTQ)
   - FASTQ format (Q-score quality lines required for BLAST conversion)
   - Minimum recommended: 500 bp reads
   - Example:
     ```
     @read_001
     ACGTACGTACGTACGT...
     +
     FFFF##FF#####FF...
     ```

---

## Usage

### ✓ Option 1: CONFIG-BASED (Recommended for most users)

**Edit the script and run it:**

1. Open `tick_species_metagenome_identifier.py` in a text editor
2. Find the `CONFIG` section at the top (lines 30-50)
3. Change these settings:
   - `input_fastq`: Your metagenome file
   - `reference_fasta`: Reference database
   - `output_dir`: Where to save results
   - `threads`: Number of CPU cores
   - Other thresholds as needed
4. Save the file
5. Run: `python tick_species_metagenome_identifier.py`

**Example:**
```python
CONFIG = {
    'reference_fasta': 'Tick_Specifier_Ticks.fasta',
    'input_fastq': 'my_sample.fastq',        # ← Change this
    'output_dir': './my_analysis',            # ← And this
    'threads': 8,
    'min_identity': 0.85,
    'min_coverage': 0.80,
}
```

### ✓ Option 2: COMMAND-LINE ARGUMENTS

**For users who prefer command-line interface:**

```bash
python tick_species_metagenome_identifier.py \
  -r Tick_Specifier_Ticks.fasta \
  -i your_metagenome_reads.fastq \
  -o my_analysis_results \
  --threads 8 \
  --min-identity 0.85 \
  --min-coverage 0.80
```

### ✓ Option 3: MIXED (Config file + command-line override)

**Best of both worlds - use config as base, override specific settings:**

1. Edit CONFIG section in script (your defaults)
2. Override with command-line arguments:

```bash
python tick_species_metagenome_identifier.py \
  --threads 16 \
  --min-identity 0.95
```

### Advanced Options

### Parameter Descriptions

Edit these values in the `CONFIG` dictionary at the top of the script:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `input_fastq` | Required | Your ONT metagenome FASTQ file |
| `reference_fasta` | Required | Tick COI reference database |
| `output_dir` | `./tick_analysis` | Where to save results |
| `threads` | 4 | Number of CPU cores to use |
| `min_query_len` | 500 | Minimum read length (bp) |
| `min_identity` | 0.85 | Minimum % identity (0-1) |
| `min_coverage` | 0.80 | Minimum query coverage (0-1) |
| `blast_evalue` | `1e-20` | BLAST E-value threshold |
| `blast_word_size` | 11 | BLAST word size (11=sensitive, 15=faster) |
| `force_reindex` | False | Force minimap2 re-indexing |

---

## Output Files

### Generated in Output Directory

1. **species_identification_report.txt**
   - Human-readable summary report
   - Species abundance table
   - Quality statistics
   - Validation metrics

2. **species_identification_results.json**
   - Machine-readable results
   - Per-read assignments with confidence scores
   - Species abundance data
   - Use this for downstream analysis

3. **alignment files**
   - `*.sam` - Minimap2 alignment (SAM format)
   - `*_blast.tsv` - BLAST results
   - `reference.mmi` - Minimap2 index
   - `reference_blast.*` - BLAST database files

### Example Report Output
```
================================================================================
TICK SPECIES IDENTIFICATION REPORT FROM ONT METAGENOME
================================================================================

Input File: sample_metagenome.fastq
Analysis Date: 2024-06-19T14:32:15.123456
Total Reads Analyzed: 1250

--------------------------------------------------------------------------------
SPECIES ABUNDANCE
--------------------------------------------------------------------------------

Species                                   Reads  %        Validated  Val. %
--------------------------------------------------------------------------------
Rhipicephalus microplus                      450  36.00%          405   90.0%
Amblyomma tholloni                           320  25.60%          298   93.1%
Rhipicephalus appendiculatus                 280  22.40%          265   94.6%
Rhipicephalus evertsi                        200  16.00%          180   90.0%

--------------------------------------------------------------------------------
QUALITY STATISTICS
--------------------------------------------------------------------------------

Total Reads with Alignments: 1250
Reads with BLAST Validation: 1148 (91.8%)
Average MAPQ: 52.3
Average Identity: 0.942
Min/Max Identity: 0.851 / 0.998
```

---

## Analysis Workflow

### Step-by-Step Pipeline

```
Input ONT FASTQ
    ↓
[1] Index Reference (minimap2)
    ↓
[2] Align Reads (minimap2 map-ont)
    ↓
[3] Parse & Filter Alignments (quality thresholds)
    ↓
[4] BLAST Validation
    ↓
[5] Combine Results & Assign Species
    ↓
[6] Estimate Abundance
    ↓
[7] Generate Reports
    ↓
Output Files (TXT + JSON)
```

### Quality Filters Applied

- **Read Length**: Minimum 500 bp (configurable)
- **Mapping Quality (MAPQ)**: ≥ 30 (high confidence)
- **Identity**: ≥ 85% by default (configurable)
- **Coverage**: ≥ 80% query coverage (configurable)
- **BLAST E-value**: < 1e-20 for validation

---

## Predefined Analysis Modes

The script includes predefined configuration templates. Uncomment one to use:

### High-Confidence Mode (Strict QC)
```python
# In the CONFIG section, uncomment this:
CONFIG.update({
    'min_identity': 0.95,
    'min_coverage': 0.90,
    'threads': 8,
})
```

### Sensitive Mode (Detect Novel Species)
```python
# In the CONFIG section, uncomment this:
CONFIG.update({
    'min_identity': 0.80,
    'min_coverage': 0.75,
    'blast_word_size': 11,
})
```

### Fast Mode (Lower Memory)
```python
# In the CONFIG section, uncomment this:
CONFIG.update({
    'threads': 2,
    'blast_word_size': 15,
})
```

---

## Utility Commands

Use `tick_analysis_utils.py` for post-analysis operations.

### How to Use Utilities

1. Open `tick_analysis_utils.py` in a text editor
2. Find the operation you want (see list below)
3. Edit its `config` settings near the top
4. Uncomment the code lines in the `main()` function for that operation
5. Run: `python tick_analysis_utils.py`

### Available Utility Operations

#### 1. Filter Results by Confidence

Edit the config near the top:
```python
OPERATION_FILTER_RESULTS = """
config = {
    'input_json': 'tick_analysis/species_identification_results.json',
    'output_file': 'tick_analysis/filtered_results.json',
    'min_mapq': 40,              # Minimum mapping quality
    'min_identity': 0.90,         # Minimum identity
    'require_blast_validation': True,  # Only BLAST-validated reads
}
"""
```

Then uncomment in `main()`:
```python
exec(OPERATION_FILTER_RESULTS)
TickAnalysisUtils.filter_by_confidence(
    config['input_json'],
    min_mapq=config['min_mapq'],
    min_identity=config['min_identity'],
    require_blast_validation=config['require_blast_validation'],
    output_file=config['output_file']
)
```

#### 2. Extract Species Statistics

Edit the config:
```python
OPERATION_EXTRACT_SPECIES = """
config = {
    'input_json': 'tick_analysis/species_identification_results.json',
    'output_file': 'tick_analysis/species_stats.json',
}
"""
```

#### 3. Generate CSV for Excel

Edit the config:
```python
OPERATION_GENERATE_CSV = """
config = {
    'input_json': 'tick_analysis/species_identification_results.json',
    'output_csv': 'results.csv',
}
"""
```

#### 4. Compare Multiple Samples

Edit the config:
```python
OPERATION_COMPARE_SAMPLES = """
config = {
    'json_files': [
        'analysis1/species_identification_results.json',
        'analysis2/species_identification_results.json',
        'analysis3/species_identification_results.json',
    ],
    'output_file': 'sample_comparison.json',
}
"""
```

#### 5. Identify Mixed/Chimeric Reads

Edit the config:
```python
OPERATION_IDENTIFY_MIXED = """
config = {
    'input_json': 'tick_analysis/species_identification_results.json',
    'output_file': 'tick_analysis/mixed_reads.json',
}
"""
```

---

## Interpreting Results

### Abundance Columns

- **Reads**: Number of reads assigned to species
- **%**: Percentage of total reads
- **Validated**: Number passing BLAST validation
- **Val. %**: Percentage validated (validation rate)

### High Validation Rate Indicates:
- Confident species assignment
- Good sequence quality
- Clear match to reference

### Low Validation Rate May Indicate:
- Novel species variant
- Degraded sequence quality
- Ambiguous species (similar COI sequences)

### MAPQ Interpretation
- MAPQ ≥ 50: Very high confidence alignment
- MAPQ 30-50: High confidence
- MAPQ < 30: Lower confidence (filtered out by default)

---

## Example Analysis

### Complete Analysis on Sample

```bash
#!/bin/bash

# Activate environment
conda activate tick_analysis

# Run main pipeline
python tick_species_metagenome_identifier.py \
  -r Tick_Specifier_Ticks.fasta \
  -i metagenome_sample.fastq \
  -o sample_analysis \
  -t 8

# Generate additional reports
python tick_analysis_utils.py csv \
  sample_analysis/species_identification_results.json \
  -o sample_analysis/results.csv

python tick_analysis_utils.py mixed \
  sample_analysis/species_identification_results.json \
  -o sample_analysis/mixed_reads.json

# View report
cat sample_analysis/species_identification_report.txt
```

### Multi-Sample Comparison

```bash
# Analyze all samples
for fastq in *.fastq; do
    python tick_species_metagenome_identifier.py \
      -r Tick_Specifier_Ticks.fasta \
      -i "$fastq" \
      -o "analysis_$(basename $fastq .fastq)"
done

# Compare results
python tick_analysis_utils.py compare \
  analysis_sample1/species_identification_results.json \
  analysis_sample2/species_identification_results.json \
  analysis_sample3/species_identification_results.json \
  -o comparison_results.json
```

---

## Troubleshooting

### Issue: "minimap2 not found"
```bash
# Reinstall minimap2
conda install -c bioconda minimap2

# Or verify PATH
which minimap2
```

### Issue: "BLAST database creation failed"
```bash
# Check write permissions in output directory
chmod -R u+w output_directory/

# Manually create BLAST database
makeblastdb -in Tick_Specifier_Ticks.fasta -dbtype nucl -out tick_ref
```

### Issue: Low species detection
- Check read quality (try `fastqc metagenome.fastq`)
- Verify reference database is appropriate
- Lower `--min-identity` threshold slightly
- Check for tick sequences in reads: `fastq-grep "ACGTACGTACGT" reads.fastq | wc -l`

### Issue: Out of memory
- Reduce `--threads` parameter
- Process reads in batches (split large FASTQ files)

---

## Performance Considerations

### Runtime Estimates
- **Reference indexing**: 1-2 seconds
- **Read alignment (minimap2)**: ~5-10 seconds per 1M reads
- **BLAST validation**: ~20-30 seconds per 1M reads
- **Total for 1M reads**: ~5-10 minutes on 8 threads

### Memory Usage
- **Base memory**: ~500 MB
- **Per thread**: ~100-200 MB
- **Total for 8 threads**: ~1.5-2 GB

### Optimization Tips
- Increase threads for faster processing
- Use SSD for temporary files
- Process multiple samples in parallel
- Keep reference database indexed

---

## Reference Database Management

### Updating Reference Database

```bash
# Combine new sequences
cat old_reference.fasta new_sequences.fasta > updated_reference.fasta

# Remove duplicates (optional)
python3 << 'EOF'
from pathlib import Path
from collections import defaultdict

seen = set()
unique_seqs = []

with open('updated_reference.fasta') as f:
    header = None
    seq = ''
    for line in f:
        if line.startswith('>'):
            if seq:
                if seq not in seen:
                    unique_seqs.append((header, seq))
                    seen.add(seq)
            header = line.strip()
            seq = ''
        else:
            seq += line.strip()
    if seq and seq not in seen:
        unique_seqs.append((header, seq))

with open('updated_reference_clean.fasta', 'w') as f:
    for header, seq in unique_seqs:
        f.write(f'{header}\n{seq}\n')
EOF
```

---

## Citation and References

- **minimap2**: Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences.
- **BLAST+**: Camacho, C., et al. (2009). BLAST+: architecture and applications.
- **COI Gene**: Folmer, O., et al. (1994). DNA primers for amplification of mitochondrial cytochrome c oxidase subunit I.

---

## Support and Contribution

For issues, suggestions, or improvements:
1. Check the troubleshooting section
2. Verify all dependencies are installed
3. Ensure input files are correctly formatted
4. Contact: [your-email@domain.com]

---

## Version History

### v1.0 (Current)
- Initial release
- Minimap2 alignment
- BLAST validation
- Species assignment
- Abundance estimation
- Utility analysis tools
