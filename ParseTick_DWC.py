#!/usr/bin/env python3
"""
Convert Darwin Core Archive (TSV) COI data to FASTA format.
Extracts sequences and creates FASTA headers from specimen metadata.
"""

import csv
import sys
from pathlib import Path

CONFIG = {
    # Input/output paths
    "INPUT_FILE": "/home/viroicbas2023/Documents/Gmoreira/TICK_DB/COI_TICKS_DWC",
    "OUTPUT_FILE": "/home/viroicbas2023/Documents/Gmoreira/Tick_Specifier_Ticks.fasta",
    
    # Column names for header construction
    "HEADER_FIELDS": [
        "materialSampleID",
        "scientificName",
        "country",
        "eventDate",
        "catalogNumber"
    ],
    
    # Sequence column name (the raw sequence data)
    "SEQUENCE_COLUMN": None,  # Will auto-detect if None
    
    # Filter options
    "MIN_SEQUENCE_LENGTH": 100,  # Skip sequences shorter than this
    "SKIP_EMPTY": True,  # Skip records with no/empty sequence
    "SKIP_DASHES_ONLY": True,  # Skip sequences that are only dashes
    "REQUIRE_SPECIES_LEVEL": True,  # Only include records with species-level identification
}


def detect_sequence_column(header_row):
    """
    Detect which column contains the sequence data.
    Sequences are long nucleotide strings (A, T, G, C, N, -, etc).
    """
    for i, col_name in enumerate(header_row):
        col_name_lower = col_name.lower()
        # Look for common sequence field names
        if any(x in col_name_lower for x in ['dna', 'sequence', 'nucleotide', 'seq']):
            return i, col_name
    
    # If no obvious sequence column, scan for the longest string columns
    # (sequences are much longer than typical metadata)
    return None, None


def has_species_level_resolution(record):
    """
    Check if record has species-level taxonomic resolution.
    Returns True if scientificName contains a binomial (genus + species epithet).
    """
    sci_name = record.get("scientificName", "").strip()
    
    if not sci_name:
        return False
    
    # Skip indeterminate/unidentified markers
    if any(marker in sci_name.lower() for marker in ['indet', 'sp.', 'cf.', 'nr.', 'aff.', '?']):
        return False
    
    # Species name should have at least 2 words (genus + epithet)
    # e.g., "Hyalomma anatolicum" or "Rhipicephalus sanguineus"
    parts = sci_name.split()
    
    # Must have at least genus + species epithet
    if len(parts) < 2:
        return False
    
    # Basic check: first word is capitalized (genus), second is lowercase (epithet)
    # This is a simple heuristic; not foolproof but catches most cases
    has_proper_format = (
        parts[0][0].isupper() and  # Genus capitalized
        parts[1][0].islower()      # Species epithet lowercase
    )
    
    return has_proper_format


def clean_sequence(seq_str):
    """
    Clean sequence string: remove whitespace, gaps-only placeholders.
    Returns None if sequence is unusable.
    """
    if not seq_str or not isinstance(seq_str, str):
        return None
    
    seq_str = seq_str.strip()
    
    # Skip if empty or only dashes/Ns
    if not seq_str or all(c in '-N?nX' for c in seq_str):
        if CONFIG["SKIP_DASHES_ONLY"]:
            return None
    
    # Remove dashes at start/end, but keep internal ones
    seq_str = seq_str.strip('-')
    
    # Check minimum length
    if len(seq_str) < CONFIG["MIN_SEQUENCE_LENGTH"]:
        return None
    
    return seq_str.upper()


def make_fasta_header(record, fields_to_use, record_number):
    """
    Build a FASTA header from Darwin Core fields.
    Format: >id_metadata1_metadata2...
    """
    header_parts = []
    
    for field in fields_to_use:
        value = record.get(field, "").strip()
        if value and value != "":
            header_parts.append(value)
    
    # If no fields collected, use record number
    if not header_parts:
        header_parts = [f"specimen_{record_number}"]
    
    # Join with underscores, sanitize
    header = "_".join(header_parts)
    header = header.replace(" ", "_").replace("/", "_").replace(":", "_")
    
    return header


def dwc_to_fasta(input_file, output_file):
    """
    Main conversion function.
    Reads DwC TSV, extracts sequences, writes FASTA.
    """
    input_path = Path(input_file)
    output_path = Path(output_file)
    
    if not input_path.exists():
        print(f"ERROR: Input file not found: {input_file}")
        sys.exit(1)
    
    # Ensure output directory exists
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    print(f"Reading DwC data from: {input_file}")
    if CONFIG["REQUIRE_SPECIES_LEVEL"]:
        print("Filter: SPECIES-LEVEL resolution required")
    
    sequences_written = 0
    records_processed = 0
    records_skipped = 0
    records_skipped_no_species = 0
    seq_col_index = None
    seq_col_name = None
    
    # Open TSV file and process
    try:
        with open(input_path, 'r', encoding='utf-8') as infile:
            # Detect delimiter and read
            reader = csv.DictReader(infile, delimiter='\t')
            
            # Get column names from header
            fieldnames = reader.fieldnames
            print(f"Found {len(fieldnames)} columns")
            
            # Find sequence column
            if CONFIG["SEQUENCE_COLUMN"] is None:
                # Auto-detect: look for the rightmost column with sequence-like data
                # In this dataset, it should be a column with long nucleotide strings
                # For now, just identify visually: it's typically after metadata, before datasetName
                # Let's look for long string columns
                seq_col_index = None
                for i, col in enumerate(fieldnames):
                    col_lower = col.lower()
                    if 'sequence' in col_lower or 'dna' in col_lower:
                        seq_col_index = i
                        seq_col_name = col
                        break
                
                # If not found by name, we'll detect from first data row
                if seq_col_index is None:
                    print("Auto-detecting sequence column from data...")
            else:
                seq_col_name = CONFIG["SEQUENCE_COLUMN"]
                if seq_col_name in fieldnames:
                    seq_col_index = fieldnames.index(seq_col_name)
            
            # Process each record
            with open(output_path, 'w', encoding='utf-8') as outfile:
                for record_num, record in enumerate(reader, start=1):
                    records_processed += 1
                    
                    # Check species-level resolution first (before sequence check)
                    if CONFIG["REQUIRE_SPECIES_LEVEL"]:
                        if not has_species_level_resolution(record):
                            records_skipped += 1
                            records_skipped_no_species += 1
                            continue
                    
                    # Find sequence column on first data row if not yet found
                    if seq_col_index is None:
                        # Find longest value (likely sequence)
                        longest_col_idx = 0
                        longest_len = 0
                        for i, val in enumerate(record.values()):
                            if val and len(val) > longest_len:
                                longest_len = len(val)
                                longest_col_idx = i
                        
                        seq_col_index = longest_col_idx
                        seq_col_name = list(record.keys())[longest_col_idx]
                        print(f"Detected sequence column: '{seq_col_name}' (column {seq_col_index})")
                    
                    # Extract sequence
                    seq = record.get(seq_col_name, "")
                    seq_clean = clean_sequence(seq)
                    
                    if not seq_clean:
                        records_skipped += 1
                        continue
                    
                    # Build header
                    header = make_fasta_header(record, CONFIG["HEADER_FIELDS"], record_num)
                    
                    # Write FASTA record
                    outfile.write(f">{header}\n")
                    outfile.write(f"{seq_clean}\n")
                    sequences_written += 1
                    
                    # Progress indicator
                    if record_num % 1000 == 0:
                        print(f"  Processed {record_num} records, written {sequences_written} sequences...")
            
            # Summary
            print(f"\n=== Conversion Summary ===")
            print(f"Total records processed: {records_processed}")
            print(f"Sequences written to FASTA: {sequences_written}")
            print(f"Records skipped (total): {records_skipped}")
            if CONFIG["REQUIRE_SPECIES_LEVEL"]:
                print(f"  - No species-level resolution: {records_skipped_no_species}")
            print(f"Output file: {output_path}")
            
            if sequences_written == 0:
                print(f"\n⚠️  WARNING: No sequences met the filtering criteria!")
                print(f"Output file NOT created.")
                output_path.unlink(missing_ok=True)
            else:
                print(f"\n✓ FASTA file ready for analysis!")
    
    except Exception as e:
        print(f"ERROR during conversion: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    dwc_to_fasta(CONFIG["INPUT_FILE"], CONFIG["OUTPUT_FILE"])