#!/usr/bin/env python3
"""
Tick Species Identification from ONT Metagenome Samples - DUAL GENE VERSION - BATCH MODE

This script identifies tick species from Oxford Nanopore Technology (ONT) 
metagenome samples with STRICT GENE-TYPE VALIDATION:
1. Dual minimap2 alignment against COI AND 18S databases
2. Gene-type validation (rejects 18S reads falsely assigned as COI and vice versa)
3. Separate species identification pipelines for COI and 18S
4. Combined reports with gene-type breakdown
5. Export COI and 18S reads as separate FASTA files for verification

BATCH MODE: Process all FASTQ files in a directory or specific files list.
"""

import subprocess
import os
import sys
import shutil
from pathlib import Path
from collections import defaultdict, Counter
import json
from typing import Dict, List, Tuple, Set, Optional, Union
import logging
import re
from datetime import datetime
import glob
import gzip

# ============================================================================
# CONFIGURATION - EDIT THESE SETTINGS
# ============================================================================

CONFIG = {
    # Input/Output paths
    'reference_coi': '/home/viroicbas2023/Documents/Gmoreira/SCOUT/Tick_Specifier_Ticks_BOLD_COI.fasta',  # Reference tick COI database
    'reference_18s': '/home/viroicbas2023/Documents/Gmoreira/SCOUT/tick_18S_high_quality.fasta',  # Reference tick 18S rRNA database
    
    # BATCH MODE: Set your input directory containing FASTQ files
    'input_dir': '/home/viroicbas2023/Documents/Gmoreira/Carraças_Metagenomica_Gui/Reads',  # Directory with FASTQ files
    'input_pattern': '*.fastq.gz*',  # Pattern to match FASTQ files (e.g., '*.fastq', '*_filtered.fastq', '*.fastq*')
    'output_base_dir': '/home/viroicbas2023/Downloads/Carraças_SCOUT_TEMP',  # Base output directory
    
    # Processing parameters
    'threads': 30,                       # Number of CPU threads to use
    'min_query_len': 500,               # Minimum read length (bp) - filter shorter reads
    'min_identity': 0.99,               # Minimum identity 0-1 (0.95=strict, 0.85=balanced, 0.80=sensitive)
    'min_coverage': 0.50,               # Minimum query coverage 0-1
    
    # BLAST parameters
    'blast_evalue': '1e-20',            # BLAST E-value threshold
    'blast_word_size': 11,              # BLAST word size (11=sensitive, 15=faster)
    'blast_num_alignments': 5,          # Number of BLAST hits to keep
    
    # Gene-type separation (CRITICAL for avoiding 18S/COI misclassification)
    'strict_gene_validation': True,     # Reject reads with conflicting COI/18S assignments
    'gene_type_identity_diff': 0.05,    # If COI hits with 95% and 18S with 94%, diff=0.01 - use COI
    
    # Analysis modes
    'force_reindex': False,             # Force minimap2 re-indexing
    'keep_temp_files': False,           # Keep temporary files (SAM, BLAST, etc.) for debugging
    'export_fasta': True,               # Export COI and 18S reads as separate FASTA files
    'generate_summary_report': True,    # Generate batch summary report
}

# ============================================================================

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class DualGeneTickIdentifier:
    """Tick species identifier with COI/18S gene-type separation."""
    
    def __init__(self, 
                 reference_coi: str,
                 reference_18s: str,
                 output_dir: str = "./tick_analysis",
                 threads: int = 4,
                 min_query_len: int = 500,
                 min_identity: float = 0.85,
                 min_coverage: float = 0.80,
                 strict_gene_validation: bool = True,
                 keep_temp_files: bool = False,
                 export_fasta: bool = True):
        """
        Initialize the dual-gene identifier.
        
        Args:
            reference_coi: Path to tick COI reference database
            reference_18s: Path to tick 18S rRNA reference database
            output_dir: Output directory for results
            threads: Number of threads for tools
            min_query_len: Minimum query length to keep
            min_identity: Minimum identity threshold (0-1)
            min_coverage: Minimum query coverage threshold (0-1)
            strict_gene_validation: Reject reads with conflicting gene assignments
            keep_temp_files: Keep temporary files for debugging
            export_fasta: Export COI and 18S reads as FASTA files
        """
        self.reference_coi = Path(reference_coi)
        self.reference_18s = Path(reference_18s)
        self.output_dir = Path(output_dir)
        self.threads = threads
        self.min_query_len = min_query_len
        self.min_identity = min_identity
        self.min_coverage = min_coverage
        self.strict_gene_validation = strict_gene_validation
        self.keep_temp_files = keep_temp_files
        self.export_fasta = export_fasta
        
        # Track temporary files for cleanup
        self.temp_files = []
        
        # Create output subdirectories
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.output_coi_dir = self.output_dir / "COI_analysis"
        self.output_18s_dir = self.output_dir / "18S_analysis"
        self.output_coi_dir.mkdir(exist_ok=True)
        self.output_18s_dir.mkdir(exist_ok=True)
        
        # Validate references exist
        if not self.reference_coi.exists():
            raise FileNotFoundError(f"COI reference database not found: {self.reference_coi}")
        if not self.reference_18s.exists():
            raise FileNotFoundError(f"18S reference database not found: {self.reference_18s}")
        
        logger.info(f"Initialized DualGeneTickIdentifier")
        logger.info(f"  COI ref: {self.reference_coi}")
        logger.info(f"  18S ref: {self.reference_18s}")
        logger.info(f"  Keep temp files: {self.keep_temp_files}")
        logger.info(f"  Export FASTA: {self.export_fasta}")
    
    def _add_temp_file(self, file_path: Path):
        """Track a temporary file for later cleanup."""
        if not self.keep_temp_files:
            self.temp_files.append(file_path)
    
    def _cleanup_temp_files(self):
        """Remove temporary files if not keeping them."""
        if self.keep_temp_files:
            logger.info("Keeping temporary files (debug mode)")
            return
        
        logger.info(f"Cleaning up {len(self.temp_files)} temporary files...")
        removed = 0
        for file_path in self.temp_files:
            try:
                if file_path.exists():
                    if file_path.is_dir():
                        shutil.rmtree(file_path)
                    else:
                        file_path.unlink()
                    removed += 1
                elif file_path.with_suffix('.mmi').exists():
                    # Handle minimap2 index files
                    file_path.with_suffix('.mmi').unlink()
                    removed += 1
            except Exception as e:
                logger.debug(f"Could not remove {file_path}: {e}")
        
        # Clean up BLAST database files
        for blast_db in self.output_dir.glob("reference_blast_*"):
            try:
                for ext in ['.nhr', '.nin', '.nsq', '.ndb', '.not']:
                    db_file = Path(str(blast_db) + ext)
                    if db_file.exists():
                        db_file.unlink()
                        removed += 1
            except Exception as e:
                logger.debug(f"Could not remove BLAST db {blast_db}: {e}")
        
        # Clean up SAM files
        for sam_file in self.output_dir.glob("*.sam"):
            try:
                sam_file.unlink()
                removed += 1
            except Exception as e:
                logger.debug(f"Could not remove SAM file {sam_file}: {e}")
        
        # Clean up BLAST output files
        for blast_file in self.output_dir.glob("*_blast.tsv"):
            try:
                blast_file.unlink()
                removed += 1
            except Exception as e:
                logger.debug(f"Could not remove BLAST file {blast_file}: {e}")
        
        # Clean up temporary FASTA files (not the exported ones)
        for fasta_file in self.output_dir.glob("*_temp.fasta"):
            try:
                fasta_file.unlink()
                removed += 1
            except Exception as e:
                logger.debug(f"Could not remove temp FASTA {fasta_file}: {e}")
        
        # Clean up minimap2 index files
        for mmi_file in self.output_dir.glob("*.mmi"):
            try:
                mmi_file.unlink()
                removed += 1
            except Exception as e:
                logger.debug(f"Could not remove mmi file {mmi_file}: {e}")
        
        logger.info(f"Removed {removed} temporary files")
    
    def check_dependencies(self) -> bool:
        """Check if required tools are installed."""
        tools = ['minimap2', 'blastn', 'makeblastdb']
        missing = []
        
        for tool in tools:
            try:
                subprocess.run([tool, '-h'], 
                             capture_output=True, 
                             timeout=5)
            except (FileNotFoundError, subprocess.TimeoutExpired):
                missing.append(tool)
        
        if missing:
            logger.error(f"Missing required tools: {', '.join(missing)}")
            logger.error("Install with: conda install -c bioconda minimap2 blast")
            return False
        
        logger.info("All required dependencies found")
        return True
    
    def index_references(self, force: bool = False) -> Tuple[Path, Path]:
        """Create minimap2 indices for both COI and 18S databases."""
        coi_index = self.output_dir / "reference_coi.mmi"
        ref_18s_index = self.output_dir / "reference_18s.mmi"
        
        # Track index files for cleanup
        self._add_temp_file(coi_index)
        self._add_temp_file(ref_18s_index)
        
        # Index COI
        if not coi_index.exists() or force:
            logger.info("Creating minimap2 index for COI...")
            cmd = ['minimap2', '-d', str(coi_index), str(self.reference_coi)]
            try:
                subprocess.run(cmd, check=True, capture_output=True)
                logger.info(f"COI index created: {coi_index}")
            except subprocess.CalledProcessError as e:
                logger.error(f"COI indexing failed: {e}")
                raise
        else:
            logger.info(f"COI index exists: {coi_index}")
        
        # Index 18S
        if not ref_18s_index.exists() or force:
            logger.info("Creating minimap2 index for 18S...")
            cmd = ['minimap2', '-d', str(ref_18s_index), str(self.reference_18s)]
            try:
                subprocess.run(cmd, check=True, capture_output=True)
                logger.info(f"18S index created: {ref_18s_index}")
            except subprocess.CalledProcessError as e:
                logger.error(f"18S indexing failed: {e}")
                raise
        else:
            logger.info(f"18S index exists: {ref_18s_index}")
        
        return coi_index, ref_18s_index
    
    def align_reads_minimap2(self, 
                            reads_fastq: str,
                            index_file: Path,
                            gene_type: str) -> Path:
        """
        Align ONT reads to reference using minimap2.
        
        Args:
            reads_fastq: Path to input FASTQ file
            index_file: Path to minimap2 index
            gene_type: 'COI' or '18S' (for output naming)
            
        Returns:
            Path to output SAM file
        """
        reads_path = Path(reads_fastq)
        if not reads_path.exists():
            raise FileNotFoundError(f"Reads file not found: {reads_fastq}")
        
        sam_file = self.output_dir / f"{reads_path.stem}_{gene_type}.sam"
        # Handle .fastq.gz files - remove .fastq from stem if present
        if reads_path.stem.endswith('.fastq'):
            sam_file = self.output_dir / f"{reads_path.stem.replace('.fastq', '')}_{gene_type}.sam"
        
        self._add_temp_file(sam_file)
        
        logger.info(f"Aligning reads to {gene_type} with minimap2...")
        cmd = [
            'minimap2',
            '-a',
            '-x', 'map-ont',
            '-t', str(self.threads),
            '--secondary=no',
            str(index_file),
            str(reads_fastq)
        ]
        
        try:
            with open(sam_file, 'w') as f:
                subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE, check=True)
            logger.info(f"{gene_type} alignment complete: {sam_file}")
            return sam_file
        except subprocess.CalledProcessError as e:
            logger.error(f"{gene_type} alignment failed: {e.stderr.decode()}")
            raise
    
    def parse_sam_alignment(self, sam_file: Path, gene_type: str) -> Dict[str, List[Dict]]:
        """Parse SAM file and extract high-quality alignments."""
        alignments = defaultdict(list)
        
        logger.info(f"Parsing {gene_type} SAM alignments...")
        with open(sam_file, 'r') as f:
            for line in f:
                if line.startswith('@'):
                    continue
                
                fields = line.strip().split('\t')
                if len(fields) < 11:
                    continue
                
                read_id = fields[0]
                ref_id = fields[2]
                mapping_quality = int(fields[4])
                seq_len = len(fields[9])
                
                # Parse CIGAR to calculate matches
                cigar = fields[5]
                # Extract M (match) and = (equal) operations
                matches = re.findall(r'(\d+)[M=]', cigar)
                matched_len = sum(int(m) for m in matches) if matches else 0
                
                identity = matched_len / seq_len if seq_len > 0 else 0
                coverage = matched_len / seq_len if seq_len > 0 else 0
                
                # Apply filters
                if (seq_len >= self.min_query_len and 
                    identity >= self.min_identity and
                    coverage >= self.min_coverage and
                    mapping_quality >= 30):
                    
                    alignments[read_id].append({
                        'read_id': read_id,
                        'ref_id': ref_id,
                        'gene_type': gene_type,
                        'mapq': mapping_quality,
                        'read_len': seq_len,
                        'matched_len': matched_len,
                        'identity': identity,
                        'coverage': coverage,
                        'sequence': fields[9]  # Store sequence for FASTA export
                    })
        
        logger.info(f"Found {len(alignments)} {gene_type} alignments")
        return alignments
    
    def determine_gene_type(self,
                          coi_alignments: Dict[str, List[Dict]],
                          ref_18s_alignments: Dict[str, List[Dict]]) -> Dict[str, Dict]:
        """
        Determine primary gene type for each read based on best alignment.
        STRICT MUTUAL EXCLUSIVITY: Each read assigned to ONLY ONE gene type.
        Uses identity as primary metric, then MAPQ as tie-breaker.
        """
        gene_assignments = {}
        
        logger.info("Determining gene type for reads (MUTUALLY EXCLUSIVE)...")
        
        # Get all reads that have alignments to either database
        all_reads = set(coi_alignments.keys()) | set(ref_18s_alignments.keys())
        
        for read_id in all_reads:
            coi_best = None
            ref_18s_best = None
            coi_identity = 0
            ref_18s_identity = 0
            
            # Get best COI alignment (by identity, then MAPQ)
            if read_id in coi_alignments and coi_alignments[read_id]:
                coi_best = max(coi_alignments[read_id], 
                              key=lambda x: (x['identity'], x['mapq']))
                coi_identity = coi_best['identity']
            
            # Get best 18S alignment (by identity, then MAPQ)
            if read_id in ref_18s_alignments and ref_18s_alignments[read_id]:
                ref_18s_best = max(ref_18s_alignments[read_id], 
                                  key=lambda x: (x['identity'], x['mapq']))
                ref_18s_identity = ref_18s_best['identity']
            
            # Determine gene type - STRICT MUTUAL EXCLUSIVITY
            gene_type = 'UNASSIGNED'
            conflict = False
            decision_reason = "No valid alignment"
            
            if coi_best and ref_18s_best:
                # Both have hits - compare identity
                identity_diff = coi_identity - ref_18s_identity
                
                # Calculate weighted scores for better decision
                coi_score = coi_identity + (coi_best['mapq'] / 10000)
                ref_18s_score = ref_18s_identity + (ref_18s_best['mapq'] / 10000)
                
                # If identity diff is small (within 3%), it's a potential conflict
                if abs(identity_diff) < 0.03:
                    if self.strict_gene_validation:
                        gene_type = 'CONFLICT'
                        conflict = True
                        decision_reason = f"Identity diff {identity_diff:.3f} < 0.03, strict mode"
                    else:
                        # Choose based on identity + MAPQ score
                        if coi_score >= ref_18s_score:
                            gene_type = 'COI'
                            decision_reason = f"COI score ({coi_score:.4f}) >= 18S score ({ref_18s_score:.4f})"
                        else:
                            gene_type = '18S'
                            decision_reason = f"18S score ({ref_18s_score:.4f}) > COI score ({coi_score:.4f})"
                else:
                    # Clear winner based on identity
                    if identity_diff > 0:
                        gene_type = 'COI'
                        decision_reason = f"COI identity ({coi_identity:.3f}) > 18S identity ({ref_18s_identity:.3f})"
                    else:
                        gene_type = '18S'
                        decision_reason = f"18S identity ({ref_18s_identity:.3f}) > COI identity ({coi_identity:.3f})"
            
            elif coi_best:
                # Only COI hit
                gene_type = 'COI'
                decision_reason = "Only COI alignment found"
            elif ref_18s_best:
                # Only 18S hit
                gene_type = '18S'
                decision_reason = "Only 18S alignment found"
            
            # Store assignment with the best alignment for the chosen gene type
            gene_assignments[read_id] = {
                'gene_type': gene_type,
                'coi_best': coi_best,
                'ref_18s_best': ref_18s_best,
                'conflict': conflict,
                'coi_identity': coi_identity,
                'ref_18s_identity': ref_18s_identity,
                'identity_diff': (coi_identity - ref_18s_identity) if (coi_best and ref_18s_best) else None,
                'decision_reason': decision_reason
            }
        
        # Stats
        stats = Counter(a['gene_type'] for a in gene_assignments.values())
        logger.info(f"Gene type breakdown: COI={stats['COI']}, 18S={stats['18S']}, "
                   f"CONFLICT={stats['CONFLICT']}, UNASSIGNED={stats['UNASSIGNED']}")
        
        # Log some examples of gene type decisions
        sample_reads = list(gene_assignments.items())[:20]
        logger.info("Sample gene type decisions:")
        for read_id, assignment in sample_reads:
            if assignment['gene_type'] in ['COI', '18S', 'CONFLICT']:
                logger.info(f"  {read_id[:30]}... -> {assignment['gene_type']} "
                           f"(COI: {assignment['coi_identity']:.3f}, "
                           f"18S: {assignment['ref_18s_identity']:.3f}) - {assignment['decision_reason']}")
        
        return gene_assignments
    
    def parse_reference_headers(self, reference_fasta: Path) -> Dict[str, str]:
        """
        Parse reference FASTA headers to extract species information.
        Improved to handle various header formats and edge cases.
        """
        species_map = {}
        
        logger.info(f"Parsing reference headers from {reference_fasta.name}...")
        
        with open(reference_fasta, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    header = line[1:].strip()
                    
                    # Extract reference ID (first word or first part before space/pipe)
                    ref_id = header.split()[0]  # First word
                    
                    # Try multiple strategies to extract species name
                    species_name = None
                    
                    # Strategy 1: Look for pattern "Genus species" (two words, first capitalized)
                    # Remove common prefixes first
                    clean_header = header
                    for prefix in ['COI_', '18S_', 'rRNA_', 'gene_', 'sp_']:
                        if clean_header.startswith(prefix):
                            clean_header = clean_header[len(prefix):]
                    
                    # Replace underscores and pipes with spaces
                    clean_header = clean_header.replace('_', ' ').replace('|', ' ')
                    
                    # Split into words and filter empty
                    words = [w for w in clean_header.split() if w.strip()]
                    
                    if len(words) >= 2:
                        # Look for Genus species pattern
                        for i in range(len(words) - 1):
                            # Check for species name pattern: Genus (capitalized) + species (lowercase)
                            if (len(words[i]) > 1 and words[i][0].isupper() and 
                                len(words[i+1]) > 1 and words[i+1][0].islower()):
                                # Make sure it's not an accession number (contains digits or special chars)
                                if (not any(c.isdigit() for c in words[i]) and 
                                    not any(c.isdigit() for c in words[i+1])):
                                    species_name = f"{words[i]} {words[i+1]}"
                                    break
                    
                    # Strategy 2: If no species found, try to extract from header using common patterns
                    if species_name is None:
                        # Look for common species name patterns
                        patterns = [
                            r'([A-Z][a-z]+)\s+([a-z]+)',  # Genus species
                            r'([A-Z][a-z]+)_([a-z]+)',   # Genus_species
                            r'([A-Z][a-z]+)\|([a-z]+)',  # Genus|species
                        ]
                        for pattern in patterns:
                            match = re.search(pattern, header)
                            if match:
                                species_name = f"{match.group(1)} {match.group(2)}"
                                break
                    
                    # Strategy 3: If still no species, use cleaned header without accession numbers
                    if species_name is None:
                        # Remove accession numbers (sequences of digits and letters with | or .)
                        cleaned = re.sub(r'[A-Z]{1,2}\d+\.?\d*', '', header)
                        cleaned = re.sub(r'[|_]\s*', ' ', cleaned)
                        words = [w for w in cleaned.split() if w.strip()]
                        if len(words) >= 2:
                            # Take first two meaningful words
                            species_name = f"{words[0]} {words[1]}"
                        else:
                            species_name = header[:50]  # Fallback
                    
                    species_map[ref_id] = species_name
        
        logger.info(f"Parsed {len(species_map)} reference sequences from {reference_fasta.name}")
        
        # Show some examples of parsed species names
        sample_items = list(species_map.items())[:5]
        if sample_items:
            logger.info("Sample parsed headers:")
            for ref_id, species in sample_items:
                logger.info(f"  {ref_id} -> {species}")
        
        return species_map
    
    def blast_validate(self, 
                      reads_fastq: str,
                      reference_fasta: Path,
                      gene_type: str) -> Path:
        """Run BLAST validation for a specific gene type."""
        reads_path = Path(reads_fastq)
        blast_db = self.output_dir / f"reference_blast_{gene_type}"
        
        # Get stem name handling .fastq.gz properly
        stem = reads_path.stem
        if stem.endswith('.fastq'):
            stem = stem.replace('.fastq', '')
        
        blast_output = self.output_dir / f"{stem}_{gene_type}_blast.tsv"
        fasta_file = self.output_dir / f"{stem}_temp.fasta"
        
        self._add_temp_file(blast_output)
        self._add_temp_file(fasta_file)
        
        # Create BLAST database if needed
        db_files = [Path(str(blast_db) + ext) for ext in ['.nhr', '.nin', '.nsq']]
        if not any(f.exists() for f in db_files):
            logger.info(f"Creating BLAST database for {gene_type}...")
            cmd = [
                'makeblastdb',
                '-in', str(reference_fasta),
                '-dbtype', 'nucl',
                '-out', str(blast_db)
            ]
            try:
                subprocess.run(cmd, check=True, capture_output=True)
                logger.info(f"BLAST database created for {gene_type}")
                # Track BLAST db files for cleanup
                for ext in ['.nhr', '.nin', '.nsq', '.ndb', '.not']:
                    self._add_temp_file(Path(str(blast_db) + ext))
            except subprocess.CalledProcessError as e:
                logger.error(f"BLAST database creation failed: {e}")
                raise
        
        # Convert FASTQ to FASTA (handling gzip if needed)
        self._fastq_to_fasta(reads_fastq, fasta_file)
        
        logger.info(f"Running BLAST validation for {gene_type}...")
        cmd = [
            'blastn',
            '-query', str(fasta_file),
            '-db', str(blast_db),
            '-out', str(blast_output),
            '-outfmt', '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore',
            '-num_threads', str(self.threads),
            '-evalue', '1e-20',
            '-word_size', '11',
            '-num_alignments', '5'
        ]
        
        try:
            subprocess.run(cmd, check=True, capture_output=True)
            logger.info(f"BLAST complete for {gene_type}: {blast_output}")
            return blast_output
        except subprocess.CalledProcessError as e:
            logger.error(f"BLAST failed for {gene_type}: {e}")
            raise
    
    def _fastq_to_fasta(self, fastq_file: str, fasta_file: Path) -> None:
        """Convert FASTQ to FASTA format, handling both compressed and uncompressed files."""
        fastq_path = Path(fastq_file)
        
        # Determine if file is gzipped based on extension
        is_gzipped = fastq_path.suffix == '.gz'
        
        try:
            if is_gzipped:
                # Handle gzipped FASTQ
                with gzip.open(fastq_path, 'rt') as infile, open(fasta_file, 'w') as outfile:
                    while True:
                        header = infile.readline()
                        if not header:
                            break
                        seq = infile.readline()
                        qual_header = infile.readline()
                        qual = infile.readline()
                        if not seq or not qual:
                            break
                        # Remove @ from header and write
                        outfile.write(f">{header[1:]}{seq}")
            else:
                # Handle uncompressed FASTQ
                with open(fastq_file, 'r') as infile, open(fasta_file, 'w') as outfile:
                    while True:
                        header = infile.readline()
                        if not header:
                            break
                        seq = infile.readline()
                        qual_header = infile.readline()
                        qual = infile.readline()
                        if not seq or not qual:
                            break
                        outfile.write(f">{header[1:]}{seq}")
            
            logger.info(f"Converted FASTQ to FASTA: {fasta_file}")
        except Exception as e:
            logger.error(f"Failed to convert {fastq_file} to FASTA: {e}")
            raise
    
    def parse_blast_results(self, blast_output: Path) -> Dict[str, List[Dict]]:
        """Parse BLAST output file."""
        blast_hits = defaultdict(list)
        
        logger.info("Parsing BLAST results...")
        with open(blast_output, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                
                fields = line.strip().split('\t')
                if len(fields) < 12:
                    continue
                
                qseqid = fields[0]
                sseqid = fields[1]
                pident = float(fields[2])
                length = int(fields[3])
                evalue = float(fields[10])
                bitscore = float(fields[11])
                
                if pident >= (self.min_identity * 100) and evalue < 1e-20:
                    blast_hits[qseqid].append({
                        'query': qseqid,
                        'subject': sseqid,
                        'pident': pident,
                        'length': length,
                        'evalue': evalue,
                        'bitscore': bitscore
                    })
        
        logger.info(f"BLAST: Found {len(blast_hits)} reads with valid hits")
        return blast_hits
    
    def assign_species(self,
                      gene_assignments: Dict,
                      coi_alignments: Dict,
                      ref_18s_alignments: Dict,
                      coi_blast_hits: Dict,
                      ref_18s_blast_hits: Dict,
                      coi_species_map: Dict,
                      ref_18s_species_map: Dict) -> Dict[str, Dict]:
        """
        Assign species separately for COI and 18S reads.
        """
        species_assignments = {}
        
        logger.info("Assigning species by gene type...")
        
        for read_id, gene_info in gene_assignments.items():
            gene_type = gene_info['gene_type']
            
            if gene_type == 'COI':
                aln = gene_info['coi_best']
                blast_hits = coi_blast_hits
                species_map = coi_species_map
            elif gene_type == '18S':
                aln = gene_info['ref_18s_best']
                blast_hits = ref_18s_blast_hits
                species_map = ref_18s_species_map
            else:
                species_assignments[read_id] = {
                    'assigned_species': gene_type,
                    'gene_type': gene_type,
                    'minimap_mapq': None,
                    'minimap_identity': None,
                    'minimap_coverage': None,
                    'blast_validated': False,
                    'blast_species': None,
                    'reference_id': None,
                    'coi_identity': gene_info.get('coi_identity', 0),
                    'ref_18s_identity': gene_info.get('ref_18s_identity', 0),
                    'identity_diff': gene_info.get('identity_diff', None),
                    'decision_reason': gene_info.get('decision_reason', ''),
                    'sequence': None
                }
                continue
            
            # Get the best reference ID (remove any trailing version numbers)
            best_ref = aln['ref_id'].split()[0]
            blast_valid = False
            blast_species = None
            
            if read_id in blast_hits and blast_hits[read_id]:
                best_blast = blast_hits[read_id][0]
                blast_ref = best_blast['subject'].split()[0]
                # Validate if BLAST hit matches minimap2 hit or is highly confident
                if blast_ref == best_ref or best_blast['pident'] >= 99.0:
                    blast_valid = True
                    blast_species = species_map.get(blast_ref, blast_ref)
            
            # Get species name, with fallback
            species = species_map.get(best_ref, best_ref)
            
            species_assignments[read_id] = {
                'assigned_species': species,
                'gene_type': gene_type,
                'minimap_mapq': aln['mapq'],
                'minimap_identity': aln['identity'],
                'minimap_coverage': aln['coverage'],
                'blast_validated': blast_valid,
                'blast_species': blast_species,
                'reference_id': best_ref,
                'coi_identity': gene_info.get('coi_identity', 0),
                'ref_18s_identity': gene_info.get('ref_18s_identity', 0),
                'identity_diff': gene_info.get('identity_diff', None),
                'decision_reason': gene_info.get('decision_reason', ''),
                'sequence': aln.get('sequence', '')
            }
        
        # Log summary
        coi_count = sum(1 for a in species_assignments.values() if a['gene_type'] == 'COI')
        ref_18s_count = sum(1 for a in species_assignments.values() if a['gene_type'] == '18S')
        conflict_count = sum(1 for a in species_assignments.values() if a['gene_type'] == 'CONFLICT')
        unassigned_count = sum(1 for a in species_assignments.values() if a['gene_type'] == 'UNASSIGNED')
        
        logger.info(f"Species assigned: COI={coi_count}, 18S={ref_18s_count}, "
                   f"CONFLICT={conflict_count}, UNASSIGNED={unassigned_count}")
        
        return species_assignments
    
    def export_fasta_files(self, 
                          species_assignments: Dict[str, Dict],
                          reads_fastq: str) -> Tuple[Path, Path]:
        """
        Export COI and 18S reads as separate FASTA files.
        
        Returns:
            Tuple of (coi_fasta_path, ref_18s_fasta_path)
        """
        if not self.export_fasta:
            logger.info("FASTA export disabled")
            return None, None
        
        logger.info("Exporting COI and 18S reads as FASTA files...")
        
        reads_path = Path(reads_fastq)
        # Get stem name handling .fastq.gz properly
        stem = reads_path.stem
        if stem.endswith('.fastq'):
            stem = stem.replace('.fastq', '')
        
        coi_fasta = self.output_dir / f"{stem}_COI_reads.fasta"
        ref_18s_fasta = self.output_dir / f"{stem}_18S_reads.fasta"
        conflict_fasta = self.output_dir / f"{stem}_CONFLICT_reads.fasta"
        
        coi_count = 0
        ref_18s_count = 0
        conflict_count = 0
        
        # Write COI reads
        with open(coi_fasta, 'w') as coi_out:
            for read_id, assignment in species_assignments.items():
                if assignment['gene_type'] == 'COI' and assignment.get('sequence'):
                    sequence = assignment['sequence']
                    # Add metadata to header
                    header = f">{read_id} gene=COI species={assignment['assigned_species']} identity={assignment['minimap_identity']:.3f}"
                    coi_out.write(f"{header}\n{sequence}\n")
                    coi_count += 1
        
        # Write 18S reads
        with open(ref_18s_fasta, 'w') as ref_18s_out:
            for read_id, assignment in species_assignments.items():
                if assignment['gene_type'] == '18S' and assignment.get('sequence'):
                    sequence = assignment['sequence']
                    header = f">{read_id} gene=18S species={assignment['assigned_species']} identity={assignment['minimap_identity']:.3f}"
                    ref_18s_out.write(f"{header}\n{sequence}\n")
                    ref_18s_count += 1
        
        # Write CONFLICT reads (if any)
        with open(conflict_fasta, 'w') as conflict_out:
            for read_id, assignment in species_assignments.items():
                if assignment['gene_type'] == 'CONFLICT' and assignment.get('sequence'):
                    sequence = assignment['sequence']
                    header = f">{read_id} gene=CONFLICT coi_identity={assignment.get('coi_identity', 0):.3f} ref_18s_identity={assignment.get('ref_18s_identity', 0):.3f}"
                    conflict_out.write(f"{header}\n{sequence}\n")
                    conflict_count += 1
        
        logger.info(f"Exported {coi_count} COI reads to: {coi_fasta}")
        logger.info(f"Exported {ref_18s_count} 18S reads to: {ref_18s_fasta}")
        if conflict_count > 0:
            logger.info(f"Exported {conflict_count} CONFLICT reads to: {conflict_fasta}")
        
        # Remove temp fasta file if it exists
        temp_fasta = self.output_dir / f"{stem}_temp.fasta"
        if temp_fasta.exists():
            temp_fasta.unlink()
        
        return coi_fasta, ref_18s_fasta
    
    def estimate_abundance(self, 
                          species_assignments: Dict[str, Dict]) -> Tuple[Dict, Dict, Dict]:
        """
        Estimate abundance separately for COI, 18S, and combined.
        
        Returns:
            Tuple of (coi_abundance, ref_18s_abundance, combined_abundance)
        """
        def _count_by_species(assignments, gene_filter: Optional[str] = None):
            species_counts = Counter()
            validated_counts = Counter()
            
            for read_id, assignment in assignments.items():
                if gene_filter and assignment['gene_type'] != gene_filter:
                    continue
                
                species = assignment['assigned_species']
                species_counts[species] += 1
                
                if assignment.get('blast_validated', False):
                    validated_counts[species] += 1
            
            total = sum(species_counts.values())
            abundance = {}
            for species in species_counts:
                abundance[species] = {
                    'read_count': species_counts[species],
                    'validated_count': validated_counts.get(species, 0),
                    'abundance_percent': (species_counts[species] / total * 100) if total > 0 else 0,
                    'validation_rate': (validated_counts.get(species, 0) / species_counts[species] * 100) if species_counts[species] > 0 else 0
                }
            return abundance
        
        coi_abundance = _count_by_species(species_assignments, 'COI')
        ref_18s_abundance = _count_by_species(species_assignments, '18S')
        combined_abundance = _count_by_species(species_assignments, None)
        
        return coi_abundance, ref_18s_abundance, combined_abundance
    
    def generate_report(self,
                       species_assignments: Dict[str, Dict],
                       gene_assignments: Dict[str, Dict],
                       coi_abundance: Dict,
                       ref_18s_abundance: Dict,
                       combined_abundance: Dict,
                       reads_fastq: str,
                       coi_fasta: Path = None,
                       ref_18s_fasta: Path = None) -> Path:
        """Generate comprehensive report with gene-type breakdown."""
        report_file = self.output_dir / "species_identification_report.txt"
        
        logger.info(f"Generating report: {report_file}")
        
        with open(report_file, 'w') as f:
            f.write("=" * 90 + "\n")
            f.write("TICK SPECIES IDENTIFICATION REPORT - DUAL GENE ANALYSIS (COI + 18S)\n")
            f.write("=" * 90 + "\n\n")
            
            f.write(f"Input File: {reads_fastq}\n")
            f.write(f"Analysis Date: {datetime.now().isoformat()}\n")
            f.write(f"Strict Gene Validation: {self.strict_gene_validation}\n")
            f.write(f"Min Identity Threshold: {self.min_identity}\n")
            f.write(f"Min Coverage Threshold: {self.min_coverage}\n")
            f.write(f"Total Reads Analyzed: {len(species_assignments)}\n")
            
            if self.export_fasta and coi_fasta and ref_18s_fasta:
                f.write(f"\nExported FASTA files:\n")
                f.write(f"  COI reads: {coi_fasta}\n")
                f.write(f"  18S reads: {ref_18s_fasta}\n")
            f.write("\n")
            
            # Gene type breakdown
            f.write("-" * 90 + "\n")
            f.write("GENE TYPE DISTRIBUTION (MUTUALLY EXCLUSIVE)\n")
            f.write("-" * 90 + "\n\n")
            
            gene_stats = Counter(a['gene_type'] for a in species_assignments.values())
            total_reads = len(species_assignments)
            
            f.write(f"{'Gene Type':<20} {'Count':>10} {'Percentage':>12}\n")
            f.write("-" * 90 + "\n")
            f.write(f"{'COI':<20} {gene_stats['COI']:>10} {(gene_stats['COI']/total_reads*100):>11.1f}%\n")
            f.write(f"{'18S rRNA':<20} {gene_stats['18S']:>10} {(gene_stats['18S']/total_reads*100):>11.1f}%\n")
            f.write(f"{'CONFLICT':<20} {gene_stats['CONFLICT']:>10} {(gene_stats['CONFLICT']/total_reads*100):>11.1f}%\n")
            f.write(f"{'UNASSIGNED':<20} {gene_stats['UNASSIGNED']:>10} {(gene_stats['UNASSIGNED']/total_reads*100):>11.1f}%\n")
            
            # Gene-type decision summary
            f.write("\n" + "-" * 90 + "\n")
            f.write("GENE-TYPE DECISION SUMMARY\n")
            f.write("-" * 90 + "\n\n")
            
            decision_reasons = Counter()
            for assignment in species_assignments.values():
                if assignment['gene_type'] in ['COI', '18S']:
                    reason = assignment.get('decision_reason', 'Unknown')
                    # Simplify reason for display
                    if 'Only COI' in reason:
                        decision_reasons['Only COI hit'] += 1
                    elif 'Only 18S' in reason:
                        decision_reasons['Only 18S hit'] += 1
                    elif 'COI score' in reason or 'COI identity' in reason:
                        decision_reasons['COI better match'] += 1
                    elif '18S score' in reason or '18S identity' in reason:
                        decision_reasons['18S better match'] += 1
                    else:
                        decision_reasons[reason] += 1
            
            f.write(f"{'Decision Type':<30} {'Count':>10}\n")
            f.write("-" * 90 + "\n")
            for reason, count in sorted(decision_reasons.items(), key=lambda x: x[1], reverse=True):
                f.write(f"{reason:<30} {count:>10}\n")
            
            # COI results
            if coi_abundance:
                f.write("\n" + "-" * 90 + "\n")
                total_coi_reads = sum(a['read_count'] for a in coi_abundance.values())
                f.write("COI SPECIES ABUNDANCE (n={} reads)\n".format(total_coi_reads))
                f.write("-" * 90 + "\n\n")
                
                sorted_species = sorted(coi_abundance.items(), 
                                       key=lambda x: x[1]['abundance_percent'], 
                                       reverse=True)
                
                f.write(f"{'Species':<45} {'Reads':>8} {'%':>8} {'Validated':>10} {'Val. %':>7}\n")
                f.write("-" * 90 + "\n")
                
                for species, stats in sorted_species:
                    f.write(f"{species:<45} {stats['read_count']:>8} "
                           f"{stats['abundance_percent']:>7.2f}% "
                           f"{stats['validated_count']:>10} "
                           f"{stats['validation_rate']:>6.1f}%\n")
            
            # 18S results
            if ref_18s_abundance:
                f.write("\n" + "-" * 90 + "\n")
                total_18s_reads = sum(a['read_count'] for a in ref_18s_abundance.values())
                f.write("18S rRNA SPECIES ABUNDANCE (n={} reads)\n".format(total_18s_reads))
                f.write("-" * 90 + "\n\n")
                
                sorted_species = sorted(ref_18s_abundance.items(), 
                                       key=lambda x: x[1]['abundance_percent'], 
                                       reverse=True)
                
                f.write(f"{'Species':<45} {'Reads':>8} {'%':>8} {'Validated':>10} {'Val. %':>7}\n")
                f.write("-" * 90 + "\n")
                
                for species, stats in sorted_species:
                    f.write(f"{species:<45} {stats['read_count']:>8} "
                           f"{stats['abundance_percent']:>7.2f}% "
                           f"{stats['validated_count']:>10} "
                           f"{stats['validation_rate']:>6.1f}%\n")
            
            # Combined overview
            if combined_abundance:
                f.write("\n" + "-" * 90 + "\n")
                f.write("COMBINED OVERVIEW (ALL SPECIES, BOTH GENE TYPES)\n")
                f.write("-" * 90 + "\n\n")
                
                sorted_species = sorted(combined_abundance.items(), 
                                       key=lambda x: x[1]['abundance_percent'], 
                                       reverse=True)
                
                f.write(f"{'Species':<45} {'Reads':>8} {'%':>8} {'Validated':>10} {'Val. %':>7}\n")
                f.write("-" * 90 + "\n")
                
                for species, stats in sorted_species:
                    f.write(f"{species:<45} {stats['read_count']:>8} "
                           f"{stats['abundance_percent']:>7.2f}% "
                           f"{stats['validated_count']:>10} "
                           f"{stats['validation_rate']:>6.1f}%\n")
            
            # Quality stats
            f.write("\n" + "-" * 90 + "\n")
            f.write("QUALITY STATISTICS\n")
            f.write("-" * 90 + "\n\n")
            
            coi_alns = [a for a in species_assignments.values() if a['gene_type'] == 'COI']
            ref_18s_alns = [a for a in species_assignments.values() if a['gene_type'] == '18S']
            
            if coi_alns:
                mapq_vals = [a['minimap_mapq'] for a in coi_alns if a['minimap_mapq'] is not None]
                identity_vals = [a['minimap_identity'] for a in coi_alns if a['minimap_identity'] is not None]
                validated = sum(1 for a in coi_alns if a['blast_validated'])
                
                f.write("COI Statistics:\n")
                f.write(f"  Total COI reads: {len(coi_alns)}\n")
                f.write(f"  BLAST validated: {validated} ({validated/len(coi_alns)*100:.1f}%)\n")
                if mapq_vals:
                    f.write(f"  Average MAPQ: {sum(mapq_vals)/len(mapq_vals):.1f}\n")
                if identity_vals:
                    f.write(f"  Average Identity: {sum(identity_vals)/len(identity_vals):.3f}\n")
                    f.write(f"  Min/Max Identity: {min(identity_vals):.3f} / {max(identity_vals):.3f}\n")
            
            if ref_18s_alns:
                mapq_vals = [a['minimap_mapq'] for a in ref_18s_alns if a['minimap_mapq'] is not None]
                identity_vals = [a['minimap_identity'] for a in ref_18s_alns if a['minimap_identity'] is not None]
                validated = sum(1 for a in ref_18s_alns if a['blast_validated'])
                
                f.write("\n18S Statistics:\n")
                f.write(f"  Total 18S reads: {len(ref_18s_alns)}\n")
                f.write(f"  BLAST validated: {validated} ({validated/len(ref_18s_alns)*100:.1f}%)\n")
                if mapq_vals:
                    f.write(f"  Average MAPQ: {sum(mapq_vals)/len(mapq_vals):.1f}\n")
                if identity_vals:
                    f.write(f"  Average Identity: {sum(identity_vals)/len(identity_vals):.3f}\n")
                    f.write(f"  Min/Max Identity: {min(identity_vals):.3f} / {max(identity_vals):.3f}\n")
            
            # Gene-type decision examples
            f.write("\n" + "-" * 90 + "\n")
            f.write("GENE-TYPE DECISION EXAMPLES\n")
            f.write("-" * 90 + "\n\n")
            
            # Show some examples of reads with both COI and 18S hits
            both_hits = [a for a in species_assignments.values() 
                        if a.get('coi_identity', 0) > 0 and a.get('ref_18s_identity', 0) > 0]
            both_hits = both_hits[:20]  # Show first 20
            
            if both_hits:
                f.write(f"{'Read ID':<25} {'Gene':<10} {'COI ID':<8} {'18S ID':<8} {'Diff':<8} {'Species':<30}\n")
                f.write("-" * 90 + "\n")
                for i, assignment in enumerate(both_hits):
                    # Find the read_id for this assignment
                    read_id = None
                    for rid, ass in species_assignments.items():
                        if ass is assignment:
                            read_id = rid
                            break
                    if read_id:
                        f.write(f"{read_id[:25]:<25} {assignment['gene_type']:<10} "
                               f"{assignment.get('coi_identity', 0):<8.3f} "
                               f"{assignment.get('ref_18s_identity', 0):<8.3f} "
                               f"{assignment.get('identity_diff', 0):<8.3f} "
                               f"{assignment['assigned_species'][:30]:<30}\n")
            
            f.write("\n" + "=" * 90 + "\n")
        
        return report_file
    
    def generate_json_results(self,
                             species_assignments: Dict[str, Dict],
                             gene_assignments: Dict[str, Dict],
                             coi_abundance: Dict,
                             ref_18s_abundance: Dict,
                             combined_abundance: Dict,
                             reads_fastq: str) -> Path:
        """Generate JSON results with gene-type separation."""
        json_file = self.output_dir / "species_identification_results.json"
        
        # Clean up assignments for JSON (remove sequence data to keep file size manageable)
        clean_assignments = {}
        for read_id, assignment in species_assignments.items():
            clean_assignment = {k: v for k, v in assignment.items() if k != 'sequence'}
            clean_assignments[read_id] = clean_assignment
        
        results = {
            'input_file': str(reads_fastq),
            'timestamp': datetime.now().isoformat(),
            'parameters': {
                'min_identity': self.min_identity,
                'min_coverage': self.min_coverage,
                'min_query_len': self.min_query_len,
                'strict_gene_validation': self.strict_gene_validation,
                'export_fasta': self.export_fasta
            },
            'gene_type_separation': {
                'COI_species': coi_abundance,
                '18S_species': ref_18s_abundance,
                'combined_overview': combined_abundance
            },
            'gene_type_summary': dict(Counter(a['gene_type'] for a in species_assignments.values())),
            'total_reads': len(species_assignments),
            'read_assignments': clean_assignments
        }
        
        with open(json_file, 'w') as f:
            json.dump(results, f, indent=2)
        
        logger.info(f"JSON results written: {json_file}")
        return json_file
    
    def run_full_pipeline(self, reads_fastq: str, sample_name: Optional[str] = None) -> Dict:
        """Run the complete dual-gene identification pipeline for a single sample."""
        
        # Use sample_name or derive from filename
        if sample_name is None:
            sample_name = Path(reads_fastq).stem
            if sample_name.endswith('.fastq'):
                sample_name = sample_name.replace('.fastq', '')
        
        logger.info("\n" + "="*60)
        logger.info(f"Processing sample: {sample_name}")
        logger.info("="*60)
        
        # Create sample-specific output subdirectory for cleaner organization
        sample_output_dir = self.output_dir / sample_name
        sample_output_dir.mkdir(parents=True, exist_ok=True)
        
        # Temporarily set output_dir to sample-specific directory
        original_output_dir = self.output_dir
        self.output_dir = sample_output_dir
        
        try:
            # Check dependencies
            if not self.check_dependencies():
                return {'status': 'failed', 'error': 'Missing dependencies'}
            
            logger.info("\n" + "="*60)
            logger.info("STEP 1: Indexing references")
            logger.info("="*60)
            coi_index, ref_18s_index = self.index_references()
            
            logger.info("\n" + "="*60)
            logger.info("STEP 2: Aligning reads to COI")
            logger.info("="*60)
            coi_sam = self.align_reads_minimap2(reads_fastq, coi_index, 'COI')
            coi_alignments = self.parse_sam_alignment(coi_sam, 'COI')
            
            logger.info("\n" + "="*60)
            logger.info("STEP 3: Aligning reads to 18S")
            logger.info("="*60)
            ref_18s_sam = self.align_reads_minimap2(reads_fastq, ref_18s_index, '18S')
            ref_18s_alignments = self.parse_sam_alignment(ref_18s_sam, '18S')
            
            logger.info("\n" + "="*60)
            logger.info("STEP 4: Determining gene type (MUTUALLY EXCLUSIVE)")
            logger.info("="*60)
            gene_assignments = self.determine_gene_type(coi_alignments, ref_18s_alignments)
            
            logger.info("\n" + "="*60)
            logger.info("STEP 5: BLAST validation for COI")
            logger.info("="*60)
            coi_blast_output = self.blast_validate(reads_fastq, self.reference_coi, 'COI')
            coi_blast_hits = self.parse_blast_results(coi_blast_output)
            
            logger.info("\n" + "="*60)
            logger.info("STEP 6: BLAST validation for 18S")
            logger.info("="*60)
            ref_18s_blast_output = self.blast_validate(reads_fastq, self.reference_18s, '18S')
            ref_18s_blast_hits = self.parse_blast_results(ref_18s_blast_output)
            
            logger.info("\n" + "="*60)
            logger.info("STEP 7: Parsing reference headers")
            logger.info("="*60)
            coi_species_map = self.parse_reference_headers(self.reference_coi)
            ref_18s_species_map = self.parse_reference_headers(self.reference_18s)
            
            logger.info("\n" + "="*60)
            logger.info("STEP 8: Assigning species")
            logger.info("="*60)
            species_assignments = self.assign_species(
                gene_assignments, coi_alignments, ref_18s_alignments,
                coi_blast_hits, ref_18s_blast_hits,
                coi_species_map, ref_18s_species_map
            )
            
            logger.info("\n" + "="*60)
            logger.info("STEP 9: Exporting FASTA files")
            logger.info("="*60)
            coi_fasta, ref_18s_fasta = self.export_fasta_files(species_assignments, reads_fastq)
            
            logger.info("\n" + "="*60)
            logger.info("STEP 10: Estimating abundance")
            logger.info("="*60)
            coi_abundance, ref_18s_abundance, combined_abundance = self.estimate_abundance(species_assignments)
            
            logger.info("\n" + "="*60)
            logger.info("STEP 11: Generating reports")
            logger.info("="*60)
            report_file = self.generate_report(
                species_assignments, gene_assignments,
                coi_abundance, ref_18s_abundance, combined_abundance,
                reads_fastq, coi_fasta, ref_18s_fasta
            )
            json_file = self.generate_json_results(
                species_assignments, gene_assignments,
                coi_abundance, ref_18s_abundance, combined_abundance,
                reads_fastq
            )
            
            # Clean up temporary files
            self._cleanup_temp_files()
            
            logger.info("\n" + "="*60)
            logger.info(f"Sample {sample_name} complete!")
            logger.info("="*60)
            logger.info(f"Report: {report_file}")
            logger.info(f"JSON Results: {json_file}")
            
            # Return results summary
            return {
                'status': 'success',
                'sample': sample_name,
                'report': str(report_file),
                'json': str(json_file),
                'coi_fasta': str(coi_fasta) if coi_fasta else None,
                'ref_18s_fasta': str(ref_18s_fasta) if ref_18s_fasta else None,
                'gene_type_summary': dict(Counter(a['gene_type'] for a in species_assignments.values())),
                'total_reads': len(species_assignments),
                'coi_abundance': coi_abundance,
                'ref_18s_abundance': ref_18s_abundance,
                'combined_abundance': combined_abundance
            }
            
        except Exception as e:
            logger.error(f"Pipeline failed for {sample_name}: {e}", exc_info=True)
            return {'status': 'failed', 'sample': sample_name, 'error': str(e)}
        finally:
            # Restore original output_dir
            self.output_dir = original_output_dir


def generate_batch_summary(results: List[Dict], output_dir: Path) -> Path:
    """Generate a summary report for batch processing."""
    summary_file = output_dir / "BATCH_SUMMARY_REPORT.txt"
    
    logger.info(f"Generating batch summary report: {summary_file}")
    
    # Filter successful runs
    successful = [r for r in results if r.get('status') == 'success']
    failed = [r for r in results if r.get('status') == 'failed']
    
    with open(summary_file, 'w') as f:
        f.write("=" * 100 + "\n")
        f.write("BATCH TICK SPECIES IDENTIFICATION - SUMMARY REPORT\n")
        f.write("=" * 100 + "\n\n")
        
        f.write(f"Analysis Date: {datetime.now().isoformat()}\n")
        f.write(f"Total Samples Processed: {len(results)}\n")
        f.write(f"Successful: {len(successful)}\n")
        f.write(f"Failed: {len(failed)}\n\n")
        
        if failed:
            f.write("-" * 100 + "\n")
            f.write("FAILED SAMPLES:\n")
            f.write("-" * 100 + "\n")
            for r in failed:
                f.write(f"  {r['sample']}: {r.get('error', 'Unknown error')}\n")
            f.write("\n")
        
        if successful:
            f.write("-" * 100 + "\n")
            f.write("SAMPLE SUMMARIES:\n")
            f.write("-" * 100 + "\n\n")
            
            # Collect all species across samples
            all_species = defaultdict(lambda: defaultdict(int))
            
            for result in successful:
                sample = result['sample']
                f.write(f"\n{'='*90}\n")
                f.write(f"SAMPLE: {sample}\n")
                f.write(f"{'='*90}\n")
                
                f.write(f"Total Reads Analyzed: {result['total_reads']}\n")
                
                # Gene type breakdown
                gene_summary = result['gene_type_summary']
                f.write("\nGene Type Breakdown:\n")
                for gene_type, count in gene_summary.items():
                    pct = (count / result['total_reads'] * 100) if result['total_reads'] > 0 else 0
                    f.write(f"  {gene_type}: {count} ({pct:.1f}%)\n")
                
                # Combined species abundance
                f.write("\nSpecies Abundance (Combined):\n")
                combined = result['combined_abundance']
                for species, stats in sorted(combined.items(), key=lambda x: x[1]['read_count'], reverse=True):
                    f.write(f"  {species}: {stats['read_count']} reads ({stats['abundance_percent']:.1f}%)\n")
                    all_species[species]['total_reads'] += stats['read_count']
                    all_species[species]['samples_present'] += 1
                
                # Output files
                f.write("\nOutput Files:\n")
                f.write(f"  Report: {result['report']}\n")
                f.write(f"  JSON: {result['json']}\n")
                if result.get('coi_fasta'):
                    f.write(f"  COI FASTA: {result['coi_fasta']}\n")
                if result.get('ref_18s_fasta'):
                    f.write(f"  18S FASTA: {result['ref_18s_fasta']}\n")
            
            # Cross-sample summary
            f.write("\n" + "=" * 100 + "\n")
            f.write("CROSS-SAMPLE SPECIES SUMMARY\n")
            f.write("=" * 100 + "\n\n")
            
            f.write(f"{'Species':<40} {'Total Reads':>15} {'Samples':>10}\n")
            f.write("-" * 100 + "\n")
            for species, stats in sorted(all_species.items(), key=lambda x: x[1]['total_reads'], reverse=True):
                f.write(f"{species:<40} {stats['total_reads']:>15} {stats['samples_present']:>10}\n")
            
            f.write("\n" + "=" * 100 + "\n")
            f.write("END OF BATCH SUMMARY\n")
            f.write("=" * 100 + "\n")
    
    logger.info(f"Batch summary report generated: {summary_file}")
    return summary_file


def find_fastq_files(input_dir: str, pattern: str = "*.fastq*") -> List[Path]:
    """Find all FASTQ files in a directory matching the pattern."""
    input_path = Path(input_dir)
    if not input_path.exists():
        logger.error(f"Input directory not found: {input_dir}")
        return []
    
    # Find files matching pattern
    files = list(input_path.glob(pattern))
    
    # Filter out common unwanted files and keep only FASTQ files
    files = [f for f in files if not f.name.startswith('.') and 
             (f.suffix in ['.fastq', '.fq'] or f.name.endswith('.fastq.gz') or f.name.endswith('.fq.gz'))]
    
    logger.info(f"Found {len(files)} FASTQ files in {input_dir}")
    for f in files[:10]:  # Show first 10
        logger.info(f"  {f.name}")
    if len(files) > 10:
        logger.info(f"  ... and {len(files) - 10} more")
    
    return sorted(files)


def main():
    """Main entry point for batch processing."""
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Identify tick species from ONT metagenome samples (COI + 18S) - BATCH MODE',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""USAGE:
  1. Edit CONFIG at top of script (RECOMMENDED)
  2. Run: python tick_species_identifier_dual_gene.py
  
BATCH MODE: Set 'input_dir' in CONFIG to your directory containing FASTQ files.
The script will process all FASTQ files in that directory.
        """
    )
    
    parser.add_argument('-r', '--reference-coi', 
                       help='Path to tick COI reference FASTA')
    parser.add_argument('-s', '--reference-18s',
                       help='Path to tick 18S reference FASTA')
    parser.add_argument('-i', '--input-dir',
                       help='Directory containing input FASTQ files')
    parser.add_argument('-p', '--pattern',
                       help='Pattern to match FASTQ files (e.g., "*.fastq", "*_filtered.fastq")')
    parser.add_argument('-o', '--output',
                       help='Base output directory for results')
    parser.add_argument('-t', '--threads', type=int,
                       help='Number of CPU threads')
    parser.add_argument('--min-identity', type=float,
                       help='Minimum identity threshold 0-1')
    parser.add_argument('--min-coverage', type=float,
                       help='Minimum query coverage 0-1')
    parser.add_argument('--no-strict-validation', action='store_true',
                       help='Disable strict gene validation (allow conflicts)')
    parser.add_argument('--keep-temp', action='store_true',
                       help='Keep temporary files for debugging')
    parser.add_argument('--no-export-fasta', action='store_true',
                       help='Disable FASTA export')
    parser.add_argument('--file-list', nargs='+',
                       help='Process specific files (space-separated list)')
    parser.add_argument('--single-file',
                       help='Process a single FASTQ file')
    
    args = parser.parse_args()
    
    try:
        # Override CONFIG with arguments if provided
        if args.reference_coi:
            CONFIG['reference_coi'] = args.reference_coi
        if args.reference_18s:
            CONFIG['reference_18s'] = args.reference_18s
        if args.input_dir:
            CONFIG['input_dir'] = args.input_dir
        if args.pattern:
            CONFIG['input_pattern'] = args.pattern
        if args.output:
            CONFIG['output_base_dir'] = args.output
        if args.threads:
            CONFIG['threads'] = args.threads
        if args.min_identity:
            CONFIG['min_identity'] = args.min_identity
        if args.min_coverage:
            CONFIG['min_coverage'] = args.min_coverage
        if args.no_strict_validation:
            CONFIG['strict_gene_validation'] = False
        if args.keep_temp:
            CONFIG['keep_temp_files'] = True
        if args.no_export_fasta:
            CONFIG['export_fasta'] = False
        
        # Create base output directory
        output_base_dir = Path(CONFIG['output_base_dir'])
        output_base_dir.mkdir(parents=True, exist_ok=True)
        
        # Determine which files to process
        fastq_files = []
        if args.single_file:
            # Process single file
            fastq_files = [Path(args.single_file)]
            if not fastq_files[0].exists():
                logger.error(f"File not found: {args.single_file}")
                sys.exit(1)
            logger.info(f"Processing single file: {fastq_files[0].name}")
        elif args.file_list:
            # Process specific files
            for f in args.file_list:
                file_path = Path(f)
                if file_path.exists():
                    fastq_files.append(file_path)
                else:
                    logger.warning(f"File not found, skipping: {f}")
            if not fastq_files:
                logger.error("No valid files found in file list")
                sys.exit(1)
        else:
            # Batch mode - find all FASTQ files
            fastq_files = find_fastq_files(CONFIG['input_dir'], CONFIG['input_pattern'])
            if not fastq_files:
                logger.error(f"No FASTQ files found in {CONFIG['input_dir']} matching {CONFIG['input_pattern']}")
                sys.exit(1)
        
        logger.info(f"\n{'='*80}")
        logger.info(f"BATCH TICK SPECIES IDENTIFICATION - PROCESSING {len(fastq_files)} SAMPLES")
        logger.info(f"{'='*80}")
        logger.info(f"COI Reference: {CONFIG['reference_coi']}")
        logger.info(f"18S Reference: {CONFIG['reference_18s']}")
        logger.info(f"Output Base: {output_base_dir}")
        logger.info(f"Strict Gene Validation: {CONFIG['strict_gene_validation']}")
        logger.info(f"Identity threshold: {CONFIG['min_identity']}")
        logger.info(f"Threads: {CONFIG['threads']}")
        logger.info(f"{'='*80}\n")
        
        # Initialize the identifier once (references and indices will be reused)
        identifier = DualGeneTickIdentifier(
            reference_coi=CONFIG['reference_coi'],
            reference_18s=CONFIG['reference_18s'],
            output_dir=str(output_base_dir),  # Base directory, will create subdirs for each sample
            threads=CONFIG['threads'],
            min_query_len=CONFIG['min_query_len'],
            min_identity=CONFIG['min_identity'],
            min_coverage=CONFIG['min_coverage'],
            strict_gene_validation=CONFIG['strict_gene_validation'],
            keep_temp_files=CONFIG['keep_temp_files'],
            export_fasta=CONFIG['export_fasta']
        )
        
        # Check dependencies once
        if not identifier.check_dependencies():
            sys.exit(1)
        
        # Process each file
        all_results = []
        total_success = 0
        total_failed = 0
        
        for i, fastq_file in enumerate(fastq_files, 1):
            sample_name = fastq_file.stem
            if sample_name.endswith('.fastq'):
                sample_name = sample_name.replace('.fastq', '')
            
            logger.info(f"\n{'='*80}")
            logger.info(f"Processing {i}/{len(fastq_files)}: {sample_name}")
            logger.info(f"{'='*80}")
            
            # Create sample-specific output directory
            sample_output_dir = output_base_dir / sample_name
            sample_output_dir.mkdir(parents=True, exist_ok=True)
            
            # Run pipeline for this sample
            result = identifier.run_full_pipeline(str(fastq_file), sample_name)
            all_results.append(result)
            
            if result.get('status') == 'success':
                total_success += 1
            else:
                total_failed += 1
            
            logger.info(f"Progress: {i}/{len(fastq_files)} samples processed")
        
        # Generate batch summary
        if CONFIG.get('generate_summary_report', True) and all_results:
            summary_file = generate_batch_summary(all_results, output_base_dir)
            logger.info(f"\nBatch summary report: {summary_file}")
        
        # Final statistics
        logger.info("\n" + "="*80)
        logger.info("BATCH PROCESSING COMPLETE")
        logger.info("="*80)
        logger.info(f"Total samples: {len(fastq_files)}")
        logger.info(f"Successful: {total_success}")
        logger.info(f"Failed: {total_failed}")
        logger.info(f"Output directory: {output_base_dir}")
        
        if total_failed > 0:
            failed_samples = [r['sample'] for r in all_results if r.get('status') == 'failed']
            logger.info(f"Failed samples: {', '.join(failed_samples)}")
            
    except Exception as e:
        logger.error(f"Batch pipeline failed: {e}", exc_info=True)
        sys.exit(1)


if __name__ == '__main__':
    main()