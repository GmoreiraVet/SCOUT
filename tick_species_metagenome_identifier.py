#!/usr/bin/env python3
"""
Tick Species Identification from ONT Metagenome Samples

This script identifies tick species from Oxford Nanopore Technology (ONT) 
metagenome samples using a two-pronged approach:
1. Minimap2 alignment for fast, sensitive mapping
2. BLAST for high-confidence species validation

Pipeline:
    1. Index reference database (if needed)
    2. Map ONT reads to reference with minimap2
    3. Validate and refine with BLAST
    4. Assign species based on best hits
    5. Generate abundance estimates and reports
"""

import subprocess
import os
import sys
from pathlib import Path
from collections import defaultdict, Counter
import json
from typing import Dict, List, Tuple, Set
import logging

# ============================================================================
# CONFIGURATION - EDIT THESE SETTINGS
# ============================================================================

CONFIG = {
    # Input/Output paths
    'reference_fasta': 'Tick_Specifier_Ticks.fasta',  # Reference tick COI database
    'input_fastq': 'your_metagenome_reads.fastq',     # Your ONT metagenome file
    'output_dir': './tick_analysis',                   # Where to save results
    
    # Processing parameters
    'threads': 8,                       # Number of CPU threads to use
    'min_query_len': 500,               # Minimum read length (bp) - filter shorter reads
    'min_identity': 0.85,               # Minimum identity 0-1 (0.95=strict, 0.85=balanced, 0.80=sensitive)
    'min_coverage': 0.80,               # Minimum query coverage 0-1
    
    # BLAST parameters
    'blast_evalue': '1e-20',            # BLAST E-value threshold
    'blast_word_size': 11,              # BLAST word size (11=sensitive, 15=faster)
    'blast_num_alignments': 5,          # Number of BLAST hits to keep
    
    # Analysis modes
    'force_reindex': False,             # Force minimap2 re-indexing
}

# ============================================================================
# PREDEFINED CONFIGURATIONS (uncomment to use)
# ============================================================================

# High-confidence mode (strict quality control)
# CONFIG.update({
#     'min_identity': 0.95,
#     'min_coverage': 0.90,
#     'threads': 8,
# })

# Sensitive mode (detect novel/divergent species)
# CONFIG.update({
#     'min_identity': 0.80,
#     'min_coverage': 0.75,
#     'blast_word_size': 11,
# })

# Fast mode (fewer threads, less memory)
# CONFIG.update({
#     'threads': 2,
#     'blast_word_size': 15,
# })

# ============================================================================

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class TickSpeciesIdentifier:
    """Main class for tick species identification from metagenome data."""
    
    def __init__(self, 
                 reference_fasta: str,
                 output_dir: str = "./tick_analysis",
                 threads: int = 4,
                 min_query_len: int = 500,
                 min_identity: float = 0.85,
                 min_coverage: float = 0.80):
        """
        Initialize the identifier.
        
        Args:
            reference_fasta: Path to tick COI reference database
            output_dir: Output directory for results
            threads: Number of threads for tools
            min_query_len: Minimum query length to keep
            min_identity: Minimum identity threshold (0-1)
            min_coverage: Minimum query coverage threshold (0-1)
        """
        self.reference_fasta = Path(reference_fasta)
        self.output_dir = Path(output_dir)
        self.threads = threads
        self.min_query_len = min_query_len
        self.min_identity = min_identity
        self.min_coverage = min_coverage
        
        # Create output directory
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Validate reference exists
        if not self.reference_fasta.exists():
            raise FileNotFoundError(f"Reference database not found: {self.reference_fasta}")
        
        logger.info(f"Initialized TickSpeciesIdentifier with reference: {self.reference_fasta}")
    
    def check_dependencies(self) -> bool:
        """Check if required tools are installed."""
        tools = ['minimap2', 'blastn']
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
    
    def index_reference(self, force: bool = False) -> Path:
        """
        Create minimap2 index for reference database.
        
        Args:
            force: Force reindexing even if index exists
            
        Returns:
            Path to the index file
        """
        index_file = self.output_dir / "reference.mmi"
        
        if index_file.exists() and not force:
            logger.info(f"Index already exists: {index_file}")
            return index_file
        
        logger.info("Creating minimap2 index...")
        cmd = [
            'minimap2',
            '-d', str(index_file),
            str(self.reference_fasta)
        ]
        
        try:
            subprocess.run(cmd, check=True, capture_output=True)
            logger.info(f"Index created: {index_file}")
            return index_file
        except subprocess.CalledProcessError as e:
            logger.error(f"Indexing failed: {e}")
            raise
    
    def align_reads_minimap2(self, 
                            reads_fastq: str,
                            index_file: Path) -> Path:
        """
        Align ONT reads to reference using minimap2.
        
        Args:
            reads_fastq: Path to input FASTQ file
            index_file: Path to minimap2 index
            
        Returns:
            Path to output PAF file
        """
        reads_path = Path(reads_fastq)
        if not reads_path.exists():
            raise FileNotFoundError(f"Reads file not found: {reads_fastq}")
        
        output_paf = self.output_dir / f"{reads_path.stem}.paf"
        
        logger.info(f"Aligning reads with minimap2: {reads_fastq}")
        # Map-ONT preset optimized for Oxford Nanopore data
        cmd = [
            'minimap2',
            '-a',  # Output SAM format instead of PAF
            '-x', 'map-ont',  # ONT preset
            '-t', str(self.threads),
            '--secondary=no',  # Skip secondary alignments
            str(self.reference_fasta),
            str(reads_fastq)
        ]
        
        sam_file = self.output_dir / f"{reads_path.stem}.sam"
        
        try:
            with open(sam_file, 'w') as f:
                subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE, check=True)
            logger.info(f"Alignment complete: {sam_file}")
            return sam_file
        except subprocess.CalledProcessError as e:
            logger.error(f"Alignment failed: {e.stderr.decode()}")
            raise
    
    def parse_sam_alignment(self, sam_file: Path) -> Dict[str, List[Dict]]:
        """
        Parse SAM file and extract high-quality alignments.
        
        Args:
            sam_file: Path to SAM alignment file
            
        Returns:
            Dictionary mapping read IDs to list of alignments
        """
        alignments = defaultdict(list)
        
        logger.info("Parsing SAM alignments...")
        with open(sam_file, 'r') as f:
            for line in f:
                if line.startswith('@'):  # Skip header
                    continue
                
                fields = line.strip().split('\t')
                if len(fields) < 11:
                    continue
                
                read_id = fields[0]
                ref_id = fields[2]
                mapping_quality = int(fields[4])
                seq_len = len(fields[9])
                
                # Parse CIGAR to get alignment stats
                cigar = fields[5]
                matched = sum(int(x) for x in __import__('re').findall(r'(\d+)M', cigar))
                alignment_len = matched
                identity = alignment_len / seq_len if seq_len > 0 else 0
                coverage = alignment_len / seq_len if seq_len > 0 else 0
                
                # Apply filters
                if (seq_len >= self.min_query_len and 
                    identity >= self.min_identity and
                    coverage >= self.min_coverage and
                    mapping_quality >= 30):  # MAPQ >= 30
                    
                    alignments[read_id].append({
                        'read_id': read_id,
                        'ref_id': ref_id,
                        'mapq': mapping_quality,
                        'read_len': seq_len,
                        'matched_len': alignment_len,
                        'identity': identity,
                        'coverage': coverage
                    })
        
        logger.info(f"Found {len(alignments)} reads with high-quality alignments")
        return alignments
    
    def parse_reference_headers(self) -> Dict[str, str]:
        """
        Parse reference FASTA headers to extract species information.
        
        Returns:
            Dictionary mapping reference IDs to species names
        """
        species_map = {}
        
        with open(self.reference_fasta, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    header = line[1:].strip()
                    # Extract species name from header
                    # Format: >ID_Species_Name_other_fields
                    parts = header.split('_')
                    if len(parts) >= 2:
                        ref_id = header.split()[0]  # First word is ID
                        # Try to extract genus and species
                        species_parts = []
                        for part in parts[1:]:
                            # Stop at likely metadata (dates, codes)
                            if len(part) > 3 or part[0].islower():
                                if not any(c.isdigit() for c in part[:2]):
                                    species_parts.append(part)
                                else:
                                    break
                        species_name = ' '.join(species_parts[:2]) if len(species_parts) >= 2 else header
                        species_map[ref_id] = species_name
        
        return species_map
    
    def blast_validate(self, reads_fastq: str) -> Path:
        """
        Run BLAST validation on reads for high-confidence species assignment.
        
        Args:
            reads_fastq: Path to input FASTQ file
            
        Returns:
            Path to BLAST output
        """
        reads_path = Path(reads_fastq)
        blast_db = self.output_dir / "reference_blast"
        blast_output = self.output_dir / f"{reads_path.stem}_blast.tsv"
        
        # Create BLAST database from reference if needed
        if not (self.output_dir / f"{blast_db}.nhr").exists():
            logger.info("Creating BLAST database...")
            cmd = [
                'makeblastdb',
                '-in', str(self.reference_fasta),
                '-dbtype', 'nucl',
                '-out', str(blast_db)
            ]
            try:
                subprocess.run(cmd, check=True, capture_output=True)
                logger.info("BLAST database created")
            except subprocess.CalledProcessError as e:
                logger.error(f"BLAST database creation failed: {e}")
                raise
        
        # Convert FASTQ to FASTA for BLAST
        fasta_file = self.output_dir / f"{reads_path.stem}.fasta"
        self._fastq_to_fasta(reads_fastq, fasta_file)
        
        logger.info("Running BLAST validation...")
        cmd = [
            'blastn',
            '-query', str(fasta_file),
            '-db', str(blast_db),
            '-out', str(blast_output),
            '-outfmt', '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore',
            '-num_threads', str(self.threads),
            '-evalue', '1e-20',
            '-word_size', '11',
            '-num_alignments', '5'  # Get top 5 hits per query
        ]
        
        try:
            subprocess.run(cmd, check=True, capture_output=True)
            logger.info(f"BLAST complete: {blast_output}")
            return blast_output
        except subprocess.CalledProcessError as e:
            logger.error(f"BLAST failed: {e}")
            raise
    
    def _fastq_to_fasta(self, fastq_file: str, fasta_file: Path) -> None:
        """Convert FASTQ to FASTA format."""
        with open(fastq_file, 'r') as infile, open(fasta_file, 'w') as outfile:
            while True:
                header = infile.readline()
                if not header:
                    break
                seq = infile.readline()
                infile.readline()  # Skip +
                infile.readline()  # Skip quality scores
                
                outfile.write(f">{header[1:]}{seq}")
    
    def parse_blast_results(self, blast_output: Path) -> Dict[str, List[Dict]]:
        """
        Parse BLAST output file.
        
        Args:
            blast_output: Path to BLAST TSV output
            
        Returns:
            Dictionary mapping read IDs to list of BLAST hits
        """
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
                      minimap_alignments: Dict,
                      blast_hits: Dict,
                      species_map: Dict[str, str]) -> Dict[str, Dict]:
        """
        Assign species based on combined minimap2 and BLAST results.
        
        Args:
            minimap_alignments: Minimap2 alignment results
            blast_hits: BLAST validation results
            species_map: Reference ID to species name mapping
            
        Returns:
            Dictionary with species assignments
        """
        species_assignments = {}
        
        logger.info("Assigning species to reads...")
        
        for read_id, minimap_alns in minimap_alignments.items():
            if not minimap_alns:
                continue
            
            # Get best minimap2 hit
            best_minimap = max(minimap_alns, 
                              key=lambda x: (x['mapq'], x['identity']))
            best_ref = best_minimap['ref_id'].split()[0]
            
            # Check BLAST validation if available
            blast_valid = False
            blast_species = None
            if read_id in blast_hits and blast_hits[read_id]:
                best_blast = blast_hits[read_id][0]
                blast_ref = best_blast['subject'].split()[0]
                # BLAST validates if same or very similar species
                if blast_ref == best_ref or best_blast['pident'] >= 99.0:
                    blast_valid = True
                    blast_species = species_map.get(blast_ref, blast_ref)
            
            species = species_map.get(best_ref, best_ref)
            
            species_assignments[read_id] = {
                'assigned_species': species,
                'minimap_mapq': best_minimap['mapq'],
                'minimap_identity': best_minimap['identity'],
                'minimap_coverage': best_minimap['coverage'],
                'blast_validated': blast_valid,
                'blast_species': blast_species,
                'reference_id': best_ref
            }
        
        return species_assignments
    
    def estimate_abundance(self, species_assignments: Dict[str, Dict]) -> Dict[str, Dict]:
        """
        Estimate species abundance from read counts.
        
        Args:
            species_assignments: Species assignment dictionary
            
        Returns:
            Dictionary with abundance statistics
        """
        species_counts = Counter()
        validated_counts = Counter()
        
        for read_id, assignment in species_assignments.items():
            species = assignment['assigned_species']
            species_counts[species] += 1
            
            if assignment['blast_validated']:
                validated_counts[species] += 1
        
        total_reads = sum(species_counts.values())
        
        abundance = {}
        for species in species_counts:
            abundance[species] = {
                'read_count': species_counts[species],
                'validated_count': validated_counts.get(species, 0),
                'abundance_percent': (species_counts[species] / total_reads * 100) if total_reads > 0 else 0,
                'validation_rate': (validated_counts.get(species, 0) / species_counts[species] * 100) if species_counts[species] > 0 else 0
            }
        
        return abundance
    
    def generate_report(self,
                       species_assignments: Dict[str, Dict],
                       abundance: Dict[str, Dict],
                       reads_fastq: str) -> Path:
        """
        Generate comprehensive analysis report.
        
        Args:
            species_assignments: Species assignments for all reads
            abundance: Abundance estimates
            reads_fastq: Input reads filename
            
        Returns:
            Path to report file
        """
        report_file = self.output_dir / "species_identification_report.txt"
        
        logger.info(f"Generating report: {report_file}")
        
        with open(report_file, 'w') as f:
            f.write("=" * 80 + "\n")
            f.write("TICK SPECIES IDENTIFICATION REPORT FROM ONT METAGENOME\n")
            f.write("=" * 80 + "\n\n")
            
            f.write(f"Input File: {reads_fastq}\n")
            f.write(f"Analysis Date: {__import__('datetime').datetime.now().isoformat()}\n")
            f.write(f"Total Reads Analyzed: {len(species_assignments)}\n\n")
            
            f.write("-" * 80 + "\n")
            f.write("SPECIES ABUNDANCE\n")
            f.write("-" * 80 + "\n\n")
            
            # Sort by abundance
            sorted_species = sorted(abundance.items(), 
                                   key=lambda x: x[1]['abundance_percent'], 
                                   reverse=True)
            
            f.write(f"{'Species':<40} {'Reads':>8} {'%':>8} {'Validated':>10} {'Val. %':>7}\n")
            f.write("-" * 80 + "\n")
            
            for species, stats in sorted_species:
                f.write(f"{species:<40} {stats['read_count']:>8} "
                       f"{stats['abundance_percent']:>7.2f}% "
                       f"{stats['validated_count']:>10} "
                       f"{stats['validation_rate']:>6.1f}%\n")
            
            f.write("\n" + "-" * 80 + "\n")
            f.write("QUALITY STATISTICS\n")
            f.write("-" * 80 + "\n\n")
            
            mapq_values = [a['minimap_mapq'] for a in species_assignments.values()]
            identity_values = [a['minimap_identity'] for a in species_assignments.values()]
            validated = sum(1 for a in species_assignments.values() if a['blast_validated'])
            
            f.write(f"Total Reads with Alignments: {len(species_assignments)}\n")
            f.write(f"Reads with BLAST Validation: {validated} ({validated/len(species_assignments)*100:.1f}%)\n")
            f.write(f"Average MAPQ: {sum(mapq_values)/len(mapq_values):.1f}\n")
            f.write(f"Average Identity: {sum(identity_values)/len(identity_values):.3f}\n")
            f.write(f"Min/Max Identity: {min(identity_values):.3f} / {max(identity_values):.3f}\n\n")
            
            f.write("=" * 80 + "\n")
        
        return report_file
    
    def generate_json_results(self,
                             species_assignments: Dict[str, Dict],
                             abundance: Dict[str, Dict],
                             reads_fastq: str) -> Path:
        """
        Generate JSON format results for downstream analysis.
        
        Args:
            species_assignments: Species assignments
            abundance: Abundance estimates
            reads_fastq: Input reads filename
            
        Returns:
            Path to JSON results file
        """
        json_file = self.output_dir / "species_identification_results.json"
        
        results = {
            'input_file': str(reads_fastq),
            'timestamp': __import__('datetime').datetime.now().isoformat(),
            'species_abundance': abundance,
            'total_reads': len(species_assignments),
            'validated_reads': sum(1 for a in species_assignments.values() if a['blast_validated']),
            'read_assignments': species_assignments
        }
        
        with open(json_file, 'w') as f:
            json.dump(results, f, indent=2)
        
        logger.info(f"JSON results written: {json_file}")
        return json_file
    
    def run_full_pipeline(self, reads_fastq: str) -> Tuple[Path, Path]:
        """
        Run the complete identification pipeline.
        
        Args:
            reads_fastq: Path to input ONT FASTQ file
            
        Returns:
            Tuple of (report_file, json_file)
        """
        # Check dependencies
        if not self.check_dependencies():
            sys.exit(1)
        
        # Create index
        index_file = self.index_reference()
        
        # Align reads
        sam_file = self.align_reads_minimap2(reads_fastq, index_file)
        minimap_alignments = self.parse_sam_alignment(sam_file)
        
        # BLAST validation
        blast_output = self.blast_validate(reads_fastq)
        blast_hits = self.parse_blast_results(blast_output)
        
        # Parse reference headers for species names
        species_map = self.parse_reference_headers()
        
        # Assign species
        species_assignments = self.assign_species(minimap_alignments, blast_hits, species_map)
        
        # Estimate abundance
        abundance = self.estimate_abundance(species_assignments)
        
        # Generate reports
        report_file = self.generate_report(species_assignments, abundance, reads_fastq)
        json_file = self.generate_json_results(species_assignments, abundance, reads_fastq)
        
        logger.info("Pipeline complete!")
        logger.info(f"Report: {report_file}")
        logger.info(f"JSON Results: {json_file}")
        
        return report_file, json_file


def main():
    """Main entry point - supports both CONFIG file and command-line arguments."""
    import argparse
    
    # Create argument parser
    parser = argparse.ArgumentParser(
        description='Identify tick species from ONT metagenome samples',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""USAGE OPTIONS:
  1. CONFIG-BASED (Edit CONFIG dict at top of script - RECOMMENDED):
     python tick_species_metagenome_identifier.py
  
  2. COMMAND-LINE (Override config with arguments):
     python tick_species_metagenome_identifier.py -r ref.fasta -i reads.fastq -t 8
  
  3. MIXED (Config file + command-line overrides):
     python tick_species_metagenome_identifier.py --min-identity 0.95 --threads 16

EXAMPLES:
  # Edit config and run (recommended for most users)
  nano tick_species_metagenome_identifier.py  # Edit CONFIG section
  python tick_species_metagenome_identifier.py
  
  # Command-line only
  python tick_species_metagenome_identifier.py \\
    -r Tick_Specifier_Ticks.fasta \\
    -i sample.fastq \\
    -o results \\
    --min-identity 0.90 \\
    -t 8
        """
    )
    
    parser.add_argument('-r', '--reference', 
                       help='Path to tick COI reference FASTA database')
    parser.add_argument('-i', '--input',
                       help='Path to input ONT FASTQ file')
    parser.add_argument('-o', '--output',
                       help='Output directory for results')
    parser.add_argument('-t', '--threads', type=int,
                       help='Number of CPU threads to use')
    parser.add_argument('--min-query-len', type=int,
                       help='Minimum query length to keep (bp)')
    parser.add_argument('--min-identity', type=float,
                       help='Minimum identity threshold 0-1')
    parser.add_argument('--min-coverage', type=float,
                       help='Minimum query coverage 0-1')
    parser.add_argument('--blast-evalue',
                       help='BLAST E-value threshold')
    parser.add_argument('--blast-word-size', type=int,
                       help='BLAST word size (11=sensitive, 15=faster)')
    parser.add_argument('--force-reindex', action='store_true',
                       help='Force minimap2 re-indexing')
    
    args = parser.parse_args()
    
    try:
        # Override CONFIG with command-line arguments if provided
        if args.reference:
            CONFIG['reference_fasta'] = args.reference
        if args.input:
            CONFIG['input_fastq'] = args.input
        if args.output:
            CONFIG['output_dir'] = args.output
        if args.threads:
            CONFIG['threads'] = args.threads
        if args.min_query_len:
            CONFIG['min_query_len'] = args.min_query_len
        if args.min_identity:
            CONFIG['min_identity'] = args.min_identity
        if args.min_coverage:
            CONFIG['min_coverage'] = args.min_coverage
        if args.blast_evalue:
            CONFIG['blast_evalue'] = args.blast_evalue
        if args.blast_word_size:
            CONFIG['blast_word_size'] = args.blast_word_size
        if args.force_reindex:
            CONFIG['force_reindex'] = args.force_reindex
        
        # Validate configuration
        if not CONFIG['reference_fasta']:
            raise ValueError("ERROR: Set 'reference_fasta' in CONFIG or use -r argument")
        if not CONFIG['input_fastq']:
            raise ValueError("ERROR: Set 'input_fastq' in CONFIG or use -i argument")
        
        logger.info(f"\n{'='*60}")
        logger.info(f"TICK SPECIES IDENTIFICATION - CONFIG SUMMARY")
        logger.info(f"{'='*60}")
        logger.info(f"Reference: {CONFIG['reference_fasta']}")
        logger.info(f"Input: {CONFIG['input_fastq']}")
        logger.info(f"Output: {CONFIG['output_dir']}")
        logger.info(f"Identity threshold: {CONFIG['min_identity']}")
        logger.info(f"Coverage threshold: {CONFIG['min_coverage']}")
        logger.info(f"Threads: {CONFIG['threads']}")
        logger.info(f"{'='*60}\n")
        
        identifier = TickSpeciesIdentifier(
            reference_fasta=CONFIG['reference_fasta'],
            output_dir=CONFIG['output_dir'],
            threads=CONFIG['threads'],
            min_query_len=CONFIG['min_query_len'],
            min_identity=CONFIG['min_identity'],
            min_coverage=CONFIG['min_coverage']
        )
        
        report_file, json_file = identifier.run_full_pipeline(CONFIG['input_fastq'])
        
        print(f"\n✓ Analysis complete!")
        print(f"  Report: {report_file}")
        print(f"  Results: {json_file}")
        print(f"  View report: cat {report_file}")
        
    except Exception as e:
        logger.error(f"Pipeline failed: {e}", exc_info=True)
        sys.exit(1)


if __name__ == '__main__':
    main()
