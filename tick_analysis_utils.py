#!/usr/bin/env python3
"""
Utility functions and analysis tools for tick species identification.
Includes filtering, visualization helpers, and downstream analysis.

USAGE: Edit the CONFIG sections below, then run the script.
"""

import json
from pathlib import Path
from collections import Counter, defaultdict
import re

# ============================================================================
# OPERATIONS - Uncomment the one you want to run
# ============================================================================\n\n\n\n# OPERATION 1: Filter results by confidence thresholds\n# Uncomment and configure to use:\nOPERATION_FILTER_RESULTS = \"\"\"\nconfig = {\n    'input_json': 'tick_analysis/species_identification_results.json',\n    'output_file': 'tick_analysis/filtered_results.json',\n    'min_mapq': 30,                     # Minimum mapping quality (30=good, 50+=excellent)\n    'min_identity': 0.85,               # Minimum identity 0-1\n    'require_blast_validation': False,  # True = only BLAST-validated reads\n}\n\"\"\"\n\n# OPERATION 2: Extract species statistics\n# Uncomment and configure to use:\nOPERATION_EXTRACT_SPECIES = \"\"\"\nconfig = {\n    'input_json': 'tick_analysis/species_identification_results.json',\n    'output_file': 'tick_analysis/species_stats.json',  # Leave empty to just print\n}\n\"\"\"\n\n# OPERATION 3: Generate CSV report for Excel\n# Uncomment and configure to use:\nOPERATION_GENERATE_CSV = \"\"\"\nconfig = {\n    'input_json': 'tick_analysis/species_identification_results.json',\n    'output_csv': 'tick_analysis/results.csv',  # REQUIRED\n}\n\"\"\"\n\n# OPERATION 4: Compare multiple samples\n# Uncomment and configure to use:\nOPERATION_COMPARE_SAMPLES = \"\"\"\nconfig = {\n    'json_files': [\n        'analysis1/species_identification_results.json',\n        'analysis2/species_identification_results.json',\n        'analysis3/species_identification_results.json',\n    ],\n    'output_file': 'sample_comparison.json',  # Leave empty to just print\n}\n\"\"\"\n\n# OPERATION 5: Identify mixed/chimeric reads\n# Uncomment and configure to use:\nOPERATION_IDENTIFY_MIXED = \"\"\"\nconfig = {\n    'input_json': 'tick_analysis/species_identification_results.json',\n    'output_file': 'tick_analysis/mixed_reads.json',  # Leave empty to just print\n}\n\"\"\""


class TickAnalysisUtils:
    """Utility functions for tick metagenome analysis."""
    
    @staticmethod
    def filter_by_confidence(json_results: str, 
                            min_mapq: int = 30,
                            min_identity: float = 0.85,
                            require_blast_validation: bool = False,
                            output_file: str = None) -> dict:
        """
        Filter species assignments by confidence thresholds.
        
        Args:
            json_results: Path to JSON results file from main pipeline
            min_mapq: Minimum mapping quality
            min_identity: Minimum identity
            require_blast_validation: Only keep BLAST-validated assignments
            output_file: Optional output file to save filtered results
            
        Returns:
            Filtered results dictionary
        """
        with open(json_results, 'r') as f:
            results = json.load(f)
        
        filtered_assignments = {}
        for read_id, assignment in results['read_assignments'].items():
            passes = True
            
            if assignment['minimap_mapq'] < min_mapq:
                passes = False
            if assignment['minimap_identity'] < min_identity:
                passes = False
            if require_blast_validation and not assignment['blast_validated']:
                passes = False
            
            if passes:
                filtered_assignments[read_id] = assignment
        
        filtered_results = {
            **results,
            'read_assignments': filtered_assignments,
            'filters_applied': {
                'min_mapq': min_mapq,
                'min_identity': min_identity,
                'require_blast_validation': require_blast_validation
            }
        }
        
        if output_file:
            with open(output_file, 'w') as f:
                json.dump(filtered_results, f, indent=2)
            print(f"Filtered results saved to: {output_file}")
        
        return filtered_results
    
    @staticmethod
    def extract_species_list(json_results: str) -> dict:
        """Extract list of detected species with confidence metrics."""
        with open(json_results, 'r') as f:
            results = json.load(f)
        
        species_stats = {}
        for read_id, assignment in results['read_assignments'].items():
            species = assignment['assigned_species']
            if species not in species_stats:
                species_stats[species] = {
                    'reads': [],
                    'mapq_scores': [],
                    'identity_scores': [],
                    'blast_validated': 0
                }
            
            species_stats[species]['reads'].append(read_id)
            species_stats[species]['mapq_scores'].append(assignment['minimap_mapq'])
            species_stats[species]['identity_scores'].append(assignment['minimap_identity'])
            if assignment['blast_validated']:
                species_stats[species]['blast_validated'] += 1
        
        # Calculate statistics
        for species in species_stats:
            scores = species_stats[species]
            mapq_list = scores['mapq_scores']
            identity_list = scores['identity_scores']
            
            scores['summary'] = {
                'num_reads': len(scores['reads']),
                'avg_mapq': sum(mapq_list) / len(mapq_list) if mapq_list else 0,
                'avg_identity': sum(identity_list) / len(identity_list) if identity_list else 0,
                'blast_validated': scores['blast_validated'],
                'validation_rate': (scores['blast_validated'] / len(scores['reads']) * 100) if scores['reads'] else 0
            }
            
            del scores['reads']
            del scores['mapq_scores']
            del scores['identity_scores']
        
        return species_stats
    
    @staticmethod
    def parse_fasta_headers(fasta_file: str) -> dict:
        """
        Parse FASTA headers to extract metadata.
        
        Returns:
            Dictionary with sequence ID to metadata mapping
        """
        headers = {}
        with open(fasta_file, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    header = line[1:].strip()
                    seq_id = header.split()[0]
                    
                    # Extract information from header
                    # Expected format: ID_Species_additional_info
                    parts = header.split('_')
                    
                    metadata = {
                        'full_header': header,
                        'id': seq_id,
                    }
                    
                    # Try to extract species (typically 2-3 parts after ID)
                    if len(parts) >= 3:
                        # Genus + species
                        if parts[1][0].isupper():
                            metadata['genus'] = parts[1]
                            metadata['species'] = parts[2] if len(parts) > 2 else ''
                            metadata['species_name'] = f"{parts[1]} {parts[2]}" if len(parts) > 2 else parts[1]
                    
                    headers[seq_id] = metadata
        
        return headers
    
    @staticmethod
    def generate_csv_report(json_results: str, output_csv: str) -> None:
        """
        Generate a CSV report from JSON results for spreadsheet analysis.
        
        Args:
            json_results: Path to JSON results
            output_csv: Path to output CSV file
        """
        import csv
        
        with open(json_results, 'r') as f:
            results = json.load(f)
        
        with open(output_csv, 'w', newline='') as f:
            writer = csv.writer(f)
            
            # Write header
            writer.writerow([
                'read_id',
                'assigned_species',
                'mapq',
                'identity',
                'coverage',
                'blast_validated',
                'blast_species',
                'reference_id'
            ])
            
            # Write data
            for read_id, assignment in results['read_assignments'].items():
                writer.writerow([
                    read_id,
                    assignment['assigned_species'],
                    assignment['minimap_mapq'],
                    f"{assignment['minimap_identity']:.4f}",
                    f"{assignment['minimap_coverage']:.4f}",
                    'Yes' if assignment['blast_validated'] else 'No',
                    assignment.get('blast_species', 'N/A'),
                    assignment['reference_id']
                ])
        
        print(f"CSV report generated: {output_csv}")
    
    @staticmethod
    def compare_samples(json_files: list) -> dict:
        """
        Compare species composition across multiple samples.
        
        Args:
            json_files: List of JSON result files
            
        Returns:
            Comparison dictionary
        """
        comparison = {}
        
        for json_file in json_files:
            sample_name = Path(json_file).stem
            
            with open(json_file, 'r') as f:
                results = json.load(f)
            
            # Extract abundance from species_abundance
            species_data = results.get('species_abundance', {})
            
            comparison[sample_name] = {
                'total_reads': results['total_reads'],
                'validated_reads': results['validated_reads'],
                'species': {}
            }
            
            for species, stats in species_data.items():
                comparison[sample_name]['species'][species] = {
                    'read_count': stats['read_count'],
                    'abundance_percent': stats['abundance_percent'],
                    'validation_rate': stats['validation_rate']
                }
        
        return comparison
    
    @staticmethod
    def identify_mixed_reads(json_results: str, 
                            output_file: str = None) -> dict:
        """
        Identify reads that might contain mixed species sequences.
        Looks for reads with unusually low identity despite high MAPQ.
        
        Args:
            json_results: Path to JSON results
            output_file: Optional file to save identified mixed reads
            
        Returns:
            Dictionary of potentially mixed reads
        """
        with open(json_results, 'r') as f:
            results = json.load(f)
        
        mixed_reads = {}
        
        # Calculate identity distribution statistics
        identities = [a['minimap_identity'] for a in results['read_assignments'].values()]
        mean_identity = sum(identities) / len(identities) if identities else 0
        std_identity = (sum((x - mean_identity) ** 2 for x in identities) / len(identities)) ** 0.5 if identities else 0
        
        # Flag reads with identity > 2 std below mean
        threshold = mean_identity - (2 * std_identity)
        
        for read_id, assignment in results['read_assignments'].items():
            if assignment['minimap_identity'] < threshold and assignment['minimap_mapq'] >= 30:
                mixed_reads[read_id] = {
                    'assigned_species': assignment['assigned_species'],
                    'identity': assignment['minimap_identity'],
                    'mapq': assignment['minimap_mapq'],
                    'reason': 'High MAPQ but unusually low identity (possible mixed sequence)'
                }
        
        if output_file:
            with open(output_file, 'w') as f:
                json.dump(mixed_reads, f, indent=2)
            print(f"Mixed read report saved to: {output_file}")
        
        return mixed_reads


# ============================================================================
# MAIN - Choose which operation to run
# ============================================================================

def main():
    """
    Configuration-based utility operations - supports both CONFIG and command-line.
    
    USAGE:
      1. CONFIG-BASED: Edit config sections above, uncomment operation, run script
      2. COMMAND-LINE: Use command-line arguments instead
      3. MIXED: Config + command-line overrides
    """
    import argparse
    import sys
    
    # Create argument parser for command-line usage
    parser = argparse.ArgumentParser(
        description='Utility tools for tick species metagenome analysis',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
USAGE OPTIONS:

1. CONFIG-BASED (Edit CONFIG sections at top, uncomment operation):
   python tick_analysis_utils.py

2. COMMAND-LINE (Use arguments):
   python tick_analysis_utils.py filter -i results.json --min-mapq 40
   python tick_analysis_utils.py csv -i results.json -o output.csv
   python tick_analysis_utils.py compare results1.json results2.json

3. MIXED (Config + override):
   python tick_analysis_utils.py filter --min-identity 0.90
        """
    )
    
    subparsers = parser.add_subparsers(dest='operation', help='Operation to perform')
    
    # Filter command
    filter_parser = subparsers.add_parser('filter', help='Filter results by confidence')
    filter_parser.add_argument('-i', '--input', required=True, help='Input JSON results file')
    filter_parser.add_argument('-o', '--output', help='Output file')
    filter_parser.add_argument('--min-mapq', type=int, default=30, help='Minimum MAPQ')
    filter_parser.add_argument('--min-identity', type=float, default=0.85, help='Minimum identity')
    filter_parser.add_argument('--require-blast', action='store_true', help='Require BLAST validation')
    
    # Species command
    species_parser = subparsers.add_parser('species', help='Extract species statistics')
    species_parser.add_argument('-i', '--input', required=True, help='Input JSON results file')
    species_parser.add_argument('-o', '--output', help='Output JSON file')
    
    # CSV command
    csv_parser = subparsers.add_parser('csv', help='Generate CSV report')
    csv_parser.add_argument('-i', '--input', required=True, help='Input JSON results file')
    csv_parser.add_argument('-o', '--output', required=True, help='Output CSV file')
    
    # Compare command
    compare_parser = subparsers.add_parser('compare', help='Compare multiple samples')
    compare_parser.add_argument('files', nargs='+', help='Input JSON results files')
    compare_parser.add_argument('-o', '--output', help='Output JSON file')
    
    # Mixed command
    mixed_parser = subparsers.add_parser('mixed', help='Identify mixed sequence reads')
    mixed_parser.add_argument('-i', '--input', required=True, help='Input JSON results file')
    mixed_parser.add_argument('-o', '--output', help='Output file')
    
    args = parser.parse_args()
    
    # Handle command-line operations
    if args.operation == 'filter':
        TickAnalysisUtils.filter_by_confidence(
            args.input,
            min_mapq=args.min_mapq,
            min_identity=args.min_identity,
            require_blast_validation=args.require_blast,
            output_file=args.output
        )
    elif args.operation == 'species':
        stats = TickAnalysisUtils.extract_species_list(args.input)
        print(json.dumps(stats, indent=2))
        if args.output:
            with open(args.output, 'w') as f:
                json.dump(stats, f, indent=2)
            print(f"✓ Species stats saved to: {args.output}")
    elif args.operation == 'csv':
        TickAnalysisUtils.generate_csv_report(args.input, args.output)
        print(f"✓ CSV report saved to: {args.output}")
    elif args.operation == 'compare':
        comparison = TickAnalysisUtils.compare_samples(args.files)
        print(json.dumps(comparison, indent=2))
        if args.output:
            with open(args.output, 'w') as f:
                json.dump(comparison, f, indent=2)
            print(f"✓ Comparison saved to: {args.output}")
    elif args.operation == 'mixed':
        mixed = TickAnalysisUtils.identify_mixed_reads(args.input, args.output)
        print(f"✓ Found {len(mixed)} potentially mixed reads")
        if not args.output:
            print(json.dumps(mixed, indent=2))
    else:
        # No command-line operation specified - show CONFIG-based instructions
        print("\n" + "="*70)
        print("TICK ANALYSIS UTILITIES")
        print("="*70)
        print("\nTWO USAGE METHODS:\n")
        
        print("✓ METHOD 1: CONFIG-BASED (Edit script, uncomment operation)")
        print("  1. Find OPERATION_* sections at top of this file")
        print("  2. Edit the 'config' settings")
        print("  3. Uncomment the operation code in main()")
        print("  4. Run: python tick_analysis_utils.py\n")
        
        print("✓ METHOD 2: COMMAND-LINE (Easier for scripting)")
        print("  # Filter results")
        print("  python tick_analysis_utils.py filter -i results.json --min-mapq 40")
        print("")
        print("  # Generate CSV")
        print("  python tick_analysis_utils.py csv -i results.json -o output.csv")
        print("")
        print("  # Compare samples")
        print("  python tick_analysis_utils.py compare sample1.json sample2.json -o comparison.json")
        print("")
        print("  # Extract species stats")
        print("  python tick_analysis_utils.py species -i results.json -o stats.json")
        print("")
        print("  # Find mixed reads")
        print("  python tick_analysis_utils.py mixed -i results.json -o mixed.json")
        print("\nFor help with any command:")
        print("  python tick_analysis_utils.py <operation> --help")
        print("="*70 + "\n")
    print("  • OPERATION_COMPARE_SAMPLES - Compare multiple samples")
    print("  • OPERATION_IDENTIFY_MIXED - Find mixed/chimeric reads")
    print("="*60 + "\n")


if __name__ == '__main__':
    main()
