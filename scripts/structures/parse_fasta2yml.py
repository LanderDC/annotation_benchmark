#!/usr/bin/env python3
"""
Create YAML configuration files for protein sequences from a multifasta file.

This script reads a multifasta file and generates individual YAML files for each
sequence with MSA file paths.
"""

import os
import sys
import argparse
import yaml
from needletail import parse_fastx_file


def create_yaml_for_sequences(fasta_file, output_dir, msa_dir, msa_extension=".a3m"):
    """
    Create YAML files for each sequence in a multifasta file.
    
    Args:
        fasta_file: Path to input multifasta file
        output_dir: Directory where YAML files will be created
        msa_dir: Directory containing MSA files
        msa_extension: File extension for MSA files (default: .a3m)
    """
    # Create the output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Track statistics
    yaml_count = 0
    
    # Read sequences from the multifasta file
    for record in parse_fastx_file(fasta_file):
        seq_id = record.id.split()[0]
        seq_sequence = str(record.seq)
        
        # Construct MSA path
        msa_filename = f"{seq_id}{msa_extension}"
        msa_path = os.path.join(msa_dir, msa_filename)
        
        # Prepare the YAML content
        yaml_content = {
            'version': 1,
            'sequences': [
                {
                    'protein': {
                        'id': 'A',
                        'sequence': seq_sequence,
                        'msa': msa_path
                    }
                }
            ]
        }
        
        # Create the YAML file
        yaml_file_path = os.path.join(output_dir, f"{seq_id}.yaml")
        with open(yaml_file_path, 'w') as yaml_file:
            yaml.dump(yaml_content, yaml_file, default_flow_style=False, sort_keys=False)
        
        yaml_count += 1
    
    print(f"Created {yaml_count} YAML files in {output_dir}")
    print(f"MSA files referenced from: {msa_dir}")


def main():
    parser = argparse.ArgumentParser(
        description="Create YAML configuration files for protein sequences from a multifasta file",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage with MSA directory
  %(prog)s input.fasta output_yamls/ -m /path/to/msas/
  
  # With custom MSA extension
  %(prog)s input.fasta output_yamls/ -m /path/to/msas/ -e .sto
  
  # Using absolute MSA path
  %(prog)s input.fasta output_yamls/ \\
    -m /user/leuven/337/vsc33750/scratch/structure_pred/boltz/msas_long_bm/ \\
    --absolute-path
  
  # Using relative path
  %(prog)s input.fasta output_yamls/ -m ./msas/
        """
    )
    
    parser.add_argument(
        'fasta_file',
        help='Input multifasta file'
    )
    parser.add_argument(
        'output_dir',
        help='Output directory for YAML files'
    )
    parser.add_argument(
        '-m', '--msa-dir',
        required=True,
        help='Directory containing MSA files (required)'
    )
    parser.add_argument(
        '-e', '--msa-extension',
        default='.a3m',
        help='File extension for MSA files (default: .a3m)'
    )
    parser.add_argument(
        '--absolute-path',
        action='store_true',
        help='Convert MSA directory to absolute path'
    )
    
    args = parser.parse_args()
    
    # Validate input file
    if not os.path.exists(args.fasta_file):
        print(f"Error: Input file '{args.fasta_file}' not found", file=sys.stderr)
        sys.exit(1)
    
    # Process MSA directory path
    msa_dir = args.msa_dir
    if args.absolute_path:
        msa_dir = os.path.abspath(msa_dir)
    
    # Create YAML files
    create_yaml_for_sequences(
        fasta_file=args.fasta_file,
        output_dir=args.output_dir,
        msa_dir=msa_dir,
        msa_extension=args.msa_extension
    )


if __name__ == "__main__":
    main()
