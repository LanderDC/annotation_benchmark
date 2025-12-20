#!/usr/bin/env python3
"""
Robust script to download all virus proteins and their metadata based on
GenBank nucleotide accessions from NCBI.

This script handles:
- Rate limiting and API throttling
- Network errors with exponential backoff retry
- Batch processing for large datasets
- Progress tracking
- Email notifications to NCBI (required)
- Proper parsing of CDS features with protein_id and translations
- Comprehensive protein metadata extraction and export to JSON
"""

import argparse
import json
import logging
import sys
import time
from datetime import datetime
from typing import Any, Dict, Iterator, List

from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[
        logging.FileHandler(
            f'protein_download_{datetime.now().strftime("%Y%m%d_%H%M%S")}.log'
        ),
        logging.StreamHandler(sys.stdout),
    ],
)
logger = logging.getLogger(__name__)

# NCBI rate limiting: max 3 requests per second without API key, 10 with key
DELAY_BETWEEN_REQUESTS = 0.4  # Conservative delay (2.5 req/s)
MAX_RETRIES = 5
BATCH_SIZE = 200  # NCBI recommends batches for large queries


def retry_with_backoff(func, max_retries=MAX_RETRIES):
    """Retry a function with exponential backoff on failure."""
    for attempt in range(max_retries):
        try:
            return func()
        except Exception as e:
            wait_time = 2**attempt
            if attempt < max_retries - 1:
                logger.warning(
                    f"Attempt {attempt + 1} failed: {e}. Retrying in {wait_time}s..."
                )
                time.sleep(wait_time)
            else:
                logger.error(f"Failed after {max_retries} attempts: {e}")
                raise


def fetch_genbank_record(accession: str, email: str) -> SeqIO.SeqRecord:
    """
    Fetch a single GenBank nucleotide record with retry logic.

    Args:
        accession: GenBank nucleotide accession number
        email: Your email address (required by NCBI)

    Returns:
        SeqRecord object
    """
    Entrez.email = email

    def _fetch():
        logger.info(f"Fetching GenBank record: {accession}")
        handle = Entrez.efetch(
            db="nucleotide", id=accession, rettype="gb", retmode="text"
        )
        record = SeqIO.read(handle, "genbank")
        handle.close()
        time.sleep(DELAY_BETWEEN_REQUESTS)
        return record

    return retry_with_backoff(_fetch)


def fetch_genbank_records_batch(
    accessions: List[str], email: str
) -> List[SeqIO.SeqRecord]:
    """
    Fetch multiple GenBank records in a single batch request.

    Args:
        accessions: List of GenBank nucleotide accession numbers
        email: Your email address (required by NCBI)

    Returns:
        List of SeqRecord objects
    """
    Entrez.email = email

    def _fetch():
        logger.info(f"Fetching {len(accessions)} GenBank records in batch")
        handle = Entrez.efetch(
            db="nucleotide", id=",".join(accessions), rettype="gb", retmode="text"
        )
        records = list(SeqIO.parse(handle, "genbank"))
        handle.close()
        time.sleep(DELAY_BETWEEN_REQUESTS)
        return records

    return retry_with_backoff(_fetch)


def extract_protein_metadata(feature, record: SeqIO.SeqRecord) -> Dict[str, Any]:
    """
    Extract comprehensive metadata from a CDS feature.

    Args:
        feature: SeqFeature object (CDS)
        record: Parent GenBank SeqRecord

    Returns:
        Dictionary containing protein metadata
    """
    metadata = {
        "protein_id": feature.qualifiers.get("protein_id", [""])[0],
        "gene_name": feature.qualifiers.get("gene", [""])[0],
        "product": feature.qualifiers.get("product", ["hypothetical protein"])[0],
        "locus_tag": feature.qualifiers.get("locus_tag", [""])[0],
        "note": feature.qualifiers.get("note", []),  # Keep as list
        "function": feature.qualifiers.get("function", [""])[0],
        "ec_number": feature.qualifiers.get("EC_number", []),  # Keep as list
        "db_xref": feature.qualifiers.get("db_xref", []),  # Keep as list
        "protein_length": 0,  # Will be filled later
        "codon_start": int(feature.qualifiers.get("codon_start", ["1"])[0]),
        "transl_table": int(feature.qualifiers.get("transl_table", ["1"])[0]),
        "nucleotide_accession": record.id,
        "nucleotide_definition": record.description,
        "organism": record.annotations.get("organism", ""),
        "taxonomy": record.annotations.get("taxonomy", []),  # Keep as list
        "source": record.annotations.get("source", ""),
        "location": str(feature.location),
        "location_start": int(feature.location.start),
        "location_end": int(feature.location.end),
        "strand": (
            "+"
            if feature.location.strand == 1
            else "-" if feature.location.strand == -1 else "."
        ),
    }

    return metadata


def extract_proteins_from_genbank(
    record: SeqIO.SeqRecord,
) -> Iterator[tuple[SeqRecord, Dict[str, Any]]]:
    """
    Extract all protein sequences and metadata from a GenBank record's CDS features.

    Args:
        record: GenBank SeqRecord object

    Yields:
        Tuple of (SeqRecord object for protein, metadata dictionary)
    """
    protein_count = 0

    for feature in record.features:
        if feature.type == "CDS":
            # Extract metadata first
            metadata = extract_protein_metadata(feature, record)

            # Extract protein_id
            protein_id = metadata["protein_id"]
            if not protein_id:
                protein_id = f"{record.id}_unknown_{protein_count}"
                metadata["protein_id"] = protein_id

            # Extract translation if available
            if "translation" in feature.qualifiers:
                protein_seq = Seq(feature.qualifiers["translation"][0])
            else:
                # If no translation provided, translate the CDS ourselves
                try:
                    cds_seq = feature.extract(record.seq)
                    # Use transl_table if specified
                    table = int(metadata["transl_table"])
                    protein_seq = cds_seq.translate(table=table, cds=True)
                except Exception as e:
                    logger.warning(f"Could not translate CDS in {record.id}: {e}")
                    continue

            # Update protein length in metadata
            metadata["protein_length"] = len(protein_seq)
            metadata["protein_sequence"] = str(protein_seq)  # Add sequence to metadata

            # Create comprehensive description for FASTA header
            description_parts = []

            # Add product (protein name)
            if metadata["product"]:
                description_parts.append(metadata["product"])

            # Add gene name if different from product
            # if metadata['gene_name'] and metadata['gene_name'] not in metadata['product']:
            #    description_parts.append(f"[gene={metadata['gene_name']}]")
            #
            ## Add organism
            # if metadata['organism']:
            #    description_parts.append(f"[organism={metadata['organism']}]")
            #
            ## Add source nucleotide accession
            # description_parts.append(f"[nucleotide={record.id}]")
            #
            ## Add location
            # description_parts.append(f"[location={metadata['location']}]")

            description = " ".join(description_parts)

            # Create protein SeqRecord
            protein_record = SeqRecord(
                protein_seq,
                id=protein_id,
                description=description,
                name=metadata["gene_name"] if metadata["gene_name"] else protein_id,
            )

            # Add annotations to SeqRecord
            protein_record.annotations["product"] = metadata["product"]
            protein_record.annotations["gene"] = metadata["gene_name"]
            protein_record.annotations["organism"] = metadata["organism"]
            protein_record.annotations["source_nucleotide"] = record.id

            protein_count += 1
            yield protein_record, metadata

    logger.info(f"Extracted {protein_count} proteins from {record.id}")


def download_virus_proteins(
    accessions: List[str],
    email: str,
    output_fasta: str,
    output_json: str = None,
    batch_mode: bool = True,
    include_sequence_in_json: bool = True,
):
    """
    Download all virus proteins and metadata from a list of GenBank nucleotide accessions.

    Args:
        accessions: List of GenBank nucleotide accession numbers
        email: Your email address (required by NCBI)
        output_fasta: Path to output FASTA file
        output_json: Path to output JSON file (optional)
        batch_mode: Use batch downloading (faster for many accessions)
        include_sequence_in_json: Include protein sequences in JSON output
    """
    logger.info(f"Starting protein download for {len(accessions)} accessions")
    logger.info(f"Output FASTA: {output_fasta}")
    if output_json:
        logger.info(f"Output JSON: {output_json}")

    total_proteins = 0
    all_metadata = []

    with open(output_fasta, "w") as fasta_handle:
        if batch_mode and len(accessions) > 1:
            # Process in batches
            for i in range(0, len(accessions), BATCH_SIZE):
                batch = accessions[i : i + BATCH_SIZE]
                logger.info(
                    f"Processing batch {i//BATCH_SIZE + 1} ({len(batch)} accessions)"
                )

                try:
                    records = fetch_genbank_records_batch(batch, email)

                    for record in records:
                        for protein_record, metadata in extract_proteins_from_genbank(
                            record
                        ):
                            SeqIO.write(protein_record, fasta_handle, "fasta")

                            # Optionally remove sequence from JSON to reduce file size
                            if output_json and not include_sequence_in_json:
                                metadata_copy = metadata.copy()
                                metadata_copy.pop("protein_sequence", None)
                                all_metadata.append(metadata_copy)
                            else:
                                all_metadata.append(metadata)

                            total_proteins += 1

                except Exception as e:
                    logger.error(f"Error processing batch: {e}")
                    # Fall back to individual processing for this batch
                    logger.info("Falling back to individual accession processing")
                    for accession in batch:
                        try:
                            record = fetch_genbank_record(accession, email)
                            for (
                                protein_record,
                                metadata,
                            ) in extract_proteins_from_genbank(record):
                                SeqIO.write(protein_record, fasta_handle, "fasta")

                                if output_json and not include_sequence_in_json:
                                    metadata_copy = metadata.copy()
                                    metadata_copy.pop("protein_sequence", None)
                                    all_metadata.append(metadata_copy)
                                else:
                                    all_metadata.append(metadata)

                                total_proteins += 1
                        except Exception as e2:
                            logger.error(f"Failed to process {accession}: {e2}")
        else:
            # Process individually
            for idx, accession in enumerate(accessions, 1):
                logger.info(f"Processing {idx}/{len(accessions)}: {accession}")
                try:
                    record = fetch_genbank_record(accession, email)
                    for protein_record, metadata in extract_proteins_from_genbank(
                        record
                    ):
                        SeqIO.write(protein_record, fasta_handle, "fasta")

                        if output_json and not include_sequence_in_json:
                            metadata_copy = metadata.copy()
                            metadata_copy.pop("protein_sequence", None)
                            all_metadata.append(metadata_copy)
                        else:
                            all_metadata.append(metadata)

                        total_proteins += 1
                except Exception as e:
                    logger.error(f"Failed to process {accession}: {e}")

    # Write JSON metadata file
    if output_json and all_metadata:
        logger.info(f"Writing JSON metadata for {len(all_metadata)} proteins...")

        with open(output_json, "w", encoding="utf-8") as json_handle:
            json.dump(all_metadata, json_handle, indent=2, ensure_ascii=False)

        logger.info(f"JSON metadata saved to: {output_json}")

        logger.info(f"Download complete! Total proteins extracted: {total_proteins}")
        logger.info(f"Proteins saved to: {output_fasta}")

    # Print summary statistics
    if all_metadata:
        logger.info("\n=== Summary Statistics ===")
        organisms = set(m["organism"] for m in all_metadata if m["organism"])
        logger.info(f"Unique organisms: {len(organisms)}")

        genes_with_names = sum(1 for m in all_metadata if m["gene_name"])
        logger.info(f"Proteins with gene names: {genes_with_names}/{total_proteins}")

        avg_length = sum(m["protein_length"] for m in all_metadata) / len(all_metadata)
        logger.info(f"Average protein length: {avg_length:.1f} amino acids")


def main():
    parser = argparse.ArgumentParser(
        description="Download virus proteins and metadata from GenBank nucleotide accessions",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Download proteins and metadata from a file
  %(prog)s -i accessions.txt -e your.email@example.com -o proteins.fasta -j metadata.json
  
  # From command line accessions
  %(prog)s -a NC_045512 MN908947 -e your.email@example.com -o covid_proteins.fasta -j covid_metadata.json
  
  # With NCBI API key (faster)
  %(prog)s -i accessions.txt -e your.email@example.com -o proteins.fasta -j metadata.json -k YOUR_API_KEY
  
  # Just FASTA without metadata file
  %(prog)s -i accessions.txt -e your.email@example.com -o proteins.fasta
  
  # JSON without protein sequences (smaller file)
  %(prog)s -i accessions.txt -e your.email@example.com -o proteins.fasta -j metadata.json --no-sequence
        """,
    )

    parser.add_argument(
        "-i", "--input", help="Input file with GenBank accessions (one per line)"
    )
    parser.add_argument(
        "-a",
        "--accessions",
        nargs="+",
        help="GenBank accession numbers (space-separated)",
    )
    parser.add_argument(
        "-e", "--email", required=True, help="Your email address (required by NCBI)"
    )
    parser.add_argument(
        "-o", "--output", required=True, help="Output FASTA file for proteins"
    )
    parser.add_argument(
        "-j", "--json", help="Output JSON file for protein metadata (optional)"
    )
    parser.add_argument(
        "-k", "--api-key", help="NCBI API key (optional, allows faster requests)"
    )
    parser.add_argument(
        "--no-batch",
        action="store_true",
        help="Disable batch mode (process accessions individually)",
    )
    parser.add_argument(
        "--no-sequence",
        action="store_true",
        help="Don't include protein sequences in JSON output (reduces file size)",
    )

    args = parser.parse_args()

    # Get accessions from file or command line
    if args.input:
        with open(args.input) as f:
            accessions = [line.strip() for line in f if line.strip()]
    elif args.accessions:
        accessions = args.accessions
    else:
        parser.error("Must provide either --input or --accessions")

    # Set API key if provided
    if args.api_key:
        Entrez.api_key = args.api_key
        global DELAY_BETWEEN_REQUESTS
        DELAY_BETWEEN_REQUESTS = 0.11  # 9 req/s with API key
        logger.info("Using NCBI API key for faster requests")

    # Run the download
    download_virus_proteins(
        accessions=accessions,
        email=args.email,
        output_fasta=args.output,
        output_json=args.json,
        batch_mode=not args.no_batch,
        include_sequence_in_json=not args.no_sequence,
    )


if __name__ == "__main__":
    main()
