#!/usr/bin/env python3
"""
Viral Protein Functional Classifier v2.0
=========================================
Classifies viral proteins into functional categories based on their names.

Author: Generated for viral protein classification
Date: December 2025
"""

import json
import re
from collections import defaultdict
from typing import Dict, List


class ViralProteinClassifier:
    """Classifier for viral proteins based on name-to-function mapping."""

    def __init__(self, categories_file="viral_protein_functional_categories.json"):
        """Initialize classifier with functional categories."""
        with open(categories_file, "r") as f:
            self.categories = json.load(f)

        # Keywords that must match as whole words even if longer than 2 chars
        self.exact_match_keywords = {"l protein", "n protein", "dna polymerase"}

        # Priority order with split categories
        self.priority_order = [
            "DNA Packaging/Capsid Maturation",
            "Tail/Baseplate Structural",
            "Genome-Linked Protein (VPg/Terminal)",
            "Reverse Transcriptase",
            "RNA-Dependent RNA Polymerase (RdRp)",
            "DNA Polymerase",
            "Viral Protease",
            "Integrase/Recombinase",
            "Primase/Primase-Polymerase",
            "mRNA Capping/Methyltransferase",
            "Nuclease (Endo/Exonuclease)",
            "Helicase",
            "NTPase/ATPase",
            "Nucleotide Metabolism Enzyme",
            "Viral Kinase/Phosphoprotein",
            "Occlusion Body Protein",
            "Nucleocapsid Protein",
            "Capsid Protein",
            "Envelope/Surface Glycoprotein",
            "Matrix/Tegument Structural",
            "Movement Protein",
            "Viroporin/Ion Channel",
            "Host Shutoff/Translation Inhibitor",
            "Innate Immune/Interferon Antagonist",
            "Transcriptional Regulator/Transactivator",
            "RNA Export/Splicing Regulator",
            "RNA-Binding Regulatory Protein",
            "Apoptosis/Cell-Cycle Modulator",
            "Cytokine/Receptor Mimic or Modulator",
            "Ubiquitin/SUMO Pathway Protein",
            "Episome Maintenance/Replication Origin Binding",
            "Replication Cofactor/Processivity Factor",
            "Membrane Remodeling/Replication Organelle",
            "Assembly/Morphogenesis Factor",
            "Virion Release/Budding Facilitator",
            "Viroplasm/Viral Factory Protein",
            "Ligase",
            "Other Viral Enzyme",
            "Polyprotein Precursor",
            "Accessory/Virulence Factor",
            "Non-Structural Protein (General)",
            "Hypothetical Protein",
        ]

    def classify_protein(
        self, protein_name: str, return_all_matches: bool = False
    ) -> List[str]:
        """Classify a single protein name into functional categories."""
        protein_lower = protein_name.lower().strip()
        matches = []

        for category in self.priority_order:
            if category not in self.categories:
                continue

            cat_info = self.categories[category]
            keywords = cat_info["keywords"]
            exclude = cat_info.get("exclude", [])

            # Check exclusions
            exclude_match = False
            for ex in exclude:
                if ex.lower() in protein_lower:
                    exclude_match = True
                    break

            if exclude_match:
                continue

            # Check keywords
            keyword_match = False
            for kw in keywords:
                kw_lower = kw.lower()

                # Short keywords need word boundary
                if len(kw_lower) <= 2 or kw_lower in [
                    "vp1",
                    "vp2",
                    "vp3",
                    "vp4",
                    "vp5",
                    "vp6",
                    "vp7",
                    "ns1",
                    "ns2",
                    "ns3",
                    "ns4",
                    "ns5",
                    "gn ",
                    "gc ",
                ]:
                    pattern = r"\b" + re.escape(kw_lower.strip()) + r"\b"
                    if re.search(pattern, protein_lower):
                        keyword_match = True
                        break
                elif kw_lower in self.exact_match_keywords:
                    pattern = r"\b" + re.escape(kw_lower) + r"\b"
                    if re.search(pattern, protein_lower):
                        keyword_match = True
                        break
                else:
                    if kw_lower in protein_lower:
                        keyword_match = True
                        break

            if keyword_match:
                matches.append(category)
                if not return_all_matches:
                    return matches

        if not matches:
            matches.append("Other function")

        return matches if return_all_matches else matches[:1]

    def classify_dataset(
        self, proteins_data: List[Dict], return_all_matches: bool = False
    ) -> Dict:
        """
        Classify all proteins in a dataset.

        Args:
            proteins_data: List of protein dictionaries from JSON
            return_all_matches: Whether to return all matching categories

        Returns:
            Dictionary mapping protein_id to classification results
        """
        classifications = {}

        for protein in proteins_data:
            protein_id = protein.get("protein_id", "unknown")
            product = protein.get("product", "hypothetical protein")

            # Product is already a string in your data
            protein_names = [product] if product else []

            all_categories = set()
            for name in protein_names:
                categories = self.classify_protein(
                    name, return_all_matches=return_all_matches
                )
                all_categories.update(categories)

            classifications[protein_id] = {
                "product": product,  # Keep as single string
                "gene_name": protein.get("gene_name", ""),
                "locus_tag": protein.get("locus_tag", ""),
                "organism": protein.get("organism", ""),
                "categories": sorted(list(all_categories)),
            }

        return classifications

    def get_category_statistics(self, classifications: Dict) -> Dict:
        """Get statistics about category distribution."""
        category_counts = defaultdict(int)
        total_proteins = len(classifications)

        for protein_id, data in classifications.items():
            for category in data["categories"]:
                category_counts[category] += 1

        statistics = {
            "total_proteins": total_proteins,
            "category_counts": dict(category_counts),
            "category_percentages": {
                cat: (count / total_proteins) * 100
                for cat, count in category_counts.items()
            },
        }

        return statistics

    def export_classifications(self, classifications: Dict, output_file: str):
        """Export classifications to JSON file."""
        with open(output_file, "w") as f:
            json.dump(classifications, f, indent=2)
        print(f"Classifications exported to {output_file}")


def main():
    """Main function demonstrating usage."""
    import argparse

    parser = argparse.ArgumentParser(
        description="Classify viral proteins into functional categories (v2.0 with split categories)"
    )
    parser.add_argument("input_file", help="Input JSON file with protein data")
    parser.add_argument(
        "-o",
        "--output",
        default="protein_classifications.json",
        help="Output file for classifications",
    )
    parser.add_argument(
        "-s",
        "--stats",
        default="classification_statistics.json",
        help="Output file for statistics",
    )
    parser.add_argument(
        "-c",
        "--categories",
        default="viral_protein_functional_categories.json",
        help="Categories JSON file",
    )
    parser.add_argument(
        "-a",
        "--all-matches",
        action="store_true",
        help="Return all matching categories",
    )

    args = parser.parse_args()

    print(f"Loading protein data from {args.input_file}...")
    with open(args.input_file, "r") as f:
        proteins_data = json.load(f)
    print(f"Loaded {len(proteins_data)} proteins")

    print(f"\nInitializing classifier...")
    classifier = ViralProteinClassifier(args.categories)
    print(f"Loaded {len(classifier.categories)} functional categories")

    print(f"\nClassifying proteins...")
    classifications = classifier.classify_dataset(
        proteins_data, return_all_matches=args.all_matches
    )

    print(f"Calculating statistics...")
    stats = classifier.get_category_statistics(classifications)

    classifier.export_classifications(classifications, args.output)

    with open(args.stats, "w") as f:
        json.dump(stats, f, indent=2)
    print(f"Statistics exported to {args.stats}")

    print(f"\n{'='*70}")
    print("CLASSIFICATION SUMMARY (v2.0)")
    print(f"{'='*70}")
    print(f"Total proteins: {stats['total_proteins']:,}")
    print(f"\nTop 15 categories:")
    sorted_cats = sorted(
        stats["category_counts"].items(), key=lambda x: x[1], reverse=True
    )
    for i, (cat, count) in enumerate(sorted_cats[:15], 1):
        pct = stats["category_percentages"][cat]
        print(f"{i:2d}. {cat:50s}: {count:7,d} ({pct:5.2f}%)")
    print(f"{'='*70}")


if __name__ == "__main__":
    main()
