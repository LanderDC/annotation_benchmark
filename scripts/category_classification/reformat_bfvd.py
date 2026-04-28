import json
import taxopy

# Load the JSON file
with open("../../data/uniparc_protein_names_eukaryotic_viruses.json", "r") as f:
    data = json.load(f)

# Initialize the taxopy TaxDb
taxdb = taxopy.TaxDb()

# Cache taxonomy lookups
taxonomy_cache = {}

def get_taxonomy_lineage(taxid):
    """Fetch full taxonomy lineage using taxopy."""
    if taxid in taxonomy_cache:
        return taxonomy_cache[taxid]

    try:
        taxon = taxopy.Taxon(taxid, taxdb)

        # name_lineage goes from the taxon up to root, so reverse it
        # to get root -> taxon order, then drop "root"
        lineage = list(reversed(taxon.name_lineage))
        if lineage and lineage[0] == "root":
            lineage = lineage[1:]

        taxonomy_cache[taxid] = lineage
        return lineage

    except taxopy.exceptions.TaxidError:
        print(f"Taxid {taxid} not found in database")
        taxonomy_cache[taxid] = []
        return []
    except Exception as e:
        print(f"Error fetching taxonomy for taxid {taxid}: {e}")
        taxonomy_cache[taxid] = []
        return []

# Process each entry
for accession, entry in data.items():
    # Remove categories
    entry.pop("categories", None)

    # Add taxonomy lineage
    taxid = entry.get("taxid")
    if taxid:
        lineage = get_taxonomy_lineage(taxid)
        entry["taxonomy"] = lineage

    print(f"Processed {accession} (taxid: {taxid})")

# Save the updated JSON
with open("../../data/bfvd_eukarytic_viruses_names_taxonomy.json", "w") as f:
    json.dump(data, f, indent=2)

print(f"\nDone. Processed {len(data)} entries.")
print(f"Unique taxids looked up: {len(taxonomy_cache)}")
