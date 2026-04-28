#!/usr/bin/env python3
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BFVD = ROOT / "data" / "bfvd_category_annotations.json"
UNIPARC = ROOT / "data" / "uniparc_protein_names_eukaryotic_viruses.json"
OUT = ROOT / "data" / "bfvd_empathi_classifications.json"

with BFVD.open("r", encoding="utf-8") as fh:
    bfvd = json.load(fh)

with UNIPARC.open("r", encoding="utf-8") as fh:
    uniparc = json.load(fh)

out = {k: v for k, v in bfvd.items() if k not in uniparc}

OUT.parent.mkdir(parents=True, exist_ok=True)
with OUT.open("w", encoding="utf-8") as fh:
    json.dump(out, fh, indent=2, ensure_ascii=False)

print(f"Wrote {len(out)} entries to {OUT}")
