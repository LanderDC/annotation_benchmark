#!/usr/bin/env python3
"""
Protein Name Functional Classifier
===================================
Classifies protein names into predefined functional categories using a local
open-source LLM via Ollama (default model: gemma3:12b).

Requirements:
    pip install ollama tqdm

Usage:
    1. Start Ollama:       ollama serve
    2. Pull the model:     ollama pull gemma3:12b
    3. Run the script:
         python classify_proteins.py \
             --proteins proteins.json \
             --categories categories.json \
             --output results.json \
             [--model gemma3:12b] \
             [--batch-size 10] \
             [--temperature 0.1] \
             [--max-retries 3] \
             [--workers 4]
"""

import argparse
import json
import logging
import re
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any

try:
    import ollama
except ImportError:
    print("ERROR: 'ollama' package not found. Install with: pip install ollama")
    sys.exit(1)

try:
    from tqdm import tqdm
except ImportError:
    print("ERROR: 'tqdm' package not found. Install with: pip install tqdm")
    sys.exit(1)

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Prompt templates
# ---------------------------------------------------------------------------

SYSTEM_PROMPT = """\
You are an expert molecular biologist specializing in protein function annotation.
Your task is to classify proteins into functional categories based on their names.
You must be precise, conservative, and only assign categories that are clearly
supported by the protein name(s). If a protein is uncharacterized or hypothetical,
include the appropriate category for that. A protein can belong to multiple categories.
Respond ONLY with valid JSON — no commentary, no markdown fences."""

SINGLE_PROTEIN_PROMPT = """\
Given the following predefined functional categories:
{categories}

Classify the protein with accession "{accession}" based on its name(s):
{protein_names}

Respond with a JSON object in exactly this format:
{{"accession": "{accession}", "categories": ["Category1", "Category2"]}}

Rules:
- Choose ONLY from the predefined categories listed above.
- Assign every category that applies; a protein may belong to multiple categories.
- If the protein is described as "uncharacterized", "hypothetical", or "putative"
  with no other functional clue, include "Hypothetical Protein" if it is one of the
  predefined categories.
- Do NOT invent new categories.
- Output ONLY the JSON object, nothing else."""

BATCH_PROTEIN_PROMPT = """\
Given the following predefined functional categories:
{categories}

Classify each of the following proteins based on their name(s).

Proteins:
{proteins_block}

Respond with a JSON array of objects, one per protein, in exactly this format:
[
  {{"accession": "XXXXX", "categories": ["Category1", "Category2"]}},
  ...
]

Rules:
- Choose ONLY from the predefined categories listed above.
- Assign every category that applies; a protein may belong to multiple categories.
- If the protein is described as "uncharacterized", "hypothetical", or "putative"
  with no other functional clue, include "Hypothetical Protein" if it is one of the
  predefined categories.
- Do NOT invent new categories.
- Output ONLY the JSON array, nothing else."""


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def load_json(path: str | Path) -> Any:
    """Load a JSON file and return its contents."""
    path = Path(path)
    if not path.exists():
        logger.error(f"File not found: {path}")
        sys.exit(1)
    with open(path, "r", encoding="utf-8") as fh:
        return json.load(fh)


def save_json(data: Any, path: str | Path) -> None:
    """Save data to a JSON file."""
    path = Path(path)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(data, fh, indent=2, ensure_ascii=False)
    logger.info(f"Results saved to {path}")


def extract_json_from_response(text: str) -> Any:
    """
    Robustly extract JSON from an LLM response that may contain markdown
    fences, preamble text, or trailing commentary.
    """
    # Strip markdown code fences if present
    text = re.sub(r"```(?:json)?", "", text).strip()

    # Try parsing the whole thing first
    try:
        return json.loads(text)
    except json.JSONDecodeError:
        pass

    # Try to find a JSON array
    match = re.search(r"\[.*\]", text, re.DOTALL)
    if match:
        try:
            return json.loads(match.group())
        except json.JSONDecodeError:
            pass

    # Try to find a JSON object
    match = re.search(r"\{.*\}", text, re.DOTALL)
    if match:
        try:
            return json.loads(match.group())
        except json.JSONDecodeError:
            pass

    return None


def validate_categories(assigned: list[str], valid_categories: set[str]) -> list[str]:
    """Filter out any categories not in the predefined set."""
    validated = [c for c in assigned if c in valid_categories]
    invalid = [c for c in assigned if c not in valid_categories]
    if invalid:
        logger.warning(f"Removed invalid categories: {invalid}")
    return validated


# ---------------------------------------------------------------------------
# LLM interaction
# ---------------------------------------------------------------------------


def query_ollama(
    model: str,
    system_prompt: str,
    user_prompt: str,
    temperature: float = 0.1,
    max_retries: int = 3,
) -> str:
    """Send a prompt to Ollama and return the response text with retries."""
    for attempt in range(1, max_retries + 1):
        try:
            response = ollama.chat(
                model=model,
                messages=[
                    {"role": "system", "content": system_prompt},
                    {"role": "user", "content": user_prompt},
                ],
                options={"temperature": temperature, "num_predict": 4096},
            )
            return response["message"]["content"].strip()
        except Exception as e:
            logger.warning(
                f"Ollama request failed (attempt {attempt}/{max_retries}): {e}"
            )
            if attempt < max_retries:
                time.sleep(2**attempt)
            else:
                raise


def classify_single_protein(
    accession: str,
    protein_names: list[str],
    categories_list: list[str],
    valid_categories: set[str],
    model: str,
    temperature: float,
    max_retries: int,
) -> dict:
    """Classify a single protein and return a result dict."""
    categories_str = "\n".join(f"  - {c}" for c in categories_list)
    names_str = "\n".join(f'  - "{n}"' for n in protein_names)

    prompt = SINGLE_PROTEIN_PROMPT.format(
        categories=categories_str,
        accession=accession,
        protein_names=names_str,
    )

    raw = query_ollama(model, SYSTEM_PROMPT, prompt, temperature, max_retries)
    parsed = extract_json_from_response(raw)

    if parsed and isinstance(parsed, dict) and "categories" in parsed:
        cats = validate_categories(parsed["categories"], valid_categories)
        return {"accession": accession, "categories": cats, "raw_response": raw}

    logger.warning(f"Failed to parse response for {accession}. Raw: {raw[:200]}...")
    return {"accession": accession, "categories": [], "raw_response": raw}


def classify_batch(
    batch: list[tuple[str, list[str]]],
    categories_list: list[str],
    valid_categories: set[str],
    model: str,
    temperature: float,
    max_retries: int,
) -> list[dict]:
    """Classify a batch of proteins in a single LLM call."""
    categories_str = "\n".join(f"  - {c}" for c in categories_list)

    proteins_lines = []
    for accession, names in batch:
        names_joined = "; ".join(f'"{n}"' for n in names)
        proteins_lines.append(f'  Accession "{accession}": [{names_joined}]')
    proteins_block = "\n".join(proteins_lines)

    prompt = BATCH_PROTEIN_PROMPT.format(
        categories=categories_str,
        proteins_block=proteins_block,
    )

    raw = query_ollama(model, SYSTEM_PROMPT, prompt, temperature, max_retries)
    parsed = extract_json_from_response(raw)

    results = []
    accession_set = {acc for acc, _ in batch}

    if parsed and isinstance(parsed, list):
        parsed_map = {}
        for item in parsed:
            if isinstance(item, dict) and "accession" in item:
                acc = item["accession"]
                cats = validate_categories(item.get("categories", []), valid_categories)
                parsed_map[acc] = cats

        for accession, names in batch:
            if accession in parsed_map:
                results.append(
                    {
                        "accession": accession,
                        "categories": parsed_map[accession],
                        "raw_response": None,  # shared batch response
                    }
                )
            else:
                logger.warning(
                    f"Accession {accession} missing from batch response; "
                    f"falling back to single classification."
                )
                result = classify_single_protein(
                    accession,
                    names,
                    categories_list,
                    valid_categories,
                    model,
                    temperature,
                    max_retries,
                )
                results.append(result)
    else:
        # Batch parse failed — fall back to individual classification
        logger.warning(
            "Batch response could not be parsed. Falling back to per-protein "
            "classification for this batch."
        )
        for accession, names in batch:
            result = classify_single_protein(
                accession,
                names,
                categories_list,
                valid_categories,
                model,
                temperature,
                max_retries,
            )
            results.append(result)

    return results


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------


def run_classification(
    proteins_path: str,
    categories_path: str,
    output_path: str,
    model: str = "gemma3:12b",
    batch_size: int = 10,
    temperature: float = 0.1,
    max_retries: int = 3,
    workers: int = 1,
) -> None:
    """Run the full classification pipeline."""

    # --- Load inputs -------------------------------------------------------
    logger.info(f"Loading proteins from {proteins_path}")
    proteins_data = load_json(proteins_path)
    logger.info(f"Loaded {len(proteins_data)} proteins")

    logger.info(f"Loading categories from {categories_path}")
    categories_raw = load_json(categories_path)

    # categories.json may be a list of strings or a dict with a key holding a
    # list.  Handle both gracefully.
    if isinstance(categories_raw, list):
        categories_list: list[str] = categories_raw
    elif isinstance(categories_raw, dict):
        # Try common key names
        for key in ("categories", "functional_categories", "classes"):
            if key in categories_raw and isinstance(categories_raw[key], list):
                categories_list = categories_raw[key]
                break
        else:
            # Use all values that are strings, or flatten first list found
            for v in categories_raw.values():
                if isinstance(v, list):
                    categories_list = v
                    break
            else:
                logger.error("Cannot determine category list from categories JSON.")
                sys.exit(1)
    else:
        logger.error("categories.json must be a list or a dict.")
        sys.exit(1)

    valid_categories = set(categories_list)
    logger.info(f"Using {len(categories_list)} predefined categories")

    # --- Verify model is available -----------------------------------------
    logger.info(f"Checking that model '{model}' is available in Ollama...")
    try:
        available = ollama.list()
        model_names = [m.model for m in available.models]
        if not any(model in name for name in model_names):
            logger.warning(f"Model '{model}' not found locally. Attempting to pull...")
            ollama.pull(model)
    except Exception as e:
        logger.error(f"Could not connect to Ollama: {e}")
        logger.error("Make sure Ollama is running ('ollama serve').")
        sys.exit(1)

    # --- Prepare batches ---------------------------------------------------
    protein_items: list[tuple[str, list[str]]] = [
        (acc, info["protein_names"]) for acc, info in proteins_data.items()
    ]

    batches: list[list[tuple[str, list[str]]]] = [
        protein_items[i : i + batch_size]
        for i in range(0, len(protein_items), batch_size)
    ]
    logger.info(
        f"Split {len(protein_items)} proteins into {len(batches)} batches "
        f"(batch_size={batch_size})"
    )

    # --- Classify ----------------------------------------------------------
    all_results: dict[str, list[str]] = {}
    errors: list[str] = []

    def process_batch(batch):
        return classify_batch(
            batch,
            categories_list,
            valid_categories,
            model,
            temperature,
            max_retries,
        )

    if workers <= 1:
        # Sequential processing
        for batch in tqdm(batches, desc="Classifying", unit="batch"):
            try:
                results = process_batch(batch)
                for r in results:
                    all_results[r["accession"]] = r["categories"]
            except Exception as e:
                for acc, _ in batch:
                    errors.append(acc)
                logger.error(f"Batch failed: {e}")
    else:
        # Parallel processing (useful if running multiple Ollama instances or
        # the model is fast enough to benefit from request pipelining)
        with ThreadPoolExecutor(max_workers=workers) as executor:
            futures = {
                executor.submit(process_batch, batch): batch for batch in batches
            }
            with tqdm(total=len(batches), desc="Classifying", unit="batch") as pbar:
                for future in as_completed(futures):
                    batch = futures[future]
                    try:
                        results = future.result()
                        for r in results:
                            all_results[r["accession"]] = r["categories"]
                    except Exception as e:
                        for acc, _ in batch:
                            errors.append(acc)
                        logger.error(f"Batch failed: {e}")
                    pbar.update(1)

    # --- Build output ------------------------------------------------------
    output_data = {}
    for accession, info in proteins_data.items():
        output_data[accession] = {
            "protein_names": info["protein_names"],
            "taxid": info.get("taxid"),
            "original_categories": info.get("categories", []),
            "predicted_categories": all_results.get(accession, []),
        }

    save_json(output_data, output_path)

    # --- Summary -----------------------------------------------------------
    total = len(proteins_data)
    classified = sum(1 for v in output_data.values() if v["predicted_categories"])
    logger.info(
        f"Classification complete: {classified}/{total} proteins received "
        f"at least one category."
    )
    if errors:
        logger.warning(f"{len(errors)} proteins encountered errors: {errors}")

    # Optional: compute agreement with existing labels if present
    if any(info.get("categories") for info in proteins_data.values()):
        compute_agreement(output_data)


def compute_agreement(output_data: dict) -> None:
    """Compute and log simple agreement metrics against original labels."""
    exact_match = 0
    total_with_labels = 0
    total_precision_sum = 0.0
    total_recall_sum = 0.0

    for accession, info in output_data.items():
        original = set(info.get("original_categories", []))
        predicted = set(info.get("predicted_categories", []))

        if not original:
            continue
        total_with_labels += 1

        if original == predicted:
            exact_match += 1

        if predicted:
            precision = len(original & predicted) / len(predicted)
            total_precision_sum += precision
        if original:
            recall = len(original & predicted) / len(original)
            total_recall_sum += recall

    if total_with_labels == 0:
        return

    avg_precision = total_precision_sum / total_with_labels
    avg_recall = total_recall_sum / total_with_labels
    f1 = (
        2 * avg_precision * avg_recall / (avg_precision + avg_recall)
        if (avg_precision + avg_recall) > 0
        else 0.0
    )

    logger.info("--- Agreement with original labels ---")
    logger.info(f"  Proteins with labels: {total_with_labels}")
    logger.info(
        f"  Exact match: {exact_match}/{total_with_labels} "
        f"({100 * exact_match / total_with_labels:.1f}%)"
    )
    logger.info(f"  Avg precision: {avg_precision:.3f}")
    logger.info(f"  Avg recall:    {avg_recall:.3f}")
    logger.info(f"  Avg F1:        {f1:.3f}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(
        description="Classify protein names into functional categories using "
        "a local LLM via Ollama.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--proteins",
        required=True,
        help="Path to JSON file with protein data (see script header for format).",
    )
    parser.add_argument(
        "--categories",
        required=True,
        help="Path to JSON file with the list of predefined functional categories.",
    )
    parser.add_argument(
        "--output",
        default="classification_results.json",
        help="Path for the output JSON file (default: classification_results.json).",
    )
    parser.add_argument(
        "--model",
        default="gemma3:12b",
        help="Ollama model to use (default: gemma3:12b). "
        "Other good options: mistral, llama3, deepseek-r1.",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=10,
        help="Number of proteins per LLM call (default: 10). "
        "Use 1 for maximum reliability, higher for speed.",
    )
    parser.add_argument(
        "--temperature",
        type=float,
        default=0.1,
        help="LLM sampling temperature (default: 0.1). Lower = more deterministic.",
    )
    parser.add_argument(
        "--max-retries",
        type=int,
        default=3,
        help="Max retries per LLM call on failure (default: 3).",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=1,
        help="Number of parallel worker threads (default: 1).",
    )

    args = parser.parse_args()

    run_classification(
        proteins_path=args.proteins,
        categories_path=args.categories,
        output_path=args.output,
        model=args.model,
        batch_size=args.batch_size,
        temperature=args.temperature,
        max_retries=args.max_retries,
        workers=args.workers,
    )


if __name__ == "__main__":
    main()
