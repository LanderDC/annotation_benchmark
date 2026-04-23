#!/usr/bin/env python3
"""
Protein Name Functional Classifier
===================================
Classifies protein names into predefined functional categories using a local
open-source LLM via HuggingFace Transformers.

Requirements:
    pip install transformers torch accelerate tqdm

Usage:
    python classify_proteins.py \
        --proteins proteins.json \
        --categories categories.json \
        --output results.json \
        [--model google/gemma-4-31b-it] \
        [--batch-size 10] \
        [--temperature 0.1] \
        [--max-new-tokens 4096] \
        [--device auto]

SLURM example:
    #!/bin/bash
    #SBATCH --job-name=protein_classify
    #SBATCH --gres=gpu:1
    #SBATCH --mem=32G
    #SBATCH --cpus-per-task=4
    #SBATCH --time=04:00:00

    module load cuda
    source ~/myenv/bin/activate
    python classify_proteins.py \\
        --proteins proteins.json \\
        --categories categories.json \\
        --output results.json \\
        --model google/gemma-4-31b-it
"""

import argparse
import json
import logging
import os
import re
import sys
from pathlib import Path
from typing import Any

try:
    import torch
except ImportError:
    print("ERROR: 'torch' not found. Install with: pip install torch")
    sys.exit(1)

try:
    from transformers import AutoModelForCausalLM, AutoTokenizer, BitsAndBytesConfig
except ImportError:
    print(
        "ERROR: 'transformers' not found. Install with: "
        "pip install transformers accelerate"
    )
    sys.exit(1)

try:
    from tqdm import tqdm
except ImportError:
    print("ERROR: 'tqdm' not found. Install with: pip install tqdm")
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
You are an expert molecular biologist in virology specializing in protein function annotation.

Your task is to classify proteins into functional categories based on their names. 
In addition to the protein name itself, you MUST also consider the virus taxonomy (e.g., family, genus, or type of virus) to infer likely protein function. 
Different virus groups often use distinct naming conventions and encode proteins with characteristic functions, so use this contextual information to guide your classification.
On the other hand, different virus groups may also use similar names for proteins with different functions.

Be precise and conservative: only assign functional categories that are clearly supported by the protein name and are biologically plausible given the virus type. 
Avoid over-interpreting ambiguous names.

If a protein is labeled as "hypothetical", "uncharacterized", or provides no clear functional indication even when considering the virus taxonomy, classify it under "Hypothetical Protein".

A protein may belong to multiple categories if justified.

Respond ONLY with valid JSON — no commentary, no markdown fences.
"""

SINGLE_PROTEIN_PROMPT = """\
Given the following predefined functional categories:
{categories}

Classify the protein with accession "{accession}" based on its name(s):
{protein_names}

Taxonomy:
{taxonomy}

Respond with a JSON object in exactly this format:
{{"accession": "{accession}", "categories": ["Category1", "Category2"]}}

Rules:
- Choose ONLY from the predefined categories listed above.
- Assign every category that applies; a protein may belong to multiple categories.
- If the protein is described as "uncharacterized", "hypothetical", or "putative"
  with no other functional clue, assign it to "Hypothetical Protein".
- Do NOT invent new categories.
- Output ONLY the JSON object, nothing else."""

BATCH_PROTEIN_PROMPT = """\
Given the following predefined functional categories:
{categories}

Classify each of the following proteins based on their name(s) and
their taxonomy.

Proteins:
{proteins_block}

Respond with a JSON array of objects, one per protein, in exactly this format:
[
  {{"accession": {"XXXXX"}, "categories": ["Category1", "Category2"]}},
  ...
]

Rules:
- Choose ONLY from the predefined categories listed above.
- Assign every category that applies; a protein may belong to multiple categories.
- If the protein is described as "uncharacterized", "hypothetical", or "putative"
  with no other functional clue, assign it to "Hypothetical Protein".
- Do NOT invent new categories.
- Output ONLY the JSON array, nothing else."""


# ---------------------------------------------------------------------------
# Model loading
# ---------------------------------------------------------------------------


def load_model(
    model_name: str,
    device: str = "auto",
    quantize: str | None = None,
    hf_token: str | None = None,
):
    """
    Load a HuggingFace causal LM and its tokenizer.

    Parameters
    ----------
    model_name : str
        HuggingFace model ID, e.g. "google/gemma-4-31b-it".
    device : str
        "auto", "cuda", "cpu", or a specific device like "cuda:0".
    quantize : str or None
        "4bit" or "8bit" for bitsandbytes quantisation (reduces VRAM usage).
    hf_token : str or None
        HuggingFace access token for gated models (e.g. Gemma).

    Returns
    -------
    model, tokenizer
    """
    logger.info(f"Loading tokenizer for {model_name}...")
    tokenizer = AutoTokenizer.from_pretrained(
        model_name,
        token=hf_token,
    )

    # --- Quantisation config -----------------------------------------------
    quant_config = None
    if quantize == "4bit":
        logger.info("Using 4-bit quantisation (bitsandbytes)")
        try:
            quant_config = BitsAndBytesConfig(
                load_in_4bit=True,
                bnb_4bit_compute_dtype=torch.float16,
                bnb_4bit_use_double_quant=True,
                bnb_4bit_quant_type="nf4",
            )
        except Exception:
            logger.error(
                "4-bit quantisation requires bitsandbytes: " "pip install bitsandbytes"
            )
            sys.exit(1)
    elif quantize == "8bit":
        logger.info("Using 8-bit quantisation (bitsandbytes)")
        try:
            quant_config = BitsAndBytesConfig(load_in_8bit=True)
        except Exception:
            logger.error(
                "8-bit quantisation requires bitsandbytes: " "pip install bitsandbytes"
            )
            sys.exit(1)

    # --- Load model --------------------------------------------------------
    logger.info(f"Loading model {model_name} (this may take a few minutes)...")
    load_kwargs: dict[str, Any] = {
        "token": hf_token,
        "dtype": torch.float16,
    }
    if quant_config:
        load_kwargs["quantization_config"] = quant_config
        load_kwargs["device_map"] = "auto"
    elif device == "auto":
        load_kwargs["device_map"] = "auto"
    else:
        load_kwargs["device_map"] = device

    model = AutoModelForCausalLM.from_pretrained(model_name, **load_kwargs)

    logger.info("Model loaded successfully")
    if hasattr(model, "hf_device_map"):
        logger.info(f"Device map: {model.hf_device_map}")

    using_gpu = False
    if hasattr(model, "hf_device_map") and model.hf_device_map:
        for mapped_device in model.hf_device_map.values():
            if isinstance(mapped_device, int):
                using_gpu = True
                break
            mapped_device_str = str(mapped_device).lower()
            if "cuda" in mapped_device_str or "mps" in mapped_device_str:
                using_gpu = True
                break
    elif hasattr(model, "device"):
        model_device = str(model.device).lower()
        using_gpu = "cuda" in model_device or "mps" in model_device

    logger.info(f"Running on {'GPU' if using_gpu else 'CPU'}")

    return model, tokenizer


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
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = path.with_suffix(path.suffix + ".tmp")
    with open(tmp_path, "w", encoding="utf-8") as fh:
        json.dump(data, fh, indent=2, ensure_ascii=False)
    os.replace(tmp_path, path)
    logger.info(f"Results saved to {path}")


def coerce_protein_fields(info: dict) -> tuple[list[str], list[str] | None]:
    """Normalize protein names and taxonomy fields for prompting/output."""
    protein_names = info.get("protein_names")
    if protein_names is None:
        protein_names = info.get("protein", [])
    if isinstance(protein_names, str):
        protein_names = [protein_names]
    if protein_names is None:
        protein_names = []

    taxonomy = info.get("taxonomy", [])
    if isinstance(taxonomy, str):
        taxonomy = [taxonomy]

    return protein_names, taxonomy


def load_existing_output(path: str | Path) -> dict[str, dict]:
    """Load an existing output JSON if present, otherwise return an empty dict."""
    path = Path(path)
    if not path.exists():
        return {}

    data = load_json(path)
    if not isinstance(data, dict):
        logger.warning(
            f"Existing output at {path} is not a JSON object; ignoring for resume."
        )
        return {}

    logger.info(f"Loaded {len(data)} existing records from {path} for resume")
    return data


def extract_json_from_response(text: str) -> Any:
    """
    Robustly extract JSON from an LLM response that may contain markdown
    fences, preamble text, or trailing commentary.
    """
    text = re.sub(r"```(?:json)?", "", text).strip()

    try:
        return json.loads(text)
    except json.JSONDecodeError:
        pass

    match = re.search(r"\[.*\]", text, re.DOTALL)
    if match:
        try:
            return json.loads(match.group())
        except json.JSONDecodeError:
            pass

    match = re.search(r"\{.*\}", text, re.DOTALL)
    if match:
        try:
            return json.loads(match.group())
        except json.JSONDecodeError:
            pass

    return None


def validate_categories(
    assigned: list[str],
    valid_categories: set[str],
    accession: str | None = None,
) -> list[str]:
    """Filter out any categories not in the predefined set."""
    validated = [c for c in assigned if c in valid_categories]
    invalid = [c for c in assigned if c not in valid_categories]
    if invalid:
        if accession:
            logger.warning(
                f"Removed invalid categories for accession {accession}: {invalid}"
            )
        else:
            logger.warning(f"Removed invalid categories: {invalid}")
    return validated


def format_taxonomy(taxonomy: list[str] | None) -> str:
    """Format taxonomy lines for the prompt."""
    if not taxonomy:
        return '  - "Unknown"'
    return "\n".join(f'  - "{rank}"' for rank in taxonomy)


def normalize_categories(categories_raw: Any) -> list[str]:
    """Normalize category inputs into a deduplicated list of category names."""
    categories_list: list[str] | None = None

    if isinstance(categories_raw, list):
        if all(isinstance(item, str) for item in categories_raw):
            categories_list = categories_raw
        elif all(isinstance(item, dict) for item in categories_raw):
            candidate_keys = (
                "Enzyme/Protein",
                "category",
                "name",
                "Activity",
                "process",
            )
            extracted: list[str] = []
            for item in categories_raw:
                for key in candidate_keys:
                    value = item.get(key)
                    if isinstance(value, str) and value.strip():
                        extracted.append(value.strip())
                        break
            if extracted:
                categories_list = extracted
    elif isinstance(categories_raw, dict):
        for key in ("categories", "functional_categories", "classes"):
            value = categories_raw.get(key)
            if isinstance(value, list):
                categories_raw = value
                break
        else:
            for value in categories_raw.values():
                if isinstance(value, list):
                    categories_raw = value
                    break
            else:
                categories_raw = None

        if isinstance(categories_raw, list):
            return normalize_categories(categories_raw)

    if not categories_list:
        logger.error(
            "Cannot determine category list from categories JSON. "
            "Expected a list of strings, a list of objects, or a dict containing a list."
        )
        sys.exit(1)

    seen: set[str] = set()
    deduped: list[str] = []
    for category in categories_list:
        if not isinstance(category, str):
            continue
        cat = category.strip()
        if cat and cat not in seen:
            seen.add(cat)
            deduped.append(cat)

    if not deduped:
        logger.error("No valid category names could be extracted from categories JSON.")
        sys.exit(1)

    return deduped


# ---------------------------------------------------------------------------
# LLM generation
# ---------------------------------------------------------------------------


def build_chat_prompt(
    tokenizer,
    system_prompt: str,
    user_prompt: str,
) -> str:
    """
    Build a prompt string using the tokenizer's chat template if available,
    otherwise fall back to a simple format.
    """
    messages = [
        {"role": "user", "content": f"{system_prompt}\n\n{user_prompt}"},
    ]

    # Many instruction-tuned models have a chat template
    if hasattr(tokenizer, "chat_template") and tokenizer.chat_template:
        try:
            return tokenizer.apply_chat_template(
                messages, tokenize=False, add_generation_prompt=True
            )
        except Exception:
            pass

    # Fallback for models without chat templates
    return (
        f"### System:\n{system_prompt}\n\n### User:\n{user_prompt}\n\n### Assistant:\n"
    )


def generate_response(
    model,
    tokenizer,
    system_prompt: str,
    user_prompt: str,
    temperature: float = 0.1,
    max_new_tokens: int = 4096,
) -> str:
    """Generate a response from the model."""
    prompt = build_chat_prompt(tokenizer, system_prompt, user_prompt)

    inputs = tokenizer(prompt, return_tensors="pt").to(model.device)

    gen_kwargs = {
        "max_new_tokens": max_new_tokens,
        "do_sample": temperature > 0,
        "pad_token_id": tokenizer.eos_token_id,
    }
    if temperature > 0:
        gen_kwargs["temperature"] = temperature
        gen_kwargs["top_p"] = 0.9

    with torch.no_grad():
        output_ids = model.generate(**inputs, **gen_kwargs)

    # Decode only the newly generated tokens
    new_tokens = output_ids[0, inputs["input_ids"].shape[1] :]
    response = tokenizer.decode(new_tokens, skip_special_tokens=True).strip()

    return response


# ---------------------------------------------------------------------------
# Classification functions
# ---------------------------------------------------------------------------


def classify_single_protein(
    accession: str,
    protein_names: list[str],
    taxonomy: list[str] | None,
    categories_list: list[str],
    valid_categories: set[str],
    model,
    tokenizer,
    temperature: float,
    max_new_tokens: int,
) -> dict:
    """Classify a single protein and return a result dict."""
    categories_str = "\n".join(f"  - {c}" for c in categories_list)
    names_str = "\n".join(f'  - "{n}"' for n in protein_names)
    taxonomy_str = format_taxonomy(taxonomy)

    prompt = SINGLE_PROTEIN_PROMPT.format(
        categories=categories_str,
        accession=accession,
        taxonomy=taxonomy_str,
        protein_names=names_str,
    )

    raw = generate_response(
        model, tokenizer, SYSTEM_PROMPT, prompt, temperature, max_new_tokens
    )
    parsed = extract_json_from_response(raw)

    if parsed and isinstance(parsed, dict) and "categories" in parsed:
        cats = validate_categories(
            parsed["categories"], valid_categories, accession=accession
        )
        return {"accession": accession, "categories": cats, "raw_response": raw}

    logger.warning(f"Failed to parse response for {accession}. Raw: {raw[:200]}...")
    return {"accession": accession, "categories": [], "raw_response": raw}


def classify_batch(
    batch: list[tuple[str, list[str], list[str] | None]],
    categories_list: list[str],
    valid_categories: set[str],
    model,
    tokenizer,
    temperature: float,
    max_new_tokens: int,
) -> list[dict]:
    """Classify a batch of proteins in a single LLM call."""
    categories_str = "\n".join(f"  - {c}" for c in categories_list)

    proteins_lines = []
    for accession, names, taxonomy in batch:
        names_joined = "; ".join(f'"{n}"' for n in names)
        taxonomy_joined = " > ".join(taxonomy) if taxonomy else "Unknown"
        proteins_lines.append(
            f'  Accession "{accession}": names=[{names_joined}] | taxonomy="{taxonomy_joined}"'
        )
    proteins_block = "\n".join(proteins_lines)

    prompt = BATCH_PROTEIN_PROMPT.format(
        categories=categories_str,
        proteins_block=proteins_block,
    )

    raw = generate_response(
        model, tokenizer, SYSTEM_PROMPT, prompt, temperature, max_new_tokens
    )
    parsed = extract_json_from_response(raw)

    results = []

    if parsed and isinstance(parsed, list):
        parsed_map = {}
        for item in parsed:
            if isinstance(item, dict) and "accession" in item:
                acc = item["accession"]
                cats = validate_categories(
                    item.get("categories", []),
                    valid_categories,
                    accession=acc,
                )
                parsed_map[acc] = cats

        for accession, names, taxonomy in batch:
            if accession in parsed_map:
                results.append(
                    {
                        "accession": accession,
                        "categories": parsed_map[accession],
                        "raw_response": None,
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
                    taxonomy,
                    categories_list,
                    valid_categories,
                    model,
                    tokenizer,
                    temperature,
                    max_new_tokens,
                )
                results.append(result)
    else:
        logger.warning(
            "Batch response could not be parsed. Falling back to per-protein "
            "classification for this batch."
        )
        for accession, names, taxonomy in batch:
            result = classify_single_protein(
                accession,
                names,
                taxonomy,
                categories_list,
                valid_categories,
                model,
                tokenizer,
                temperature,
                max_new_tokens,
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
    model_name: str = "google/gemma-4-31b-it",
    batch_size: int = 10,
    save_every: int = 1000,
    temperature: float = 0.1,
    max_new_tokens: int = 4096,
    device: str = "auto",
    quantize: str | None = None,
    hf_token: str | None = None,
) -> None:
    """Run the full classification pipeline."""

    # --- Load inputs -------------------------------------------------------
    logger.info(f"Loading proteins from {proteins_path}")
    proteins_data = load_json(proteins_path)
    logger.info(f"Loaded {len(proteins_data)} proteins")

    logger.info(f"Loading categories from {categories_path}")
    categories_raw = load_json(categories_path)
    categories_list = normalize_categories(categories_raw)

    valid_categories = set(categories_list)
    logger.info(f"Using {len(categories_list)} predefined categories")

    if save_every < 1:
        logger.error("--save-every must be >= 1")
        sys.exit(1)

    if batch_size != 1:
        logger.info(
            "Single-prompt mode is enabled for maximum accuracy; "
            f"ignoring --batch-size={batch_size} and using per-protein inference."
        )

    # --- Load model --------------------------------------------------------
    model, tokenizer = load_model(model_name, device, quantize, hf_token)

    # --- Resume state ------------------------------------------------------
    output_data = load_existing_output(output_path)
    processed_accessions = {
        acc for acc, info in output_data.items() if isinstance(info, dict)
    }
    if processed_accessions:
        logger.info(
            f"Resuming run: skipping {len(processed_accessions)} already-saved proteins"
        )

    # --- Classify ----------------------------------------------------------
    errors: list[str] = []
    newly_processed = 0
    total = len(proteins_data)

    for accession, info in tqdm(
        proteins_data.items(), desc="Classifying", unit="protein", total=total
    ):
        if accession in processed_accessions:
            continue

        protein_names, taxonomy = coerce_protein_fields(info)
        try:
            result = classify_single_protein(
                accession,
                protein_names,
                taxonomy,
                categories_list,
                valid_categories,
                model,
                tokenizer,
                temperature,
                max_new_tokens,
            )

            output_data[result["accession"]] = {
                "protein_names": protein_names,
                "taxid": info.get("taxid"),
                "taxonomy": taxonomy if taxonomy is not None else [],
                "original_categories": info.get("categories", []),
                "predicted_categories": result["categories"],
            }
            newly_processed += 1

            if newly_processed % save_every == 0:
                logger.info(
                    f"Checkpoint: saving after {newly_processed} newly processed proteins"
                )
                save_json(output_data, output_path)
        except Exception as e:
            errors.append(accession)
            logger.error(f"Classification failed for {accession}: {e}")

    # --- Build complete output snapshot -----------------------------------
    for accession, info in proteins_data.items():
        if accession in output_data:
            continue
        protein_names, taxonomy = coerce_protein_fields(info)
        output_data[accession] = {
            "protein_names": protein_names,
            "taxid": info.get("taxid"),
            "taxonomy": taxonomy if taxonomy is not None else [],
            "original_categories": info.get("categories", []),
            "predicted_categories": [],
        }

    save_json(output_data, output_path)

    # --- Summary -----------------------------------------------------------
    classified = sum(1 for v in output_data.values() if v["predicted_categories"])
    logger.info(
        f"Classification complete: {classified}/{total} proteins received "
        f"at least one category."
    )
    logger.info(f"Newly processed in this run: {newly_processed}")
    if errors:
        logger.warning(f"{len(errors)} proteins encountered errors: {errors}")

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
        "a HuggingFace Transformers LLM.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--proteins",
        required=True,
        help="Path to JSON file with protein data.",
    )
    parser.add_argument(
        "--categories",
        required=True,
        help="Path to JSON file with predefined functional categories.",
    )
    parser.add_argument(
        "--output",
        default="classification_results.json",
        help="Path for the output JSON (default: classification_results.json).",
    )
    parser.add_argument(
        "--model",
        default="google/gemma-4-31b-it",
        help="HuggingFace model ID (default: google/gemma-4-31b-it). "
        "Other good options: mistralai/Mistral-7B-Instruct-v0.3, "
        "meta-llama/Llama-3.1-8B-Instruct, "
        "deepseek-ai/DeepSeek-R1-Distill-Llama-8B.",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=10,
        help="Proteins per LLM call (default: 10). Use 1 for max reliability.",
    )
    parser.add_argument(
        "--save-every",
        type=int,
        default=1000,
        help=(
            "Save progress to --output after every N newly processed proteins "
            "(default: 1000)."
        ),
    )
    parser.add_argument(
        "--temperature",
        type=float,
        default=0.1,
        help="Sampling temperature (default: 0.1).",
    )
    parser.add_argument(
        "--max-new-tokens",
        type=int,
        default=4096,
        help="Max new tokens to generate per call (default: 4096).",
    )
    parser.add_argument(
        "--device",
        default="auto",
        help="Device: 'auto', 'cuda', 'cuda:0', 'cpu' (default: auto).",
    )
    parser.add_argument(
        "--quantize",
        choices=["4bit", "8bit"],
        default=None,
        help="Quantise the model to reduce VRAM. Requires bitsandbytes.",
    )
    parser.add_argument(
        "--hf-token",
        default=None,
        help="HuggingFace access token for gated models (e.g. Gemma). "
        "Can also set HF_TOKEN env var.",
    )

    args = parser.parse_args()

    # Allow token from environment
    hf_token = args.hf_token
    if hf_token is None:
        hf_token = os.environ.get("HF_TOKEN")

    run_classification(
        proteins_path=args.proteins,
        categories_path=args.categories,
        output_path=args.output,
        model_name=args.model,
        batch_size=args.batch_size,
        save_every=args.save_every,
        temperature=args.temperature,
        max_new_tokens=args.max_new_tokens,
        device=args.device,
        quantize=args.quantize,
        hf_token=hf_token,
    )


if __name__ == "__main__":
    main()
