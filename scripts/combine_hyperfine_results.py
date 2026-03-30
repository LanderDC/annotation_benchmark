#!/usr/bin/env python3
import argparse
import json
from pathlib import Path


def read_gpu_max_mem_mb(path: Path) -> float | None:
    if not path.exists():
        return None
    value = path.read_text(encoding="utf-8").strip()
    if not value:
        return None
    try:
        return float(value)
    except ValueError:
        return None


def parse_hyperfine_json(json_path: Path) -> list[dict]:
    data = json.loads(json_path.read_text(encoding="utf-8"))
    results = data.get("results", [])

    script_name = json_path.parent.name
    step_name = json_path.stem
    gpu_max_mem_mb = read_gpu_max_mem_mb(json_path.with_suffix(".gpu_max_mem_mb.txt"))
    gpu_trace_path = json_path.with_suffix(".gpu_mem_mb.tsv")

    parsed = []
    for result in results:
        parsed.append(
            {
                "script": script_name,
                "step": step_name,
                "command_name": result.get("command"),
                "timing": {
                    "mean_seconds": result.get("mean"),
                    "stddev_seconds": result.get("stddev"),
                    "min_seconds": result.get("min"),
                    "max_seconds": result.get("max"),
                    "median_seconds": result.get("median"),
                    "user_seconds": result.get("user"),
                    "system_seconds": result.get("system"),
                    "times_seconds": result.get("times", []),
                },
                "memory": {"memory_usage_bytes": result.get("memory_usage_bytes")},
                "gpu_memory": {
                    "max_memory_mb": gpu_max_mem_mb,
                },
                "source_files": {
                    "hyperfine_json": str(json_path),
                    "gpu_max_mem_mb_txt": str(
                        json_path.with_suffix(".gpu_max_mem_mb.txt")
                    ),
                    "gpu_memory_trace_tsv": str(gpu_trace_path),
                },
            }
        )

    return parsed


def combine_results(input_dir: Path, output_path: Path) -> dict:
    rows = []
    for json_path in sorted(input_dir.rglob("*.json")):
        if json_path.resolve() == output_path.resolve():
            continue
        rows.extend(parse_hyperfine_json(json_path))

    combined = {
        "n_records": len(rows),
        "input_dir": str(input_dir),
        "records": rows,
    }

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(combined, indent=2), encoding="utf-8")
    return combined


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Combine hyperfine benchmark outputs into one JSON file."
    )
    parser.add_argument("--input-dir", default="results/hyperfine", type=Path)
    parser.add_argument(
        "--output", default="results/hyperfine/hyperfine_combined.json", type=Path
    )
    args = parser.parse_args()

    combine_results(args.input_dir, args.output)


if __name__ == "__main__":
    main()
