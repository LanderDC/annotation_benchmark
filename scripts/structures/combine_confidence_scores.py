import argparse
import json
import shutil
from pathlib import Path


def copy_confidence_files(
    project_root: Path, source_glob: str, destination_directory: Path
) -> list[Path]:
    destination_directory.mkdir(parents=True, exist_ok=True)

    copied_files: list[Path] = []
    for source_path in sorted(project_root.glob(source_glob)):
        if not source_path.is_file():
            continue

        # Avoid filename collisions by prefixing with the parent folder name
        # (folder typically identifies the protein/sample in Boltz output).
        destination_name = f"{source_path.parent.name}_{source_path.name}"
        destination_path = destination_directory / destination_name

        shutil.copy2(source_path, destination_path)
        copied_files.append(destination_path)

    return copied_files


def combine_json_files(input_directory: Path, output_file: Path) -> dict:
    combined_data: dict = {}

    for file_path in sorted(input_directory.glob("*.json")):
        if not file_path.is_file():
            continue

        with file_path.open("r", encoding="utf-8") as file:
            content = json.load(file)
            key = file_path.stem
            combined_data[key] = content

    output_file.parent.mkdir(parents=True, exist_ok=True)
    output_file.write_text(json.dumps(combined_data, indent=4), encoding="utf-8")
    return combined_data


def build_default_paths() -> tuple[Path, Path, str]:
    project_root = Path(__file__).resolve().parents[2]
    confidence_dir = project_root / "results" / "boltz" / "confidence_scores"
    output_json = project_root / "results" / "boltz" / "combined_plddt_scores.json"
    source_glob = (
        "results/boltz/predictions/boltz_results_inputs/predictions/*/confidence_*.json"
    )
    return confidence_dir, output_json, source_glob


def main() -> dict:
    default_confidence_dir, default_output_json, default_source_glob = (
        build_default_paths()
    )

    parser = argparse.ArgumentParser(
        description="Copy Boltz confidence JSON files and combine them into one JSON file."
    )
    parser.add_argument(
        "--source-glob",
        default=default_source_glob,
        help="Glob for source confidence files, relative to --project-root.",
    )
    parser.add_argument(
        "--confidence-dir",
        default=str(default_confidence_dir),
        type=Path,
        help="Directory where confidence JSON files will be copied.",
    )
    parser.add_argument(
        "--output",
        default=str(default_output_json),
        type=Path,
        help="Path for the combined JSON output.",
    )
    parser.add_argument(
        "--project-root",
        default=str(Path(__file__).resolve().parents[2]),
        type=Path,
        help="Project root used to resolve --source-glob.",
    )
    args = parser.parse_args()

    copied_files = copy_confidence_files(
        project_root=args.project_root,
        source_glob=args.source_glob,
        destination_directory=args.confidence_dir,
    )
    combined = combine_json_files(
        input_directory=args.confidence_dir, output_file=args.output
    )

    return {
        "copied_files": [str(path) for path in copied_files],
        "combined_output": str(args.output),
        "n_combined": len(combined),
    }


if __name__ == "__main__":
    main()
