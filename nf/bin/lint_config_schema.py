#!/usr/bin/env python3
"""Lightweight lint checks for nextflow.config and nextflow_schema.json."""

import json
import re
import sys
from pathlib import Path


def _read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def _extract_config_value(config_text: str, key: str):
    pattern = rf"^\s*{re.escape(key)}\s*=\s*['\"]([^'\"]+)['\"]"
    match = re.search(pattern, config_text, re.MULTILINE)
    return match.group(1) if match else None


def main() -> int:
    repo_nf = Path(__file__).resolve().parents[1]
    config_path = repo_nf / "nextflow.config"
    schema_path = repo_nf / "nextflow_schema.json"

    errors = []
    config_text = _read_text(config_path)
    schema = json.loads(_read_text(schema_path))

    config_version = _extract_config_value(config_text, "version")
    config_chemistry = _extract_config_value(config_text, "chemistry")
    schema_chemistry = schema.get("properties", {}).get("chemistry", {}).get("default")

    if not config_version:
        errors.append("manifest.version not found in nextflow.config")
    if not config_chemistry:
        errors.append("params.chemistry not found in nextflow.config")
    if config_chemistry != schema_chemistry:
        errors.append(
            f"chemistry default mismatch: nextflow.config={config_chemistry!r}, "
            f"nextflow_schema.json={schema_chemistry!r}"
        )

    for workflow_path in [
        repo_nf / "subworkflows" / "rna_met" / "main.nf",
        repo_nf / "subworkflows" / "met_only" / "main.nf",
        repo_nf / "subworkflows" / "forcecell" / "main.nf",
    ]:
        text = _read_text(workflow_path)
        if config_version and config_version not in text:
            errors.append(f"{workflow_path.relative_to(repo_nf)} does not mention version {config_version}")

    if errors:
        print("Lint failed:", file=sys.stderr)
        for error in errors:
            print(f"- {error}", file=sys.stderr)
        return 1

    print("Lint passed: config/schema chemistry and workflow version references are consistent.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
