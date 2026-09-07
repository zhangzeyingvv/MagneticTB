#!/bin/sh

# Canonicalize every editable documentation source into a generated notebook.

set -eu

script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
project_root=$(dirname "$(dirname "$script_dir")")
codex_root=${CODEX_HOME:-"$HOME/.codex"}
guard=${MATHEMATICA_NOTEBOOK_GUARD:-"$codex_root/skills/mathematica-notebook/scripts/notebook_guard.py"}

if [ ! -f "$guard" ]; then
  echo "Notebook guard not found: $guard" >&2
  echo "Set MATHEMATICA_NOTEBOOK_GUARD to notebook_guard.py." >&2
  exit 2
fi

build_language() {
  source_dir=$1
  target_dir=$2

  find "$source_dir" -type f -name '*.source.wl' -print | sort |
    while IFS= read -r source_path; do
      relative_path=${source_path#"$source_dir"/}
      notebook_relative=${relative_path%.source.wl}.nb
      notebook_path="$target_dir/$notebook_relative"

      mkdir -p "$(dirname "$notebook_path")"
      echo "Canonicalizing $notebook_path"
      python3 "$guard" canonicalize "$source_path" "$notebook_path"
    done
}

build_language \
  "$script_dir/Sources" \
  "$project_root/Documentation/English"

build_language \
  "$script_dir/SourcesChineseSimplified" \
  "$project_root/Documentation/ChineseSimplified"

echo "All bilingual MagneticTB notebooks were generated."
