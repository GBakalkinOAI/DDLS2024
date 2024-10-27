#!/bin/bash

show_help() {
  echo "Usage: combine_py_md.sh combined_file SOURCE_DIR1 [SOURCE_DIR2] ..."
  echo "combined_file.py file is created, where each .py file from SOURCE_DIR(s) is a section."
  echo "combined_file.md file is created, where each .md file from SOURCE_DIR(s) is a section."
}

if [[ "$1" == "--help" ]]; then
  show_help
  exit 0
fi

if [[ $# -lt 2 ]]; then
  show_help
  exit 1
fi

# Create or empty the combined files
combined_file_base=$(realpath "$1" | sed 's/\.md$//' | sed 's/\.py$//')
combined_file_py="${combined_file_base}.py"
> "$combined_file_py"
combined_file_md="${combined_file_base}.md"
> "$combined_file_md"
shift

combine_files() {
  local extension="$1"
  local combined_file="$2"
  shift 2

  for SOURCE_DIR in "$@"; do
    SOURCE_DIR=$(realpath "$SOURCE_DIR")
    find "$SOURCE_DIR" -name "*.$extension" ! -path "$combined_file" -print | while read -r file; do
      echo -e "\n# Section: $(basename "$file" .$extension)\n" >> "$combined_file"
      cat "$file" >> "$combined_file"
      echo -e "\n" >> "$combined_file" # Add extra newline for separation
    done
  done
}

# Combine Python and Markdown files
combine_files "py" "$combined_file_py" "$@"
combine_files "md" "$combined_file_md" "$@"

echo "Combined python file created at: $combined_file_py"
echo "Combined markdown file created at: $combined_file_md"