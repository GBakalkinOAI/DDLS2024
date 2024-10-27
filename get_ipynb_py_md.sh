#!/bin/bash

show_help() {
  echo "Usage: get_ipynb_py_md.sh OUTPUT_DIR SOURCE_DIR1 [SOURCE_DIR2] ..."
  echo "Collect all .ipynb, .py and .md files in SOURCE_DIR(s) in OUTPUT_DIR"
  echo "Each .ipynb file is converted into .py to save space and remove figures, which needs 'mamba activate ipynb2py'"
  echo "Next step: Manually remove unnecessary files and run combine_py_md.sh"
}

if [[ "$1" == "--help" ]]; then
  show_help
  exit 0
fi

if [[ $# -lt 2 ]]; then
  show_help
  exit 1
fi

OUTPUT_DIR=$(realpath "$1")
mkdir -p "$OUTPUT_DIR"
shift

# Iterate over all source directories
for SOURCE_DIR in "$@"; do
  SOURCE_DIR=$(realpath "$SOURCE_DIR")
  echo "Getting *.ipynb *.py *.md from ${SOURCE_DIR} to ${OUTPUT_DIR} sluggified..."

  # get all .ipynb files and convert them to sluggified .py files
  find "$SOURCE_DIR" -name "*.ipynb" -print | while read -r file; do
    file=$(realpath --relative-to="$SOURCE_DIR" "$file")
    slug_name=$(echo "$file" | sed 's/ /-/g; s/\//-/g; s/\.ipynb//g')
    jupyter nbconvert --to script "$SOURCE_DIR/$file" --output "${OUTPUT_DIR}/${slug_name}"
  done

  # get all .py and .md files and convert them to sluggified .py files
  find "$SOURCE_DIR" \( -name "*.md" -o -name "*.py" \) -print | while read -r file; do
    file=$(realpath --relative-to="$SOURCE_DIR" "$file")
    slug_name=$(echo "$file" | sed 's/ /-/g; s/\//-/g; s/\.ipynb//g')
    # echo "cp  ${SOURCE_DIR}/${file}  ${OUTPUT_DIR}/${slug_name}"
    cp "$SOURCE_DIR/$file" "${OUTPUT_DIR}/${slug_name}"
  done
done

# Remove duplicate files with identical content, keeping the one with the shortest name
find "$OUTPUT_DIR" -type f \( -name "*.md" -o -name "*.py" \) | while read -r file1; do
  find "$OUTPUT_DIR" -type f \( -name "*.md" -o -name "*.py" \) | while read -r file2; do
    if [[ "$file1" != "$file2" && -e "$file1" && -e "$file2" ]]; then
      if cmp -s "$file1" "$file2"; then
        if [[ "${#file1}" -le "${#file2}" ]]; then
          echo "Removing duplicate file: $file2"
          rm "$file2"
        else
          echo "Removing duplicate file: $file1"
          rm "$file1"
          break
        fi
      fi
    fi
  done
done