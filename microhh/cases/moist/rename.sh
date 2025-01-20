#!/bin/bash
for file in drycbl*; do
  # Replace 'drycbl' with 'moistcbl' in the filename
  new_file=$(echo "$file" | sed 's/^drycbl/moistcbl/')
  mv "$file" "$new_file"
done

