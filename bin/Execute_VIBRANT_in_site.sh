# This script has to be executed from a folder containing .fna files
# it executes VIBRANT using 10 cores

ls -1 | grep "fna" | while IFS= read -r file; do echo VIBRANT_run.py -i "$file" -t 10; VIBRANT_run.py -i "$file" -t 10; done
