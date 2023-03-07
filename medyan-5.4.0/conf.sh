#!/bin/sh -e

# System requirements for this script:
#   - git >= 2.7.0
#   - cmake >= 3.12.0
#   - Latest C++ compiler

# Set directories and files
export medyan_root_dir="$(X= cd -- "$(dirname -- "$0")" && pwd -P)"

# Run bootstrapping
echo "Bootstrapping..."
. "$medyan_root_dir/scripts/bootstrap.sh"
