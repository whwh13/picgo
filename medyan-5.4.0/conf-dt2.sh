#!/bin/sh -e

# The configuration script for Deepthought2 HPC cluster
#
# Before using this script, the following modules need to be loaded
#   gcc/9.1.0
#   boost
#   cmake/3.13.2

medyan_root_dir="$(X= cd -- "$(dirname -- "$0")" && pwd -P)"

# Install git for vcpkg (requires >= 2.7.0)
export PATH="~/bin:$PATH"
if [ "$(git version | cut -d"." -f2)" -lt 7 ]; then
    echo "git version is too old. Building a new git version..."
    (
        mkdir -p "$medyan_root_dir/scripts/.build" &&
        cd "$medyan_root_dir/scripts/.build" &&
        wget https://github.com/git/git/archive/v2.25.0.zip -O git.zip &&
        unzip -qq git.zip &&
        cd git-2.25.0 &&
        make -j10 &&
        make install
    )
fi

# Set up variables
export MEDYAN_NO_GUI="true"
export MEDYAN_ADDITIONAL_LINK_DIRS="$GCC_ROOTDIR/lib64"
export MEDYAN_RPATH="$GCC_ROOTDIR/lib64"

export MEDYAN_SPECIAL_ENVIRONMENT="Deepthought2"

# Run the script
. "$medyan_root_dir/conf.sh"
