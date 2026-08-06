#!/bin/sh
# Install the exact open-source AmberClassic/NAB fallback used by Phase 2A.
set -eu

repository=https://github.com/Amber-MD/AmberClassic.git
commit=bdb3e0dee5b90f2be2950e26cfad1ae5a7440cae

if [ "$#" -ne 1 ]; then
    echo "usage: $0 TARGET_DIRECTORY" >&2
    exit 2
fi

target=$1
if [ -e "$target" ]; then
    echo "refusing to replace existing path: $target" >&2
    exit 2
fi

git clone --no-checkout "$repository" "$target"
git -C "$target" checkout --detach "$commit"
(
    cd "$target"
    ./configure --noboost
    make -j1 LEX=flex install
)

test -x "$target/bin/nab"
test -x "$target/bin/teLeap"
actual=$(git -C "$target" rev-parse HEAD)
test "$actual" = "$commit"

echo "AmberClassic NAB installed at $target"
echo "export NUCLEIC_BUILDER_AMBERCLASSIC_HOME=$target"
