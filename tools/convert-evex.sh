#!/bin/bash

# Convert EVEX standoff data to PubAnnotation JSON.

# Assumes .tar.gz packages as input, creates .zip

set -u
set -e
set -o pipefail

SCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
OUTDIR="$( pwd )"
SO2PA="$SCRIPTDIR/../standoff2pa.py"
TYPEMAP="$SCRIPTDIR/../typemap-evex.txt"

if [ "$#" -lt 1 ]; then
    echo "Usage: $0 FILE [FILE [...]]" >&2
    exit 1
fi

TMPDIR=`mktemp -d tmp.XXXXXXXXXX`
function cleanup {
    rm -rf "$TMPDIR"
}
trap cleanup EXIT

for f in "$@"; do
    echo "Unpacking $f ..." >&2
    tar xzf "$f" -C "$TMPDIR"

    echo "Converting $f ..." >&2
    find "$TMPDIR" -name '*.ann' | \
	xargs python "$SO2PA" -m "$TYPEMAP" -o "$TMPDIR"
    
    o=`basename "$f" .tar.gz`.tar.gz
    echo "Packaging into $o ..." >&2
    cd "$TMPDIR"
    ls | egrep '\.json$' > filelist
    tar czf "$OUTDIR/$o" -T filelist
    cd ..

    echo "Cleaning up for $f ..." >&2
    rm -rf "$TMPDIR"
    mkdir "$TMPDIR"
done
