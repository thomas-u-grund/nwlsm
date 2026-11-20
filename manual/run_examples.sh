#!/bin/zsh
# Re-runs the example do-files of the manual with Stata 17 (-e) and rebuilds the listings in manual/ex/out.
# Usage: NWCORE=/path/to/nwcommands manual/run_examples.sh [ex01_karate ex02_monks ...]
# NWLSM_ROOT is the directory of the package files to test (default: the repository that holds manual/; use a checkout of the tag to reproduce a release).
# NWPLUS defaults to a fresh scratch directory, so that no installed copy shadows the repository build. The data sets (Zachary, Sampson, Lazega)
# are loaded with nwwebuse (internet access needed); the simulated data of ex04_influence are made in the do-file.
here="$(cd "$(dirname "$0")" && pwd)"
root=${NWLSM_ROOT:-$(dirname "$here")}
export NWPLUS=${NWPLUS:-$(mktemp -d)}
if [[ -z $NWCORE || ! -f $NWCORE/unw_core.do ]]; then echo "set NWCORE to the nwcommands (core) directory"; exit 1; fi
STATA=${STATA:-/Applications/Stata 17.0/StataBE.app/Contents/MacOS/StataBE}
cd "$here/ex" || exit 1
files=("$@"); [[ ${#files} -eq 0 ]] && files=(ex[0-9]*.do) || files=(${^files}.do)
mkdir -p out ../fig
for f in $files; do
  echo "running $f"
  { echo "capture program drop nwinstall"; echo "sysdir set PLUS \"$NWPLUS\""; echo "adopath ++ \"$NWCORE\""; echo "adopath ++ \"$NWCORE/lib\""
    echo "adopath ++ \"$root\""; echo "adopath ++ \"$root/lib\""; echo "adopath ++ \"$root/lib/plugins/macos\""; echo "mata: mata mlib index"
    echo "do $f"; } > run_$f
  perl -e 'alarm 3000; exec @ARGV' "$STATA" -e do run_$f > /dev/null 2>&1
  rm -f run_$f
  if grep -qE "^r\([0-9]+\)" run_${f%.do}.log; then echo "ERROR in $f (log kept: ex/run_${f%.do}.log)"; else rm -f run_${f%.do}.log; fi
done
python3 mklistings.py
