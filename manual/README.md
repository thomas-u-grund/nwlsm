# The nwlsm manual

`nwlsm_manual.pdf` is the user manual and reference of nwlsm (title "Manual", no version number in the citation). The LaTeX source is `nwlsm_manual.tex` with one file per chapter (`ch*.tex`); the listings of the examples are `ex/out/*.txt`, the figures `fig/*.pdf`.

## Bumping the version (one-line change)

The manual's title and citation carry no version number; the title page states the series it describes (`\manualseries`, valid for 1.0.x). The exact version and date of the release that produced the listings are two more macros in `nwlsm_manual.tex`, and the tag appears in the installation sentence of `ch01_intro.tex`:

    \newcommand{\manualseries}{1.0}
    \newcommand{\pkgversion}{1.0.0}
    \newcommand{\pkgdate}{20 November 2026}

Change `\manualseries` only for a new minor or major release whose documented behaviour changed. For a patch release rebuild (`make -C manual pdf`) and rerun the examples if their output changed.

## The DOI

`\manualdoi` is empty until the Zenodo record exists; an empty DOI prints nothing in the citation. `ZENODO.md` holds the metadata of the record and the steps to do.

## Building

    make -C manual pdf        # pdflatex, 3 passes
    make -C manual check      # undefined references, overfull boxes
    NWCORE=/path/to/nwcommands make -C manual examples   # reruns the do-files of ex/ in Stata 17 (several minutes; the Zachary, Sampson and Lazega data come from the web)

`NWLSM_ROOT` is the directory of the package files that the examples use (default: this repository; set it to a checkout of a release tag to reproduce that release). The do-files `ex/ex01_karate.do` to `ex/ex05_compare.do` write one text log per listing into `ex/out` with `log using`; `ex/mklistings.py` removes the log header and footer lines. The listings use fixed seeds; the timings in `ex/out/lazega_time.txt` depend on the computer and its load. The data of `ex04_influence.do` are simulated in the do-file.
