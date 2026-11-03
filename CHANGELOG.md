# Changelog

## 2.1.0

- The `R9` preset now uses the R10 patterns (the canonical repeat and every 1 bp substitution of it) as well as the Guppy miscall patterns. The Guppy patterns alone miss variant repeats common in the proximal telomere, such as TGAGGG, so the R9 preset could call the boundary up to about 500 bp short on G-strand reads. **R9 results differ from 2.0.0.** `R10` and `R9_archive` are unchanged.

## 2.0.0

- TeloNP now has presets for nanopore chemistry, chosen with `chemistry` in `getTeloNPBoundary` or `--chemistry` in `teloBPCmd.py`:
  - `R10` (new default): canonical repeats plus 1 bp variants, and a 2 kb limit on stepping over gaps in the telomere, so an interstitial telomeric sequence (ITS) far from the telomere is no longer counted as telomere.
  - `R9`: the Guppy miscall patterns with the same 2 kb gap limit.
  - `R9_archive`: the TeloNP settings from v1, unchanged. **To reproduce v1 results, pass `chemistry="R9_archive"`.**
- `getTeloBoundary` has a new `maxDiscontinuityGap` argument for the bounded gap search. Its default (0) keeps the old behaviour.
- `teloBPCmd.py` writes `teloBP_run_metadata.json` with the version and full settings of each run.
- New `getTeloNPSettings()` and `getTeloBPSettings()` return the settings a call uses.

## 1.0.0

- First release.
