# Changelog

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
