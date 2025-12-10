Memory safety notes for fit_veg_mixture_mesma.R

Problem:
- A full Cartesian expansion of per-vegetation candidate variants can create an enormous number of combinations and exhaust available RAM ("cannot allocate vector of size ... GB").

Mitigations implemented in code:
- evaluate_all_combinations now avoids materializing the entire expand.grid when the number of combos exceeds COMBO_SAFE_EXPAND_LIMIT (default 1e6). It generates combinations lazily using a mixed-radix counter and processes them in chunks to keep memory usage low.
- If the combination space is extremely large (greater than COMBO_ABORT_LIMIT, default 5e7), the code will abort with an informative message to reduce search space.
- New config variables can be tuned in your environment before running the script: COMBO_SAFE_EXPAND_LIMIT and COMBO_ABORT_LIMIT.

Practical recommendations:
- Reduce TOPK_VARIANTS (default TOPK_VARIANTS) to a smaller number (e.g., 2-5) before running when you have many variant candidates per vegetation.
- Reduce MAX_VEG_COMPONENTS if you have many vegetation classes in the library.
- If you need to search very large spaces, consider sampling strategies, heuristic search (greedy/local search), or approximate algorithms — exhaustive search beyond millions to tens of millions of combos is generally impractical on a workstation.

If you want, I can add a configurable sampling-based fallback (random sampling of combinations) or implement a heuristic greedy search to keep compute bounded while still producing good results.

New preprocessing step (soil-first MESMA):
- The pipeline now supports constructing a soil prototype from geojson-labeled 'no soil' points and subtracting the estimated soil fraction from each observation's spectral signature before vegetation library construction.
- This is controlled by ENABLE_SOIL_PREPROCESS (default TRUE). Parameters SOIL_PURE_THRESHOLD and SOIL_MIN_SAMPLES let you tune how pure the soil set is selected.
- After subtraction, the original values are preserved under `raw_<index>` columns and the estimated per-location-year soil fraction is stored in `soil_frac` (constant for all rows within a location-year).

Important: there is no per-row soil-subtraction fallback. If the provisional PCA cannot be built (insufficient data) or the per-location-year alpha estimator fails, the script will skip soil subtraction rather than try a per-row subtraction.

This helps removing soil contamination from vegetation endmember extraction and preserves the time-dimension reduction (PCA) downstream.

Plotting note: Plots that overlay vegetation and barren/soil fractions (e.g. the global coverage plots) now compute an adaptive scaling factor so that barren and vegetation series are displayed on the same figure but with separate primary and secondary y-axes. This avoids hard-coded scale values and makes the visualization readable when barren and vegetation are on very different ranges.
