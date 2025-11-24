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
