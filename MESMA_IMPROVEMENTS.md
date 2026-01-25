# Potential Improvements for MESMA Methodology

Based on an analysis of `fit_veg_mixture_mesma.R`, here are three key areas where the MESMA methodology can be improved (excluding the spline-based approach).

## 1. Endmember Selection: Shift from SAM to EAR
**Current State:** The function `sam_extract_endmembers` uses a greedy Spectral Angle Mapper (SAM) approach. It selects the pixel closest to the mean, and then iteratively adds pixels that have the largest spectral angle (maximum difference) from the current set.
**Limitation:** SAM tends to select "extreme" pixels (vertices of the convex hull). While distinct, these "purest" pixels may not be the most *representative* of the class variability for unmixing purposes, potentially leading to higher residual errors when modeling "average" pixels.
**Improvement:** Implement **Endmember Average RMSE (EAR)** selection.
*   **Method:** For a given class and size $k$, select the subset of $k$ endmembers that minimizes the average Root Mean Square Error (RMSE) when unmixing all other samples in that class.
*   **Benefit:** Produces a library that is statistically more representative of the intra-class variance, often improving unmixing accuracy for non-extreme samples.

## 2. Global Optimization: Heuristic Search instead of Grid Search
**Current State:** The `optimize_library` function performs a "Global Combinatorial Optimization" using `expand.grid` to test all combinations of cluster counts ($k$) across vegetation types.
**Limitation:** This approach suffers from the "curse of dimensionality". Adding more vegetation classes or increasing the range of $k$ (e.g., 1-20 instead of 1-11) makes the grid size explode, forcing the use of random sampling (`sample(seq_len(n_combos), 500)`), which guarantees sub-optimality.
**Improvement:** Replace Grid Search with a **Genetic Algorithm (GA)** or **Sequential Forward Floating Selection (SFFS)**.
*   **Method:** Use a GA to evolve the vector of $k$ values (e.g., `[3, 5, 2]`) by maximizing the OOB Kappa or Accuracy.
*   **Benefit:** Efficiently searches a much larger hyperparameter space and converges to better optima than random sampling, without the computational cost of a full grid search.

## 3. Solver Efficiency: Quadprog/OSQP instead of Augmented NNLS
**Current State:** The `solve_weights_fcls` function enforces the sum-to-one constraint by "augmenting" the matrix (adding a row of 1s weighted by a large lambda) and using `nnls::nnls`.
**Limitation:** Augmented NNLS is robust but can be slower and numerically less stable than dedicated Quadratic Programming (QP) solvers for large problems. The code defines `ENABLE_QP_SOLVER` in the config but `solve_weights_fcls` currently hardcodes `nnls::nnls`.
**Improvement:** Implement a dedicated FCLS solver using `quadprog` or `osqp`.
*   **Method:** Formulate the unmixing as a QP problem: minimize $\frac{1}{2}x^T Q x + c^T x$ subject to $Ax = b$ (sum-to-one) and $x \ge 0$.
*   **Benefit:** faster execution for the millions of unmixing operations required in the bootstrap/inference loops, with potentially better precision for the sum-to-one constraint.
