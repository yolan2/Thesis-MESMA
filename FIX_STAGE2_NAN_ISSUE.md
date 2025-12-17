# Stage 2 Signal Loss - Root Cause and Fix

## Problem
Stage 2 is being skipped with message: "insufficient Stage2 signal (norm=0)"

## Root Cause
The vegetation signal (`y_s2_masked`) has norm=0 after Stage 1 subtraction, normalization, masking, and weighting.

Looking at the log output, most locations show very low vegetation fractions:
- `veg_frac=0.0000` (pure barren)  
- `veg_frac=0.0698` (mostly barren)
- `veg_frac=0.2020` (mostly barren)

When Stage 1 estimates very little vegetation, the vegetation signal for Stage 2 is near-zero.

## Key Question: Is This Correct or a Bug?

Two possibilities:
1. **Correct**: There really is minimal vegetation, so Stage 2 has nothing to unmix
2. **Bug**: The masking/weighting is incorrectly zeroing out the vegetation signal

## Investigation: Check What's Happening

The signal flow is:
1. `y_veg_raw = y_raw - barren_frac * RAW_BARREN_PROTO` 
2. Normalize: `y_s2 = (y_veg_raw - mu) / sigma`
3. **Mask**: `y_s2_masked = y_s2[valid_mask]` 
4. **Weight**: `y_s2_masked = y_s2_masked * sqrt(weights_s2_masked)`
5. Compute norm: `sqrt(sum(y_s2_masked^2))` → This is 0!

## The Most Likely Culprit: Masking + Weighting Interaction

Your suspicion was correct! The issue is likely the combination of:
1. **Sparse observations** (only 7-14 valid pentads out of 73)
2. **Stage 2 weights** applied to masked positions

When you have:
- Few valid observations (sparse mask)
- Low weights on those specific masked positions
- The weighted masked signal becomes near-zero even if there IS vegetation signal

## The Fix: Don't Apply Position-Specific Weights to Masked Stage 2 Data

Change line 6402 from:
```r
y_s2_masked <- y_s2_masked * sqrt(pmax(weights_s2_masked, 0))
```

To use a global average weight instead:
```r
# Use global mean weight instead of position-specific weights for masked data
mean_weight_s2 <- mean(STAGE2_PARAMS$weights, na.rm = TRUE)
y_s2_masked <- y_s2_masked * sqrt(mean_weight_s2)
```

This way:
- Stage 2 still benefits from the PCA-LDA discriminative information (via library weighting)
- But the observation vector isn't unfairly penalized at sparse masked positions

