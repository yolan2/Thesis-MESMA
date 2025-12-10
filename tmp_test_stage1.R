# quick test script for stage1 unmix
src <- 'c:/Users/yolan/OneDrive/Documenten/UGENT/Master/R_MESMA/fit_veg_mixture_mesma.R'
source(src)

m1 <- c(0.1,0.2,0.3,0.4,0.5)
m2 <- c(0.5,0.4,0.3,0.2,0.1)
m1_large <- m1 * 100
y_scaled <- 0.6 * m1_large + 0.4 * m2
res <- geometric_stage1_unmix(y_scaled, m1_large, m2)
print(res)
