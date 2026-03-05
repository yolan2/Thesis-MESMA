source('c:/Users/yolan/OneDrive/Documenten/UGENT/Master/R_MESMA/january_averages.R')
print('Testing compute_global_index_snr with normal input')
df <- data.frame(location_id=c('A','A','B'), pheno_year=c(2000,2000,2001), PPI=c(0.1,0.2,0.3))
print(compute_global_index_snr(df, c('PPI')))

print('Testing compute_global_index_snr with weird pheno_years')
df2 <- data.frame(location_id=c('A','A'), pheno_year=I(list(environment(),environment())), PPI=c(0.1,0.2))
print(compute_global_index_snr(df2, c('PPI')))
