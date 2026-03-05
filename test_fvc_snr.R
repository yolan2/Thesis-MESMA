source('c:/Users/yolan/OneDrive/Documenten/UGENT/Master/R_MESMA/january_averages.R')
print('Testing compute_fvc_snr structure')
df_signal <- data.frame(location_id=c('A','A','B'), idx1=c(1,2,3))
df_noise <- df_signal
res <- compute_fvc_snr(df_signal, df_noise, c('idx1'))
print(res)
str(res)
