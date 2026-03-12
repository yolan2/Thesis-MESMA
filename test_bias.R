library(data.table)
library(ggplot2)

INPUT_CSV <- 'C:/Users/yolan/Downloads/Landsat_Harmonized_Bands_1985_2025_low (3).csv'
df <- fread(INPUT_CSV)
setnames(df, names(df), trimws(names(df)))
df[,date:=as.IDate(date)]
ls_old <- df[satellite=='LANDSAT_457']
ls_new <- df[satellite=='LANDSAT_89']
pairs<-ls_old[ls_new,on=.(location_id,date),nomatch=0L,allow.cartesian=TRUE]
FEATURES<-c('Blue','Green','Red','NIR','SWIR1','SWIR2','NDVI','EVI','PPI','MSAVI')
for(f in FEATURES){ if(f %in% names(pairs)) setnames(pairs,f,paste0(f,'_457')); if(paste0('i.',f) %in% names(pairs)) setnames(pairs,paste0('i.',f),paste0(f,'_89')); }
drop_cols<-grep('^i\\.',names(pairs),value=TRUE)
pairs[, (drop_cols):=NULL]
pairs_long <- melt(pairs,id.vars=c('location_id','date','year','lat','lon','vegetation'), measure.vars=patterns('_457$','_89$'), variable.name='band_idx',value.name=c('val_457','val_89'))
FEATURES<-c('Blue','Green','Red','NIR','SWIR1','SWIR2','NDVI','EVI','PPI','MSAVI')
pairs_long[, band := FEATURES[band_idx]]
pairs_long[, diff := val_457 - val_89]

for(f in c('SWIR1','SWIR2')){
  dif <- pairs_long[band==f,diff]
  cat(f,': n=',length(dif[!is.na(dif)]),' unique diffs=',length(unique(dif)),' range=',range(dif,na.rm=TRUE),'\n')
}

p_bias_dist <- ggplot(pairs_long, aes(x = band, y = diff, fill = band)) +
  geom_violin(trim = TRUE, alpha = 0.6, colour = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, colour = "grey10", fill = "white") +
  geom_hline(yintercept = 0, colour = "red", linewidth = 0.7) +
  scale_fill_brewer(palette = "Set2", guide = "none") +
  labs(title = "Bias by Feature (L457 - L89)", x = "Feature", y = "Difference") +
  theme_bw(base_size = 11)

ggsave('c:/tmp/bias_test.png', p_bias_dist, width=10, height=6, dpi=150)
cat('Saved plot to c:/tmp/bias_test.png\n')
