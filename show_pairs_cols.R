library(data.table)
df <- fread('C:/Users/yolan/Downloads/Landsat_Harmonized_Bands_1985_2025_low (3).csv')
setnames(df, names(df), trimws(names(df)))
df[,date:=as.IDate(date)]
ls_old <- df[satellite=='LANDSAT_457']
ls_new <- df[satellite=='LANDSAT_89']
pairs <- ls_old[ls_new, on=.(location_id,date), nomatch=0L, allow.cartesian=TRUE]
FEATURES <- c('Blue','Green','Red','NIR','SWIR1','SWIR2','NDVI','EVI','PPI','MSAVI')
for(f in FEATURES){ if(f %in% names(pairs)) setnames(pairs,f,paste0(f,'_457')); if(paste0('i.',f) %in% names(pairs)) setnames(pairs,paste0('i.',f),paste0(f,'_89')); }
drop_cols<-grep('^i\\.',names(pairs),value=TRUE)
pairs[,(drop_cols):=NULL]
cat('columns in pairs after renaming:\n')
cat(paste(seq_along(names(pairs)), names(pairs), sep=': '), sep='\n')
