p <- 'C:/Users/yolan/Downloads/Landsat_Harmonized_Bands_1985_2025_train (3).csv'
df <- read.csv(p, stringsAsFactors=FALSE)
str(df)
df$date <- as.Date(df$date)
cat('date range:', as.character(min(df$date,na.rm=TRUE)), as.character(max(df$date,na.rm=TRUE)), '\n')
print(table(format(df$date,'%m')))
