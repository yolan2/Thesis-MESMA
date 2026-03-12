library(readxl); library(dplyr); library(lubridate)
gw_xlsx <- "C:/Users/yolan/OneDrive/Documenten/UGENT/Master/masterproef/GIS/groundwater_depths (1).xlsx"
gw_raw <- read_excel(gw_xlsx, sheet = "Yingsu")
cat("=== GW COLUMNS ===\n"); cat(paste(names(gw_raw), collapse = " | "), "\n")
cat("=== FIRST 4 ROWS ===\n"); print(head(gw_raw, 4))
cat("=== GW YEARLY ===\n")
names_low <- tolower(names(gw_raw))
date_col  <- names(gw_raw)[grepl("date|day|time", names_low)][1]
depth_col <- names(gw_raw)[grepl("depth|gwd|groundwater|waterlevel|wl|gw|level", names_low)][1]
cat("date_col:", date_col, "\ndepth_col:", depth_col, "\n")
cat("depth sample:", paste(head(gw_raw[[depth_col]], 5), collapse=", "), "\n")
