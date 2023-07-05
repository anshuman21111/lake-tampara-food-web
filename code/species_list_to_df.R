#####################################################################
### combine species data into a single dataframe
#####################################################################


library(dplyr)

# getting info about all excel sheets
fname = "data/Species list.xlsx"
sheets <- readxl::excel_sheets(fname)
tibble <- lapply(sheets, function(x) readxl::read_excel(fname, sheet = x))
data_frame <- lapply(tibble, as.data.frame)

# assigning names to data frames
names(data_frame) <- sheets

big_data_frame <- bind_rows(data_frame, .id = "guild")

write.csv2(big_data_frame, file = "data/species_data_frame.csv", quote = FALSE, row.names = FALSE)
