# Export just the strat names for GBDB cross-reference
library(readr)
library(dplyr)

readr::locale(encoding = "UTF-8")

extracts <- read_csv(paste(getwd(), "data", "results.csv", sep="/"))

# rather than remove unresolved strat_name_id hits, try to compare them
# extracts <- extracts[!grepl(pattern = "\\~", extracts$strat_name_id),]

extracts <- extracts[,c('strat_phrase_root', 'strat_flag')]

distinct <- distinct(as_tibble(extracts))

write_csv(distinct, 'strat_names.csv')
