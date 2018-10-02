options(stringsAsFactors = FALSE)

pcks <- c("tidyverse", "openair")
easypackages::libraries(pcks)

# names for columns ------------------------------------------------------------
var_names <- readLines("ETRHEQ_Model.m")[89:102] %>%
  str_split("-") %>%
  unlist() %>%
  magrittr::extract(c(FALSE, TRUE)) %>%
  str_split(",") %>%
  unlist

short_vnames <- var_names %>%
  magrittr::extract(c(TRUE, FALSE)) %>% 
  str_trim()

long_names <- var_names %>%
  magrittr::extract(c(FALSE, TRUE)) %>% 
  str_trim()



# import sample data------------------------------------------------------------
sample_data <- R.matlab::readMat("sample_data.mat")
str(sample_data)

# site name
sample_data$sample.data.site.name[[1]]
# year
YEAR <- sample_data$sample.data.year[[1]]
# data array
## one day
sample_data$sample.data[, 1, ]
dim(sample_data$sample.data)
dimnames(sample_data$sample.data) <- list(NULL, NULL, short_vnames)
head(sample_data$sample.data[, 1, ])
# sample_data$sample.data[half-hours, days, variables]

# save sample data array in RDS, columns (dim = 2) with names set
saveRDS(
  object = sample_data$sample.data,
  file = "data/sample_data_array.RDS"
  )



# sample data in a data frame -------------------------------------------------
sample_data_df <- sample_data$sample.data %>%
  reshape2::melt(.) %>%
  spread(Var3, value) %>%
  as_tibble() %>% 
  arrange(Var2)
#View(sample_data_df)
names(sample_data_df)[1:2] <- c("hhour", "doy")


#  <- plyr::ldply(
#   1:ncol(sample_data$sample.data), 
#   function(iday) {
#     sample_data$sample.data[, iday, ]
#   }
# )

sample_data_df


names(sample_data_df)[-c(1:2)] <- short_vnames
head(sample_data_df)
dim(sample_data_df)

# date column
sample_data_df <- 
mutate(sample_data_df, 
       dhh = hhour/2 * lubridate::dhours(1),
       date = lubridate::as_date(paste0(YEAR, "-", doy), format = "%Y-%j", tz = "UTC"),
       date = date + dhh,
       dhh = NULL,
       hhour = NULL,
       doy = NULL
       ) %>%
  select(date, everything())

sample_data_df
long_names
sample_data_df %>%
timePlot(., names(.)[-1], scales = list(y = list(relation = "free")))

# save sample data array in RDS, columns (dim = 2) with names set
saveRDS(
  object = sample_data_df,
  file = "data/sample_data_df.RDS"
)
