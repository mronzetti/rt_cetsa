# #  888888ba  d888888P           a88888b.  88888888b d888888P .d88888b   .d888888
# #  88    `8b    88             d8'   `88  88           88    88.    "' d8'    88
# #  a88aaaa8P'   88             88        a88aaaa       88    `Y88888b. 88aaaaa88a
# #  88   `8b.    88    88888888 88         88           88          `8b 88     88
# #  88     88    88             Y8.   .88  88           88    d8'   .8P 88     88
# #  dP     dP    dP              Y88888P'  88888888P    dP     Y88888P  88     88
# # oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
# #
# # Preparation of microplate file from MatLab readout
# # 
# # Michael Ronzetti
library(tidyverse)
library(readxl)
library(janitor)
grid_96w <- expand.grid(row = LETTERS[1:8], col = c(1:12)) %>%
  arrange(row) %>%
  mutate(address = paste(row, col, sep = '')) %>%
  select(-c('row', 'col'))
# Read in matlab file, drop well/col values, and prepare to assign wells
df <- read_excel(
  "matlab_530.xlsx",
  sheet = "Sheet1",
  col_names = FALSE,
  .name_repair = "unique"
) %>%
  select(-c('...1', '...2')) %>%
  rownames_to_column() %>% rename('well' = 'rowname')
# Construct r-friendly column names from temperature
tracker <- 1
for (val in 2:ncol(df) - 1) {
  names(df)[val + 1] <- paste('t_', val, sep = '')
  tracker <- tracker + 1
}
# Pivot around data to tidy
df <- df %>%
  pivot_longer(., cols = 2:ncol(df)) %>%
  pivot_wider(names_from = well) %>%
  rename(., 'Temperature' = 'name') %>%
  mutate(., Temperature = as.integer(gsub("[^0-9.]", "", Temperature))) %>%
  mutate(., Temperature = Temperature + 36)
# Assemble the data frames for moltenprot analysis with proper column identifiers.
q1 <- df %>%
  select(., 1, 2:97)
tracker <- 1
for (val in 1:nrow(grid_96w)) {
  colnames(q1)[val + 1] <- grid_96w$address[val]
  tracker <- tracker + 1
}
q2 <- df %>%
  select(., 1, 98:193)
tracker <- 1
for (val in 1:nrow(grid_96w)) {
  colnames(q2)[val + 1] <- grid_96w$address[val]
  tracker <- tracker + 1
}
q3 <- df %>%
  select(., 1, 194:289)
tracker <- 1
for (val in 1:nrow(grid_96w)) {
  colnames(q3)[val + 1] <- grid_96w$address[val]
  tracker <- tracker + 1
}
q4 <- df %>%
  select(., 1, 290:385)
tracker <- 1
for (val in 1:nrow(grid_96w)) {
  colnames(q4)[val + 1] <- grid_96w$address[val]
  tracker <- tracker + 1
}
write.csv(q1, 'cleaned_expt1.csv', row.names = FALSE)
write.csv(q2, 'cleaned_expt2.csv', row.names = FALSE)
write.csv(q3, 'cleaned_expt3.csv', row.names = FALSE)
write.csv(q4, 'cleaned_expt4.csv', row.names = FALSE)