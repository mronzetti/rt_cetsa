# #  888888ba  d888888P           a88888b.  88888888b d888888P .d88888b   .d888888
# #  88    `8b    88             d8'   `88  88           88    88.    "' d8'    88
# #  a88aaaa8P'   88             88        a88aaaa       88    `Y88888b. 88aaaaa88a
# #  88   `8b.    88    88888888 88         88           88          `8b 88     88
# #  88     88    88             Y8.   .88  88           88    d8'   .8P 88     88
# #  dP     dP    dP              Y88888P'  88888888P    dP     Y88888P  88     88
# # oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
# #
# # Process moltenprot file and prepare for readout and analysis
# #
# # Michael Ronzetti
library(tidyverse)
library(stringr)
library(readxl)
library(PharmacoGx)
# Data columns to create well assignments
letter <- LETTERS[1:12]
number <- c(1:24)
number <- str_pad(number, 2, pad = '0')
col_by_row <-
  expand.grid(row = sprintf('%.2d', 1:16), col = sprintf('%.2d', 1:24)) %>%
  arrange(., row)
# Gather deltaG and curve fit parameters for 384-well plate
#
# Read in 4 experimental files from MoltenProt readout
exp1_param <-
  read_excel('./data/cleaned_expt1/Signal_resources/Signal_results.xlsx',
             sheet = 'Fit parameters') %>%
  select(-c('Condition'))
exp2_param <-
  read_excel('./data/cleaned_expt2/Signal_resources/Signal_results.xlsx',
             sheet = 'Fit parameters') %>%
  select(-c('Condition'))
exp3_param <-
  read_excel('./data/cleaned_expt3/Signal_resources/Signal_results.xlsx',
             sheet = 'Fit parameters') %>%
  select(-c('Condition'))
exp4_param <-
  read_excel('./data/cleaned_expt4/Signal_resources/Signal_results.xlsx',
             sheet = 'Fit parameters') %>%
  select(-c('Condition'))
# Reformat ID column in each exp from MoltenProt format (A1, not A01) to arrange
exp1_param$ID <-
  gsub('([A-Z])(\\d)(?!\\d)', '\\10\\2\\3', exp1_param$ID, perl = TRUE)
exp1_param <- exp1_param %>% arrange(ID)
exp2_param$ID <-
  gsub('([A-Z])(\\d)(?!\\d)', '\\10\\2\\3', exp2_param$ID, perl = TRUE)
exp2_param <- exp2_param %>% arrange(ID)
exp3_param$ID <-
  gsub('([A-Z])(\\d)(?!\\d)', '\\10\\2\\3', exp3_param$ID, perl = TRUE)
exp3_param <- exp3_param %>% arrange(ID)
exp4_param$ID <-
  gsub('([A-Z])(\\d)(?!\\d)', '\\10\\2\\3', exp4_param$ID, perl = TRUE)
exp4_param <- exp4_param %>% arrange(ID)
# Combine all experiments and add identifiers
exp_param_full <-
  exp1_param %>% rbind(., exp2_param, exp3_param, exp4_param) %>%
  rownames_to_column() %>% rename('well' = 'rowname') %>%
  select(
    -c(
      'ID',
      'kN_init',
      'bN_init',
      'kU_init',
      'bU_init',
      'dHm_init',
      'Tm_init',
      'kN_fit',
      'bN_fit',
      'kU_fit',
      'bU_fit',
      'S',
      'dCp_component'
    )
  ) %>%
  bind_cols(col_by_row) %>%
  relocate(c('row', 'col'), .after = well)
# Construct well assignment
tracker <- 1
for (val in letter) {
  for (num in number) {
    exp_param_full$well[tracker] <- paste(val, num, sep = '')
    print(paste(val, num, sep = ''))
    tracker <- tracker + 1
  }
}
write.csv(exp_param_full, './data/exp_parameters.csv', row.names = FALSE)
# Gather baseline-corrected fit curves for 384-well plate and pivot plate.
#
# Read in 4 experimental files from MoltenProt readout
exp1_curve <-
  read_excel('./data/cleaned_expt1/Signal_resources/Signal_results.xlsx',
             sheet = 'Baseline-corrected')
exp2_curve <-
  read_excel('./data/cleaned_expt2/Signal_resources/Signal_results.xlsx',
             sheet = 'Baseline-corrected') %>%
  select(-c('Temperature'))
exp3_curve <-
  read_excel('./data/cleaned_expt3/Signal_resources/Signal_results.xlsx',
             sheet = 'Baseline-corrected') %>%
  select(-c('Temperature'))
exp4_curve <-
  read_excel('./data/cleaned_expt4/Signal_resources/Signal_results.xlsx',
             sheet = 'Baseline-corrected') %>%
  select(-c('Temperature'))
exp_curve_all <-
  cbind(
    xp1 = exp1_curve,
    xp2 = exp2_curve,
    xp3 = exp3_curve,
    xp4 = exp4_curve
  ) %>%
  rename(., Temperature = xp1.Temperature) %>%
  mutate(., Temperature = paste('val_t_', Temperature, sep = ''))
exp_curve_all <- exp_curve_all %>%
  pivot_longer(cols = 2:ncol(exp_curve_all)) %>%
  pivot_wider(names_from = Temperature) %>%
  rownames_to_column() %>% rename('well' = 'rowname') %>%
  bind_cols(col_by_row) %>%
  select(-'name')
tracker <- 1
for (val in letter) {
  for (num in number) {
    exp_curve_all$well[tracker] <- paste(val, num, sep = '')
    print(paste(val, num, sep = ''))
    tracker <- tracker + 1
  }
}
write.csv(exp_curve_all, './data/exp_curvefit.csv', row.names = FALSE)
# Attach platemap to assign ID, conc, and pivot
#
# Read in platemap ID sheet and replace empty with DMSO, pivot both conc and id long
id_df <- read_excel('./data/platemap.xlsx', sheet = 'sample') %>%
  select(-1) %>%
  pivot_longer(., cols = 1:ncol(.))
id_df$value <- gsub('empty', 'DMSO', id_df$value)
conc_df <- read_excel('./data/platemap.xlsx', sheet = 'conc') %>%
  select(-1) %>%
  pivot_longer(., cols = 1:ncol(.))
# Construct the full experimental details file and export to csv
exp_curve_all <- exp_curve_all %>%
  select(-c('well', 'row', 'col'))
full_df <- exp_param_full[, 1] %>%
  cbind(id_df$value, conc_df$value) %>%
  cbind(exp_param_full[, 2:ncol(exp_param_full)]) %>%
  cbind(exp_curve_all) %>%
  rename(compound_id = 'id_df$value', compound_conc = 'conc_df$value')
#
# Find AUC for each well and attach for full data file
# 
auc_df <- select(full_df,starts_with('val')) %>%
  rownames_to_column() %>%
  pivot_longer(c(2:ncol(.))) %>%
  pivot_wider(names_from='rowname')
auc_df$name <- gsub('val_t_','',auc_df$name)
auc_df <- auc_df %>%
  mutate(name=as.numeric(name)-273.15)
auc_vals <- tibble(vals = 1:384)
for(i in 2:ncol(auc_df)){
  temp_df <- auc_df[,1]
  temp_df[,2] <- auc_df[,(i)]
  auc_vals$vals[(i-1)] <- computeAUC(temp_df[,1],temp_df[,2],viability_as_pct = FALSE,verbose=FALSE)
}
full_df <- full_df %>% cbind(.,auc_vals) %>% rename(auc=vals)
write.csv(full_df, './data/fullexpdata.csv', row.names = FALSE)