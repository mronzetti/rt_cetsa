# #  888888ba  d888888P           a88888b.  88888888b d888888P .d88888b   .d888888
# #  88    `8b    88             d8'   `88  88           88    88.    "' d8'    88
# #  a88aaaa8P    88             88        a88aaaa       88    `Y88888b. 88aaaaa88
# #  88   `8b.    88    88888888 88         88           88          `8b 88     88
# #  88     88    88             Y8.   .88  88           88    d8'   .8P 88     88
# #  dP     dP    dP              Y88888P'  88888888P    dP     Y88888P  88     88
# # ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
# #
# Analysis Comparison Script
# Comparing the parameters for: Tagg, dGu, AUC, Signal difference along the curve.
# -Michael Ronzetti
library(tidyverse)
library(readxl)
library(janitor)
library(ggbeeswarm)
library(rlang)
library(cowplot)
library(ggrepel)
#
# Read in full experimental data and prepare tibble for statistical tests
raw_df <- read.csv('./exp_fulldata.csv')
analysis_df <-
  tibble(
    measurement = colnames(raw_df[6:length(colnames(raw_df))]),
    avg_dmso = 0,
    sd_dmso = 0,
    avg_530 = 0,
    sd_530 = 0,
    zprime = 0
  )
#Perform descriptive statistics and Zprime calculation for plate
for (i in 6:ncol(raw_df)) {
  print(colnames(raw_df)[i])
  avg_dmso <- mean(raw_df[, i][raw_df$compound_id == 'DMSO'])
  sd_dmso <- sd(raw_df[, i][raw_df$compound_id == 'DMSO'])
  avg_530 <- mean(raw_df[, i][raw_df$compound_id == 'NCGC00372530'])
  sd_530 <- sd(raw_df[, i][raw_df$compound_id == 'NCGC00372530'])
  analysis_df[(i - 5), 2] <- avg_dmso
  analysis_df[(i - 5), 3] <- sd_dmso
  analysis_df[(i - 5), 4] <- avg_530
  analysis_df[(i - 5), 5] <- sd_530
  analysis_df[(i - 5), 6] <-
    1 - ((3 * sd_530 + 3 * sd_dmso) / abs(avg_530 - avg_dmso))
}
#Add bool for measurements that are above Z' of 0.5
analysis_df <- analysis_df %>%
  mutate(zprime_check = ifelse(zprime >= 0.5, TRUE, FALSE))
write.csv(analysis_df,'./analysis.csv')
#
# Plot of the Z' values against each analysis type
#
analysis_plot <-
  ggplot(filter(analysis_df, zprime > 0),
         aes(x = measurement, y = zprime, fill = zprime_check)) +
  geom_point(size = 3,
             shape = 21,
             color = 'black') +
  geom_hline(
    yintercept = 0.5,
    linetype = 'longdash',
    size = 1.5,
    alpha = 0.7
  ) +
  ylim(0, 1) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1,
      face = 'bold'
    ),
    legend.background = element_blank(),
    legend.box.background = element_rect(colour = "black")
  ) +
  labs(title = 'RT-CETSA Z\' Calculations', fill = 'Acceptable Z\'') +
  scale_fill_manual(values = c('#DDAA33', '#004488'))
print(analysis_plot)
png(
  filename = './models/analysis_plot.png',
  width = 3200,
  height = 1800,
  res = 300
)
print(analysis_plot)
dev.off()
#
# Loop through analysis tests and create ggplots of the point distributions and statistics
#
loop_df <- raw_df %>%
  select(-c('well', 'col', 'row', 'compound_conc')) %>%
  arrange(., compound_id)
plot_df <-
  tibble(measurement = colnames(raw_df[6:length(colnames(raw_df))]))
for (i in 2:(ncol(loop_df))) {
  plt <-
    ggplot(loop_df, aes(x = loop_df[, 1], y = loop_df[, (i)], fill = loop_df[, 1])) +
    geom_beeswarm(size = 2, shape = 21) +
    geom_boxplot(alpha = 0.35, size = 1.25) +
    scale_fill_manual(values = c('#DDAA33', '#004488')) +
    labs(
      title = paste('RT-CETSA Data Analysis by ', colnames(loop_df)[(i)], sep =
                      ''),
      subtitle = paste('Z\': ', signif(analysis_df$zprime[(i) - 1], 3), sep =
                         ''),
      y = colnames(loop_df)[(i)]
    ) +
    theme_cowplot() +
    theme(legend.position = 'None',
          axis.title.x = element_blank())
  print(plt)
  png(
    filename = paste('./models/', colnames(loop_df)[(i)], '.png', sep = ''),
    width = 3200,
    height = 1800,
    res = 300
  )
  print(plt)
  dev.off()
  print(paste(i, colnames(loop_df)[(i)], sep = ' '))
}
#
# Pivot around data sheet for response-temp curves
#
dr_df <- raw_df %>%
  select(c(2, 11:ncol(.))) %>%
  pivot_longer(., cols = starts_with('val'), names_to = 'temp')
dr_df$temp <- as.numeric(gsub('val_t_', '', dr_df$temp))
dr_df <- dr_df %>%
  mutate(temp = temp - 273.15)
dr_plot <-
  ggplot(dr_df, aes(
    x = temp,
    y = value,
    group = temp,
    color = compound_id
  )) +
  geom_jitter(
    data = filter(dr_df, compound_id == 'DMSO'),
    aes(x = temp, y =
          value),
    size = 0.1,
    alpha = 0.25
  ) +
  geom_jitter(
    data = filter(dr_df, compound_id == 'NCGC00372530'),
    aes(x = temp, y =
          value),
    size = 0.1,
    alpha = 0.25
  ) +
  geom_boxplot(
    data = filter(dr_df, compound_id == 'DMSO'),
    aes(x = temp, y =
          value),
    outlier.shape = NA,
    alpha = 0.8
  ) +
  geom_boxplot(
    data = filter(dr_df, compound_id != 'DMSO'),
    aes(x = temp, y =
          value),
    outlier.shape = NA,
    alpha = 0.8
  ) +
  scale_color_manual(values = c('#DDAA33', '#004488')) +
  theme_cowplot() +
  theme(
    legend.background = element_blank(),
    legend.box.background = element_rect(colour = "black")
  ) +
  labs(title = 'RT-CETSA Temperature-Response Aggregate Curves',
       x = 'Temperature (C)',
       y = 'Normalized Curve') +
  scale_x_continuous(breaks = seq(35, 95, 10))
print(dr_plot)
png(
  filename = './models/aggregate_curves.png',
  width = 3200,
  height = 1800,
  res = 300
)
print(dr_plot)
dev.off()
