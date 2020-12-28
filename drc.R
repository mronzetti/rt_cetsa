# #  888888ba  d888888P           a88888b.  88888888b d888888P .d88888b   .d888888
# #  88    `8b    88             d8'   `88  88           88    88.    "' d8'    88
# #  a88aaaa8P    88             88        a88aaaa       88    `Y88888b. 88aaaaa88
# #  88   `8b.    88    88888888 88         88           88          `8b 88     88
# #  88     88    88             Y8.   .88  88           88    d8'   .8P 88     88
# #  dP     dP    dP              Y88888P'  88888888P    dP     Y88888P  88     88
# # ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
# #
# # Dose-Response Curves w/ RT-CETSA technique
# #
# # Michael Ronzetti
library(tidyverse)
library(drc)
library(ggbeeswarm)
library(cowplot)
library(ggrepel)
dir.create('./data/dr_curves/')
#Read in experimental file
raw_df <- read.csv('./data/fullexpdata.csv') %>%
  dplyr::select(-c('well', 'row', 'col'))
model_df <-
  tibble(compound = (unique(
    filter(raw_df, compound_id != 'DMSO')$compound_id
  )))
for (i in 3:ncol(raw_df)) {
  col.nm <- colnames(raw_df)[i]
  model_df[, col.nm] <- NA
  print(col.nm)
}
modelfit_df <- tibble(colnames(model_df)[2:ncol(model_df)])
names(modelfit_df)[1] <- 'analysis'
modelfit_df <- modelfit_df %>%
  filter(analysis != 'BS_factor')
for (i in 1:nrow(model_df)) {
  temp_df <-
    filter(raw_df, raw_df$compound_id == model_df$compound[(i)]) %>%
    dplyr::select(!'BS_factor')
  print(paste('Analyzing: ', model_df$compound[i]), sep = '')
  temp_modelfit_df <- modelfit_df[1] %>%
    mutate(ic50 = 0, noEffect = 0)
  for (n in 3:ncol(temp_df)) {
    dr_df <- temp_df %>% dplyr::select(c(2, n))
    colnames(dr_df)[1] <- 'conc'
    colnames(dr_df)[2] <- 'resp'
    temp.model <-
      drm(
        resp ~ conc,
        data = dr_df,
        fct = LL.4(),
        control = drmc(
          errorm = FALSE,
          maxIt = 500,
          noMessage = TRUE
        )
      )
    pred.fit <-
      expand.grid(pr.x = exp(seq(log(min(
        dr_df[1]
      )), log(max(
        dr_df[1]
      )), length = 10000)))
    if ("convergence" %in% names(temp.model) == FALSE) {
      pm <-
        predict(object = temp.model,
                newdata = pred.fit,
                interval = 'confidence')
      pred.fit$p <- pm[, 1]
      pred.fit$pmin <- pm[, 2]
      pred.fit$pmax <- pm[, 3]
      dr_plot <- ggplot(dr_df, aes(x = conc, y = resp)) +
        geom_point() +
        geom_ribbon(data = pred.fit,
                    aes(
                      x = pr.x,
                      y = p,
                      ymin = pmin,
                      ymax = pmax
                    ),
                    alpha = 0.2) +
        geom_line(data = pred.fit, aes(x = pr.x, y = p)) +
        scale_x_log10() +
        theme_cowplot() +
        labs(
          title = paste(
            'Analysis of ',
            model_df$compound[i],
            ' by ',
            colnames(temp_df)[n],
            sep = ''
          ),
          subtitle = paste(
            'EC50: ',
            signif(temp.model$coefficients[4], 3),
            ' nM',
            '\n',
            'Significance of noEffect Test: ',
            signif(noEffect(temp.model)[3], 3),
            sep = ''
          ),
          x = 'Concentration'
        )
      print(dr_plot)
      png(
        filename = paste(
          './data/dr_curves/',
          model_df$compound[i],
          colnames(temp_df)[n],
          '.png',
          sep = ''
        ),
        width = 3200,
        height = 1800,
        res = 300
      )
      print(dr_plot)
      dev.off()
      print(n)
      temp_modelfit_df$ic50[(n-2)] <-
        signif(temp.model$coefficients[4], 3)
      temp_modelfit_df$noEffect[(n-2)] <-
        signif(noEffect(temp.model)[3], 3)
      print(paste('******  ',colnames(temp_df)[n],'  ******',sep=''))
      print(paste('IC50: ', signif(temp.model$coefficients[4], 3), sep =
                    ''))
      if (signif(noEffect(temp.model)[3], 3) < 0.01) {
        print(paste('noEffect: ', signif(noEffect(temp.model)[3], 3), ' (PASS)', sep = ''))
      } 
      else {
        print(paste('noEffect: ', signif(noEffect(temp.model)[3], 3), ' (FAIL)', sep = ''))
      }
    }
  }
  modelfit_df <- modelfit_df %>%
    cbind(., temp_modelfit_df[2:3])
  names(modelfit_df)[names(modelfit_df) == 'ic50'] <-
    paste('ic50_', model_df$compound[i], sep = '')
  names(modelfit_df)[names(modelfit_df) == 'noEffect'] <-
    paste('noEffect_', model_df$compound[i], sep = '')
}