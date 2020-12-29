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
library(viridis)
library(ggrepel)
library(ggpubr)
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
    mutate(ic50 = 0, noEffect = 0) %>%
    filter(.,analysis != 'BS_factor')
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
      )), length = 1000)))
    if ("convergence" %in% names(temp.model) == FALSE) {
      pm <-
        predict(object = temp.model,
                newdata = pred.fit,
                interval = 'confidence')
      pred.fit$p <- pm[, 1]
      pred.fit$pmin <- pm[, 2]
      pred.fit$pmax <- pm[, 3]
      dr_plot <- ggplot(dr_df, aes(x = conc, y = resp)) +
        geom_line(data = pred.fit, aes(x = pr.x, y = p),size=1.5,color='black') +
        geom_point(size=4,shape=21,fill='orange',color='black') +
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
write.csv(modelfit_df,'./data/modelfit_params.csv',row.names = FALSE)

# Assay wide plots
ic50_df <- modelfit_df %>%
  dplyr::select(starts_with('analysis') | starts_with('ic50'))
names(ic50_df) <- gsub(pattern='ic50_',replacement='',x=names(ic50_df))
ic50_df <- ic50_df %>%
  filter(str_detect(analysis,'val')) %>%
  pivot_longer(!analysis) %>%
  filter(value<10000) %>%
  mutate(value=log(value,10))
ic50_df$name <- gsub('[[:punct:]]','-',ic50_df$name)
ic50_distribution <- ggplot(ic50_df,aes(x=analysis,y=name,fill=value)) +
  geom_tile() +
  scale_fill_viridis(discrete=FALSE,direction = -1)
print(ic50_distribution)
png(
  filename = './models/ic50_dist.png',
  width = 3200,
  height = 1800,
  res = 300
)
print(ic50_distribution)
dev.off()

noEffect_df <- read.csv('./data/modelfit_params.csv') %>%
  dplyr::select(starts_with('analysis') | starts_with('noEffect'))
names(noEffect_df) <- gsub(pattern='noEffect_',replacement='',x=names(noEffect_df))
noEffect_df <- noEffect_df %>%
  pivot_longer(!analysis) %>%
  mutate(value=log(value,10))
noEffect_df$name <- gsub('[[:punct:]]','-',noEffect_df$name)
noEff_distribution <- ggplot(noEffect_df,aes(x=analysis,y=name,fill=value)) +
  geom_tile(size=1) +
  geom_vline(xintercept = 5.5) +
  theme_pubclean() +
  theme(
    axis.text.y = element_text(size=8, face='bold'),
    axis.text.x = element_text(size=8,angle=90,hjust=1,vjust=0.25),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    legend.position = 'right',
    legend.background = element_rect(size=1, linetype="solid", 
                                     colour ="black")
  ) +
  labs(
    title='NoEffect Log10(P) By Analysis Method',
    fill='Log10 P-val'
  ) +
  scale_fill_gradient2(
    midpoint=-2,
    high='#BBBBBB',
    mid='#BBBBBB',
    low='#EE3377'
  )
print(noEff_distribution)
png(
  filename = './models/noEff_dist.png',
  width = 3200,
  height = 1800,
  res = 300
)
print(noEff_distribution)
dev.off()
