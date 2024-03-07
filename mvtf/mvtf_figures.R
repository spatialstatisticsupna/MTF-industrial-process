rm(list=ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# Load packages
library(tidyverse)
library(viridis)

source('theme_mtf.R')



#################################
# Load the multivariate results #
#################################

mv_forecast = readRDS('results/MTF_mv_forecast.Rds')
# load('results/MTF_mv_forecast.RData'); mv_forecast=mtf_mv_forecast; rm(mtf_mv_forecast)



########################################################################################################
# Compute average prediction error, MAE and RMSE for each model, each variable and for each work shift #
########################################################################################################
yn      = mv_forecast$var %>% unique()
shifts  = mv_forecast$shift %>% unique()
models  = mv_forecast$model %>% unique()
metrics = c('pr.err','mae','rmse')

MVTF.metrics = tibble(model=rep(models,length(shifts)),
                      shift=rep(shifts,each=length(models)))
MVTF.metrics[,do.call('paste',c(expand.grid(metrics,yn),sep='_'))] = NA

for (r in 1:nrow(MVTF.metrics)) {
  mod = MVTF.metrics$model[r]
  shf = MVTF.metrics$shift[r]
  for (v in 1:length(yn)) {
    X = mv_forecast %>%
      filter(model==mod & shift==shf & var==yn[[v]])
    pr   = X %>% pull(pr.err) %>% mean() %>% round(.,3)
    mae  = X %>% pull(abs.err) %>% mean() %>% round(.,3)
    rmse = X %>% pull(sq.err) %>% mean() %>% sqrt() %>% round(.,3)
    MVTF.metrics[r,3*(v-1) + 3] = pr
    MVTF.metrics[r,3*(v-1) + 4] = mae
    MVTF.metrics[r,3*(v-1) + 5] = rmse
  }
}



################################################################################################################
# Figure 5.1: Boxplot of average prediction error, MAE and RMSE of each model and each variable by work shifts #
################################################################################################################

bp.metr_by_sh = MVTF.metrics %>% 
  reshape2::melt(id.vars=1:2, measure.vars=c(4,5,7,8,10,11)) %>% 
  ggplot(aes(x=model, y=value)) +
  geom_boxplot(aes(fill=model)) + 
  coord_flip() + 
  scale_fill_viridis(discrete=TRUE, option="H") +
  labs(y=NULL, fill=NULL, x=NULL) +
  guides(fill='none') + 
  scale_x_discrete(labels=c("q0"="no lags", "q1"="1 lag", "q2"="2 lags",
                            "q3"="3 lags", "q4"="4 lags", "q5"="5 lags",
                            'persistence'='persistence')) +
  facet_wrap(vars(variable), nrow=3, scales='free_x') +
  theme_mtf 


# Folder to save results
if(!file.exists("figures")) dir.create("figures")

# Save the results
ggsave(filename='figures/Figure5_1.pdf', plot=bp.metr_by_sh)




#################
# Other figures #
#################

# Evolution of metrics across models for each workshift
MVTF.metrics %>% 
  reshape2::melt(id.vars=1:2,measure.vars=c(4,5,7,8,10,11)) %>% 
  ggplot(aes(x=model,y=value)) +
  geom_point(aes(color=variable)) + geom_line(aes(group=variable,color=variable))  + 
  scale_color_viridis(discrete=TRUE, option="H") +
  labs(y=NULL,color=NULL)  + 
  facet_wrap(vars(shift),nrow=5,scale='free_y') + 
  theme(legend.position='bottom',
        legend.direction = 'horizontal',
        legend.background = element_rect(fill = NA, colour = NA),
        axis.text.x=element_text(family='sans',face='italic',size=7),
        axis.title.x=element_text(family='sans',face='bold',size=8),
        plot.title=element_text(family='sans',face='bold',size=10,hjust=0.5),
        axis.title.y=element_text(family='sans',face='bold',size=8),
        axis.text.y=element_text(family='sans',size=7), 
        axis.ticks.x=element_blank(),
        legend.text=element_text(family='sans',size=8),
        legend.title=element_text(family='sans',face='bold',size=8),
        legend.key.size=unit(0.4,'cm'),
        strip.text=element_text(family='sans',face='bold',size=8))



# Empirical Cumulative Distribution of absolute and squared for each model
# (no distinction between shifts)
mv_forecast %>% 
  reshape2::melt(id.vars='model', measure.vars=c('abs.err','sq.err')) %>% 
  ggplot(aes(x=value, color=model)) +
  stat_ecdf(geom="step", pad=FALSE, linewidth=1) +
  geom_hline(yintercept=0.95, linetype='dashed') + 
  labs(x=NULL, y='Empirical Cumulative Distribution Function', color=NULL) +
  scale_colour_brewer(palette = "Set1",
                      labels=c("q0" = "no lags", "q1" = "1 lag", "q2" = "2 lags",
                               "q3" = "3 lags", "q4" = "4 lags", "q5" = "5 lags",
                               'persistence'='persistence')) +
  theme_mtf + 
  facet_wrap(vars(variable),nrow=1,scales='free_x')


# Empirical Cumulative Distribution of absolute error for each model
# (distinction between shifts)
mv_forecast %>% 
  ggplot(aes(x=abs.err, color=model)) +
  stat_ecdf(geom="step", pad=FALSE, linewidth=.85) +
  geom_hline(yintercept=0.95, linetype='dashed') + 
  labs(x=NULL, y='Empirical Cumulative Distribution Function MAE', color=NULL) +
  scale_colour_brewer(palette = "Set1",
                      labels=c("q0" = "no lags", "q1" = "1 lag", "q2" = "2 lags",
                               "q3" = "3 lags", "q4" = "4 lags", "q5" = "5 lags",
                               'persistence'='persistence')) +
  theme_mtf + 
  facet_wrap(vars(shift),nrow=3,scales='free_x') 


# Empirical Cumulative Distribution of squared error for each model
# (distinction between shifts)
mv_forecast %>% 
  ggplot(aes(x=sq.err, color=model)) +
  stat_ecdf(geom="step", pad=FALSE) +
  geom_hline(yintercept=0.95, linetype='dashed') + 
  labs(x=NULL, y='Empirical Cumulative Distribution Function RMSE', color=NULL) +
  scale_colour_brewer(palette = "Set1",
                      labels=c("q0" = "no lags", "q1" = "1 lag", "q2" = "2 lags",
                               "q3" = "3 lags", "q4" = "4 lags", "q5" = "5 lags",
                               'persistence'='persistence')) +
  theme_mtf + 
  facet_wrap(vars(shift),nrow=3,scales='free_x') 




#####################################################################
# Figure 5.2: Comparison between multivariate and univariate models #
#####################################################################

###############################
# Load the univariate results #
###############################

uv_forecast = readRDS('results/MTF_uv_forecast.Rds')
# load('results/MTF_uv_forecast.RData'); uv_forecast=mtf_uv_forecast; rm(mtf_uv_forecast)


# gather multivariate and univariate results
mv_forecast = mv_forecast %>% add_column(vers='multiv', .after=2)
uv_forecast = uv_forecast %>% add_column(vers='univ', .after=2)

all_forecast = rbind(mv_forecast, uv_forecast)



########################################################################################################
# Compute average prediction error, MAE and RMSE for each model, each variable and for each work shift #
########################################################################################################
# yn      = uv_forecast$var %>% unique()
# shifts  = uv_forecast$shift %>% unique()
# models  = uv_forecast$model %>% unique()
# metrics = c('pr.err','mae','rmse')

UVTF.metrics = tibble(model=rep(models,length(shifts)),
                      shift=rep(shifts,each=length(models)))
UVTF.metrics[,do.call('paste',c(expand.grid(metrics,yn),sep='_'))] = NA

for (r in 1:nrow(UVTF.metrics)) {
  mod = UVTF.metrics$model[r]
  shf = UVTF.metrics$shift[r]
  for (v in 1:length(yn)) {
    X = uv_forecast %>%
      filter(model==mod & shift==shf & var==yn[[v]])
    pr   = X %>% pull(pr.err) %>% mean() %>% round(.,3)
    mae  = X %>% pull(abs.err) %>% mean() %>% round(.,3)
    rmse = X %>% pull(sq.err) %>% mean() %>% sqrt() %>% round(.,3)
    UVTF.metrics[r,3*(v-1) + 3] = pr
    UVTF.metrics[r,3*(v-1) + 4] = mae
    UVTF.metrics[r,3*(v-1) + 5] = rmse
  }
}


# gather multivariate and univariate metrics
MVTF.metrics = MVTF.metrics %>% add_column(vers='multiv', .after=1)
UVTF.metrics = UVTF.metrics %>% add_column(vers='univ', .after=1)

all.metrics = rbind(MVTF.metrics, UVTF.metrics) %>% filter(model!='persistence')



bp.metr_mv.uv = all.metrics %>% 
  reshape2::melt(id.vars=1:2, measure.vars=4:12) %>%
  ggplot(aes(x=model,y=value)) + 
  geom_boxplot(aes(linetype=vers,fill=model)) +
  scale_fill_viridis(discrete=TRUE, option="H") +
  scale_linetype_manual(values=c('solid','11')) +
  labs(y=NULL,x=NULL,fill=NULL,linetype=NULL) + 
  guides(fill='none') + 
  scale_x_discrete(labels=c("q0" = "no lags", "q1" = "1 lag", "q2" = "2 lags",
                            "q3" = "3 lags", "q4" = "4 lags", "q5" = "5 lags")) + 
  coord_flip() + 
  facet_wrap(vars(variable), nrow=3, scales='free_x') + 
  theme_mtf +
  theme(legend.position='bottom',
        legend.direction = 'horizontal',
        legend.box.spacing = unit(-0.2,'cm'))

# Folder to save results
if(!file.exists("figures")) dir.create("figures")

# Save the results
ggsave(filename='figures/Figure5_2.pdf', plot=bp.metr_mv.uv)


