rm(list=ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# Libraries ---------------------------------------------------------------

library(tidyverse)

source('functions.R')
source('update_model.R')

#################
# Read the data #
#################

# Rates (%)
# lo  : loading
# av  : availability
# pf  : performance
# qu  : quality
# oee : overall equipment effectiveness (av x pf x qu)
# -------------------------------------

# ics : ideal cycle speed (units/min)
# rcs : real cycle speed (= units/LT) (new: 29/12/2023)

# Times (minutes)
# OT      : opening time
# NWT     : non working time (lack of demand, stock out, strikes, security drillings, ...)
# SBT     : stand by time (planned maintenance, starving, loading, handling, ...)
# LT      : loading time (LT = OT - NWT - SBT)
# lo      : LT/OT
# -------------------------------------
# DT      : down time (failures, reactive maintenance, set-up, recalibration, ...)
# OpT     : operating time (OpT = LT - DT)
# av      : OpT/LT
# -------------------------------------
# PLT     : performance losses time (micro-stoppages, start up, shut down, heat up, reduced speed, ...)
# NOpt    : net operating time (NOpT = OpT - PLT)
# pf      : NOpT/OpT 
# -------------------------------------
# QLT     : quality rate time (defects, scraps, reworks, ...)
# VT      : valuable time (VT = NOpT - QRT)
# qu      : VT/NOpT
# 
# 
# ---------------------------- OT ----------------------------
# - NWT -|- SBT -|-------------------- LT --------------------
#                 --- DT ---|-------------- OpT --------------
#                            -- PLT --|--------- NOpT --------
#                                      - QLT -|------ VT -----


OEEbw = readRDS('../data/oee_data_4weeks.Rds')


######################
# Variable selection #
######################

rates = c('lo','av','pf','qu','oee')

responses = c('OpT','NOpT','VT')                # response variables
tn = c(rates, 'OT', 'rcs', 'TU')                # classification variables
xn = c('new.of','new.sh','ics')                 # co-variates
ini.st.var = 'new.sh'

r  = length(rates)
p  = length(tn)
m2 = length(xn)


##################################################################
#  Leave-One-Out estimation: for each week wk, select weeks[-wk] #
#     as training set, weeks[wk] as test set. Run the model      #
#     and save aggregated results.                               #
##################################################################

weeks = OEEbw$week.id %>% unique()
LAGS  = 0:5



########################## 
# Object to save results #
##########################
ForecastResults = ForecastList(responses, nrow(OEEbw), c(paste0('q',LAGS),'persistence'))
for (mod in names(ForecastResults)) {
  for (v in responses) {
    ForecastResults[[mod]][[v]]$shift = interaction(OEEbw$wday,OEEbw$tday,sep='')
    ForecastResults[[mod]][[v]]$model = mod
    ForecastResults[[mod]][[v]]$var   = v
    ForecastResults[[mod]][[v]]$value = OEEbw %>% pull(v)
  }
}



############################
# Start loop for responses #
############################
for (h in 1:length(responses)) {
  yn = responses[h]
  m1 = 1
  
  Fi = 0
  
  
  ########################
  # Start loop for weeks #
  ########################
  for (wk in 1:length(weeks)) {
    
    # train and test datasets
    OEEtr = OEEbw %>% filter(week.id %in% weeks[-wk])
    OEEte = OEEbw %>% filter(week.id == weeks[wk])
    N0    = nrow(OEEtr)
    
    message('\n**********************************')
    message(paste('****   Leave-One-Out method    ***\n**** leaving out week',weeks[wk],'***'))
    message('**********************************\n')
    
    
    set.seed(1234)
    
    # K-means clustering in the training set using tn (goodness-of-fit threshold: 90%)
    KM.Clustering = PerformKMClustering(OEEtr[,tn], threshold=0.9)
    C = KM.Clustering$size %>% length()                        
    
    # Centroids 
    Centers = ClusteringCentroids(data         = OEEtr,
                                  classif.vars = tn,
                                  add.vars     = yn,
                                  clustering   = KM.Clustering)
    
    # Save state labels in the training set
    OEEtr$class = factor(KM.Clustering$cluster)
    
    
    
    # Some labels and parameters
    wnames = OEEtr$wday %>% unique()
    snames = OEEtr$tday %>% unique()
    fnames = paste('class', 1:C)
    W      = length(wnames)
    K      = length(snames)
    reini  = FALSE
    ffx    = 0.99; ffr = 0.95; ff = c(ffx,ffr) 
    
    
    
    # Pull out data for learning step 
    ClassTrain = data.frame('lsp'   = OEEtr$lsp,
                            'class' = OEEtr$class,
                            'wday'  = OEEtr$wday,
                            'tday'  = OEEtr$tday,
                            OEEtr[,union(yn,tn)],
                            OEEtr[,union(xn,ini.st.var)]
    )
    
    # Min-max normalization for further use in the knn assignment
    maxTrain = apply(ClassTrain[,tn],2,max) %>% as.numeric()                 
    minTrain = apply(ClassTrain[,tn],2,min) %>% as.numeric()  
    minmax = function(X, minx, maxx) (X-minx)/(maxx-minx)     
    minmaxTrain = t(apply(ClassTrain[,tn],1,minmax,minx=minTrain,maxx=maxTrain)) 
    
    Parameters = list()
    
    
    
    #################
    # Learning Step #
    #################
    for (ii in 1:length(LAGS)) {   # loop for models
      q = LAGS[ii]
      qname = paste0('q',q)
      Parameters[[qname]] = Initialization(
        W, K, m1, m2, C, p, lags=q, wnames,
        snames, xn, fnames, vnames=yn, pnames=tn)
      
      message(paste('  *** Learning parameters in model with',q,'lags, variable',yn,'... ***'))
      Parameters[[qname]] = UpdateModel(Parameters[[qname]], ClassTrain, ff=ff, responses=yn,
                                        others=xn, ini.st.var=ini.st.var, clus.vars=tn,
                                        class.var='class', group.vars=c('wday','tday'), 
                                        lags=q, reini=reini)
    }
    
    ####################
    # Forecasting Step #
    ####################
    L = nrow(OEEte)
    
    # (method for cluster assignment: knn / mahalanobis)
    cluster.assign = 'knn'
    
    for (ii in 1:length(LAGS)) {   # loop for models
      q = LAGS[ii]
      qname = paste0('q',q)
      Params = Parameters[[qname]]
      
      first = ifelse(q==0,N0,N0-q+1)
      ClassTest = rbind(ClassTrain[first:N0,],
                        data.frame('lsp'   = OEEte$lsp,
                                   'class' = NA,
                                   'wday'  = OEEte$wday,
                                   'tday'  = OEEte$tday,
                                   OEEte[,union(yn,tn)],
                                   OEEte[,union(xn,ini.st.var)]),
                        make.row.names = FALSE)
      
      
      # Prediction 
      Centroids = Centers
      cat('\n  *** Prediction step in model with',q,'lags, variable',yn,'... ***')
      for (i in 1:L) {    # loop for prediction-update
        j = i + Fi
        
        U = predictors(ClassTest, i, yn, q, xn, ini.st.var, 'class')
        forecast = Prediction(Params, U, 'class', c('wday','tday'))
        
        # compute absolute and root squared errors of last prediction
        w    = U[,'wday']
        s    = U[,'tday']
        prev = U[,'class']
        
        if (q == 0) {  # compute only once for persistence model
          
          ForecastResults$persistence[[yn]][j, 'pred']       = ClassTest[i,yn]
          ForecastResults$persistence[[yn]][j, 'abs.err']    = abs(ClassTest[i,yn] - ClassTest[i+1,yn])
          ForecastResults$persistence[[yn]][j, 'sq.err']     = (ClassTest[i,yn] - ClassTest[i+1,yn])**2
          
        }
        
        ind = ifelse(q==0,i+1,i+q)
        
        yhat = forecast$prediction[1]
        shat = sqrt(forecast$error[1,1])
        y    = ClassTest[ind, yn]
        ForecastResults[[qname]][[yn]][j, 'pred']       = yhat
        ForecastResults[[qname]][[yn]][j, 'pr.err']     = shat
        ForecastResults[[qname]][[yn]][j, 'lower']      = yhat-1.96*shat
        ForecastResults[[qname]][[yn]][j, 'upper']      = yhat+1.96*shat
        ForecastResults[[qname]][[yn]][j, 'abs.err']    = abs(yhat - y)
        ForecastResults[[qname]][[yn]][j, 'sq.err']     = (yhat - y)**2
      
        
        # last observation assignment and centroids updating 
        #if (i%%100==0) cat('*** Updating centroids and model parameters ... ***\n')
        sizes = Centroids$n
        
        if (cluster.assign == 'knn') {
          z.new = OEEte[i,tn]
          nearest_class = class::knn(train = minmaxTrain,
                                     test = minmax(z.new, minTrain, maxTrain),
                                     cl = ClassTrain$class, k=sqrt(nrow(ClassTrain)))
        }
        if (cluster.assign == 'mahala') {
          P = ResponsibilityBasedOnMahalDist(ClassTest[ind,],Params,tn)
          nearest_class = which.max(P)
        }
        
        ClassTest[ind,'class'] = nearest_class
        z.new = OEEte[i,union(tn, yn)] 
        nv = ncol(Centroids)-5
        Centroids = update_centroids(centr.mat=Centroids,new_data=z.new,nc=nearest_class %>% as.integer(),vars=1:nv+1)
        
        # update parameters
        first = ifelse(q==0,ind-1,ind-q)
        Params = UpdateModel(Params, ClassTest[first:ind,], ff=ff, responses=yn,
                             others=xn, ini.st.var=ini.st.var, clus.vars=tn,
                             class.var='class', group.vars=c('wday','tday'),
                             lags=q, reini=reini)
      }
      cat('\n  *** End of prediction step in model with',q,'lags, variable',yn,'... ***')
    }
    Fi = Fi + L
  }
}



# gather prediction results in a supermatrix
models = names(ForecastResults)
mtf_uv_forecast.aux = vector('list', length=length(models))
names(mtf_uv_forecast.aux) = models
for (model in models) {
  mtf_uv_forecast.aux[[model]] = do.call(rbind,ForecastResults[[model]])
}

mtf_uv_forecast = do.call(rbind,mtf_uv_forecast.aux)
rm(mtf_uv_forecast.aux)



# Folder to save results
if(!file.exists("../results")) dir.create("../results")

# Save the results
saveRDS(mtf_uv_forecast, file="../results/MTF_uv_forecast.Rds")
#save(mtf_uv_forecast, file='results/MTF_uv_forecast.RData')


