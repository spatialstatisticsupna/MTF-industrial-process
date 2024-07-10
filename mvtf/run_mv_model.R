rm(list=ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# Libraries 

library(tidyverse)
library(MTS)        # for VARX models

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
# rcs : real cycle speed (= units/LT) 

# Times (minutes)
# OT      : opening time
# SBT     : stand by time (planned maintenance, starving, loading, handling, ...)
# LT      : loading time (LT = OT - SBT)
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
# TIME SCHEDULE:
# ---------------------------- OT ----------------------------
# ----- SBT -----|-------------------- LT --------------------
#                 --- DT ---|-------------- OpT --------------
#                            -- PLT --|--------- NOpT --------
#                                      - QLT -|------ VT -----


OEEbw = readRDS('../data/oee_data_4weeks.Rds')


######################
# Variable selection #
######################

rates = c('av','pf','oee')

yn = c('OpT','NOpT')                            # response variables
tn = c(rates, 'OT', 'rcs', 'TU')                # classification variables
xn = c('new.of','new.sh','ics')                 # co-variates
ini.st.var = 'new.sh'

r  = length(rates)
m1 = length(yn)
p  = length(tn)
m2 = length(xn)


#################################################################
#  Leave-One-Out estimation: for each week wk select weeks[-wk] #
#     as training set, weeks[wk] as test set. Run the model     #
#     and save aggregated results.                              #
#################################################################

weeks = OEEbw$week.id %>% unique()
LAGS  = 0:5


########################## 
# Object to save results #
##########################
ForecastResults = ForecastList(yn, nrow(OEEbw), c(paste0('q',LAGS),'persistence',paste0('VARX_',LAGS)))
for (mod in names(ForecastResults)) {
  for (v in yn) {
    ForecastResults[[mod]][[v]]$shift = interaction(OEEbw$wday,OEEbw$tday,sep='')
    ForecastResults[[mod]][[v]]$model = mod
    ForecastResults[[mod]][[v]]$var   = v
    ForecastResults[[mod]][[v]]$value = OEEbw %>% pull(v)
  }
}

Fi = 0

for (wk in 1:length(weeks)) {
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
  
  # Pull out data for VARX model
  ytr = as.matrix(OEEtr[,yn],rownames.force=T)
  xtr = as.matrix(OEEtr[,xn],rownames.force=T)
  xtr = cbind(xtr,dummy.code(OEEtr$wday))
  xtr = cbind(xtr,dummy.code(OEEtr$tday))
  # Set Monday Morning as baseline shift
  bs  = c('Mo','M')
  xtr = xtr[,setdiff(colnames(xtr),bs)]
  
  
  # Min-max normalization for further use in the knn assignment
  maxTrain = apply(ClassTrain[,tn],2,max) %>% as.numeric()                 
  minTrain = apply(ClassTrain[,tn],2,min) %>% as.numeric()  
  minmax = function(X, minx, maxx) (X-minx)/(maxx-minx)     
  minmaxTrain = t(apply(ClassTrain[,tn],1,minmax,minx=minTrain,maxx=maxTrain)) 
  
  Parameters = list()
  
  #################
  # Learning Step #
  #################
  for (ii in 1:length(LAGS)) {
    q = LAGS[ii]
    qname = paste0('q',q)
    Parameters[[qname]] = Initialization(
      W, K, m1, m2, C, p, lags=q, wnames,
      snames, xn, fnames, vnames=yn, pnames=tn)
    
    message(paste('  *** Learning parameters in multivariate model with',q,'lags ... ***'))
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
  
  for (ii in 1:length(LAGS)) {
    q = LAGS[ii]
    qname  = paste0('q',q)
    vxname = paste0('VARX_',q)
    Params = Parameters[[qname]]
    
    # VARX model:    > STATIC!!, meaning that coefficients do not update (Eqs. [3.3a - 3.3d]) with new information
    #                > Do not include the classes in the model (hidden states in Fig. [3.1], or Eq. [3.2b] in the paper)
    varx = MTS::VARX(zt=ytr,p=q,xt=xtr,m=0,output=FALSE)
    varx = MTS::refVARX(varx)
    
    # Warning!! We have to make the next trick to obtain the predictions in the VARX model with 0 lags 
    #           (that is, a regular regression)
    if (q==0) {
      varx$aror = 1
      varx$Phi=diag(c(0,0))
    }
    
    
    first = ifelse(q==0,N0,N0-q+1)
    ClassTest = rbind(ClassTrain[first:N0,],
                      data.frame('lsp'   = OEEte$lsp,
                                 'class' = NA,
                                 'wday'  = OEEte$wday,
                                 'tday'  = OEEte$tday,
                                 OEEte[,union(yn,tn)],
                                 OEEte[,union(xn,ini.st.var)]),
                      make.row.names = FALSE)
    
    # Pull out data for VARX prediction
    yte = as.matrix(OEEte[,yn],rownames.force=T)
    xte = as.matrix(OEEte[,xn],rownames.force=T)
    xte = cbind(xte,dummy.code(OEEte$wday))
    xte = cbind(xte,dummy.code(OEEte$tday))
    xte = xte[,colnames(xtr)]
    
    # Prediction
    Centroids = Centers
    cat('  *** Prediction step in model with',q,'lags ... ***\n')
    for (i in 1:L) {
      j = i + Fi
      
      # MVTF-Model prediction 
      U = predictors(ClassTest, i, yn, q, xn, ini.st.var, 'class')
      forecast = Prediction(Params, U, 'class', c('wday','tday'))
      
      # VARX model prediction
      newxt = xte[i,] %>% as.matrix() %>% t()
      varx.pred = MTS::VARXpred(m1=varx,newxt=newxt)

      # Add last data for the next prediction
      varx$data = rbind(varx$data,yte[i,])
      varx$xt   = rbind(varx$xt,xte[i,])
      
      
      w    = U[,'wday']
      s    = U[,'tday']
      prev = U[,'class']
      
      
      # compute absolute and squared errors
      for (jj in 1:m1) {
        v    = yn[jj]
        
        if (q == 0) {
          
          ForecastResults$persistence[[v]][j, 'pred']       = ClassTest[i,v]
          ForecastResults$persistence[[v]][j, 'abs.err']    = abs(ClassTest[i,v] - ClassTest[i+1, v])
          ForecastResults$persistence[[v]][j, 'sq.err']     = (ClassTest[i,v] - ClassTest[i+1, v])**2
          
        }
        
        # MVTF forecast
        yhat = forecast$prediction[jj]
        shat = sqrt(forecast$error[jj,jj])
        y    = ForecastResults[[qname]][[v]][j, 'value']
        lower = yhat-1.96*shat; upper = yhat+1.96*shat
        ForecastResults[[qname]][[v]][j, 'pred']       = yhat
        ForecastResults[[qname]][[v]][j, 'pr.err']     = shat
        ForecastResults[[qname]][[v]][j, 'lower']      = lower
        ForecastResults[[qname]][[v]][j, 'upper']      = upper
        ForecastResults[[qname]][[v]][j, 'abs.err']    = abs(yhat - y)
        ForecastResults[[qname]][[v]][j, 'sq.err']     = (yhat - y)**2
        ForecastResults[[qname]][[v]][j, 'coverage']   = as.integer(y>=lower & y<=upper)
        
        # VARX forecast
        yhat = varx.pred$pred[jj]
        shat = varx.pred$se[jj]
        lower = yhat-1.96*shat; upper = yhat+1.96*shat
        ForecastResults[[vxname]][[v]][j, 'pred']       = yhat
        ForecastResults[[vxname]][[v]][j, 'pr.err']     = shat
        ForecastResults[[vxname]][[v]][j, 'lower']      = lower
        ForecastResults[[vxname]][[v]][j, 'upper']      = upper
        ForecastResults[[vxname]][[v]][j, 'abs.err']    = abs(yhat - y)
        ForecastResults[[vxname]][[v]][j, 'sq.err']     = (yhat - y)**2
        ForecastResults[[vxname]][[v]][j, 'coverage']   = as.integer(y>=lower & y<=upper)
      }
      
      # update centroids 
      ind   = ifelse(q==0,i+1,i+q)
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
      Centroids = update_centroids(centr.mat=Centroids,
                                   new_data=z.new,
                                   nc=nearest_class %>% as.integer(),
                                   vars=1:nv+1)
      
      # update parameters
      first = ifelse(q==0,ind-1,ind-q)
      Params = UpdateModel(Params, ClassTest[first:ind,], ff=ff, responses=yn,
                           others=xn, ini.st.var=ini.st.var, clus.vars=tn,
                           class.var='class', group.vars=c('wday','tday'),
                           lags=q, reini=reini)
    }
    cat('  *** End of prediction step in model with',q,'lags ... ***\n')
  }
  Fi = Fi + L
}



# Gather prediction results in a supermatrix

models = names(ForecastResults)
mtf_mv_forecast.aux = vector('list', length=length(models))
names(mtf_mv_forecast.aux) = models
for (model in models) {
  mtf_mv_forecast.aux[[model]] = do.call(rbind,ForecastResults[[model]])
}

mtf_mv_forecast = do.call(rbind,mtf_mv_forecast.aux)
rm(mtf_mv_forecast.aux)



# Folder to save results

if(!file.exists("../results")) dir.create("../results")



# Save the results

saveRDS(mtf_mv_forecast, file="../results/MTF_mv_forecast.Rds")


