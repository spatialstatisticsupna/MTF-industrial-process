library(tidyverse)
library(psych)

InitPriorDirichlet = function(C,K,W,alpha,wnames=NULL,snames=NULL) {

  # Sets an uninformative (symmetric) Dirichlet as a priori distribution

  # INPUT
  #  > C     : number of classes
  #  > K     : number of states
  #  > W     : number of work shifts
  #  > alpha : initial common pseudo-count (scalar)
  #  > wnames: working days names
  #  > snames: working shifts names

  # OUTPUT
  #  > Pi: parameters for between-classes transition probabilities distribution
  #       -  Pi$pi: initial state probabilities
  #       -  Pi$A : transition probabilities

  if (is.null(wnames)) wnames = paste('shift',1:W)
  if (is.null(snames)) snames = paste('state',1:K)
  cnames = paste('class',1:C)

  Pi = list()
  Pi$pi = array(alpha,dim=c(W,K,1,C), dimnames=list('workshift'= wnames,
                                                    'state' = snames,
                                                    NULL,
                                                    'to' = cnames))
  Pi$A  = array(alpha,dim=c(W,K,C,C), dimnames=list('workshift'= wnames,
                                                 'state' = snames,
                                                 'from' = cnames,
                                                 'to' = cnames))
  return(Pi)
}


Initialization = function(W,K,m,R1,R2,p,lags=1,wnames=NULL,snames=NULL,onames=NULL,fnames=NULL,vnames=NULL,pnames=NULL){

  # INPUT
  #  > W     : number of working days
  #  > K     : number of working shifts
  #  > m     : number of responses (continuous variables to predict)
  #  > R1,R2 : number of predictors, for the AR model and for the class model
  #  > p     : number of classification variables used in the clustering step
  #  > lags  : number of responses lags
  #  > wnames: working days names
  #  > snames: working shifts names
  #  > onames: predictor names (AR model)
  #  > fnames: predictor names (class model)
  #  > pnames: classification variables names

  # OUTPUT
  #  > Parameters: model parameters initialization

  if (is.null(wnames)) wnames=paste('shift',1:W)
  if (is.null(snames)) snames=paste('state',1:K)
  if (is.null(onames)) fnames=paste('other',1:R1)
  if (is.null(fnames)) fnames=paste('feat',1:R2)
  if (is.null(vnames)) vnames=paste('var',1:m)
  if (is.null(pnames)) pnames=paste('clus.var',1:p)
  
  vnames2 = c('intercept',paste0(vnames,'_',rep(1:lags, each=m)),onames)
  if (lags==0) vnames2 = c('intercept',onames)

  Pi     = InitPriorDirichlet(R2,K,W,alpha=0.5,wnames=wnames,snames=snames)
  Psi    = list()
  Omega  = list()

  Psi$Hx = array(0, dim=c(W,K,lags*m+R1+1,m),
                 dimnames=list('wday' = wnames,
                               'shift'= snames,
                               'var'  = vnames2,
                               'on'   = vnames
                               )
                 )
  Psi$Vx = array(0, dim=c(W,K,m,m),
                 dimnames=list('wday'  = wnames,
                               'shift' = snames,
                               'var.i' = vnames,
                               'var.j' = vnames
                               )
                 )

  Psi$Hr = array(0, dim=c(W,K,R2,m),
                 dimnames=list('wday'   = wnames,
                               'shift'  = snames,
                               'feature'= fnames,
                               'var'    = vnames
                               )
                 )
  Psi$Vr = array(0, dim=c(W,K,m,m),
                 dimnames=list('wday'  = wnames,
                               'shift' = snames,
                               'var.i' = vnames,
                               'var.j' = vnames
                               )
                 )
  Psi$Hw = array(0, dim=c(1,p,R2),
                 dimnames=list(NULL,
                               'clsf.var' = pnames,
                               'feature'  = fnames)
                 )
  Psi$Vw = array(0, dim=c(p,p,R2),
                 dimnames=list('clsf.var' = pnames,
                               'clsf.var' = pnames,
                               'feature'  = fnames
                 )
  )
  Omega$Px = array(0, dim=c(W,K,lags*m+R1+1,lags*m+R1+1),
                   dimnames=list('wday'  = wnames,
                                 'shift' = snames,
                                 'var.i' = vnames2,
                                 'var.j' = vnames2
                                 )
                   )
  Omega$gx = array(1, dim=c(W,K),
                        dimnames=list('wday' = wnames,
                                      'shift'= snames)
                   )
  Omega$Pr = array(0, dim=c(W,K,R2,R2),
                   dimnames=list('wday'  = wnames,
                                 'shift'  = snames,
                                 'feat.i' = fnames,
                                 'feat.j' = fnames
                                 )
                   )
  Omega$gr = array(1, dim=c(W,K),
                   dimnames=list('wday' = wnames,
                                 'shift'= snames)
  )
  Omega$Pw = rep(1e23, R2); names(Omega$Pw) = fnames
  Omega$gw = rep(1, R2); names(Omega$gw) = fnames

  for (w in 1:W) {
    for (k in 1:K) {
      Omega$Px[w,k,,] = diag(1,lags*m+R1+1)
      Omega$Pr[w,k,,] = diag(1,R2)
    }
  }


  return(list('Pi'=Pi, 'Psi'=Psi, 'Omega'=Omega))
}


DirExpVal = function(alpha) {

  # INPUT
  #  > alpha: Dirichlet distribution parameters
  # 
  # OUTPUT
  #  > Expected value of the Dirichlet distribution

  return(alpha/sum(alpha))
}

ClassSelection = function(probs) {

  # Random class selection based on a probability vector
  #
  # INPUT
  #  > probs: probability vector
  #
  # OUTPUT
  #  > i : selected class
  
  ap = cumsum(probs)
  u  = runif(1)
  i  = which(u <= ap)[1]
  return(i)
}

ResponsibilityBasedOnMahalDist = function(X,ParamSet,clsf.vars) {
  
  # Returns the membership probability of a set of observations as the normalized inverse of
  # the Mahalanobis distance
  #
  # INPUT:
  #  > X        : new observations
  #  > ParamSet : set of parameters
  #  > clsf.vars: classification variables
  #  
  # OUTPUT:
  #  > P : matrix of probabilities, P(i,j) = P[Xi is in class j]
  
  N      = nrow(X)
  C      = length(ParamSet$Omega$gw)
  cnames = names(ParamSet$Omega$gw)
  
  mu    = ParamSet$Psi$Hw
  sigma = ParamSet$Psi$Vw
  
  D = matrix(0, N, C, dimnames=list(paste0('obs.',1:N),cnames))
  for (c in 1:C) {
    m = mu[,,c][clsf.vars]
    S = corpcor::pseudoinverse(sigma[,,c][clsf.vars,clsf.vars]) 
    D[,c] = mahalanobis(X[,clsf.vars], m, S, inverted=T)
  }
  
  P = 1/D
  P = P/rowSums(P)
  return(P)
}

DistribPrediction = function(params,w,s,prec) {

  # Predictive Dirichlet distribution when the previous class is c(n-1)=prec and 
  # current covariates are zn=(w,s).
  
  # INPUT
  #  > params: set of Dirichlet parameters
  #  > w     : working day of the actual sample
  #  > s     : shift of the actual sample
  #  > prec  : precedent class of the actual sample 

  # OUPUT
  #  > Expected value of multinomial probabilities for the next jump

  alpha = params[w,s,prec,]
  return(DirExpVal(alpha))
}


###############################
#  Equations [3.3a] - [3.3d]  #
###############################

UpdateCovMatrix = function(V, H, P, g, y, u, ff) {
  
  # Equation [3.3c]
  # INPUT:
  #  > V : prior covarianve matrix
  #  > H : prior coefficient matrix
  #  > P : prior state matrix
  #  > g : prior discount factor
  #  > y : response
  #  > u : covariates
  #  > ff: forgetting factor (scalar)
  #
  # OUTPUT:
  #  > Updated V
  
  if (!is.null(dim(P))) {
    HTu = t(H) %*% u
    d   = as.numeric(ff + (t(u) %*% P %*% u))
    e   = y - HTu
    
    V = V - ( V - ( ff * e %*% t(e) / d ) ) / g
  }
  
  if (is.null(dim(P))) {
    HTu = H*u
    d = as.numeric(ff + u*P*u)
    e = as.matrix(y-HTu)
    V = V - (V - (ff*e%*%t(e)/d))/g
  }
  
  return(V)
}


UpdateCoeffMatrix = function(H, P, y, u, ff) {
  
  # Equation [3.3b]
  # INPUT:
  #  > H : prior coefficient matrix
  #  > P : prior state matrix
  #  > y : response
  #  > u : covariates
  #  > ff: forgetting factor (scalar)
  #
  # OUTPUT:
  #  > Updated H

  if (!is.null(dim(P))) {
    HTu = t(H) %*% u
    Pu  = P %*% u
    d   = as.numeric(ff + (t(u) %*% Pu))
    e   = y - HTu
    
    H = H + ( Pu %*% t(e) / d )
  }
  
  if (is.null(dim(P))) {
    Hu = H*u
    Pu = P*u
    d = as.numeric(ff + u*Pu)
    e = y-Hu
    H = H + (Pu*e/d)
  }
  
  return(H)
}

UpdateStateMatrix = function(P, u, ff, reini=TRUE, max.kappa=1000) {
  
  # Equation [3.3d]
  # INPUT:
  #  > P         : prior state matrix
  #  > u         : covariates
  #  > ff        : forgetting factor (scalar)
  #  > reini     : reinitiate to identity in case of bad conditioning
  #  > max.kappa : maximum condition number allowed before reinitiate
  #
  # OUTPUT:
  #  > Updated P

  if (!is.null(dim(P))) {
    Pu  = P %*% u
    d   = as.numeric(ff + (t(u) %*% Pu))
    
    P = ( P - Pu %*% t(Pu) / d ) / ff
    
    # reinitialize in case of bad conditioning
    if (reini) {
      if (kappa(P, exact=TRUE) > max.kappa) { P = diag(length(diag(P))) }
    }
  }
  
  if (is.null(dim(P))) {
    P = 1/(u*u + ff/P)
  }

  return(P)
}

UpdateWeight = function(g, ff) { 
  
  # Equation [3.3a]
  # INPUT:
  #  > g : prior discount factor
  #  > ff: forgetting factor (scalar)
  #
  # OUTPUT:
  #  > Updated g
  
  return(1 + g*ff) 
}




UpdateModel = function(ModelParams,data,ff,responses,others,ini.st.var,clus.vars,class.var,group.vars,lags=1,reini=FALSE,max.kappa=100,predictors='probs') {

  # Update model parameters 
  
  # INPUT
  #  > ModelParams: set of Distribution (coeff H and covariance V matrices),
  #                 State Parameters (state matrix P and weights g) and Dirichlet
  #                 Distribution parameters (Pi)
  #  > data       : complete data
  #  > ff         : forgetting factor (vector of 2)
  #  > responses  
  #  > others     : covariates names
  #  > ini.st.var : which covariate marks the beginning of the sequence
  #  > clus.vars  : classification variables
  #  > class.var  : which variable contains the class label
  #  > group.vars : weekday and shift
  #  > lags       : number of responses lags
  #  > reini      : allows reinitialization of state matrix in case of bad conditioning
  #  > max.kappa  : maximum condition number allowed before reinitiate
  #  > predictors : which predictors are used in the class model, probabilities ('probs') or dummies ('dummies')

  # OUTPUT
  #  > Returns updated sets of Dirichlet, distribution and state parameters.

  # ** Possible numerical instability of state matrix P is addressed by its
  #    re-initialization in case of high condition number.

  # number of responses and observations
  m <- length(responses)
  N <- nrow(data)

  # forgetting factors
  ffx <- ff[1]
  ffr <- ff[2]

  # split parameters
  Pi    <- ModelParams$Pi
  Psi   <- ModelParams$Psi
  Omega <- ModelParams$Omega

  first = ifelse(lags==0,2,1+lags)
  for (n in first:N) { 
    y  <- data[1:n,]
    yn <- y[n, responses] %>% as.numeric()
    tn <- y[n, clus.vars] %>% as.numeric()

    # a priori info and current class
    w    <- (y %>% pull(group.vars[1]))[n] %>% as.character()
    #                                                     ^
    #                                 /\_________________/                     
    # HE TENIDO QUE AÑADIR %>% as.character() PORQUE SINO EN ALGUNOS CASOS NO LEE BIEN Pi[w,s,c,] NO ME DIGAS POR QUÉ!!!!!
    
    s      <- (y %>% pull(group.vars[2]))[n]
    prec   <- (y %>% pull(class.var))[n-1]
    c      <- (y %>% pull(class.var))[n]
    ini.st <- (y %>% pull(ini.st.var))[n]

    # covariate vectors
    ux <- c(1)
    if (lags > 0) {
      for (j in 1:lags) ux = append(ux, y[n-j,responses] %>% as.numeric())
    }
    ux = append(ux,y[n,others]) %>% as.numeric()
    if (ini.st) probs = DirExpVal(Pi$pi[w,s,1,]) %>% as.numeric()
    if (!ini.st) probs = DirExpVal(Pi$A[w,s,prec,]) %>% as.numeric()
    
    C = length(probs)
    i = ClassSelection(probs)
    dummies = diag(C)[i,]
    ur <- get(predictors)
    
    # update r parameters
    Psi$Vr[w,s,,]   <- UpdateCovMatrix(Psi$Vr[w,s,,],Psi$Hr[w,s,,],Omega$Pr[w,s,,],Omega$gr[w,s],yn,ur,ffr)
    Psi$Hr[w,s,,]   <- UpdateCoeffMatrix(Psi$Hr[w,s,,],Omega$Pr[w,s,,],yn,ur,ffr)
    Omega$Pr[w,s,,] <- UpdateStateMatrix(Omega$Pr[w,s,,],ur,ffr)
    Omega$gr[w,s]   <- UpdateWeight(Omega$gr[w,s],ffr)

    # update x parameters
    Psi$Vx[w,s,,]   <- UpdateCovMatrix(Psi$Vx[w,s,,],Psi$Hx[w,s,,],Omega$Px[w,s,,],Omega$gx[w,s],yn,ux,ffx)
    Psi$Hx[w,s,,]   <- UpdateCoeffMatrix(Psi$Hx[w,s,,],Omega$Px[w,s,,],yn,ux,ffx)
    Omega$Px[w,s,,] <- UpdateStateMatrix(Omega$Px[w,s,,],ux,ffx,reini=reini,max.kappa=max.kappa)
    Omega$gx[w,s]   <- UpdateWeight(Omega$gx[w,s],ffx)
    
    # update centroid parameters
    Psi$Vw[,,c] <- UpdateCovMatrix(Psi$Vw[,,c],Psi$Hw[,,c],Omega$Pw[c],Omega$gw[c],tn,1,1)
    Psi$Hw[,,c] <- UpdateCoeffMatrix(Psi$Hw[,,c],Omega$Pw[c],tn,1,1)
    Omega$Pw[c] <- UpdateStateMatrix(Omega$Pw[c],1,1,reini=reini,max.kappa=max.kappa)
    Omega$gw[c] <- UpdateWeight(Omega$gw[c],1)

    # update dirichlet parameters
    if (ini.st) Pi$pi[w,s,1,c] = Pi$pi[w,s,1,c] + 1
    if (!ini.st) Pi$A[w,s,prec,c] = Pi$A[w,s,prec,c] + 1
  }

  return(list('Pi'=Pi, 'Psi'=Psi, 'Omega'=Omega))
}



predictors = function(test.data,row,responses,lags,others,ini.st.var,class.var) {
  
  # Collects all of the covariates, including lagged responses, for prediction
  #
  # INPUT:
  #  > test.data : test set
  #  > row       : predicted row
  #  > responses
  #  > lags
  #  > others    : xn
  #  > ini.st.var: which covariate marks the beginning of the sequence
  #  > class.var : which variable contains the class label
  #
  # OUTPUT: 
  #  > U: previous class, current workday and shift, lagged responses, xn
  
  m = length(responses)
  
  rrow          = ifelse(lags==0,row+1,row+lags)
  U             = test.data[rrow,] %>% select(!where(is.numeric))
  U[,class.var] = test.data[rrow-1,class.var]
  
  if (lags>0) {
    for (j in 1:lags) U = cbind(U,test.data[rrow-j,] %>% select(all_of(responses)))
    colnames(U)[1:(lags*m) + 4] = paste0(colnames(U)[1:(lags*m) + 4],'_',rep(1:lags, each=m))
  }
  U = cbind(U,test.data[rrow,] %>% select(all_of(union(others,ini.st.var))))
  attributes(U)$is.init.st.var = ini.st.var
  attributes(U)$predictors     = others
  return(U)
}



Prediction = function(Params,U,class.var,group.vars,preds='probs') {
  
  # Computes prediction [3.4a] and accuracy [3.4b] using weights [3.5]

  # INPUT
  #  > Params    : set of parameters
  #  > U         : full covariates vector
  #  > class.var : which variable contains the class label
  #  > group.vars: weekday and shift
  #  > preds     : which predictors are used in the class model, probabilities ('probs') or dummies ('dummies')
  #
  # OUTPUT
  # > Prediction of relevant variables and errors 

  c = U %>% pull(class.var)
  w = U %>% pull(group.vars[1]) %>% as.character()
  s = U %>% pull(group.vars[2])
  ini.st.var = attr(U,'is.init.st.var')
  ini.st     = U %>% pull(ini.st.var)
  wn         = attr(U,'predictors')

  if (!ini.st) probs = DistribPrediction(Params$Pi$A,w,s,c) %>% as.numeric()
  if (ini.st) probs = DistribPrediction(Params$Pi$pi,w,s,1) %>% as.numeric()
  C = length(probs)
  cl = ClassSelection(probs)
  dummies = diag(C)[cl,]
  ur <- get(preds)
  
  ux = U[5:length(U)]
  if (length(intersect(ini.st.var,wn))==0) ux = ux %>% select(-ini.st.var)
  ux = c(1,ux) %>% as.numeric()

  Hr = Params$Psi$Hr[w,s,,]
  Hx = Params$Psi$Hx[w,s,,]
  Vr = Params$Psi$Vr[w,s,,]
  Vx = Params$Psi$Vx[w,s,,]

  m  = ifelse(is.null(dim(Vr)), 1,dim(Vr)[1])
  
  if (m==1) {
    W = Vr/(Vr+Vx)
    I = 1
    pred.x = sum(ux*Hx)
  }
  
  if (m>1) {
    W = diag(diag(Vr/(Vr+Vx)))
    I = diag(m)
    pred.x = t(Hx) %*% ux
  }

  pred.r = t(Hr) %*% ur
  
  pred  = W %*% pred.x + (I-W) %*% pred.r
  error = W %*% Vx %*% t(W) + (I-W) %*% Vr %*% t((I-W))

  dimnames(pred) = dimnames(error) = NULL

  return(list('pred.1'=pred.x, 'pred.2'=pred.r, 'prediction'=pred, 'error'=error))
}
