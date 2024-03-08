


###########################
# Objects to save results #
###########################

ForecastList = function(vars, nrow, models=c('model')) {
  
  # List of lists of matrices, first level is 'model', second level is the response variable.
  # Each matrix contains, for each observation of a given variable in a given model, the following values:
  #   - workshift of the observation
  #   - model (i.e., number of lags or persistence)
  #   - response variable name
  #   - real value
  #   - predicted value
  #   - prediction error (i.e., squared root of the correspondent diagonal element of the covariance matrix)
  #   - prediction interval lower and upper values computed as (pred +- 1.96 sq.err)
  #   - absolute error
  #   - squared error
  #
  # INPUT:
  #  > vars  : response variables
  #  > nrow  : number of predicted observations
  #  > models: vector with model names 
  #
  # OUTPUT: 
  #  > FR: list of lists of matrices
  
  FR = list()
  cnames = c('shift','model','var','value','pred',
             'pr.err','lower','upper','abs.err','sq.err')
  nc = length(cnames)
  for (m in models) { 
    FR[[m]] = list()
    for (v in vars) {
      FR[[m]][[v]] = matrix(NA, nrow=nrow, ncol=nc, 
                            dimnames=list(seq(1,nrow), cnames)) %>% 
        as_tibble()
    }
  }
  if (length(models)==1) FR = FR[[1]]
  return(FR)
}




##################
#   Clustering   #
##################

PerformKMClustering = function(data, threshold=0.7) {
  
  # Given a data set and a threshold, performs a K-Means clustering selecting the smallest
  # number of clusters such that the ratio between the within SS and the total SS surpasses
  # the threshold.
  #
  # INPUT:
  #  > data
  #  > threshold: a number between 0 and 1 
  #
  # OUTPUT:
  #  > km : K-Means clustering as returned by stats::kmeans 
  #
  # The data set must contain only the variables that will be used for performing the clustering. 
  
  n_clusters = 1
  gof        = 0
  while (gof < threshold) {
    n_clusters = n_clusters + 1
    km  = kmeans(centers=n_clusters, x=data, iter.max=1000, nstart=50)
    gof = km$betweenss/km$totss
  }
  message(sprintf('Extracted %i clusters. Ratio BSS/TSS = %.2f\n', n_clusters, gof*100))
  return(km)
}


# ++++++++++++++++++++++++++
# Compute cluster centroids 
# ++++++++++++++++++++++++++
ClusteringCentroids = function(data, classif.vars, add.vars, clustering) {
  
  # Given a data set, the classification variables and a kmeans model returns 
  # a data frame with the numerical description of the centroids.
  #
  # INPUT:
  #  > data
  #  > classif.vars  : vector of names or indexes of classification variables
  #  > add.vars      : vector of names or indexes of additional variables 
  #                    used to describe the centroids
  #  > clustering    : kmeans model
  #
  # OUTPUT:
  #  > CENTROIDS: a data frame whose rows are the clusters and columns are the average of 
  #               each variable, the number and proportion, the within SS and the average 
  #               silhouette within each cluster
  #
  
  if (!is.null(add.vars) & class(add.vars) == 'character') add.vars=match(add.vars, colnames(data))
  if (!is.null(classif.vars) & class(classif.vars) == 'character') classif.vars=match(classif.vars, colnames(data))
  N0 = dim(data)[1]
  distances    = dist(data[classif.vars])
  silhouettes  = cluster::silhouette(clustering$cluster,distances)
  summary.silh = summary(silhouettes)
  KM.count = data %>% 
    mutate(KM = clustering$cluster) %>%
    group_by(KM) %>%
    count()
  CENTROIDS = data %>% 
    mutate(KM = clustering$cluster) %>%
    group_by(KM) %>%
    select(all_of(c(classif.vars, add.vars))) %>%
    summarise_all(.funs=mean)
  CENTROIDS$n    = KM.count$n
  CENTROIDS$prop = round(100*CENTROIDS$n/N0,2)
  CENTROIDS$WSS  = clustering$withinss
  CENTROIDS      = arrange(CENTROIDS, KM)
  CENTROIDS$silh = round(summary.silh$clus.avg.widths[CENTROIDS$KM],2)
  rm(distances)
  return(CENTROIDS)
}


# ++++++++++++++++
# Update centroids (naive)
# ++++++++++++++++
update_centroids = function(centr.mat, new_data, nc, vars) {
  
  # Naive update of the centroids given a new_data that belongs to cluster nc
  #
  # INPUTS: 
  #  > centr.mat: data frame with centroids
  #  > new_data  
  #  > nc       : updated centroid
  #  > vars     : indexes of the variables to be updated
  #
  # OUTPUTS:
  #  > centr.mat: updated centroids
  
  old     = centr.mat[nc, vars]                                # centroid to be updated   
  old.n   = pull(centr.mat[nc,], n)                            # size to be updated
  old.wss = pull(centr.mat[nc,], WSS)                          # wss to be updated
  
  new.n   = old.n + 1                                          # new size
  new     = (old*old.n + new_data)/new.n                       # new centroid
  new.wss = old.wss + Rfast::dista(new_data, new, square=TRUE) # new wss
  
  new     -> centr.mat[nc ,vars]                               # .....
  new.n   -> centr.mat[nc,'n']                                 # .....
  new.wss -> centr.mat[nc,'WSS']                               # updating centroid matrix values
  centr.mat$prop = centr.mat$n/sum(centr.mat$n)*100            # ......
  
  return(centr.mat)
}
