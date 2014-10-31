source_url('https://github.com/PieceMaker/rscripts/blob/master/waterData.r', local = F)
lnWaterData <- log(waterData)
names(lnWaterData) <- c('lnQ',
                        'lnX1',
                        'lnX2',
                        'lnX3',
                        'lnX4',
                        'lnX5',
                        'lnX6',
                        'lnX7',
                        'lnX8',
                        'lnX9')

lnWaterDataModels <- function(seed = NULL) {
  set.seed(seed)

  library(leaps)

  lnq <- lnWaterData$lnQ
  lnx <- lnWaterData[,-1]

  rleaps <- regsubsets(lnx, lnq, int = T, nbest = 500, nvmax = 25, really.big = T, method = c('ex'))
  cleaps <- summary(rleaps, matrix=T)
  Models<-cleaps$which
  Models <- rbind(c(T, rep(F, dim(lnx)[2])), Models)
  colnames(Models)[1] <- 'Intercept'

  K <- 10
  sample.i <- sample(nrow(lnx), nrow(lnx))
  foldsize <- floor(nrow(lnx)/K)
  foldsize <- rep(foldsize, K)
  unalloc <- nrow(lnx)-K*unique(foldsize)
  if(unalloc > 0) {
    foldsize[1:unalloc] <- foldsize[1:unalloc]+1
  }

  prederrors <- matrix(0, 2^(ncol(lnWaterData)-1), K)
  sample.index <- 0
  X <- as.matrix(cbind(1, lnx))
  for(k in 1:K) {
    testset <- sample.i[(sample.index+1):(sample.index+foldsize[k])]
    trainset <- setdiff(sample.i, testset)
    sample.index <- sample.index+foldsize[k]
    for(m in 1:2^(ncol(X)-1)) {
      betahat <- solve(t(X[trainset,Models[m,]])%*%X[trainset,Models[m,]])%*%t(X[trainset,Models[m,]])%*%lnq[trainset]
      ypred <- X[testset,Models[m,]]%*%betahat
      prederrors[m,k] <- sum((lnq[testset]-ypred)^2)
    }
  }
  PE <- apply(prederrors, 1, sum)/nrow(lnWaterData)

  topmodels <- Models[sort.list(PE)[1:5],]
  #Prevent warning when using as.data.frame
  rownames(topmodels) <- c(1:5)
  return(cbind(as.data.frame(topmodels), PredError = sort(PE)[1:5]))
}

