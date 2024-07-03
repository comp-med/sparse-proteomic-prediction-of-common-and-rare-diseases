###################################################
#### Script Name: coxnet_optim
#### Author: Julia Carrasco Zanini Sanchez
#### Description: Function to run coxnet optimization and testing  
#### Date : 17/07/2022

coxnet.optim.r <- function(Train.data, pred.vec, Train.surv.data, times, Test.data, Test.surv.data, boot.samples){
  ## Train.data : data set containing predictor variables used for model optimization 
  ## pred.vec : vector containing names of predictor variable 
  ## Train.surv.data : survival object from optimization set
  ## times : number of boostrap samples to compute confidence intervals around the C-index
  ## Test.data : data set containing predictor variables for the test set
  ## Test.surv.data : survival object from test set
  ## boot.samples : list of bootstrap samples 
  
  if(sum(!complete.cases(Train.data[,pred.vec]))!=0)
    stop("Missing training values for predictors")
  if(sum(!complete.cases(Test.data[,pred.vec]))!=0)
    stop("Missing testing values for predictors")
  if(sum(is.na(Train.surv.data))!=0)
    stop("Missing training survival values")
  if(sum(is.na(Test.surv.data))!=0)
    stop("Missing testing survival values")
  
  require(caret)
  require(glmnet)
  require(ROSE)
  
  ## optimize
  set.seed(354)
  las.morb <- cv.glmnet(as.matrix(Train.data[,pred.vec]),
                        as.matrix(Train.surv.data),
                        family="cox",
                        alpha=0,
                        lambda = 10^-seq(10,.25,-.25),
                        nfolds = 5,
                        sampling = "rose")
  
  ## parameters for the best model
  coxnet.opt <- las.morb$glmnet.fit
  s.opt <- las.morb$lambda.min
  opt.coef <- coef(coxnet.opt, s = s.opt)
  
  ## Testing - C-index by boostraping
  boot.cindex <- NULL
  jj <- boot.samples
  tmp.pred <- list()
  tmp.index <- list()
  for (i in 1:times) {
    # jj[[i]]  <- sample(nrow(Test.data), round(nrow(Test.data)*1), replace = T)
    tmp.pred[[i]] <- as.numeric(predict(coxnet.opt, 
                                        type = "response", 
                                        newx = as.matrix(Test.data[jj[[i]], pred.vec]),
                                        s = s.opt))
    tmp.index[[i]] <- glmnet::Cindex(tmp.pred[[i]], as.matrix(Test.surv.data[jj[[i]]]))
    boot.cindex <- c(boot.cindex, tmp.index[[i]])
  }
  
  ## linear predictor for the entire test set 
  pred.lp <- as.numeric(predict(coxnet.opt, 
                        type = "link", 
                        newx = as.matrix(Test.data[, pred.vec]),
                        s = s.opt))
  ## relative risk estimates for the test set 
  pred.risk <- as.numeric(predict(coxnet.opt, 
                          type = "response", 
                          newx = as.matrix(Test.data[, pred.vec]),
                          s = s.opt))
  ## aggregate all results
  coxnet.res <- list(
    glmnet.opt = coxnet.opt,
    lambda.opt = s.opt,
    opt.coefficients = opt.coef,
    Cindex.vec = boot.cindex,
    Cindex.mean = mean(boot.cindex,na.rm=T),
    ci.low = quantile(boot.cindex, 0.025,na.rm=T),
    ci.upp = quantile(boot.cindex, 0.975,na.rm=T),
    linear.predictor = pred.lp,
    relative.risk = pred.risk
  )
  return(coxnet.res)
  }
