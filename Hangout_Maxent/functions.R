#evModel: evaluate maxent distribution model using kfold partitioning
#Arguments:
# fold: vector of fold indicator
# pres.covs: matrix with presence covariates
# bkg.covs: matrix with background covariates
# factor.ind: vector of column indices of categorical covariates (if any)
# mxnt.args: arguments for maxent run
# path: path to save results
#Returns
# A list with auc and maximum tss values for each fold

evModel <-function(fold, pres.covs, bkg.covs, factor.ind=NULL,mxnt.args,path){
  auc<-rep(NA,max(fold))
  max.tss<-rep(NA,max(fold))
  for (i in 1:max(fold)){
    occtest <- pres.covs[fold == i, ]
    occtrain <- pres.covs[fold != i, ]
    env.values<-data.frame(rbind(pres.covs, bkg.covs))

    if(!is.null(factor.ind)){
      env.values[, factor.ind]<-as.factor(env.values[, factor.ind])
    }
    me <- maxent(env.values, y, args=mxnt.args, path="./outputR")
    e<-evaluate(me, p=data.frame(occtest), a=data.frame(bkg.covs))
    auc[i]<-e@auc
    tss<-e@TPR+e@TNR-1
    max.tss[i]<-e@t[which.max(tss)]
  }
 return(list(auc=auc, max.tss=max.tss))
}

#getLambdaTable
#Convert the lambda object in maxent to a data frame
#Arguments:
# lambdas(character vector): character vector returned by maxent function. This is usually obtained
#  by accessing the lambda slot in a maxent object, e.g. mxnt.obj@lambdas
#Returns:
# A data frame of lambda values, with columns feature name, coefficient, min value, max value.

getLambdaTable<-function(lambdas){
  lambdas.list <- strsplit(lambdas,",")
  nparams = length(lambdas) - 4
  varnames=rep("NA",nparams)
  result<-data.frame(lambdas=rep(0,nparams))
  for (i in 1:nparams){
    varnames[i]<-lambdas.list[[i]][1]
    result[i,1]<-as.numeric(lambdas.list[[i]][2])
  }
  result<-data.frame(varnames,result,stringsAsFactors=F)
  return(result)
}

#evModel2: evaluate maxent distribution model using kfold partitioning
#Arguments:
# fold: vector of fold indicator
# pres.covs: matrix with presence covariates
# bkg.covs: matrix with background covariates
# factor.ind: vector of column indices of categorical covariates (if any)
# mxnt.args: arguments for maxent run
# path: path to save results
#Returns
# A list with number of parameters, auc and maximum tss values for each fold

#Función de evaluación
evModel2 <-function(fold, pres.covs, bkg.covs, factor.ind=NULL,mxnt.args,path){
  auc<-rep(NA,max(fold))
  max.tss<-rep(NA,max(fold))
  nparams<-rep(NA,max(fold))
  for (i in 1:max(fold)){
    occtest <- pres.covs[fold == i, ]
    occtrain <- pres.covs[fold != i, ]
    env.values<-data.frame(rbind(pres.covs, bkg.covs))
    
    if(!is.null(factor.ind)){
      env.values[, factor.ind]<-as.factor(env.values[, factor.ind])
    }
    me <- maxent(env.values, y, args=mxnt.args, path="./outputR")
    nparams[i]<-length(which(getLambdaTable(me@lambdas)$lambdas!=0))
    e<-evaluate(me, p=data.frame(occtest), a=data.frame(bkg.covs))
    auc[i]<-e@auc
    tss<-e@TPR+e@TNR-1
    max.tss[i]<-e@t[which.max(tss)]
    
  }
  return(list(nparams=nparams,auc=auc, max.tss=max.tss))
}