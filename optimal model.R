rm(list = ls())
library(survival)
library(caret)
library(plotrix)
library(survminer)
library(glmnet)
library(gbm)
library(randomForestSRC)
###############################################################################
## data pre-processing
###############################################################################

dat = read.csv("D://study//Semester 3//sarcoma.csv", stringsAsFactors = TRUE)
dat = subset(dat, select = -c(X,dfstim,dfsind) )
dat = na.omit(dat)
rownames(dat) <- 1:nrow(dat)
dat$sex =( as.numeric(as.factor(dat$sex)))
dat$grade = (as.numeric(as.factor(dat$grade)))
dat$tumor.subtype =( as.numeric(as.factor(dat$tumor.subtype)))

data_subset = subset(dat,select=c(suv.mean, grade, max.grad_HIST,H0,
                                  sum.var_GLCM,Reg.grad.min,energy_GLCM,
                                  surind,surtim))
repeats = 10
folds = 5

prod = folds*repeats
threshold = (60*prod)/100
original_cols = ncol(dat)-2
subset_cols = ncol(data_subset)-2
subset_rows = nrow(data_subset)
subset_nms = colnames(subset(data_subset, select = -c(surtim, surind)))

####################################################################
# VAriable declaration
####################################################################
stp_pred_matrx = matrix(0,nrow=prod,ncol=subset_rows)
full_step_model = vector("list",folds)
# stp_train_conc = vector("list",folds)
# stp_test_conc = vector("list",folds)

final_stp_train_conc  = vector("list",prod)
final_stp_test_conc = vector("list",prod)
final_stp_models = vector("list",prod)

lasso_train_conc = lasso_test_conc = numeric(prod)
lasso_pred_matrix  = matrix(0, nrow=prod, ncol=nrow(data_subset))
lasso.feat.selected = matrix(0, nrow=prod, ncol=subset_cols)
colnames(lasso.feat.selected) = subset_nms

ridge_train_conc = ridge_test_conc =  numeric(prod)
ridge_pred_matrix  = matrix(0, nrow=prod, ncol=nrow(data_subset))

gbm_summary = vector("list",prod) 
gbm.train.conc = gbm.test.conc = vector("double",prod) 
gbm_pred_matrix = matrix(0,nrow=prod,ncol=subset_rows)

rsf_model = rsf_pred = vector("list",prod)
rsf_conc_train = rsf_bs = vector("list",prod)
rsf_conc_test = vector("list",prod)
rsf_conc_test = vector("list",prod)

#####################################################################
## KM curve
####################################################################
split.kms <- function(zos,os,oevn,NQ=100,zqmin=.05,zqmax=.95){
  qmin = quantile(zos,zqmin,na.rm=TRUE)
  qmax = quantile(zos,zqmax,na.rm=TRUE)
  p0q = numeric(NQ)
  med0s = seq(qmin,qmax,length=NQ)
  izs = matrix(NA,ncol=NQ,nrow=length(zos))
  for(iq in 1:NQ){
    IZ = as.integer(zos<med0s[iq])
    p0q[iq] = 1-pchisq(survdiff(Surv(os,oevn)~factor(IZ))$chisq,1)
    izs[,iq] = IZ
  }
  best.med0 = med0s[which.min(p0q)]
  IZ = as.integer(zos<best.med0)
  return(list(iz=IZ,izs=izs,best.t=best.med0,
              ps=p0q,best.p=p0q[which.min(p0q)]))
}

#######################################################################
# model fitting
#######################################################################
counter = 1
set.seed(10)
for(a in 1:repeats){
  shuffled_dat = data_subset[sample(1:nrow(data_subset),nrow(data_subset)), ]
  # print(dim(shuffled_dat))
  cvIndex = createFolds(factor(shuffled_dat$surind), folds,returnTrain = T)
  for(b in 1:length(cvIndex)){
    train_data = data_subset[cvIndex[[b]],]
    test_data = data_subset[-cvIndex[[b]],]
    
    stp_cox_fit = coxph(Surv(surtim,surind)~.,data=train_data)
    stp_cox_pred = predict(stp_cox_fit,newdata = test_data,type="lp")
    
    for(sw_pred_op in 1:length(stp_cox_pred)){
      stp_var = stp_cox_pred[sw_pred_op]
      stp_indx = as.numeric(names(stp_var))
      stp_pred_matrx[counter,stp_indx] = as.numeric(stp_var)[1]
    }
    
    final_stp_train_conc[counter]<-concordance(stp_cox_fit, newdata = train_data)$concordance
    final_stp_test_conc[counter] <- concordance(stp_cox_fit, newdata = test_data)$concordance
    
    x.train = subset(train_data, select = -c(surtim, surind))
    y.train = subset(train_data, select = c(surtim, surind))
    surv.train = Surv(y.train$surtim,y.train$surind)
    x.train.m = model.matrix(surv.train~.+0,data=x.train)
    y.test = subset(test_data, select = c(surtim, surind))
    surv.test = Surv(y.test$surtim,y.test$surind)
    x.test = subset(test_data, select = -c(surtim, surind))
    x.test.m = model.matrix(surv.test ~.+0,data=x.test)
    
    lasso.cvfit = cv.glmnet(x.train.m,  surv.train, family = "cox",type.measure ="C", folds=folds,maxit=100000)
    
    lasso.glmnet.fit = glmnet(x.train.m, surv.train,family ="cox", lambda = lasso.cvfit$lambda.min)
    lasso.predict.train = predict(lasso.glmnet.fit,newx=x.train.m,type="link")[,1]
    lasso.predict.test = predict(lasso.glmnet.fit,newx=x.test.m,type="link")[,1]
    for(lass_pred_op in 1:length(lasso.predict.test)){
      lasso_var = lasso.predict.test[lass_pred_op]
      lasso_indx = as.numeric(names(lasso_var))
      lasso_pred_matrix[counter,lasso_indx] = as.numeric(lasso_var)[1]
    }
    lasso_train_conc[counter] =  Cindex(lasso.predict.train,y= surv.train)
    lasso_test_conc[counter] =  Cindex(lasso.predict.test,y= surv.test)
    lasso.feat.selected[counter,] = as.numeric(coef(lasso.glmnet.fit)!=0)
    
    ridge.cvfit = cv.glmnet(x.train.m,  surv.train, family = "cox",alpha=0,type.measure ="C", folds=folds,maxit=100000)
    
    ridge.glmnet.fit = glmnet(x.train.m, surv.train,family ="cox", alpha=0,lambda = ridge.cvfit$lambda.min)
    ridge.predict.train = predict(ridge.glmnet.fit,newx=x.train.m,type="link")[,1]
    ridge.predict.test = predict(ridge.glmnet.fit,newx=x.test.m,type="link")[,1] 
    for(ridg_pred_op in 1:length(ridge.predict.test)){
      ridge_var = ridge.predict.test[ridg_pred_op]
      ridge_indx = as.numeric(names(ridge_var))
      ridge_pred_matrix[counter,ridge_indx] = as.numeric(ridge_var)[1]
    }
    ridge_train_conc[counter] =  Cindex(ridge.predict.train,y= surv.train)
    ridge_test_conc[counter] =  Cindex(ridge.predict.test,y= surv.test)
    
    gbm_model = gbm(surv.train~., data=x.train, distribution="coxph")
    gbm_summary[[counter]] = summary(gbm_model)
    gbm_pred_train = predict(gbm_model, newdata=train_data, type="link")
    gbm_pred_test = predict(gbm_model, newdata=test_data, type="link")
    gbm.train.conc[counter] = Cindex(gbm_pred_train, y=surv.train)
    gbm.test.conc[counter] = Cindex(gbm_pred_test, y=surv.test)
    gbm_indxs = rownames(test_data)
    for(indx_pos in 1:length(gbm_indxs)){
      gbm.indx.pos = as.numeric(gbm_indxs[indx_pos])
      gbm_pred_matrix[counter,gbm.indx.pos] = gbm_pred_test[indx_pos]
    }
    
    temp_rsf_model = rfsrc(Surv(surtim, surind) ~ .,data = train_data, importance = TRUE)
    temp_rsf_pred = predict.rfsrc(temp_rsf_model,newdata = test_data,type="lp")
    rsf_conc_train[counter] = Cindex(temp_rsf_model$predicted.oob,y=surv.train)
    rsf_conc_test[counter] = Cindex(temp_rsf_pred$predicted,y=surv.test)
    counter = counter + 1
  }
}


round(median(unlist(final_stp_train_conc)),4)
round(std.error(unlist(final_stp_train_conc)),4)
round(median(unlist(final_stp_test_conc)),4)
round(std.error(unlist(final_stp_test_conc)),4)
round(median(lasso_train_conc),4)
round(std.error(lasso_train_conc),4)
round(median(lasso_test_conc),4)
round(std.error(lasso_test_conc),4)

round(median(ridge_train_conc),4)
round(std.error(ridge_train_conc),4)
round(median(ridge_test_conc),4)
round(std.error(ridge_test_conc),4)

round(median(gbm.train.conc),4)
round(std.error(gbm.train.conc),4)
round(median(gbm.test.conc),4)
round(std.error(gbm.test.conc),4)

round(median(unlist(rsf_conc_train)),4)
round(std.error(unlist(rsf_conc_train)),4)
round(median(unlist(rsf_conc_test)),4)
round(std.error(unlist(rsf_conc_test)),4)
# 

########################################
# boxplot
########################################
boxplot(unlist(final_stp_train_conc),unlist(final_stp_test_conc), lasso_train_conc, lasso_test_conc,
        ridge_train_conc,ridge_test_conc,gbm.train.conc,gbm.test.conc, unlist(rsf_conc_train),
        unlist(rsf_conc_test), col= c(0,5,0,5,0,5,0,5,0,5),
        names=c("Stepwise","","Lasso","","Ridge","", "GBM","","RSF",""), 
        main="Train Test concordances for optimal features")
legend("topright",c("train concordance","test concordance"),fill=c(0,5),cex=0.8)


##################################
## KM plots
#################################

stp_mean_pred = as.array(apply(stp_pred_matrx,2,mean))
stp_kmo = split.kms(stp_mean_pred, os=data_subset$surtim, oevn=data_subset$surind)
stp_temp_var = stp_kmo$iz
stp_km.split = survfit(Surv(data_subset$surtim, data_subset$surind)~stp_temp_var, data=data_subset)
ggsurvplot(stp_km.split,title="Overall Survival for Cox model",pval = TRUE,conf.int = TRUE,risk.table = TRUE,legend.labs =
             c("high risk", "low risk"))


lasso_mean_pred = as.array(apply(lasso_pred_matrix,2,mean))
lasso_kmo = split.kms(lasso_mean_pred, os=data_subset$surtim, oevn=data_subset$surind)
lasso_temp_var = lasso_kmo$iz
lasso_km.split = survfit(Surv(data_subset$surtim, data_subset$surind)~lasso_temp_var, data=data_subset)
ggsurvplot(lasso_km.split,title="Overall Survival for Penalised Cox model(lasso)",pval = TRUE,risk.table = TRUE,conf.int = TRUE,legend.labs =
             c("high risk", "low risk"))


ridge_mean_pred = as.array(apply(ridge_pred_matrix,2,mean))
ridge_kmo = split.kms(ridge_mean_pred, os=data_subset$surtim, oevn=data_subset$surind)
ridge_temp_var = ridge_kmo$iz
ridge_km.split = survfit(Surv(data_subset$surtim, data_subset$surind)~ridge_temp_var, data=data_subset)
ggsurvplot(ridge_km.split,title="Overall Survival for Penalised Cox model(ridge)",pval = TRUE,risk.table = TRUE,conf.int = TRUE,legend.labs =
             c("high risk", "low risk"))

gbm_mean_pred = as.array(apply(gbm_pred_matrix,2,mean))
gbm_kmo = split.kms(gbm_mean_pred, os=data_subset$surtim, oevn=data_subset$surind)
gbm_temp_var = gbm_kmo$iz
gbm_km.split = survfit(Surv(data_subset$surtim, data_subset$surind)~gbm_temp_var, data=data_subset)
ggsurvplot(gbm_km.split,title="Overall Survival for Boosting Cox model",pval = TRUE,risk.table = TRUE,conf.int = TRUE,legend.labs =
             c("high risk", "low risk"))
