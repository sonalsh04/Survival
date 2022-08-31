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

table(dat$surind)
# hist(as.numeric(dat$surtim))
#convert to numeric
dat$sex =( as.numeric(as.factor(dat$sex)))
dat$grade = (as.numeric(as.factor(dat$grade)))
dat$tumor.subtype =( as.numeric(as.factor(dat$tumor.subtype)))

###############################################################################
## removing correlated variables
###############################################################################

cor_matrix <- cor(dat)
cor_matrix_rm <- cor_matrix
cor_matrix_rm[upper.tri(cor_matrix_rm)] <- 0
diag(cor_matrix_rm) <- 0
data_subset <- dat[ , !apply(cor_matrix_rm,
                             2,
                             function(x) any(abs(x) > 0.90))]

###############################################################################
data_subset = dat

###############################################################################
## Distribution with Data
###############################################################################
hist(dat$surtim, main="Distribution of Survival Time", xlab="Survival Time",col=5)

count_sex = table(dat$sex)
barplot(count_sex, names.arg = c("Female","Male"), col=c(2,5), main="Gender")

count_grade = table(dat$grade)
barplot(count_grade, names.arg = c("High","Intermediate","Low"), col=c(2,4,5), main="Grade")

count_subtype = table(dat$tumor.subtype)
barplot(count_grade, names.arg = c("STS","Bone","Cartilage"), col=c(2,4,5), main="Tumour Subtype")

count_surind = table(dat$surind)
barplot(count_surind,col=c(2,4), main="Survival Indicator")

###############################################################################
## Function for plotting KM
###############################################################################

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

###############################################################################
# Common
###############################################################################
repeats = 10
folds = 5

prod = folds*repeats
threshold = (60*prod)/100
original_cols = ncol(dat)-2
subset_cols = ncol(data_subset)-2
subset_rows = nrow(data_subset)
subset_nms = colnames(subset(data_subset, select = -c(surtim, surind)))

stp_feat_selected = matrix(0,nrow=prod,ncol=subset_cols)
colnames(stp_feat_selected) = subset_nms
stp_pred_matrx = matrix(0,nrow=prod,ncol=subset_rows)
full_step_model = vector("list",folds)
# stp_train_conc = vector("list",folds)
# stp_test_conc = vector("list",folds)

final_stp_train_conc  = vector("list",prod)
final_stp_test_conc = vector("list",prod)
final_stp_models = vector("list",prod)

lasso_train_conc = lasso_test_conc = numeric(prod)
lasso.feat.selected = matrix(0, nrow=prod, ncol=subset_cols)
colnames(lasso.feat.selected) = subset_nms
lasso_pred_matrix  = matrix(0, nrow=prod, ncol=nrow(data_subset))

ridge_cvIndex = createFolds(factor(data_subset$surind), folds,returnTrain = T)
ridge_train_conc = ridge_test_conc =  numeric(prod)
ridge.coef.matrix = matrix(0, nrow=prod, ncol=subset_cols)
colnames(ridge.coef.matrix) = subset_nms
ridge_pred_matrix  = matrix(0, nrow=prod, ncol=nrow(data_subset))

gbm_summary = vector("list",prod) 
gbm.train.conc = gbm.test.conc = vector("double",prod) 
gbm.coef.matrix = matrix(0, nrow=prod, ncol= subset_cols)
colnames(gbm.coef.matrix) = subset_nms
gbm_pred_matrix = matrix(0,nrow=prod,ncol=subset_rows)

rsf_model = rsf_pred = vector("list",prod)
rsf_conc_train = rsf_bs = vector("list",prod)
rsf_conc_test = vector("list",prod)
rsf_conc_test = vector("list",prod)
rsf_feat_sel = matrix(0,nrow=prod,ncol=subset_cols)
colnames(rsf_feat_sel) = subset_nms


##################################
## Fitting all models
##################################
n=nrow(data_subset)
counter = 1
set.seed(4060)
for(outer_loop in 1:repeats){
  shuffled_dat = data_subset[sample(1:nrow(data_subset),nrow(data_subset)), ]
  cvIndex = createFolds(factor(shuffled_dat$surind), folds,returnTrain = T)
   for(inner_loop in 1:length(cvIndex)){
     train_data = data_subset[cvIndex[[inner_loop]],]
     test_data = data_subset[-cvIndex[[inner_loop]],]

    start_cox = coxph(Surv(surtim,surind) ~ 1, data = train_data)
    full_cox = coxph(Surv(surtim,surind) ~ ., data = train_data)
    fit_step = step(start_cox, direction = "both", scope = full_cox$formula)
    stp_names = names(fit_step$coefficients)
    for(c in 1:length(stp_names)){
      stp_feat_selected[counter,stp_names[c]] = 1
    }
    full_form = fit_step$formula[3]
    full_form_1 = formula(paste("Surv(surtim,surind)~",full_form))

    #fitting cox model on features selected by stepwise
    stp_cox_fit = coxph(full_form_1,data=train_data)
    stp_cox_pred = predict(stp_cox_fit,newdata = test_data,type="lp")

    for(sw_pred_op in 1:length(stp_cox_pred)){
      stp_var = stp_cox_pred[sw_pred_op]
      stp_indx = as.numeric(names(stp_var))
      stp_pred_matrx[counter,stp_indx] = as.numeric(stp_var)[1]
    }

    final_stp_models[counter] <- full_form
    final_stp_train_conc[counter]<-concordance(stp_cox_fit, newdata = train_data)$concordance
    final_stp_test_conc[counter] <- concordance(stp_cox_fit, newdata = test_data)$concordance

  #   ######
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
    lasso.predict.test = predict(lasso.glmnet.fit,newx=x.test.m,type="link")[,1] # the fitted relative-risk for "cox";
    for(lass_pred_op in 1:length(lasso.predict.test)){
      lasso_var = lasso.predict.test[lass_pred_op]
      lasso_indx = as.numeric(names(lasso_var))
      lasso_pred_matrix[counter,lasso_indx] = as.numeric(lasso_var)[1]
    }
    lasso_train_conc[counter] =  Cindex(lasso.predict.train,y= surv.train)
    lasso_test_conc[counter] =  Cindex(lasso.predict.test,y= surv.test)
    lasso.feat.selected[counter,] = as.numeric(coef(lasso.glmnet.fit)!=0)


    ############
    ridge.cvfit = cv.glmnet(x.train.m,  surv.train, family = "cox",alpha=0,type.measure ="C", folds=folds,maxit=100000)

    ridge.glmnet.fit = glmnet(x.train.m, surv.train,family ="cox", alpha=0,lambda = ridge.cvfit$lambda.min)
    ridge.predict.train = predict(ridge.glmnet.fit,newx=x.train.m,type="link")[,1] # the fitted relative-risk for "cox";
    ridge.predict.test = predict(ridge.glmnet.fit,newx=x.test.m,type="link")[,1] # the fitted relative-risk for "cox";
    for(ridg_pred_op in 1:length(ridge.predict.test)){
      ridge_var = ridge.predict.test[ridg_pred_op]
      ridge_indx = as.numeric(names(ridge_var))
      ridge_pred_matrix[counter,ridge_indx] = as.numeric(ridge_var)[1]
    }
    ridge_train_conc[counter] =  Cindex(ridge.predict.train,y= surv.train)
    ridge_test_conc[counter] =  Cindex(ridge.predict.test,y= surv.test)

    ridge.coefs = coef(ridge.glmnet.fit)[,1]
    for(indx in 1:length(ridge.coefs)){
      ridge.coef.matrix[counter,indx] = abs(as.numeric(ridge.coefs[indx]))
    }


    #################
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

    for(i in 1:prod ){
      for(j in 1:subset_cols){
        gbm_index_var=which(subset_nms==gbm_summary[[i]][j,]$var)
        gbm.coef.matrix[i,gbm_index_var]=gbm_summary[[i]][j,]$rel.inf
      }
    }

    #####
    temp_rsf_model = rfsrc(Surv(surtim, surind) ~ .,data = train_data, importance = TRUE)
    temp_rsf_pred = predict.rfsrc(temp_rsf_model,newdata = test_data,type="lp")
    rsf_model[[counter]] = temp_rsf_model
    rsf_pred[[counter]] = temp_rsf_pred
    #
    rsf_conc_train[counter] = Cindex(temp_rsf_model$predicted.oob,y=surv.train)
    rsf_conc_test[counter] = Cindex(temp_rsf_pred$predicted,y=surv.test)
    #
    rsf_var_imp = subsample(temp_rsf_model)
    rsf_feats = head(subset_nms[order(temp_rsf_model$importance, decreasing=TRUE)],10)
    for(q in 1:length(rsf_feats)){
      rsf_feat_sel[counter,rsf_feats[q]] = 1
    }
    ######
    ######
    counter = counter + 1
   }
}



######################################################################
## RFE + RFSRC
######################################################################


set.seed(4060)
rfe_train_conc = vector("integer",subset_cols)
rfe_test_conc = vector("integer",subset_cols)
variables_picked = list()

model_under_consideration = data_subset
rfe_counter = 1
while(length(model_under_consideration)>3){
  shuffled_orig_dat = model_under_consideration[sample(1:nrow(model_under_consideration),nrow(model_under_consideration)), ]
  cvIndex1 = createFolds(factor(shuffled_orig_dat$surind), folds,returnTrain = T)
  orig_train_data = shuffled_orig_dat[cvIndex1[[1]],]
  orig_test_data = shuffled_orig_dat[-cvIndex1[[1]],]
  orig.surv.train = Surv(orig_train_data$surtim,orig_train_data$surind)
  orig.surv.test = Surv(orig_test_data$surtim,orig_test_data$surind)
  model_rfsrc=rfsrc(Surv(surtim,surind) ~ .,data = orig_train_data, importance = TRUE)
  pred = predict.rfsrc(model_rfsrc,newdata = orig_test_data,type="lp")
  rfe_train_conc[rfe_counter] = Cindex(model_rfsrc$predicted.oob, y = orig.surv.train)
  rfe_test_conc[rfe_counter] = Cindex(pred$predicted, y = orig.surv.test)
  response_subset=subset(model_under_consideration, select = c(surtim,surind))
  var_imp=vimp.rfsrc(model_rfsrc)$importance
  sorted_vars= sort(var_imp,decreasing=TRUE)
  selected_vars = head(sorted_vars,-1)
  variables_picked[[rfe_counter]] = toString(names(selected_vars))
  new_subset=subset(model_under_consideration,select=c(names(selected_vars)))
  new_new_subset=cbind(new_subset,response_subset)
  model_under_consideration=new_new_subset
  rfe_counter = rfe_counter + 1
}
##################################
## Finding important variables
##################################
variables_picked = unlist(variables_picked)
indx_with_max_train_conc = which.max(rfe_train_conc)
indx_with_max_test_conc = which.max(rfe_test_conc)
rfe_feat_sel = variables_picked[indx_with_max_test_conc]

sw_feat_sel = which(colSums(stp_feat_selected)>threshold) == TRUE

lasso_feat_sel = which(colSums(lasso.feat.selected)>threshold) == TRUE

ridge_coef_sum = colSums(ridge.coef.matrix)
ridge_feat_sel = names(head(ridge_coef_sum[order(-ridge_coef_sum)],6))

gbm_coef_sum = colSums(gbm.coef.matrix)
gbm_feat_sel = names(head(gbm_coef_sum[order(-gbm_coef_sum)],6))

rsfc_feat_sel = which(colSums(rsf_feat_sel)>threshold) == TRUE

##################################
## finding concordances
##################################

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


round(median(rfe_train_conc),4)
round(std.error(rfe_train_conc),4)
round(median(rfe_test_conc),4)
round(std.error(rfe_test_conc),4)
##################################
##boxplots without pre-filtering
##################################

boxplot(unlist(final_stp_train_conc),lasso_train_conc,
        ridge_train_conc,gbm.train.conc,unlist(rsf_conc_train),
        names=c("Stepwise","Lasso","Ridge", "GBM","RSF"),
        col=c(0,2,7,4,5),main="Train Concordances")

boxplot(unlist(final_stp_test_conc),lasso_test_conc,
        ridge_test_conc,gbm.test.conc,unlist(rsf_conc_test),
        names=c("Stepwise","Lasso","Ridge", "GBM","RSF"),
        col=c(0,2,7,4,5),main="Test Concordances")

boxplot(unlist(final_stp_train_conc),unlist(final_stp_test_conc), lasso_train_conc, lasso_test_conc,
        ridge_train_conc,ridge_test_conc,gbm.train.conc,gbm.test.conc, unlist(rsf_conc_train),
        unlist(rsf_conc_test),rfe_train_conc,  col= c(0,5,0,5,0,5,0,5,0,5) ,
        names=c("Stepwise","","Lasso","","Ridge","", "GBM","","RSF",""), 
        main="Train - test concordances")
legend("topright",c("train concordance","test concordance"),fill=c(0,5))

##################################
##boxplots with pre-filterin
##################################

boxplot(unlist(final_stp_train_conc),lasso_train_conc,
  ridge_train_conc,gbm.train.conc,unlist(rsf_conc_train),rfe_train_conc,
  names=c("Stepwise","Lasso","Ridge", "GBM","RSF","RFE"),
  col=c(0,2,7,4,5,6),main="Train Concordances")

boxplot(unlist(final_stp_test_conc),lasso_test_conc,
        ridge_test_conc,gbm.test.conc,unlist(rsf_conc_test),rfe_test_conc,
        names=c("Stepwise","Lasso","Ridge", "GBM","RSF","RFE"),
        col=c(0,2,7,4,5,6),main="Test Concordances")

boxplot(unlist(final_stp_train_conc),unlist(final_stp_test_conc), lasso_train_conc, lasso_test_conc,
        ridge_train_conc,ridge_test_conc,gbm.train.conc,gbm.test.conc, unlist(rsf_conc_train),
        unlist(rsf_conc_test),rfe_train_conc,rfe_test_conc ,  col= c(0,5,0,5,0,5,0,5,0,5,0,5) ,
        names=c("Stepwise","","Lasso","","Ridge","", "GBM","","RSF","","RFE+RSF",""), 
        main="Train - Test concordances for pre-filtered Dataset")
legend("topright",c("train concordance","test concordance"),fill=c(0,5),cex=0.8)

##################################
## KM plots
#################################

stp_mean_pred = as.array(apply(stp_pred_matrx,2,mean))
stp_kmo = split.kms(stp_mean_pred, os=data_subset$surtim, oevn=data_subset$surind)
stp_temp_var = stp_kmo$iz
stp_km.split = survfit(Surv(data_subset$surtim, data_subset$surind)~stp_temp_var, data=data_subset)
ggsurvplot(stp_km.split,title="Overall Survival for Cox model",pval = TRUE,risk.table = TRUE,conf.int = TRUE,legend.labs =
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

