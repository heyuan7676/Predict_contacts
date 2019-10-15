library(BayesTree)
library(neuralnet)
library(AUC)
library(randomForest)
library(e1071)

chr = "chr22"

### read in train data and test data
base_dir = "/scratch1/battle-fs1/heyuan/HiC/chromovar3d.stanford.edu/project"
data = read.table(paste0(base_dir,"/histone/chr22_60smp_train.txt"),header=T)
test_data1 = read.table(paste0(base_dir,"/histone/chr22_60smp_test.txt"),header=T)
chr21_1  = read.table(paste0(base_dir,"/histone/chr21_pos.csv"),header=T)
chr21_1 = cbind(chr21_1,"y"=rep("1",nrow(chr21_1)))
chr21_2  = read.table(paste0(base_dir,"/histone/chr21_neg.csv"),header=T)
chr21_2 = cbind(chr21_2,"y"=rep("0",nrow(chr21_2)))
test_data2 = rbind(chr21_1,chr21_2) 





#################################################################################
## Feature set 1: 9 features
#################################################################################


###########################
## logistic model
###########################
logit_model = glm(y~.,data=data,family="binomial")
print("Print the odds ratios and 95% CIs for the coefficients")
print(exp(cbind(OR = coef(logit_model), confint(logit_model))))
summary(logit_model)

logit_function <- function(model,new_data){
    predictions = predict(model,newdata=new_data,type="response")
    ys = data.frame("predictions"=predictions,"response"=factor(new_data$y,levels=c(0,1),labels=c(0,1)))
    roc_object = roc(ys$predictions,ys$response)
    auc_area = auc(roc_object)
    print(paste0("The AUC for the logistic model is ",auc_area))
    return(roc_object)
}


## predict
roc1_logit = logit_function(logit_model,test_data1)   ### 0.800
roc2_logit = logit_function(logit_model,test_data2)   ### 0.801



#########################
## random forest
#########################

rf_function <- function(model,new_data){
    predictions = predict(model,newdata=new_data,type="prob")[,2]
    ys = data.frame("predictions"=predictions,"response"=factor(new_data$y,levels=c(0,1),labels=c(0,1)))
    roc_object = roc(ys$predictions,ys$response)
    auc_area = auc(roc_object)
    print(paste0("The AUC for the logistic model is ",auc_area))
    return(roc_object)
}


    ntree = 200
    rf = randomForest(factor(y) ~ ., data, importance=T, ntree=ntree, norm.votes=FALSE)
    roc1_rf = rf_function(rf,test_data1)
    roc2_rf = rf_function(rf,test_data2)


#########################
## neural network
#########################

nn_function <- function(model,new_data){
    predictions = compute(model,new_data[,1:9])$net.result
    ys = data.frame("predictions" = predictions,"response" = factor(new_data$y,levels=c(0,1),labels=c("0","1")))
    roc_object = roc(ys$predictions,ys$response)
    auc_area = auc(roc_object)
    print(paste0("The AUC for the logistic model is ",auc_area))
    return(roc_object)
}

    i=5
    nn = neuralnet(y~H3K4ME3_1+H3K4ME3_2+H3K4ME3_3+H3K4ME1_1+H3K4ME1_2+H3K4ME1_3+H3K27AC_1+H3K27AC_2+H3K27AC_3,data=data,hidden=i,linear.output=FALSE)
    roc1_nn = nn_function(nn,test_data1)
    roc2_nn = nn_function(nn,test_data2)




#########################
## SVM
#########################


#### read in SVM manullay
roc1_svm = read.table("/home/qliu24/ML_project/SVM/ROC_60smp.txt",header=T)
roc2_svm = read.table("/home/qliu24/ML_project/SVM/chr21/ROC_chr21.txt",header=T)
colnames(roc1_svm) = c("fpr","tpr")
colnames(roc2_svm) = c("fpr","tpr")



#########################
## naiveBayes
#########################

nBayes = naiveBayes(y~H3K4ME3_1+H3K4ME3_2+H3K4ME3_3+H3K4ME1_1+H3K4ME1_2+H3K4ME1_3+H3K27AC_1+H3K27AC_2+H3K27AC_3,data=data)

nBayes_function <- function(model,new_data){
    predictions = predict(model,newdata=new_data,type="raw")[,2]
    ys = data.frame("predictions"=predictions,"response"=factor(new_data$y,levels=c(0,1),labels=c(0,1)))
    roc_object = roc(ys$predictions,ys$response)
    auc_area = auc(roc_object)
    print(paste0("The AUC for the logistic model is ",auc_area))
    return(roc_object)
}

roc1_nB = nBayes_function(nBayes,test_data1)   ### 0.7639
roc2_nB = nBayes_function(nBayes,test_data2)   ### 0.7523


#########################
## BART model
#########################

bart_function <- function(train_data,test_data){
    bart_model = bart(train_data[,c(1:9)],train_data$y,x.test = test_data[,c(1:9)],verbose=F)
    ### variable importance
    times = apply(bart_model$varcount,2,sum)
    prop = times/sum(times)
    yhat_train_mean = apply(pnorm(bart_model$yhat.train),2,mean) 
    yhat_test_mean = apply(pnorm(bart_model$yhat.test),2,mean)
    ys = data.frame("predictions"=yhat_test_mean,"response"=factor(test_data$y,levels=c(0,1),labels=c(0,1)))
    roc_object = roc(ys$predictions,ys$response) 
    auc_area = auc(roc_object)
    print(paste0("The AUC for the logistic model is ",auc_area))
    return(roc_object) 
}


roc1_bart = bart_function(data,test_data1)   ### 0.9294
roc2_bart = bart_function(data,test_data2)   ### 0.9352



#########################
## plot
#########################
library(ggplot2)

roc_plot <- function(model,roc_object){
    roc_for_plot = data.frame("model"=rep(model,length(roc_object$fpr)),"fpr" = roc_object$fpr, "tpr" = roc_object$tpr ) 
    return(roc_for_plot)
}

logit_plot1 = roc_plot("logistic(0.80)",roc1_logit)
rf_plot1 = roc_plot("random_forest(0.98)",roc1_rf)
nn_plot1 = roc_plot("neural_network(0.90)",roc1_nn)
svm_plot1 = roc_plot("SVM (0.98)",roc1_svm)
nB_plot1 = roc_plot("Naives_Bayes(0.76)",roc1_nB)
bart_plot1 = roc_plot("BART(0.93)",roc1_bart)

plot_df1 = rbind(logit_plot1,rf_plot1,nn_plot1,svm_plot1,nB_plot1,bart_plot1)

g1 = ggplot(plot_df1,aes(x=fpr,y=tpr,col=model))+
  geom_line(aes(y=tpr))+
  labs(title="Fig 2a. ROC curve for the test sample 1(Feature set 1)",x="False Positive Rate",y="True Positive Rate") +
  theme(plot.title=element_text(size=7),axis.text.x=element_text(angle=45, hjust=1, size=5),axis.text.y=element_text(vjust=0, size=5),legend.text=element_text(size=5)) +
  theme(legend.title=element_text(size=5),axis.title.x=element_text(size=6),axis.title.y=element_text(size=6))



logit_plot2 = roc_plot("logistic(0.69)",roc2_logit)
rf_plot2 = roc_plot("random_forest(0.81)",roc2_rf)
nn_plot2 = roc_plot("neural_network(0.90)",roc2_nn)
svm_plot2 = roc_plot("SVM (0.80)",roc2_svm)
nB_plot2 = roc_plot("Naives_Bayes(0.75)",roc2_nB)
bart_plot2 = roc_plot("BART(0.83)",roc2_bart)
          
plot_df2 = rbind(logit_plot2,rf_plot2,nn_plot2,svm_plot2,nB_plot2,bart_plot2)


g2 = ggplot(plot_df2,aes(x=fpr,col=model))+
  geom_line(aes(y=tpr))+
  labs(title="Fig 2b. ROC curve for the test sample 2(Feature set 1)",x="False Positive Rate",y="True Positive Rate") + 
  theme(plot.title=element_text(size=7),axis.text.x=element_text(angle=45, hjust=1, size=5),axis.text.y=element_text(vjust=0, size=5),legend.text=element_text(size=5)) +
  theme(legend.title=element_text(size=5),axis.title.x=element_text(size=6),axis.title.y=element_text(size=6))








#################################################################################
## Feature set 2: 3 sums
#################################################################################



### read in train data and test data
base_dir = "/scratch1/battle-fs1/heyuan/HiC/chromovar3d.stanford.edu/project"
data = read.table(paste0(base_dir,"/histone/chr22_60smp_train.txt"),header=T)
test_data1 = read.table(paste0(base_dir,"/histone/chr22_60smp_test.txt"),header=T)
chr21_1  = read.table(paste0(base_dir,"/histone/chr21_pos.csv"),header=T)
chr21_1 = cbind(chr21_1,"y"=rep("1",nrow(chr21_1)))
chr21_2  = read.table(paste0(base_dir,"/histone/chr21_neg.csv"),header=T)
chr21_2 = cbind(chr21_2,"y"=rep("0",nrow(chr21_2)))
test_data2 = rbind(chr21_1,chr21_2)

computeSum <- function(data){
    data = data.frame("H3K4ME3"=(data[,1]+data[,2]+data[,3]),"H3K4ME1"=(data[,4]+data[,5]+data[,6]),"H3K27AC"=(data[,7]+data[,8]+data[,9]),"y"=data[,10])
    return(data)
}

data = computeSum(data)
test_data1 = computeSum(test_data1)
test_data2 = computeSum(test_data2)


###########################
## logistic model
###########################
logit_model = glm(y~.,data=data,family="binomial")
print("Print the odds ratios and 95% CIs for the coefficients")
print(exp(cbind(OR = coef(logit_model), confint(logit_model))))
summary(logit_model)

logit_function <- function(model,new_data){
    predictions = predict(model,newdata=new_data,type="response")
    ys = data.frame("predictions"=predictions,"response"=factor(new_data$y,levels=c(0,1),labels=c(0,1)))
    roc_object = roc(ys$predictions,ys$response)
    auc_area = auc(roc_object)
    print(paste0("The AUC for the logistic model is ",auc_area))
    return(roc_object)
}


## predict
roc1_logit = logit_function(logit_model,test_data1)   ### 0.8151
roc2_logit = logit_function(logit_model,test_data2)   ### 0.7997



#########################
## random forest
#########################

rf_function <- function(model,new_data){
    predictions = predict(model,newdata=new_data,type="prob")[,2]
    ys = data.frame("predictions"=predictions,"response"=factor(new_data$y,levels=c(0,1),labels=c(0,1)))
    roc_object = roc(ys$predictions,ys$response)
    auc_area = auc(roc_object)
    print(paste0("The AUC for the logistic model is ",auc_area))
    return(roc_object)
}

    i=10
    ntree = 20*i
    rf = randomForest(factor(y) ~ ., data, importance=T, ntree=ntree, norm.votes=FALSE)
    roc1_rf = rf_function(rf,test_data1)   ### 0.97
    roc2_rf = rf_function(rf,test_data2)   ### 0.94








#########################
## neural network
#########################

nn_function <- function(model,new_data){
    predictions = compute(model,new_data[,1:3])$net.result
    ys = data.frame("predictions" = predictions,"response" = factor(new_data$y,levels=c(0,1),labels=c("0","1")))
    roc_object = roc(ys$predictions,ys$response)
    auc_area = auc(roc_object)
    print(paste0("The AUC for the logistic model is ",auc_area))
    return(roc_object)
}

    nn = neuralnet(y~H3K4ME3+H3K4ME1+H3K27AC,data=data,hidden=2,linear.output=FALSE)
    roc1_nn = nn_function(nn,test_data1) ### 0.86
    roc2_nn = nn_function(nn,test_data2)   ### 0.88





#########################
## SVM
#########################

#### read in SVM manullay
roc1_svm = read.table("/home/qliu24/ML_project/SVM/3binsumup/ROC_60smp.txt",header=T)
roc2_svm = read.table("/home/qliu24/ML_project/SVM/chr21/ROC_chr21_b1.txt",header=T)
colnames(roc1_svm) = c("fpr","tpr")
colnames(roc2_svm) = c("fpr","tpr")



#########################
## naiveBayes
#########################

nBayes = naiveBayes(y~H3K4ME3 + H3K4ME1 + H3K27AC,data=data)

nBayes_function <- function(model,new_data){
    predictions = predict(model,newdata=new_data,type="raw")[,2]
    ys = data.frame("predictions"=predictions,"response"=factor(new_data$y,levels=c(0,1),labels=c(0,1)))
    roc_object = roc(ys$predictions,ys$response)
    auc_area = auc(roc_object)
    print(paste0("The AUC for the logistic model is ",auc_area))
    return(roc_object)
}

roc1_nB = nBayes_function(nBayes,test_data1)   ### 0.79
roc2_nB = nBayes_function(nBayes,test_data2)   ### 0.77


#########################
## BART model
#########################

bart_function <- function(train_data,test_data){
    bart_model = bart(train_data[,c(1:3)],train_data$y,x.test = test_data[,c(1:3)],verbose=F)
    ### variable importance
    times = apply(bart_model$varcount,2,sum)
    prop = times/sum(times)
    yhat_train_mean = apply(pnorm(bart_model$yhat.train),2,mean) 
    yhat_test_mean = apply(pnorm(bart_model$yhat.test),2,mean)
    ys = data.frame("predictions"=yhat_test_mean,"response"=factor(test_data$y,levels=c(0,1),labels=c(0,1)))
    roc_object = roc(ys$predictions,ys$response) 
    auc_area = auc(roc_object)
    print(paste0("The AUC for the logistic model is ",auc_area))
    return(roc_object) 
}


roc1_bart = bart_function(data,test_data1)   ### 0.92
roc2_bart = bart_function(data,test_data2)   ### 0.91




#########################
## plot
#########################
library(ggplot2)

roc_plot <- function(model,roc_object){
    roc_for_plot = data.frame("model"=rep(model,length(roc_object$fpr)),"fpr" = roc_object$fpr, "tpr" = roc_object$tpr ) 
    return(roc_for_plot)
}

logit_plot1 = roc_plot("logistic(0.82)",roc1_logit)
rf_plot1 = roc_plot("random_forest(0.97)",roc1_rf)
nn_plot1 = roc_plot("neural_network(0.85)",roc1_nn)
svm_plot1 = roc_plot("SVM(0.95)",roc1_svm)
nB_plot1 = roc_plot("Naives_Bayes(0.79)",roc1_nB)
bart_plot1 = roc_plot("BART(0.92)",roc1_bart)

plot_df1 = rbind(logit_plot1,rf_plot1,nn_plot1,svm_plot1,nB_plot1,bart_plot1)

g3 = ggplot(plot_df1,aes(x=fpr,y=tpr,col=model))+
  geom_line(aes(y=tpr))+
  labs(title="Fig 2c. ROC curve for the test sample 1 (Feature set 2)",x="False Positive Rate",y="True Positive Rate") +
  theme(plot.title=element_text(size=7),axis.text.x=element_text(angle=45, hjust=1, size=5),axis.text.y=element_text(vjust=0, size=5),legend.text=element_text(size=5)) +
  theme(legend.title=element_text(size=5),axis.title.x=element_text(size=6),axis.title.y=element_text(size=6))



logit_plot2 = roc_plot("logistic(0.62)",roc2_logit)
rf_plot2 = roc_plot("random_forest(0.76)",roc2_rf)
nn_plot2 = roc_plot("neural_network(0.85)",roc2_nn)
svm_plot2 = roc_plot("SVM(0.79)",roc2_svm)
nB_plot2 = roc_plot("Naives_Bayes(0.52)",roc2_nB)
bart_plot2 = roc_plot("BART(0.81)",roc2_bart)
          
plot_df2 = rbind(logit_plot2,rf_plot2,nn_plot2,svm_plot2,nB_plot2,bart_plot2)


g4 = ggplot(plot_df2,aes(x=fpr,col=model))+
  geom_line(aes(y=tpr))+
  labs(title="Fig 2d. ROC curve for the test sample 2(Feature set2)",x="False Positive Rate",y="True Positive Rate") + 
  theme(plot.title=element_text(size=7),axis.text.x=element_text(angle=45, hjust=1, size=5),axis.text.y=element_text(vjust=0, size=5),legend.text=element_text(size=5)) +
  theme(legend.title=element_text(size=5),axis.title.x=element_text(size=6),axis.title.y=element_text(size=6))





#################################################################################
## Feature set 2: 3 best
#################################################################################



### read in train data and test data
base_dir = "/scratch1/battle-fs1/heyuan/HiC/chromovar3d.stanford.edu/project"
data = read.table(paste0(base_dir,"/histone/chr22_60smp_train.txt"),header=T)
test_data1 = read.table(paste0(base_dir,"/histone/chr22_60smp_test.txt"),header=T)
chr21_1  = read.table(paste0(base_dir,"/histone/chr21_pos.csv"),header=T)
chr21_1 = cbind(chr21_1,"y"=rep("1",nrow(chr21_1)))
chr21_2  = read.table(paste0(base_dir,"/histone/chr21_neg.csv"),header=T)
chr21_2 = cbind(chr21_2,"y"=rep("0",nrow(chr21_2)))
test_data2 = rbind(chr21_1,chr21_2)

pickout = c("H3K4ME3_1","H3K4ME3_2","H3K4ME1_1")
data = data[,c(pickout,"y")]
test_data1 = test_data1[,c(pickout,"y")]
test_data2 = test_data2[,c(pickout,"y")]


###########################
## logistic model
###########################
logit_model = glm(y~.,data=data,family="binomial")
print("Print the odds ratios and 95% CIs for the coefficients")
print(exp(cbind(OR = coef(logit_model), confint(logit_model))))
summary(logit_model)

logit_function <- function(model,new_data){
    predictions = predict(model,newdata=new_data,type="response")
    ys = data.frame("predictions"=predictions,"response"=factor(new_data$y,levels=c(0,1),labels=c(0,1)))
    roc_object = roc(ys$predictions,ys$response)
    auc_area = auc(roc_object)
    print(paste0("The AUC for the logistic model is ",auc_area))
    return(roc_object)
}


## predict
roc1_logit = logit_function(logit_model,test_data1)   ### 0.8151
roc2_logit = logit_function(logit_model,test_data2)   ### 0.7997



#########################
## random forest
#########################

rf_function <- function(model,new_data){
    predictions = predict(model,newdata=new_data,type="prob")[,2]
    ys = data.frame("predictions"=predictions,"response"=factor(new_data$y,levels=c(0,1),labels=c(0,1)))
    roc_object = roc(ys$predictions,ys$response)
    auc_area = auc(roc_object)
    print(paste0("The AUC for the logistic model is ",auc_area))
    return(roc_object)
}

    ntree=200
    rf = randomForest(factor(y) ~ ., data, importance=T, ntree=ntree, norm.votes=FALSE)
    roc1_rf = rf_function(rf,test_data1)   ### 0.97
    roc2_rf = rf_function(rf,test_data2)   ### 0.94


#########################
## neural network
#########################

nn_function <- function(model,new_data){
    predictions = compute(model,new_data[,1:3])$net.result
    ys = data.frame("predictions" = predictions,"response" = factor(new_data$y,levels=c(0,1),labels=c("0","1")))
    roc_object = roc(ys$predictions,ys$response)
    auc_area = auc(roc_object)
    print(paste0("The AUC for the logistic model is ",auc_area))
    return(roc_object)
}

    nn = neuralnet(y~H3K4ME3_1+H3K4ME3_2+H3K4ME1_1,data=data,hidden=2,linear.output=FALSE)
    roc1_nn = nn_function(nn,test_data1) ### 0.86
    roc2_nn = nn_function(nn,test_data2)   ### 0.88




#########################
## SVM
#########################


#### read in SVM manullay
roc1_svm = read.table("/home/qliu24/ML_project/SVM/3binbest/ROC_60smp.txt",header=T)
roc2_svm = read.table("/home/qliu24/ML_project/SVM/chr21/ROC_chr21_b2.txt",header=T)
colnames(roc1_svm) = c("fpr","tpr")
colnames(roc2_svm) = c("fpr","tpr")




#########################
## naiveBayes
#########################

nBayes = naiveBayes(y~H3K4ME3_1+H3K4ME3_2+H3K4ME1_1,data=data)

nBayes_function <- function(model,new_data){
    predictions = predict(model,newdata=new_data,type="raw")[,2]
    ys = data.frame("predictions"=predictions,"response"=factor(new_data$y,levels=c(0,1),labels=c(0,1)))
    roc_object = roc(ys$predictions,ys$response)
    auc_area = auc(roc_object)
    print(paste0("The AUC for the logistic model is ",auc_area))
    return(roc_object)
}

roc1_nB = nBayes_function(nBayes,test_data1)   ### 0.79
roc2_nB = nBayes_function(nBayes,test_data2)   ### 0.77


#########################
## BART model
#########################

bart_function <- function(train_data,test_data){
    bart_model = bart(train_data[,c(1:3)],train_data$y,x.test = test_data[,c(1:3)],verbose=F)
    ### variable importance
    times = apply(bart_model$varcount,2,sum)
    prop = times/sum(times)
    yhat_train_mean = apply(pnorm(bart_model$yhat.train),2,mean) 
    yhat_test_mean = apply(pnorm(bart_model$yhat.test),2,mean)
    ys = data.frame("predictions"=yhat_test_mean,"response"=factor(test_data$y,levels=c(0,1),labels=c(0,1)))
    roc_object = roc(ys$predictions,ys$response) 
    auc_area = auc(roc_object)
    print(paste0("The AUC for the logistic model is ",auc_area))
    return(roc_object) 
}


roc1_bart = bart_function(data,test_data1)   ### 0.92
roc2_bart = bart_function(data,test_data2)   ### 0.91




#########################
## plot
#########################
library(ggplot2)

roc_plot <- function(model,roc_object){
    roc_for_plot = data.frame("model"=rep(model,length(roc_object$fpr)),"fpr" = roc_object$fpr, "tpr" = roc_object$tpr ) 
    return(roc_for_plot)
}

logit_plot1 = roc_plot("logistic(0.79)",roc1_logit)
rf_plot1 = roc_plot("random_forest(0.95)",roc1_rf)
nn_plot1 = roc_plot("neural_network(0.80)",roc1_nn)
svm_plot1 = roc_plot("SVM(0.89)",roc1_svm)
nB_plot1 = roc_plot("Naives_Bayes(0.76)",roc1_nB)
bart_plot1 = roc_plot("BART(0.83)",roc1_bart)

plot_df1 = rbind(logit_plot1,rf_plot1,nn_plot1,svm_plot1,nB_plot1,bart_plot1)

g5 = ggplot(plot_df1,aes(x=fpr,y=tpr,col=model))+
  geom_line(aes(y=tpr))+
  labs(title="Fig 2e.  ROC curve for the test sample 1(Feature set3)",x="False Positive Rate",y="True Positive Rate") +
  theme(plot.title=element_text(size=7),axis.text.x=element_text(angle=45, hjust=1, size=5),axis.text.y=element_text(vjust=0, size=5),legend.text=element_text(size=5)) +
  theme(legend.title=element_text(size=5),axis.title.x=element_text(size=6),axis.title.y=element_text(size=6))



logit_plot2 = roc_plot("logistic(0.49)",roc2_logit)
rf_plot2 = roc_plot("random_forest(0.63)",roc2_rf)
nn_plot2 = roc_plot("neural_network(0.69)",roc2_nn)
svm_plot2 = roc_plot("SVM(0.58)",roc2_svm)
nB_plot2 = roc_plot("Naives_Bayes(0.52)",roc2_nB)
bart_plot2 = roc_plot("BART(0.72)",roc2_bart)
          

plot_df2 = rbind(logit_plot2,rf_plot2,nn_plot2,svm_plot2,nB_plot2,bart_plot2)


g6 = ggplot(plot_df2,aes(x=fpr,col=model))+
  geom_line(aes(y=tpr))+
  labs(title="Fig 2f. ROC curve for the test sample 2(Feature set3)",x="False Positive Rate",y="True Positive Rate") + 
  theme(plot.title=element_text(size=7),axis.text.x=element_text(angle=45, hjust=1, size=5),axis.text.y=element_text(vjust=0, size=5),legend.text=element_text(size=5)) +
  theme(legend.title=element_text(size=5),axis.title.x=element_text(size=6),axis.title.y=element_text(size=6))






####### multiplot


multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}





pdf("Fig 2.pdf")
multiplot(g1,g3,g5,g2,g4,g6,cols=2)
dev.off()

