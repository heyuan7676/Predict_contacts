function [fpr_test1, tpr_test1, fpr_test2, tpr_test2] = train_and_test(x, x_test1, x_test2, y, y_test1, y_test2, algo)
if strcmp(algo, 'LogisticRegression')    % Logistic regression
    [b,dev,stats] = glmfit(x, y, 'binomial');
    disp('Logistic regression coefficient p-values:');
    disp(stats.p);
    
    yhat1 = glmval(b,x_test1,'logit');
    yhat2 = glmval(b,x_test2,'logit');

    [fpr_test1,tpr_test1,T,AUC_test1] = perfcurve(y_test1,yhat1,1);
    [fpr_test2,tpr_test2,T,AUC_test2] = perfcurve(y_test2,yhat2,1);
    disp('Logistic regression test1 AUC:');
    disp(AUC_test1);
    disp('Logistic regression test2 AUC:');
    disp(AUC_test2);
elseif strcmp(algo, 'RandomForests')    %random forest
    B = TreeBagger(100,x,y,'OOBVarImp','on');
    plot(oobError(B))
    xlabel('number of grown trees')
    ylabel('out-of-bag classification error')
    
    disp('RandomForests feature importance:');
    disp(B.OOBPermutedVarDeltaError);
    
    [yhat1 yscore1] = predict(B,x_test1);
    [yhat2 yscore2] = predict(B,x_test2);
    
    [fpr_test1,tpr_test1,T,AUC_test1] = perfcurve(y_test1,yscore1(:,2),1);
    [fpr_test2,tpr_test2,T,AUC_test2] = perfcurve(y_test2,yscore2(:,2),1);
    disp('Random forests test1 AUC:');
    disp(AUC_test1);
    disp('Random forests test2 AUC:');
    disp(AUC_test2);
elseif strcmp(algo, 'SVM')    %SVM
    SVMModel = fitcsvm(x,y,'KernelFunction','rbf');
    [yhat1 yscore1]= predict(SVMModel,x_test1);
    [yhat2 yscore2]= predict(SVMModel,x_test2);
    [fpr_test1,tpr_test1,T,AUC_test1] = perfcurve(y_test1,yscore1(:,2),1);
    [fpr_test2,tpr_test2,T,AUC_test2] = perfcurve(y_test2,yscore2(:,2),1);
    disp('SVM test1 AUC:');
    disp(AUC_test1);
    disp('SVM test2 AUC:');
    disp(AUC_test2);
else
    disp(strcat('Algorithm ',algo,' not support.'));
end


