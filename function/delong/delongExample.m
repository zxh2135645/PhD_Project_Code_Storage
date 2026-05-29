% In this demonstration, we will be training two SVM models and compare
% their AUCs. We will be using the fisheriris dataset. Here, one SVM model
% will be trained using all features of the dataset, while the other will
% only be trained using a single feature. 

clc; clear all;

% Loading the fisherirs dataset
load fisheriris

% We will only utilize two of the enumerated species: setosa and
% versicolor. Each class has 50 samples each.
X = meas(1:100,:);

% We encode our dataset in such a way that 0 - setosa; 1 - versicolor
Y = [zeros(50,1);ones(50,1)];

% Creation of the training set
n = randperm(100,70); % randomizing 70 samples for training.
X_tr = X(n,:);
Y_tr = Y(n,:);

% Creation of the test set
X_ts = X;
X_ts(n,:) = [];
Y_ts = Y;
Y_ts(n,:) = [];

% Model 1 training and evaluation using all 4 features
model1 = fitcsvm(X_tr,Y_tr);
[~,prob] = predict(model1,X_ts);
a = logsig(prob(:,2));

% Model 2 training and evaluation using the only the first feature
model2 = fitcsvm(X_tr(:,1),Y_tr);
[~,prob] = predict(model2,X_ts(:,1));
b = logsig(prob(:,2));

% Computing for the AUC and delong statistic
[auc, p] = delong(a,b,Y_ts)
