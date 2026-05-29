function [auc, p] = delong(a,b,y)
% NOTE: This implementation of delong's test is highly non-vectorized. It
% is primarily intended for education purposes; demystifying the 
% complicated math behind delong's paper regarding the comparison of ROC 
% AUCs. 
%
% delong computes for the AUC of two models trained using the
% same training and test set but of which differs in the number of
% features considered. This function quantifies the AUC difference
% between the models by calculating the DeLong's p statistic. For more
% information about DeLong's test, you may read about the paper here:
% https://www.jstor.org/stable/2531595.
%
% 'a' - Vector of values containing the prediction probability of the
% first model for each sample in the test set.
%
% 'b' - Vector of values containing the prediction probability of the
% second model for each sample in the test set. Variables 'a' and 'b' are
% interchangable. However, 'a' is customarily the probability vector
% resulting to the training and evaluation of a model using L features,
% while 'b' is the result of training a model using L-k features, where
% L-k denotes the reduction of the feature space.
%
% 'y' - Vector of values encoding the true class label of the test
% samples. Classes should be encoded as: 1 - positive class, and 0 - 
% negative class.
%
% 'auc' - A vector containing the AUCs of the two models being compared.
% The vector is arranged in such a way: [AUC_model1, AUC_model2]
%
% 'p' - DeLong's p-value statistic.
%
% EXAMPLES: 
% 1. Compute for the p-value between the probability vectors 'a' and
% 'b' resulting from two models, given the true class label vector 'y'.
%
%       a = [0.2 0.3 0.6 0.9 0.8];
%       b = [0.1 0.15 0.8 0.7 0.4];
%       y = [0 0 1 1 1];
%       [~,p] = delong(a,b,y);
%
% 2. Compute for the p-value and AUC of two SVM models trained using the
% fisheriris dataset.
%
%       % Loading the fisherirs dataset
%       load fisheriris
% 
%       % We will only utilize two of the enumerated species: setosa and
%       % versicolor. Each class has 50 samples each.
%       X = meas(1:100,:);
% 
%       % We encode our dataset in such a way that 0 - setosa; 1 - versicolor
%       Y = [zeros(50,1);ones(50,1)];
% 
%       % Creation of the training set
%       n = randperm(100,70); % randomizing 70 samples for training.
%       X_tr = X(n,:);
%       Y_tr = Y(n,:);
% 
%       % Creation of the test set
%       X_ts = X;
%       X_ts(n,:) = [];
%       Y_ts = Y;
%       Y_ts(n,:) = [];
% 
%       % Model 1 training and evaluation using all 4 features
%       model1 = fitcsvm(X_tr,Y_tr);
%       [~,prob] = predict(model1,X_ts);
%       a = logsig(prob(:,2));
% 
%       % Model 2 training and evaluation using the only the first feature
%       model2 = fitcsvm(X_tr(:,1),Y_tr);
%       [~,prob] = predict(model2,X_ts(:,1));
%       b = logsig(prob(:,2));
% 
%       % Computing for the AUC and delong statistic
%       [auc, p] = delong(a,b,Y_ts)

    n = find(y == 0);
    m = find(y == 1);

    N = length(n);
    M = length(m);

    % computation of emperical AUCs
    % You can also use the built-in matlab function "perfcurve" to estimate
    % the AUCs, a_hat and b_hat. You will have to code it as such:
    
    % [~,~,AUC] = perfcurve(a,y,c)
    
    % where c is the value in y which encodes the positive class. 

    a_hat = 0;
    for i = 1:N
        for j = 1:M
            a_hat = a_hat + h(a(m(j)),a(n(i)));
        end
    end
    a_hat = a_hat/(M*N);

    b_hat = 0;
    for i = 1:N
        for j = 1:M
            b_hat = b_hat + h(b(m(j)),b(n(i)));
        end
    end
    b_hat = b_hat/(M*N);
    
    auc = [a_hat, b_hat];

    % computations of structural components for covarriance estimation
    v_a_10 = zeros(M,1);
    for i = 1:M
        for j = 1:N
            v_a_10(i) = v_a_10(i) + h(a(m(i)),a(n(j)));
        end
        v_a_10(i) = v_a_10(i)/N;
    end

    v_a_01 = zeros(N,1);
    for i = 1:N
        for j = 1:M
            v_a_01(i) = v_a_01(i) + h(a(m(j)),a(n(i)));
        end
        v_a_01(i) = v_a_01(i)/M;
    end

    v_b_10 = zeros(M,1);
    for i = 1:M
        for j = 1:N
            v_b_10(i) = v_b_10(i) + h(b(m(i)),b(n(j)));
        end
        v_b_10(i) = v_b_10(i)/N;
    end

    v_b_01 = zeros(N,1);
    for i = 1:N
        for j = 1:M
            v_b_01(i) = v_b_01(i) + h(b(m(j)),b(n(i)));
        end
        v_b_01(i) = v_b_01(i)/M;
    end

    % calculation of covariance matrixes
    s_10 = [];
    s_10(1,1) = s(v_a_10,v_a_10,a_hat,a_hat);
    s_10(2,2) = s(v_b_10,v_b_10,b_hat,b_hat);
    s_10(1,2) = s(v_a_10,v_b_10,a_hat,b_hat);
    s_10(2,1) = s_10(1,2);

    s_01 = [];
    s_01(1,1) = s(v_a_01,v_a_01,a_hat,a_hat);
    s_01(2,2) = s(v_b_01,v_b_01,b_hat,b_hat);
    s_01(1,2) = s(v_a_01,v_b_01,a_hat,b_hat);
    s_01(2,1) = s_01(1,2);
    
    S = s_10/M + s_01/N;
    
    % computation of z-score
    z = (a_hat-b_hat)/sqrt(S(1,1)+S(2,2)-2*S(1,2));
    z(isnan(z)) = 0;

    % translation of z-score to p-value (2 tailed)
    p = 2*normcdf(-abs(z));

    % psi function
    function out = h(x,y)
        out = 0.5;
        if(y<x)
            out = 1;
        elseif(y>x)
            out = 0;
        end
    end
    
    % covarriance matrix block function
    function out = s(x,y,x_hat,y_hat)
        out = 0;
        for k = 1:length(x)
            out = out + (x(k)-x_hat)*(y(k)-y_hat);
        end
        out = out/(length(x)-1);
    end

end

