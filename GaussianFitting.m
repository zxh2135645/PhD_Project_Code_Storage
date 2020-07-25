

x = [15:65];

%% 
%nel = unique(forhist1_remote);
count_remote = zeros(1, numel(x));
for n = 1:length(nel)
count_remote(n) = sum(forhist1_remote == x(n));
end
figure(); plot(x, count_remote); hold on; yyaxis right; plot(x_remote, y_remote1); 
%%
%nel = unique(forhist1);
count = zeros(1, numel(x));
for n = 1:length(x)
count(n) = sum(forhist1 == x(n));
end
figure(); plot(x, count); hold on; yyaxis right; plot(x, y_mi1); 

count_normalize = count / max(count);
count_remote_normalize = count_remote / max(count_remote);
%%
lambda = [0:0.1:10];
lambda2 = [0:0.1:10];
mu = [10:1:50];
sigma = [1:1:100];
f_hat = zeros(numel(x), numel(mu), numel(sigma));
%F_hat = zeros(numel(x), numel(lambda), numel(mu), numel(sigma),numel(lambda2));
sumLsq = zeros(numel(lambda), numel(mu), numel(sigma),numel(lambda2));
for j = 1:numel(mu)
    for k = 1:numel(sigma)
        f_hat(:,j,k) = normpdf(x, mu(j), sigma(k));
        for i = 1:numel(lambda)
            for l = 1:numel(lambda2)
                sumLsq(i,j,k,l) = sum((lambda(i) * count_remote_normalize' + lambda2(l) * f_hat(:,j,k) - count_normalize').^2);
            end
        end
    end
end

%%
temp = sumLsq - min(sumLsq(:));
ind = find(temp == 0);
[I,J,K,L] = ind2sub(size(sumLsq), ind);
f_op = normpdf(x, mu(J), sigma(K));

figure(); 
plot(x, f_op); 
hold on; 
plot(x, count_remote_normalize); 
plot(x, lambda2(L)*f_op+lambda(I)*count_remote_normalize);
plot(x, count_normalize)
