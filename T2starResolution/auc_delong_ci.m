function [auc, ci95, se] = auc_delong_ci(y, s, alpha)
% auc_delong_ci  DeLong (1988) 95% CI for ROC AUC.
% y: binary labels (1=positive, 0=negative) as Nx1
% s: scores (higher = more positive) as Nx1
% alpha: e.g., 0.05 for 95% CI (optional, default 0.05)

if nargin < 3, alpha = 0.05; end

y = y(:); s = s(:);
assert(numel(y)==numel(s), 'y and s must have same length.');
assert(all(ismember(unique(y), [0 1])), 'y must be binary 0/1.');

pos = s(y==1);
neg = s(y==0);
m = numel(pos);
n = numel(neg);
assert(m>0 && n>0, 'Need at least one positive and one negative.');

% Midrank-style kernel: psi(a,b) = 1 if a>b, 0.5 if a==b, 0 if a<b
% Compute pairwise comparisons efficiently
P = pos(:);
N = neg(:)';

cmp = (P > N) + 0.5*(P == N);  % m x n matrix

% AUC
auc = mean(cmp(:));

% DeLong components
V10 = mean(cmp, 2);  % per-positive average over negatives (m x 1)
V01 = mean(cmp, 1)'; % per-negative average over positives (n x 1)

% Sample variances
s10 = var(V10, 1);   % population variance (divide by m)
s01 = var(V01, 1);   % population variance (divide by n)

% AUC variance estimate
var_auc = s10/m + s01/n;
se = sqrt(var_auc);

% Normal approximation CI (standard with DeLong SE)
z = norminv(1 - alpha/2);
ci95 = [auc - z*se, auc + z*se];
ci95 = max(min(ci95, 1), 0);   % clamp to [0,1]
end