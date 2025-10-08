function model = fit_gamma_like(x, y, doPlot)
%FIT_GAMMA_LIKE Fit a gamma-shaped curve (not necessarily Gamma).
%   model = fit_gamma_like(x, y)
%   - Tries scaled Gamma, Weibull, Lognormal pdf-like shapes.
%   - Chooses best by AIC; also returns a robust LOESS smooth.
%   - Returns struct with best model, parameters, AIC, predict function.
%
%   model fields: .bestName, .params, .predict, .AIC, .SSE, .n, .k
%                 .all (table of candidates), .xgrid, .yhat_best, .yhat_loess
%
%   Example:
%       x = linspace(0.05,10,60)'; y = 5*gampdf(x,3,1.2) + 0.02*randn(size(x));
%       M = fit_gamma_like(x,y,true);

if nargin < 3, doPlot = true; end

x = x(:); y = y(:);
% Keep finite, positive x (for lognormal/gamma-like shapes)
valid = isfinite(x) & isfinite(y) & x > 0;
x = x(valid); y = y(valid);

% Sort by x
[x, ix] = sort(x); y = y(ix);

n = numel(x);

% Simple robust scaling for numeric stability (optional)
xscale = median(x);
if xscale <= 0, xscale = 1; end
xs = x / xscale;

% -------- Candidate model definitions (scaled PDFs) ----------
% y = A * pdf_like(x; params)
% Params are constrained positive where appropriate.

% 1) Scaled Gamma-like: y = A * x^(k-1) * exp(-x/theta) / (theta^k * gamma(k))
g_fun = @(p, xx) max(0, p(1) * (xx.^(p(2)-1) .* exp(-xx./p(3)) ./ (p(3).^p(2) .* gamma_safe(p(2)))));
% p = [A, k, theta], k>0, theta>0

% 2) Scaled Weibull pdf: y = A * (k/lambda) * (x/lambda)^(k-1) * exp(-(x/lambda)^k)
w_fun = @(p, xx) max(0, p(1) * (p(2)./p(3)) .* (xx./p(3)).^(p(2)-1) .* exp(- (xx./p(3)).^p(2)));
% p = [A, k, lambda], k>0, lambda>0

% 3) Scaled Lognormal pdf: y = A * [1/(x*sigma*sqrt(2*pi))] * exp(-(ln x - mu)^2/(2 sigma^2))
ln_fun = @(p, xx) max(0, p(1) * (1./(xx * p(3) * sqrt(2*pi))) .* exp(-(log(xx)-p(2)).^2/(2*p(3)^2)));
% p = [A, mu, sigma], sigma>0

% Initial guesses (from data moments / peak)
[~, imax] = max(y);
xmode = x(imax);
ymode = y(imax);
xmean = sum(x.*max(y,0))/sum(max(y,0));
if ~isfinite(xmean), xmean = median(x); end

% Gamma-ish guesses
A0_g   = max(ymode, eps);
k0_g   = max( (xmean/xmode), 1.2 );     % crude
theta0 = max(xmean / max(k0_g,1e-3), 1e-3);

% Weibull guesses
A0_w     = max(ymode, eps);
k0_w     = 1.5;                 % mild skew
lambda0  = max(xmode, 1e-3);

% Lognormal guesses
A0_ln   = max(ymode, eps);
mu0_ln  = log(max(xmode, 1e-3));
sig0_ln = 0.5;

% Bounds
LB_g  = [0,  1e-3, 1e-6];
UB_g  = [Inf, 50,   Inf];
LB_w  = [0,  1e-3, 1e-6];
UB_w  = [Inf, 50,   Inf];
LB_ln = [0, -Inf,  1e-6];
UB_ln = [Inf, Inf,  Inf];

% Use lsqcurvefit if available; otherwise fminsearch wrapper
useLS = exist('lsqcurvefit','file') == 2;

opts = optimset('Display','off');
if useLS
    opts2 = optimoptions('lsqcurvefit','Display','off');
end

% Rescaled-domain wrappers (improve conditioning)
wrap = @(f) @(p) f(p, xs) - y;

% Fit Gamma 
% ---------- Fit Gamma ----------
res_g = @(p) g_fun(p, xs) - y;                 % residuals in scaled x
p0_g  = [A0_g, k0_g, theta0];
if useLS
    [pg,~,~] = lsqcurvefit(@(p,xx) g_fun(p,xx), p0_g, xs, y, LB_g, UB_g, opts2);
else
    pg = fminsearch(@(p) sum(res_g(p).^2), p0_g, opts);
    pg = max(pg, LB_g); pg = min(pg, UB_g);
end
SSEg = sum((y - g_fun(pg, xs)).^2); kg = 3; AICg = n*log(SSEg/n) + 2*kg;

% ---------- Fit Weibull (lambda fitted in scaled domain, then unscaled) ----------
LB_w_s = [LB_w(1), LB_w(2), LB_w(3)/xscale];
UB_w_s = [UB_w(1), UB_w(2), UB_w(3)/xscale];
p0_w_s = [A0_w, k0_w, max(lambda0/xscale, 1e-6)];

res_w_s = @(p) w_fun([p(1), p(2), p(3)], xs) - y;

if useLS
    [pw_s,~,~] = lsqcurvefit(@(p,xx) w_fun([p(1),p(2),p(3)], xx), ...
                              p0_w_s, xs, y, LB_w_s, UB_w_s, opts2);
else
    pw_s = fminsearch(@(p) sum(res_w_s(p).^2), p0_w_s, opts);
    pw_s = max(pw_s, LB_w_s); pw_s = min(pw_s, UB_w_s);
end
pw = [pw_s(1), pw_s(2), pw_s(3)*xscale];   % unscale lambda
SSEw = sum((y - w_fun(pw, x)).^2); kw = 3; AICw = n*log(SSEw/n) + 2*kw;

% ---------- Fit Lognormal (mu fitted in log(xs), then shifted back) ----------
p0_ln_s = [A0_ln, mu0_ln - log(xscale), sig0_ln];
LB_ln_s = LB_ln;               % same bounds; sigma > 0 already ensured
UB_ln_s = UB_ln;

res_ln_s = @(p) ln_fun(p, xs) - y;

if useLS
    [pln_s,~,~] = lsqcurvefit(@(p,xx) ln_fun(p,xx), p0_ln_s, xs, y, LB_ln_s, UB_ln_s, opts2);
else
    pln_s = fminsearch(@(p) sum(res_ln_s(p).^2), p0_ln_s, opts);
    pln_s = max(pln_s, LB_ln_s); pln_s = min(pln_s, UB_ln_s);
end
pln = [pln_s(1), pln_s(2) + log(xscale), pln_s(3)];  % shift mu back
SSEln = sum((y - ln_fun(pln, x)).^2); kln = 3; AICln = n*log(SSEln/n) + 2*kln;


% Nonparametric robust LOESS (nice for visualization/fallback)
y_loess = smooth(x, y, 0.15, 'rloess');  % adjust span as needed

% Compare and pick best
cands = table(["Gamma";"Weibull";"Lognormal"], ...
              [AICg; AICw; AICln], ...
              [SSEg; SSEw; SSEln], ...
              'VariableNames', {'Model','AIC','SSE'});

[~, ibest] = min(cands.AIC);
bestName = cands.Model(ibest);

switch bestName
    case "Gamma"
        params = struct('A',pg(1),'k',pg(2),'theta',pg(3));
        predict = @(xx) g_fun([params.A, params.k, params.theta], xx);
        SSE=SSEg; kpar=3; AIC=AICg;
    case "Weibull"
        params = struct('A',pw(1),'k',pw(2),'lambda',pw(3));
        predict = @(xx) w_fun([params.A, params.k, params.lambda], xx);
        SSE=SSEw; kpar=3; AIC=AICw;
    case "Lognormal"
        params = struct('A',pln(1),'mu',pln(2),'sigma',pln(3));
        predict = @(xx) ln_fun([params.A, params.mu, params.sigma], xx);
        SSE=SSEln; kpar=3; AIC=AICln;
end

xgrid = linspace(min(x), max(x), 400)';
yhat = predict(xgrid);

model = struct('bestName',bestName,'params',params,'predict',predict, ...
               'AIC',AIC,'SSE',SSE,'n',n,'k',kpar, ...
               'all',cands,'xgrid',xgrid,'yhat_best',yhat,'yhat_loess',y_loess);

if doPlot
    figure; hold on
    scatter(x,y,18,'filled','MarkerFaceAlpha',0.35); 
    plot(xgrid, yhat, 'LineWidth',2);
    plot(x, y_loess, '--', 'LineWidth',1.5);
    legend('Data', sprintf('Best: %s', bestName), 'LOESS (fallback)','Location','best');
    xlabel('x'); ylabel('y'); title('Gamma-like Curve Fit');
    grid on
end
end

function g = gamma_safe(k)
% Safe gamma for vector k>0 (fallback to exp(gammaln))
g = exp(gammaln(k));
end