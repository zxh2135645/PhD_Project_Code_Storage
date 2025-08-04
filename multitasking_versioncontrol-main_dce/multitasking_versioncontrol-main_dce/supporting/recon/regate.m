function [X,mask] = regate(X_rt,Segidx,Hidx,Ridx,wallClock,echoIdx)

if nargin < 6
    echoIdx = ones(size(Segidx));
end
if nargin < 5
    wallClock = ones(size(Segidx));
end
if nargin < 4
    Ridx = ones(size(Segidx));
end
if nargin < 4
    Hidx = ones(size(Segidx));
end

T = size(X_rt,1);
segments = max(Segidx);
Hframes  = max(Hidx);
Rframes  = max(Ridx);
Necho    = max(echoIdx);

X    = zeros(size(X_rt,2),segments,Hframes,Rframes,max(wallClock),Necho);
mask = zeros(size(X_rt,2),segments,Hframes,Rframes,max(wallClock),Necho);
for t = 1:T
    if Ridx(t) > 0 && wallClock(t) > 0
        X(:,Segidx(t),Hidx(t),Ridx(t),wallClock(t),echoIdx(t))    = X(:,Segidx(t),Hidx(t),Ridx(t),wallClock(t),echoIdx(t)) + X_rt(t,:).';
        mask(:,Segidx(t),Hidx(t),Ridx(t),wallClock(t),echoIdx(t)) = mask(:,Segidx(t),Hidx(t),Ridx(t),wallClock(t),echoIdx(t)) + 1;
    end
end
if max(mask(:))>1
    X = X./mask;
    X(isnan(X)) = 0;
    mask = sqrt(mask);
else
    mask = logical(mask);
end

