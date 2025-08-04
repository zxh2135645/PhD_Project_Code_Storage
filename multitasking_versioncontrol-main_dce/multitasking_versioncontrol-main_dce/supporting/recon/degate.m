function X_rt = degate(X,sizes,Segidx,Hidx,Ridx,wallClock,echoIdx)

if nargin < 7
    echoIdx = ones(size(Segidx));
end

T    = numel(Segidx);
dims = numel(sizes);

collapse  = @(x,dim) reshape(permute(reshape(x,sizes),[1:(dim-1), (dim+1):dims, dim]),[],sizes(dim));
icollapse = @(x,dim) ipermute(reshape(x,sizes([1:(dim-1), (dim+1):dims, dim])),[1:(dim-1), (dim+1):dims, dim]);
vec = @(x) x(:);

X    = icollapse(X,1);
X_rt = complex(zeros(sizes(1),T));
for t = 1:T
    if Ridx(t) > 0 
        X_rt(:,t) = X(:,Segidx(t),Hidx(t),Ridx(t),wallClock(t),echoIdx(t));
    else
        Ridx(t) = Ridx(t-1);
        X_rt(:,t) = X(:,Segidx(t),Hidx(t),Ridx(t),wallClock(t),echoIdx(t));
    end
end

return