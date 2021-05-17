function X_rt = degate_5D(X,sizes,Segidx,Hidx,Ridx,wall_clock,seg_multiplier)

T=numel(Hidx);
dims=numel(sizes);

collapse=@(x,dim)reshape(permute(reshape(x,sizes),[1:(dim-1), (dim+1):dims, dim]),[],sizes(dim));
icollapse=@(x,dim)ipermute(reshape(x,sizes([1:(dim-1), (dim+1):dims, dim])),[1:(dim-1), (dim+1):dims, dim]);
vec=@(x)x(:);

X = icollapse(X,1);
X_rt=complex(zeros(sizes(1),T));
for t=1:T
   X_rt(:,t)=X(:,Segidx(t),Hidx(t),Ridx(t),wall_clock(t),seg_multiplier(t));
end

return