function [X,mask] = regate(X_rt,Segidx,Hidx,Ridx,wall_clock)

N=numel(Segidx);
segments=max(Segidx);
Hframes=max(Hidx);
Rframes=max(Ridx);

X=zeros(size(X_rt,2),segments,Hframes,Rframes,max(wall_clock));
mask=zeros(size(X_rt,2),segments,Hframes,Rframes,max(wall_clock));
for t=1:N
  X(:,Segidx(t),Hidx(t),Ridx(t),wall_clock(t))=X(:,Segidx(t),Hidx(t),Ridx(t),wall_clock(t))+X_rt(t,:).';
  mask(:,Segidx(t),Hidx(t),Ridx(t),wall_clock(t))=mask(:,Segidx(t),Hidx(t),Ridx(t),wall_clock(t))+1;
end
if max(mask(:))>1
  X=X./mask;
  X(isnan(X))=0;
else
  mask=logical(mask);
end

return