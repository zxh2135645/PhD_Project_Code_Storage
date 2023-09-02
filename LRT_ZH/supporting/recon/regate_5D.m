function [X,mask] = regate_5D(X_rt,Segidx,Hidx,Ridx,wall_clock,seg_multiplier)

N=numel(Segidx);
segments=max(Segidx);
Hframes=max(Hidx);
Rframes=max(Ridx);
Aframes=max(seg_multiplier);

X=zeros(size(X_rt,2),segments,Hframes,Rframes,max(wall_clock), Aframes);
mask=zeros(size(X_rt,2),segments,Hframes,Rframes,max(wall_clock), Aframes);
for t=1:N
  X(:,Segidx(t),Hidx(t),Ridx(t),wall_clock(t),seg_multiplier(t))=X(:,Segidx(t),Hidx(t),Ridx(t),wall_clock(t),seg_multiplier(t))+X_rt(t,:).';
  mask(:,Segidx(t),Hidx(t),Ridx(t),wall_clock(t),seg_multiplier(t))=mask(:,Segidx(t),Hidx(t),Ridx(t),wall_clock(t),seg_multiplier(t))+1;
end
if max(mask(:))>1
  X=X./mask;
  X(isnan(X))=0;
else
  mask=logical(mask);
end

return