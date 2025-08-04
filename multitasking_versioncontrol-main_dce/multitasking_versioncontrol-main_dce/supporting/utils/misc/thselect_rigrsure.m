function thr = thselect_rigrsure(x)

if isrow(x)
    x = x(:);
end
[n,m] = size(x);

sx = sort(abs(x),1);
sx2 = sx.^2;
N1 = repmat((n-2*(1:n))',1,m);
N2 = repmat((n-1:-1:0)',1,m);
CS1 = cumsum(sx2,1);
%risks = (N1+CS1+N2.*sx2)./n;
risks = (N1+CS1)./n;
[~,best] = min(risks,[],1);
% thr will be row vector
thr = sx(best);


