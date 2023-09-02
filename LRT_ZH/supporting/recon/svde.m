function [U,S,V] = svde(X)
% misses one side of nullspace

[M, N] = size(X);

if M <= N %if square or wide
  [U,Lambda] = eig(X*X','vector');
else
  [V,Lambda] = eig(X'*X,'vector');
end

[S, idx] = sort(sqrt(abs(Lambda)),'descend');

if nargout == 1
  U = S; %actually S
else
  if M <= N %if square or wide
    U = U(:,idx);
    V = X'*(U.*repmat(1./reshape(S,1,[]),[M 1]));
    V(isnan(V))=0; %missing nullspace for now
  else
    V = V(:,idx);
    U = X*(V.*repmat(1./reshape(S,1,[]),[N 1]));
    U(isnan(U))=0; %missing nullspace
  end
  S = diag(S);
end