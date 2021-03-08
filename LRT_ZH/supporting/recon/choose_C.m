function [C, UU, ranks] = choose_C(X,ranks,Xorig,mask)
% [C, UU] = choose_C(Navdata,ranks) (individual truncation)
% [C, UU] = choose_C(Navdata,rank) (truncate only C)

orig_flag = (nargin > 2);

if numel(ranks)==1
  rank=ranks;
  ranks=size(X);
  ranks(1)=rank;
end

sizes=size(X);
dims=numel(sizes);

collapse=@(x,dim)reshape(permute(reshape(x,sizes),[1:(dim-1), (dim+1):dims, dim]),[],sizes(dim));
icollapse=@(x,dim)ipermute(reshape(x,sizes([1:(dim-1), (dim+1):dims, dim])),[1:(dim-1), (dim+1):dims, dim]);
vec=@(x)x(:);

UU=1;
for j=(dims):-1:2
  [~,S,V]=svd(collapse(X,j),'econ');
  figure,plot(20*log10(diag(S)/S(1)),'.-'),grid on
%   temp=input(sprintf('Rank [%d]: ',ranks(j)));
%   if ~isempty(temp)
%     ranks(j)=temp;
%   end
  V=V(:,1:ranks(j));
  UU=kron(UU,V);
end

UC=collapse(X,1).'*UU;
if orig_flag
  mat = @(x)reshape(x,sizes(1),[]);

  M2 = collapse(mask.^2,1).';

  AhA = @(UC)vec(((mat(UC)*UU').*M2)*UU);
  Ahb = vec((collapse(Xorig,1).'.*M2)*UU);
  UC = mat(pcg(AhA,Ahb,[],50,[],[],UC(:)));
end
if size(UC,1)>size(UC,2)
  [~,C1,C2]=svd(UC,'econ');
else
  [C2,C1,]=svd(UC','econ');
end
C=C1*C2';

curve=sqrt(real(diag(C*C')));
figure,plot(20*log10(curve(1:min(200,end))/curve(1)),'.-'),grid on
ranks(1)=min(ranks(1),prod(ranks(2:end)));
temp=input(sprintf('Rank [%d]: ',ranks(1)));
% temp = 32;

if ~isempty(temp)
  ranks(1)=temp;
end

C=C(1:ranks(1),:);

return
