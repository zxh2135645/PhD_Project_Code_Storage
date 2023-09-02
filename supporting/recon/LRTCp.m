function recovered=LRTCp(undersampled,lambda,mask,iterations,init_guess,curveU)

lowmem=true;

orig_size=size(undersampled);
undersampled=squeeze(undersampled);
mask=squeeze(mask);

if nargin < 2
  lambda=norm(undersampled(:))*1e-4
elseif isempty(lambda)
  lambda=norm(undersampled(:))*1e-4
end
if nargin >= 3
  if isempty(mask)
    if sum(isnan(undersampled))>0
      mask = ~isnan(undersampled);
    else
      mask = (undersampled~=0);
    end
  end
else
  if sum(isnan(undersampled))>0
    mask = ~isnan(undersampled);
  else
    mask = (undersampled~=0);
  end
end
if nargin < 4
  iterations=25;
end
if nargin > 5
  fixed_flag = ~isempty(curveU);
  if fixed_flag
    [curveU,~,~]=svd(curveU,'econ'); %orthonormalize
  end
else
  fixed_flag = false;
end

sizes=size(undersampled);
dims=numel(sizes);

if nargin < 5
  command = 'ndgrid(1:sizes(1)';
  outputs2 = '[xq1';
  for j=2:dims
    command = strcat(command,sprintf(',1:sizes(%d)',j));
    outputs2 = strcat(outputs2,sprintf(',xq%d',j));
  end
  eval(strcat(outputs2,'] = ',command,');'));
  
  outputs1='x1';
  x1=xq1(logical(mask));
  for j=2:dims
    outputs1 = strcat(outputs1,sprintf(',x%d',j));
    eval(sprintf('x%d=xq%d(logical(mask));',j,j))
  end
  
  disp('Interpolating missing data for initial guess...')
  tic; eval(strcat('recovered=griddata(',outputs1,',undersampled(logical(mask)),',outputs2(2:end),',''nearest'');')); toc;
else
  recovered=init_guess;
  clear init_guess;
end
recovered=recovered(:);

undersampled(isnan(undersampled))=0; %in case of nans
recovered(isnan(recovered))=0;

tensor=@(x)reshape(x,sizes);
collapse=@(x,dim)reshape(permute(tensor(x),[1:(dim-1), (dim+1):dims, dim]),[],sizes(dim));
icollapse=@(x,dim)ipermute(reshape(x,sizes([1:(dim-1), (dim+1):dims, dim])),[1:(dim-1), (dim+1):dims, dim]);
vec=@(x)x(:);

mask2=mask(:).^2;
Ahb=mask2.*undersampled(:);
for j=1:dims
  So{j}=svde(collapse(recovered,j));
  secondSVs(j)=So{j}(2);
  if ~lowmem %himem
    Z{j}=zeros(size(undersampled(:)));
  end
  Y{j}=zeros(size(undersampled(:)));
end
alpha=secondSVs;
if fixed_flag
  alpha(2) = alpha(2)/1e2; %1e3;
end

rho=lambda./alpha;
rho(isinf(rho) | isnan(rho))=0;

drho=1./rho;
drho(isinf(drho) | isnan(drho))=0;

step=zeros(1,dims);
for i=1:(iterations+1)
  dcp=norm((recovered(:)-undersampled(:)).*mask(:))^2;
  nn=0;
  for j=1:dims
    fprintf('SVD #%d ... ',j)
    tic;
    Sx{j}=svde(collapse(recovered,j));
    nn=nn+sum(Sx{j}(2:end));
    
    if (j == 2) && fixed_flag
      Znew{j}=vec(icollapse((collapse(recovered+Y{j}*drho(j),j)*curveU)*curveU',j));
    elseif lambda~=0
%       [U,S,V]=svd(collapse(recovered+Y{j}*drho(j),j),'econ');
      [U,S,V]=svde(collapse(recovered+Y{j}*drho(j),j));
      S=diag(S);
      Znew{j}=max(S-alpha(j),0);
      Znew{j}(1)=S(1);
      Znew{j}=vec(icollapse(U*diag(Znew{j})*V',j));
    else
      Znew{j}=Z{j};
    end
    
    if lowmem
      Z{j}=Znew{j};
      Znew{j}=[];
      pr=1;
      dr = pr/1.25;
    else %himem
      dr=rho(j)*norm(Z{j}-Znew{j});
      Z{j}=Znew{j};
      pr=norm(recovered-Z{j});
    end
    Y{j}=Y{j}+rho(j)*(recovered-Z{j}); %ADMM
    
    step(j)=dr/pr;
    step(j)=median([1/1.25, step(j), 1]);

    fprintf('%0.0f seconds\n',toc)
  end
  step(isnan(step))=1/1.25;

  [dcp, lambda*nn, dcp+lambda*nn]
% dcp
  
  alpha=alpha.*step;
  rho=lambda./alpha;
  rho(isinf(rho) | isnan(rho))=0;
  
  drho=1./rho;
  drho(isinf(drho) | isnan(drho))=0;
  
%   for j=1:dims
%     figure(j),plot(log10([So{j}, Sx{j}, svde(collapse(Z{j}-Y{j}*drho(j),j))]),'linewidth',2),hold on
%     plot(xlim,log10(alpha(j)/step(j))*[1 1],'k--','linewidth',2),hold off,drawnow
%     %figure(j),semilogy([So{j}, S{j}, Z-Y{j}*drho(j)],'linewidth',2),drawnow
%   end
  
  if i==(iterations+1)
    break
  end
  
  Ahbtemp=Ahb;
  for j=1:dims
    Ahbtemp=Ahbtemp+rho(j)/2*vec(Z{j}-Y{j}*drho(j));
    if lowmem
      Z{j}=[];
    end
  end
  recovered = Ahbtemp./(mask2+sum(rho)/2);
  clear Ahbtemp; 
  
  % profile viewer
  % profile resume
end
recovered=reshape(recovered,orig_size);

return

function xl=signedlog10(x)
xl=zeros(size(x));
xl(x>0)=log10(x(x>0));
xl(x==0)=min(xl(x>0));
xl(x<0)=min(xl(x>0))-max(-log10(-x(x<0)))-log10(-x(x<0));
return