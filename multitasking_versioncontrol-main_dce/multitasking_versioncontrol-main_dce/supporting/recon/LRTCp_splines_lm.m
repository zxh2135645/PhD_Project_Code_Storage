function recovered = LRTCp_splines_lm(undersampled,lambda,mus,mask,iterations,init_guess,curveU,orders,circs,ts_proj,flagCommandLine)

lowmem = true;
function_method = true;

undersampled = double(squeeze(undersampled));
mask = squeeze(mask);

if nargin < 2
    lambda = norm(undersampled(:))*1e-4;
elseif isempty(lambda)
    lambda = norm(undersampled(:))*1e-4;
end
if nargin >= 4
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
sizes_undersampled = size(undersampled);
sizes = sizes_undersampled;
dims  = numel(sizes_undersampled);

collapse  = @(x,dim,sizes) reshape(permute(reshape(x,sizes),[1:(dim-1), (dim+1):dims, dim]),[],sizes(dim));
icollapse = @(x,dim,sizes) ipermute(reshape(x,sizes([1:(dim-1), (dim+1):dims, dim])),[1:(dim-1), (dim+1):dims, dim]);
vec = @(x) x(:);

if nargin < 5
    iterations = 25;
end

if nargin > 6
    fixed_flag = ~isempty(curveU);
    if fixed_flag
        curveU = double(curveU);
        curveU_orig = curveU;
        [curveU,~,~] = svd(curveU_orig,'econ'); %orthonormalize
        cL = size(curveU,2);
        sizes(2) = cL;
    end
else
    fixed_flag = false;
end

if nargin < 6
    command = 'ndgrid(1:sizes(1)';
    outputs2 = '[xq1';
    for j = 2:dims
        command = strcat(command,sprintf(',1:sizes(%d)',j));
        outputs2 = strcat(outputs2,sprintf(',xq%d',j));
    end
    eval(strcat(outputs2,'] = ',command,');'));

    outputs1='x1';
    x1 = xq1(logical(mask));
    for j = 2:dims
        outputs1 = strcat(outputs1,sprintf(',x%d',j));
        eval(sprintf('x%d=xq%d(logical(mask));',j,j))
    end

    disp('Interpolating missing data for initial guess...')
    tic; eval(strcat('recovered=griddata(',outputs1,',undersampled(logical(mask)),',outputs2(2:end),',''nearest'');')); toc;
else
    recovered = icollapse(collapse(init_guess,2,sizes) * curveU_orig.' * pinv(curveU.'),2,sizes);
    clear init_guess;
end
recovered = double(recovered(:));

if nargin < 8
    orders = ones(size(mus)); %3 for piecewise quadratic
end

if nargin < 9
    circs = false(size(mus));
end

sizes_big = sizes;
if nargin < 10
    ts_proj = 1;
elseif isempty(ts_proj)
    ts_proj = 1;
else
    function_method = true; %force function_method
    sizes_big(1) = size(ts_proj,1);
    ts_proj = double(ts_proj);
end

if nargin < 11
    flagCommandLine = true;
end

smdims = (dims-numel(mus)+1):dims;

undersampled(isnan(undersampled)) = 0; %in case of nans
recovered(isnan(recovered)) = 0;

if ~flagCommandLine
    % Progress bar
    progress = waitbar(0,'','Name', sprintf('Progress'),...
            'CreateCancelBtn',...
            'setappdata(gcbf,''canceling'',1)');
    setappdata(progress,'canceling',0)
    waitbar(0, progress, sprintf('Calculating tensor subspace...'));
end

try
    
%% build gradient functions
if function_method
    for dim = 1:numel(smdims)
        G{dim} = speye(sizes(smdims(dim)));
        if circs(dim)
            G{dim} = G{dim}-G{dim}(:,[end, 1:end-1]);
            G{dim} = G{dim}^orders(dim);
        else
            for j = 1:orders(dim)
                G{dim} = G{dim}(:,2:end)-G{dim}(:,1:end-1);
            end
        end
        GGh{dim} = G{dim}*G{dim}.';
    end
else
    for dim = 1:numel(smdims)
        G{dim} = speye(prod(sizes));
        for j = 1:orders(dim)
            sz1 = sizes;
            sz1(smdims(dim)) = 1;
            sz2 = sizes;
            sz2(smdims(dim)) = sz2(smdims(dim))-j;
            indices1 = cat(smdims(dim),false(sz1),true(sz2));
            indices2 = cat(smdims(dim),true(sz2),false(sz1));
            if circs(dim)
                G{dim} = cat(1,G{dim}(indices1(:),:)-G{dim}(indices2(:),:),G{dim}(~indices1(:),:)-G{dim}(~indices2(:),:));
            else
                G{dim} = G{dim}(indices1(:),:)-G{dim}(indices2(:),:);
            end
        end
        GhG{dim} = G{dim}.'*G{dim};
    end
end
%%
mask2 = mask(:).^2;
Ahb   = icollapse(collapse(mask2 .* undersampled(:),2,sizes_undersampled) * conj(curveU),2,sizes);
Ahb   = Ahb(:);
for j = 1:dims
    So{j} = svde(collapse(recovered,j,sizes)); %undersampled?
    secondSVs(j) = So{j}(2);
    if ~lowmem %himem
        Z{j} = zeros(size(recovered(:)));
    end
    Y{j} = zeros(size(recovered(:)));
end
alpha = secondSVs;
if fixed_flag
    alpha(2) = alpha(2)/1e2; %1e3;
end
% alpha([1 3:end])=inf; %TEMPORARILY cancel other LR dimensions

rho = lambda./alpha;
rho(isinf(rho) | isnan(rho)) = 0;

drho = 1./rho;
drho(isinf(drho) | isnan(drho)) = 0;

for dim = 1:numel(smdims)
    if function_method
        Gx{dim} = collapse(double(recovered),smdims(dim),sizes)*G{dim};
        Gx{dim} = reshape(conj(ts_proj)*reshape(Gx{dim},sizes(1),[]),[],size(G{dim},2));
    else
        Gx{dim} = G{dim}*recovered(:);
    end
  
    Z2{dim} = zeros(size(Gx{dim}));
    Y2{dim} = zeros(size(Gx{dim}));
    alpha2{dim} = max(abs(Gx{dim}(:)));
    rho2{dim} = mus(dim)/alpha2{dim};
    rho2{dim}(isinf(rho2{dim}) | isnan(rho2{dim})) = 0;

    drho2{dim} = 1/rho2{dim};
    drho2{dim}(isinf(drho2{dim}) | isnan(drho2{dim})) = 0;
end

% [dcp, mus*sum(abs(Gx(:))), lambda*nn, dcp+lambda*nn] %will show later

% profile on
fprintf('\n%d iterations: ',iterations);
for i = 1:(iterations+1)
    for j = 1:dims
        tic;
    
%         if lambda ~= 0
%             Sx{j} = svde(collapse(recovered,j));
%         end
    
        if lambda ~= 0
      %       [U,S,V]=svd(collapse(recovered+Y{j}*drho(j),j),'econ');
            [U,S,V] = svde(collapse(recovered+Y{j}*drho(j),j,sizes));
            S = diag(S);
            Znew{j} = max(S-alpha(j),0);
            Znew{j}(1) = S(1);
            Znew{j} = vec(icollapse(U*diag(Znew{j})*V',j,sizes));
        else
            Znew{j} = 0; %Z{j};
        end
    
        if lowmem
            Z{j} = Znew{j};
            Znew{j} = [];
            pr = 1;
            dr = pr/1.25;
        else %himem
            dr = rho(j) * norm(Z{j}-Znew{j});
            Z{j} = Znew{j};
            pr = norm(recovered-Z{j});
        end
        Y{j} = Y{j} + rho(j) * (recovered-Z{j});    %ADMM
    
        step(j) = dr/pr;
        step(j) = median([1/1.25, step(j), 1]);
    end
    step(isnan(step)) = 1/1.25;
  
    tic;
%     if debug_flag
%         d2 = 0;
%     end
    for dim = 1:numel(smdims)
        if function_method
            Gx{dim} = collapse(double(recovered),smdims(dim),sizes) * G{dim};
            Gx{dim} = reshape(conj(ts_proj) * reshape(Gx{dim},sizes(1),[]),[],size(G{dim},2));
        else
            Gx{dim} = G{dim} * recovered(:);
        end
        Z2new{dim} = sign(Gx{dim}+Y2{dim}*drho2{dim}).*max(abs(Gx{dim}+Y2{dim}*drho2{dim})-alpha2{dim},0);
        if function_method
            dr2{dim} = rho2{dim}*norm(vec(reshape(ts_proj.'*reshape(Z2{dim}-Z2new{dim},sizes_big(1),[]),[],size(G{dim},2))*G{dim}.'));
        else
            dr2{dim}=rho2{dim}*norm((Z2{dim}-Z2new{dim}).'*G{dim});
        end
        Z2{dim}  = Z2new{dim};
        pr2{dim} = norm(Gx{dim}-Z2{dim},'fro');
        Y2{dim}  = Y2{dim} + rho2{dim}*(Gx{dim}-Z2{dim}); %ADMM
%         if debug_flag
%             d2 = d2+mus(dim)*sum(abs(Gx{dim}(:)));
%         end
    end
%     if debug_flag
%         fprintf('%0.0f seconds\n',toc)
%         [dcp, d2, lambda*nn, dcp+d2+lambda*nn]
%     end
  
    for dim=1:numel(smdims)
        step2{dim} = dr2{dim}/pr2{dim};
        step2{dim} = median([1/1.25, step2{dim}, 1]);
        step2{dim}(isnan(step2{dim})) = 1/1.25;
    end
  
    alpha = alpha.*step;
    rho   = lambda./alpha;
    rho(isinf(rho) | isnan(rho)) = 0;

    drho = 1./rho;
    drho(isinf(drho) | isnan(drho)) = 0;

    for dim = 1:numel(smdims)
        alpha2{dim} = alpha2{dim} * step2{dim};
        rho2{dim}   = mus(dim)/alpha2{dim};
        rho2{dim}(isinf(rho2{dim}) | isnan(rho2{dim})) = 0;
        drho2{dim} = 1/rho2{dim};
        drho2{dim}(isinf(drho2{dim}) | isnan(drho2{dim})) = 0;
    end
  %   for j=1:dims
  %     figure(j),plot(log10([So{j}, Sx{j}, svde(collapse(Z{j}-Y{j}*drho(j),j))]),'linewidth',2),hold on
  %     plot(xlim,log10(alpha(j)/step(j))*[1 1],'k--','linewidth',2),hold off,drawnow
  %     %figure(j),semilogy([So{j}, S{j}, Z-Y{j}*drho(j)],'linewidth',2),drawnow
  %   end
  
    if i == (iterations+1)
        break
    end
    fprintf('\n -- %d/%d ...', i, iterations);
    
    Ahbtemp = double(Ahb);
    if function_method
        newbench  = 0;
        benchmark = inf;
        AhA = @(x) double(vec(icollapse((collapse(x,2,sizes)*curveU.' .* collapse(mask2,2,sizes_undersampled)) * conj(curveU),2,sizes)) + sum(rho)/2.*x(:));
        for dim = 1:numel(smdims)
            Ahbtemp = Ahbtemp + rho2{dim}/2*vec(icollapse(reshape(ts_proj.'*reshape(Z2{dim}-Y2{dim}*drho2{dim},sizes_big(1),[]),[],size(G{dim},2))*G{dim}.',smdims(dim),sizes));
            AhA = @(x) double((AhA(x) + rho2{dim}/2*vec(icollapse(collapse(x,smdims(dim),sizes)*GGh{dim},smdims(dim),sizes))));
        end
    else
        newbench  = inf;
        benchmark = 0;
        AhA = spdiags(mask2(:)+sum(rho)/2,0,numel(mask),numel(mask));
        for dim = 1:numel(smdims)
            Ahbtemp = Ahbtemp+rho2{dim}/2*((Z2{dim}-Y2{dim}*drho2{dim}).'*G{dim}).';
            AhA = AhA + rho2{dim}/2 * GhG{dim};
        end
    end
    
    for j = 1:dims
        Ahbtemp = Ahbtemp + rho(j)/2 * vec(Z{j}-Y{j}*drho(j));
        if lowmem
            Z{j} = [];
        end
    end
  
    tic;
    if (mod(i,2)==0) || (newbench<benchmark) %every other, try PCG
        [recovered,flag] = pcg2(AhA,double(Ahbtemp),[],median([5 i 10]),[],[],double(recovered));
        newbench  = toc;
    else
        try
          recovered = AhA\Ahbtemp; %sparse solver
          benchmark = toc;
        catch
          [recovered,flag] = pcg2(AhA,Ahbtemp,[],median([5 i 10]),[],[],recovered);
          benchmark = inf;
          newbench  = toc;
        end
    end
  
    if ~flagCommandLine
        if getappdata(progress,'canceling')
            break
        else
            waitbar(i/(iterations+1), progress);
        end
    end
end

catch errormsg
    if isfield(errormsg,'message')
        fprintf(2, '\n%s\n', errormsg.message);
    end
    fprintf(2, 'tensor subspace calculation: LRTCp_splines_lm failed. \n');
end

recovered = icollapse(collapse(recovered,2,sizes) * (curveU.' * pinv(curveU_orig.')),2,sizes);

if ~flagCommandLine
    delete(progress);
end


% return

% function xl=signedlog10(x)
% xl=zeros(size(x));
% xl(x>0)=log10(x(x>0));
% xl(x==0)=min(xl(x>0));
% xl(x<0)=min(xl(x>0))-max(-log10(-x(x<0)))-log10(-x(x<0));
% return


