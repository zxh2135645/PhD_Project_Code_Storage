function [lr,sms] = autolambdaTensor(navData_bloch,ts_proj,msdev,flagCommandLine)

if nargin < 4
    flagCommandLine = true;
end

if ~flagCommandLine
    % Progress bar
    progress = waitbar(0,'','Name', sprintf('Progress'),...
                'CreateCancelBtn',...
                'setappdata(gcbf,''canceling'',1)');
    setappdata(progress,'canceling',0)
    waitbar(0, progress, sprintf('Estimating lambda...'));
end

fprintf('Estimating lambda... ')

try
    
    %% Only trying first unfolding, but others are relevant. Maybe eventually analyze all unfoldings (not counting 2nd Bloch-modeled unfolding)?
    temp1 = svde(navData_bloch(:,:));    % first unfolding
  
    %figure,plot(temp);
    [temp2, temp]  = kmeans(temp1,2);                   % k-means of zerofilled
    laplace_b = mean(temp)/(lambertw(-2*exp(-2)) + 2);  % estimate laplacian parameter b from k-means threshold
    lr = 2*msdev^2/laplace_b;                           % calculate regularization parameter from k-space noise variance and laplacian b

    temp2 = sqrt(mean(temp1.^2));
    laplace_b = temp2/(lambertw(-2*exp(-2)) + 2); 
    lr = 2*msdev^2/laplace_b; 

    fprintf('lr done. ')
    if ~flagCommandLine
        waitbar(0.2, progress);
    end
    
    %% Calculate sms in dimension 3~5 of the tensor

    dimIdx = 3:5;
    for n = numel(dimIdx):-1:1
        if size(navData_bloch,dimIdx(n)) == 1
            dimIdx(n) = [];
        end
    end
    if numel(dimIdx) > 0
        sms = zeros(1,numel(dimIdx));
        if ~isempty(dimIdx)
            for n = 1:numel(dimIdx)
                temp1 = diff(navData_bloch,1,dimIdx(n));
                temp1 = conj(ts_proj)*temp1(:,:);
                temp1 = abs(temp1(:));
                % figure;plot(sort(temp1(:)),'.');
                [a,temp]  = kmeans(temp1,2,'Start','cluster');
                laplace_b = mean(temp)/(lambertw(-2*exp(-2)) + 2);
                % laplace_b = sqrt(mean(temp.^2))/(lambertw(-2*exp(-2)) + 2);
                sms(n) = 2*msdev^2/laplace_b;
                
                [~,idx] = max(temp);
                sms(n) = min(temp1(a==idx));
                fprintf('sms (%d/%d) done. ',n,numel(dimIdx));
                if ~flagCommandLine
                    waitbar(0.6, progress);
                end
            end
        end
    else
        sms = [0 0];
    end
catch
    fprintf('lambda estimation failed. \n')
end

if ~flagCommandLine
    delete(progress); 
end

% other future improvements: deterministic initialization of k-means? Faster computation?


%% this is where the conversion from k-means threshold to b comes from:
% syms x a positive
% 
% f1=exp(-x)/int(exp(-x),x,0,a)
% f2=exp(-x)/int(exp(-x),x,a,inf)
% 
% EX_1 = int(f1*x,x,0,a)
% EX_2 = int(f2*x,x,a,inf)
% 
% EX2_1 = int(f1*x^2,x,0,a)
% EX2_2 = int(f2*x^2,x,a,inf)
% 
% var1 = EX2_1-EX_1^2
% var2 = EX2_2-EX_2^2
% 
% cost = var1*int(exp(-x),x,0,a) + var2*int(exp(-x),x,a,inf)
% 
% dcost=diff(simplify(cost),a)
% solve(dcost==0,a)
% 
% costfun=matlabFunction(simplify(cost))
% figure,plot(linspace(0,2,101),costfun(linspace(0,2,101)))



