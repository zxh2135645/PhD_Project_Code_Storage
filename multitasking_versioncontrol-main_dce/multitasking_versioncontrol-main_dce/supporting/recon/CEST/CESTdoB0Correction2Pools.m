function [imgCEST_corr, B0map, zSpectFitRef] = CESTdoB0Correction2Pools(imgCEST, df_ppm, mask ,x0, xlb, xub)
%   Do B0 correction with parfor loop
        
    size3 = @(x) [size(x,1),size(x,2),size(x,3)];

    if nargin < 3
        mask = ones(size3(imgCEST));
    end
    if nargout < 3
        flagSpectFitRef = false;
    else
        flagSpectFitRef = true;
    end

    mask(mean(imgCEST,3)>1) = 0;
    imgCEST = reshape(imgCEST,[],numel(df_ppm));
    imgCEST_corr = zeros(size(imgCEST));
    zSpectFitRef = zeros(size(imgCEST));
    B0map = zeros(size(imgCEST,1),1);
    zSpect  = imgCEST(mask(:),:);
    zSpect_corr = zeros(size(zSpect));
    zSpect_fit  = zeros(size(zSpect));
    tempB0map = zeros(size(zSpect,1),1);

    x0  = [x0(1:6),  0.05];
    xlb = [xlb(1:6), 0.0];
    xub = [xub(1:6), 0.5];

    zfun = @(x,xdata) (1 - x(1).*(x(3)^2/4)./((x(3)^2/4)+(xdata-x(2)).^2)).*(1 - x(4).*(x(6)^2/4)./((x(6)^2/4)+(xdata-x(5)).^2)) - x(7);
    
%     B0_fit_range = abs(df_ppm) <= 1 || abs(df_ppm) >= 6.5;
%     df_ppm_cen = df_ppm(B0_fit_range);
    df_ppm_orig = df_ppm(:).';
    df_ppm_full = logspace(log10(0.1),log10(max(abs(df_ppm))),100);
    df_ppm_full = unique([-df_ppm_full 0 df_ppm_full df_ppm_orig]);
    B0_fit_range = abs(df_ppm_full) <= 1 | abs(df_ppm_full) >= 6.5;
    df_ppm_cen = df_ppm_full(B0_fit_range);
    zSpect_temp = zeros(length(tempB0map),numel(df_ppm_full));

    fprintf('\nStart B0 correction: ')
    parfor j = 1:length(tempB0map)
        zSpecttemp  = interp1(df_ppm_orig,zSpect(j,:),df_ppm_full);

        % 2-pool fitting (water/MT)
        [deltaB0,~] = estimateB0(df_ppm_cen, zSpecttemp(B0_fit_range), x0, xlb, xub);
        
        % Save B0 map and correct B0
        tempB0map(j) = deltaB0;
        zSpect_temp(j,:) = interp1(df_ppm_orig-deltaB0, zSpect(j,:), df_ppm_full, 'spline');
    end

    % PCA denoising
    [components, score, latent] = pcaogc(zSpect_corr(:));        
    cutoff = 10;
    temp1 = mean(temp) + score(:,1:cutoff)*components(1:cutoff).';
    [~,idx] = intersect(df_ppm_full,df_ppm_orig);
    zSpect_corr = temp1(idx);

    parfor j = 1:length(tempB0map)
        % 2-pool fitting (water/MT)
        if flagSpectFitRef
            flagNodeltaB0 = true;
            [~,xfit] = estimateB0(df_ppm_cen, temp1(B0_fit_range), x0, xlb, xub, flagNodeltaB0);
            xfit(5) = xfit(5) - xfit(2);
            xfit(2) = 0;
            zSpect_fit(j,:) = zfun(xfit,df_ppm_orig);
        end
    end
    fprintf('done.\n')

    imgCEST_corr(mask(:),:) = zSpect_corr;
    imgCEST_corr = reshape(imgCEST_corr,size(imgCEST));
    B0map(mask(:),:) = tempB0map;
    B0map = reshape(B0map,size(mask));
    if flagSpectFitRef
        zSpectFitRef(mask(:),:) = zSpect_fit;
        zSpectFitRef = reshape(zSpectFitRef,size(imgCEST));
    end
end


function [deltaB0,xfit] = estimateB0(dB0,Zspec,x0,xlb,xub,flagNodeltaB0)
    if nargin < 6
        flagNodeltaB0 = false;
    end
    Zspec = Zspec/max(Zspec);
    if flagNodeltaB0
        x0(2)  = 0;
        xlb(2) = 0;
        xub(2) = 0;
    else
        [~,idx] = min(Zspec);
        x0(2)  = dB0(idx);
        xlb(2) = dB0(idx)-1;
        xub(2) = dB0(idx)+1;
    end
    
    [deltaB0,xfit] = WASSR(dB0, Zspec, x0, xlb, xub);

    % x(1): Amplitude water
    % x(2): dB0 water // here in ppm
    % x(3): FWHM(\omega_1?) water // here in ppm
    % x(4): Amplitude MT
    % x(5): dB0 MT // here in ppm
    % x(6): FWHM(\omega_1?) MT // here in ppm
    % x(7): Zbase
end


function [B0, xfit] = WASSR(db0_sampled, zspec_sampled, x0, xlb, xub)
    zfun = @(x,xdata) (1 - x(1).*(x(3)^2/4)./((x(3)^2/4)+(xdata-x(2)).^2)).*(1 - x(4).*(x(6)^2/4)./((x(6)^2/4)+(xdata-x(5)).^2)) - x(7);   
    opts = optimset('Display','off');
    [xfit,~] = lsqcurvefit(@(x,xdata) zfun(x,xdata), double(x0), double(db0_sampled), double(zspec_sampled), double(xlb), double(xub), opts);
    B0 = xfit(2);
end


function [coeff, score, latent] = pcaogc(X)
    % De-mean
    X = bsxfun(@minus,X,mean(X));
    % Calculate eigenvalues and eigenvectors of the covariance matrix
    covarianceMatrix = cov(X);
    [V,D] = eig(covarianceMatrix);
    % "coeff" are the principal component vectors. These are the eigenvectors of the covariance matrix. Compare ...
    coeff = V(:,end:-1:1);
    % Multiply the original data by the principal component vectors to get the projections of the original data on the
    % principal component vector space. This is also the output "score". Compare ...
    score = X*coeff;
    % The columns of X*coeff are orthogonal to each other. This is shown with ...
    % The variances of these vectors are the eigenvalues of the covariance matrix, and are also the output "latent". Compare
    % these three outputs
    latent = var(score)';
end
