function [imgCEST_corr, B0map] = CESTdoB0Correction(imgCEST, df_ppm, mask)
%   Do B0 correction with parfor loop
        
    size3 = @(x) [size(x,1),size(x,2),size(x,3)];

    if nargin < 3
        mask = ones(size3(imgCEST));
    end
    
    imgCEST = reshape(imgCEST,[],numel(df_ppm));
    imgCEST_corr = zeros(size(imgCEST));
    B0map = zeros(size(imgCEST,1),1);
    zSpect  = imgCEST(mask(:),:);
    zSpect_corr = zeros(size(zSpect));
    tempB0map = zeros(size(zSpect,1),1);

    B0_fit_range = abs(df_ppm) <= 1;
    df_ppm_cen = df_ppm(B0_fit_range);
    fprintf('\nStart B0 correction: ')
    parfor j = 1:length(tempB0map)
        zSpecttemp = zSpect(j,:);
        deltaB0 = estimateB0(df_ppm_cen, zSpecttemp(B0_fit_range), 'lsq');
        
        % Save B0 map and correct B0
        tempB0map(j) = deltaB0;
        zSpect_corr(j,:) = interp1(df_ppm-deltaB0, zSpecttemp, df_ppm, 'spline');
    end
    fprintf('done.\n')

    imgCEST_corr(mask,:) = zSpect_corr;
    imgCEST_corr = reshape(imgCEST_corr,size(imgCEST));
    B0map(mask,:) = tempB0map;
    B0map = reshape(B0map,size(mask));
end


function deltaB0 = estimateB0(dB0, Zspec, optMethod)
    optMethods = {'lsq', 'min'};
    
    Zspec = Zspec/max(Zspec);
    minSearchB0 = linspace(min(dB0), max(dB0), 100);
    
    if strcmp(optMethod, optMethods(1))      % lsq
        minB0 = WASSR(dB0, minSearchB0, Zspec);
    elseif strcmp(optMethod, optMethods(2))  % min
        Zspec_interp = spline(dB0,Zspec);
        [~, tmpIndex] = min(ppval(Zspec_interp, minSearchB0));
        minB0 = minSearchB0(tmpIndex);
    end
    
    deltaB0 = minB0;
end


function B0 = WASSR(db0_sampled, minSearchB0, zspec_sampled)
    myfun = @(x,xdata) (1 - x(1).*(x(3)^2/4)./((x(3)^2/4)+(xdata-x(2)).^2));
    x0 = [1, 0, 10];
    lb = [0, -2, 0.1];
    ub = [1, 2, 15];
    
    ydata1 = spline(db0_sampled, zspec_sampled, minSearchB0);
    [~,b] = min(ydata1);
    x0(2) = minSearchB0(b);
    
    opts = optimset('Display','off');
    
    [x,~,~,~] = lsqcurvefit(@(x,xdata) myfun(x,xdata), x0, db0_sampled, zspec_sampled, lb, ub, opts);
    B0 = x(2);
    
    % x(1): Amplitude
    % x(2): dB0 // here in ppm
    % x(3): FWHM(\omega_1?) // here in ppm
end