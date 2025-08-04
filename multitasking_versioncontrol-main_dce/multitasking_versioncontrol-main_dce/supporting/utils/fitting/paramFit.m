function fitResult = paramFit(fitParams)

if ~isfield(fitParams,'flag2stepFitting')
    fitParams.flag2stepFitting = false;
end

try 
    rbins = size(fitParams.Phi,4);
    if fitParams.rphase > rbins || fitParams.rphase < 1
        fitParams.rphase = 1;
    end
catch
    fitParams.rphase = 1;
end

try
    cbins = size(fitParams.Phi,3);
    for n = numel(fitParams.cphases):-1:1
        if fitParams.cphases(n) > cbins || fitParams.cphases(n) < 1
            fitParams.cphases(n) = [];
        end
    end
    fitParams.cphases = unique(fitParams.cphases);
    if isempty(fitParams.cphases)
        fitParams.cphases = 1:cbins;
    end
catch
    cbins = size(fitParams.Phi,3);
    fitParams.cphases = 1:cbins;
end

try
    for n = numel(fitParams.fitSlice):-1:1
        if fitParams.fitSlice(n) > fitParams.Nz || fitParams.fitSlice(n) < 1
            fitParams.fitSlice(n) = [];
        end
    end
    fitParams.fitSlice = unique(fitParams.fitSlice);
    if isempty(fitParams.fitSlice)
        fitParams.fitSlice = 1:fitParams.Nz;
    end
catch
    fitParams.fitSlice = 1:fitParams.Nz;
end

switch fitParams.ScanType
    case {'CEST'}
        fitResult = fit_CESTLorentzian(fitParams);
    case {'IR','SR'}
        fitResult = fit_T1_1FA(fitParams);
    case {'IR_VFA','SR_VFA'}
        fitResult = fit_T1_VFA(fitParams);
    case {'Cine'}
        fprintf('Parameter fitting: scanType Cine not supported\n');
    otherwise
        if fitParams.flag2stepFitting
            fitResult = fit_T1rhoT2IRVFA_2step(fitParams);
        else
            fitResult = fit_T1rhoT2IRVFA(fitParams);
        end
end

if size(fitParams.Phi,6) > 1
    fitParams.mask = fitResult.mask;
    temp = fit_T2star(fitParams);
    fitResult.T2starmap = temp.T2starmap;
    fitResult.RSST2starmap = temp.RSSmap;
%     fitResult.dFreqmap = temp.dFreqmap;

%     temp = fit_ME_deltaB0(fitParams);
%     fitResult.delta_B0  = temp.delta_B0;
% %     fitResult.R2starmap = temp.R2star;
%     fitResult.water = temp.water;
%     fitResult.fat = temp.fat;
%     fitResult.ff = temp.ff;
end


% switch fitParams.ScanType
%     case {'IR','SR'}
%         fitResult = fitT1_1FA(fitParams);
%     case 'T2IR'
%         fitResult = fitT1T2_1FA(fitParams);
%     case 'T2IR_VFA'
%         fitResult = fitT1T2_2FA(fitParams);
% end
