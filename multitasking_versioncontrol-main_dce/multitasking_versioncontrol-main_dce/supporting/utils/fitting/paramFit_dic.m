function fitResult = paramFit_dic(fitParams)

switch fitParams.ScanType
    case {'IR','SR'}
    otherwise
        fitResult = fit_T1rhoT2IRVFA_dic(fitParams);
end

% if strcmp(fitParams.ScanType,'T2IR')
%     fitResult = fit_T1T2_1FA_dic(fitParams);
% elseif strcmp(fitParams.ScanType,'T2IR_VFA')
%     fitResult = fit_T1T2_2FA_dic(fitParams);
% elseif strcmp(fitParams.ScanType,'IR')
% 
% end
