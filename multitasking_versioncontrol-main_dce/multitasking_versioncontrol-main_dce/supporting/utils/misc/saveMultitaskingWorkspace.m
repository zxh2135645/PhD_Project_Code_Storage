fprintf('Saving multitasking recon workspace... ')

% Remove redundant variables

if isfield(DataArray,'kspaceData_orig')
    DataArray = rmfield(DataArray,'kspaceData_orig');
end
if isfield(DataArray,'navData_orig')
    DataArray = rmfield(DataArray,'navData_orig');
end
if isfield(DataArray,'binsResp')
    DataArray = rmfield(DataArray,'binsResp');
end

vars = fieldnames(DataArray);
for n = 1:length(vars)
    if contains(vars{n},'recon')
        DataArray = rmfield(DataArray,vars{n});
    end
end

% Collect parameters needed for parametric fitting
try
    FitParams = createFitParams(Params,ReconOptions,DataArray,TemporalBasis,SpatialCoeff,FitParams);
catch
end

% Save workspace
try
    save([Params.filePath '/reconstructedImages'],'ReconstructedImages','-v7.3'); 
    save([Params.filePath '/fitParams'],'FitParams','-v7.3'); 
    save(Params.fileString,'Params','ReconOptions','DataArray','TemporalBasis','SpatialCoeff','TwixObj','-v7.3'); 
catch
    save('reconstructedImages','ReconstructedImages','-v7.3'); 
    save('fitParams','FitParams','-v7.3'); 
    save('multitasking.mat','Params','ReconOptions','DataArray','TemporalBasis','SpatialCoeff','TwixObj','-v7.3'); 
end

fprintf('done.\n')
diary off
