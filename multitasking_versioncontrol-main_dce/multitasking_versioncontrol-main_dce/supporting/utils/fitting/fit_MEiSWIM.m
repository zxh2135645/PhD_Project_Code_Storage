function fitResult = fit_MEiSWIM(fitParams)

vec = @(x) x(:);
row = @(x) x(:).';

cphases = 1;
rphase  = 1;

fitSlice = 1:size(fitParams.U,3);

combinedEchoes = 1:size(fitParams.Phi,6);       % which echoes to combine (if one of the echoes does not
                                                % have good quality you can exclude it from the final averaging)
                        
extractVarFromStruct(fitParams);

%% Generate ME images
Necho = size(Phi,6);
L     = size(Gr,1);
tempNecho = size(Phi,1)/L;

Phi = reshape(permute(Phi,[1 2 5 3 4 6]),L*tempNecho,Nseg*moduleLength,size(Phi,3),size(Phi,4),Necho);

for n = numel(fitSlice):-1:1
    if fitSlice(n) > Nzdisp || fitSlice(n) < 1
        fitSlice(n) = [];
    end
end

[Ny,Nx,Nz,~,~] = size(U);
dispim = @(x) x(floor(Ny/2)-floor(Nydisp/2) + (1:Nydisp), floor(Nx/2)-floor(Nxdisp/2) + (1:Nxdisp), floor(Nz/2)-floor(Nzdisp/2) + (1:Nzdisp), :, :);
Utemp = reshape(dispim(U),Nydisp,Nxdisp,Nzdisp,tempNecho,L);
Utemp = Utemp(:,:,fitSlice,:,:);

for echo = 1:numel(combinedEchoes)
    if combinedEchoes(echo) > Necho || combinedEchoes(echo) < 1
        combinedEchoes(echo) = [];
    end
end

recon = zeros(Nydisp,Nxdisp,Nzdisp,numel(combinedEchoes));
for echo = 1:numel(combinedEchoes)
    if tempNecho == Necho
        PhiT1 = Gr\reshape(Phi(echo:tempNecho:end,Nseg,cphases(1),rphase,echo),L,[]);
        temp = reshape(Utemp(:,:,:,echo,:),[],L);
    else
        PhiT1 = Gr\reshape(Phi(:,Nseg,cphases(1),rphase,echo),L,[]);
        temp = reshape(Utemp(:,:,:,1,:),[],L);
    end
    recon(:,:,:,echo) = reshape(temp*PhiT1,Nydisp,Nxdisp,[],1);
end

rawVoxelSpacing = params.rawVoxelSpacing;
newVoxelSpacing = rawVoxelSpacing;
% [~,minSpacing] = min(rawVoxelSpacing);
% temp = 1:3; temp(minSpacing) = [];
% tempIdx1 = temp(1);
% tempIdx2 = temp(2);
% newRatio1 = rawVoxelSpacing(temp(1))/rawVoxelSpacing(minSpacing);
% newRatio2 = rawVoxelSpacing(temp(2))/rawVoxelSpacing(minSpacing);
% [~,Idx] = sort([[tempIdx1 tempIdx2 minSpacing]]);
% newRatio = [newRatio1 newRatio2 1];
% newVoxelSpacing = row(rawVoxelSpacing)./row(newRatio(Idx));
% newVoxelSpacing(3) = rawVoxelSpacing(3);
% newRatio = newRatio(Idx(1:2));
% size2 = @(x) [size(x,1) size(x,2)];

recon = recon./max(abs(recon(:)));
reconMag = abs(recon);
reconMag = reconMag/max(reconMag(:))*2048;
reconPhs = angle(recon);
reconPhs = (reconPhs/2/pi + 0.5)*4095;

[Nydisp,Nxdisp,Nzdisp,~] = size(recon);
Ny = 2.^(ceil(log2(Nydisp)));
Nx = 2.^(ceil(log2(Nxdisp)));
Nz = 2.^(ceil(log2(Nzdisp)));
newDim = [Ny Nx Nz];

vecCol  = params.RotMatQ*[1;0;0];
vecRow  = params.RotMatQ*[0;1;0];
ImageOrientationPatient = round([vecRow; vecCol], 6);

B0 = params.lResonanceFrequency/42.5764e6;

TEarray = TEarray*1e3;  % convert from sec to msec

if isunix
    sep = '/';
else
    sep = '\';
end

filePath = [pwd sep 'QSMoutput' sep];

% MEQSM parameters
correctPhaseShift = 0;
swimsign = 1;
thBETmask = 0.35; % BET threhsold

disp('MEQSM processing is starting...');
% [fitResult.MEQSM,fitResult.MEQSM_filled] = runMEQSM_forTianle(reconMag,reconPhs,size(reconMag),newDim,newVoxelSpacing,ImageOrientationPatient,B0,TEarray,combinedEchoes,filePath,correctPhaseShift,swimsign);
[fitResult.MEQSM,fitResult.MEQSM_filled] = runMEQSM_forTianle_v2(reconMag,reconPhs,size(reconMag),newDim,newVoxelSpacing,ImageOrientationPatient,B0,TEarray,combinedEchoes,filePath,correctPhaseShift,swimsign,thBETmask);
disp(['MEQSM processing is done and results are stored at: ' filePath]);

