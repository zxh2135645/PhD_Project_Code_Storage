function [dataArray,reconOptions] = binningResp_DCE(params,reconOptions,dataArray,temporalBasis,spatialCoeff)
%v0.2

% License check
if ~checkLicense(reconOptions)
    dlg = errordlg('Multitasking license check failed');
    waitfor(dlg);
    return;
end

if nargin < 6
    selectROI = false;
end

extractVarFromStruct(params);
extractVarFromStruct(reconOptions);
extractVarFromStruct(dataArray);
extractVarFromStruct(temporalBasis);
extractVarFromStruct(spatialCoeff);

L_init = size(Phi_rt_small_init,1);
cL = size(curvePhi,2);
Segidx = mod(navIndices-1,linesPerShot*moduleLength) + 1;

vec = @(x) x(:);
prep = @(x,st) reshape(x,st.Nd(1),st.Nd(2),st.Nz,[]);

% Relative slice thickness
imageDim = [Ny Nx Nz];
pNx = imageDim(abs(xDir));
pNy = imageDim(abs(yDir));
pNz = imageDim(abs(zDir));
slthick  = voxelSpacing(abs(zDir))*Norig/volumeFOV(abs(xDir));
newNz    = round(slthick*pNz);

if rbins > 1
    % liver dome tracking
    binningRespDomeTrackDCE;
else
    Ridx = ones(1,size(Phi_rt_small_init,2));
end

dataArray.Ridx   = Ridx(:).';
dataArray.Segidx = Segidx(:).';

temp = fft(Ridx);
temp(1) = 0;
[~,idx] = max(abs(temp));
dataArray.meanRPeriod = (lEchoSpacing*SGBlock*numel(Ridx))/idx;

reconOptions.rbins = max(Ridx);

%% Images for individual bins

dispim = @(x) x(floor(Ny/2)-floor(Nydisp/2) + (1:Nydisp), floor(Nx/2)-floor(Nxdisp/2) + (1:Nxdisp), :, :);

imageDim = [Nydisp Nxdisp Nz];
pNx   = imageDim(abs(xDir));
pNy   = imageDim(abs(yDir));
pNz   = imageDim(abs(zDir));
newNz = round(slthick*pNz);

Utemp = reshape(spatialCoeff.U_init,Ny,Nx,Nz,[]);
if MBfactor == 1
    Utemp = fftshift(Utemp,3);
end
if isCartesian
    Utemp = fftshift(Utemp,1);
end
Utemp = permute(dispim(Utemp),abs([zDir xDir yDir 4]));
Utemp = Utemp(:,:,floor(pNy/2) + 1,:);
if sign(zDir) > 0
    Utemp = flip(Utemp,1);
end

temp = zeros(pNz,pNx,rbins);
for j = 1:rbins
    temp(:,:,j) = abs(reshape(reshape(Utemp,[],L_init)...
                    *mean(Phi_rt_small_init(:,Ridx==j),2),pNz,pNx,[]));
end
cw = prctile(temp(:),99.5);
dataArray.binsRespMean = imresize(abs(temp)/cw,[newNz pNx]);

% dataArray.binsResp = cell(rbins,1);
% for j = 1:rbins
%     temp = flipud(abs(reshape(reshape(U_init,[],L_init)...
%                     *Phi_rt_small_init(:,Ridx==j),pNz,pNx,[])));
%     cw = prctile(temp(:),99);
%     dataArray.binsResp{j} = abs(temp)/cw;
% end

h = findall(groot,'Type','figure','Name','Binning Results');
if isempty(h)
    figure('Name','Binning Results','units','normalized','OuterPosition',[0.1 0.4 0.8 0.5]);
else
    figure(h);
end
subplot(2,2,3),plot(Ridx,'.-');axis([-inf inf 0 rbins+1]);title(sprintf('Ridx, mean cycle %.2f sec',dataArray.meanRPeriod));
subplot(2,2,4),plot(Segidx(:), Ridx(:),'.');axis([0 linesPerShot*moduleLength 0 rbins+1]);title('Ridx / Segment Index');


%% License check
function isValid = checkLicense(reconOptions)

% get machine ID
[~,strAdd] = genID;

md = java.security.MessageDigest.getInstance('MD5');    
try
    % load license file
    licenseID = loadLicenseFile(reconOptions);
    [hashes,dateNr] = getHash(licenseID);
    
    % check hash
    for n = 1:length(strAdd)
        ID   = double(strAdd{n})*sum(dateNr);
        hash = dec2hex(uint8(double(md.digest(ID))+128));
        hash = hash(:).';
        isValid = contains(hashes,hash);
        if isValid; break; end
    end
    if ~isValid
        fprintf(2,'License check error: invalid hash.\n');
        disp(['Current date: ' date]);
        disp('Current system info:');
        for n = 1:length(strAdd)
            disp(['    ' strAdd{n}]);
        end
    end
catch errormsg
    fprintf(2,'License check error: %s\n', errormsg.message);
    isValid = false;
end
