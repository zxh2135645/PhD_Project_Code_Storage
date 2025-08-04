function [dataArray,temporalBasis] = tensorSubspace(tempStructTensor,params,reconOptions,dataArray,temporalBasis)

% License check
% if ~checkLicense(reconOptions)
%     dlg = errordlg('Multitasking license check failed');
%     waitfor(dlg);
%     return;
% end

vec = @(x) x(:);
row = @(x) x(:).';

collapse  = @(x,dim,sizes) reshape(permute(reshape(x,sizes),[1:(dim-1), (dim+1):numel(sizes), dim]),[],sizes(dim));
icollapse = @(x,dim,sizes) ipermute(reshape(x,sizes([1:(dim-1), (dim+1):numel(sizes), dim])),[1:(dim-1), (dim+1):numel(sizes), dim]);

extractVarFromStruct(params);
extractVarFromStruct(reconOptions);
extractVarFromStruct(dataArray);
extractVarFromStruct(temporalBasis);
extractVarFromStruct(tempStructTensor);

lr  = tensor.lr;
sms = tensor.sms;
smorders = tensor.smorders;
circs = tensor.circs;

% curvePhi_orig = temporalBasis.curvePhi;
% curvePhi      = tempStructTensor.curvePhi;

%%
% initial guess
if ~doSpline && isUndersampled && ~doBloch
    tempdim = 3; %find(sizes~=1,1,'last');
    tempsizes = ones(size(sizes));
    tempsizes(tempdim) = sizes(tempdim);
    navData_sm = repmat(sum(navData_tensor.*mask,tempdim)./sum(mask,tempdim),tempsizes);
    navData_sm(logical(mask)) = navData_tensor(logical(mask));
    navData_sm(isnan(navData_sm)) = 0;
end

fprintf('Calculating tensor subspace...');

%%
% curvePhi = [];
% doBloch  = ~isempty(curvePhi);
% doSpline = sum(abs(sms)) ~= 0;
% isUndersampled = ~exp(sum(log(double(mask(:))))); 

if doBloch && doSpline
    if lowmem == 2
        navData_sm = LRTCp_splines_lm(navData_bloch.*logical(mask),0,sms,mask,5,navData_bloch_cL,curvePhi,smorders,circs,ts_proj);
        navData_sm = LRTCp_splines_lm(navData_tensor,lr,sms,mask,20,navData_sm,curvePhi,smorders,circs,ts_proj);
    else
        if strcmp(ScanType,'T2prep')
            navData_sm = repmat(sum(navData_bloch.*mask,5)./sum(mask,5),[1 1 1 1 sizes(5)]);
            navData_sm(logical(mask)) = navData_tensor(logical(mask));
        else
            navData_sm = LRTCp_splines(navData_bloch.*logical(mask),0,sms,mask,5,navData_bloch.*logical(mask),[],smorders,circs,ts_proj);
        end
        navData_sm = LRTCp_splines(navData_tensor,lr,sms,mask,20,navData_sm,curvePhi,smorders,circs,ts_proj);
        %     morozov_new
    end
elseif doSpline
    if lowmem == 2
        navData_sm = LRTCp_splines_lm(navData_bloch.*logical(mask),0,sms,mask,5,navData_bloch_cL,curvePhi,smorders,circs,ts_proj);
        navData_sm = LRTCp_splines_lm(navData_tensor,lr,sms,mask,20,navData_sm,curvePhi,smorders,circs,ts_proj);
    else
        navData_sm = LRTCp_splines(navData_tensor,0,sms,mask,5,navData_tensor.*logical(mask),[],smorders,circs,ts_proj);
        navData_sm = LRTCp_splines(navData_tensor,lr,sms,mask,20,navData_sm,[],smorders,circs,ts_proj);
    end
elseif doBloch && isUndersampled
    navData_sm = LRTCp(navData_tensor,lr,mask,20,navData_bloch.*logical(mask),curvePhi);
elseif doBloch && ~isUndersampled
    if lowmem == 2
        navData_sm = navData_bloch_cL;
    else
        navData_sm = navData_bloch;
    end
elseif isUndersampled
    navData_sm = LRTCp(navData_tensor,lr,mask,20,navData_sm);
else
    navData_sm = navData_tensor;
end

sizes = size(navData_sm);

switch ScanType
    case 'Cine'
        ranks=[L_tensor, cL, sizes(3:end)];
    case 'SR'
        ranks=[L_tensor, min(10,cL), sizes(3:end)];
    case {'IR', 'IR_VFA'}
        ranks = [L_tensor, min(10,cL), sizes(3:end)];
    case 'T2prep'
        ranks=[L_tensor, cL, sizes(3:end-1), 32];
    otherwise
        ranks = [L_tensor, cL, sizes(3:end)];
end

for j = numel(ranks):-1:2
    if (ranks(j) == 1) && (size(navData_sm,j) == 1)
        ranks(j) = [];
    end
end

[C,UU,ranks] = choose_C(squeeze(navData_sm),ranks);

% if T2_IR, make sure the dimension of UU, Vincent X. Mao
% 02/11/20
% if doBloch && (strcmp(ScanType,'T2IR') || strcmp(ScanType,'T2IR_VFA')) && lowmem~=2  %reinforce dictionary subspace
%     UU = reshape(curvePhi*(pinv(curvePhi)*reshape(UU,Nseg,[])),size(UU));
% end
% if strcmp(ScanType,'SR')
%     fc = lcm(Nseg,SGBlock)/Nseg;
%     Segidx = mod((1:Ntpoint/(SGBlock/fc)).'-1,Nseg/(SGBlock/fc))+1;
%     Segidx = Segidx(1:fc:end);
%     if doBloch
%         UU = reshape(temporalBasis.curvePhi*pinv(curvePhi)*reshape(UU,ceil(Nseg/(SGBlock/fc)),[]),[], size(UU,2));
%     end
% end
%

Phi = C*UU';
clear UU;

if lowmem == 2
    Phi = reshape(Phi, [size(C,1) sizes_lm(2:end)]);
else
    Phi = reshape(Phi, [size(C,1) sizes(2:end)]);
end
sizes = size(Phi);
sizes_orig(1) = sizes(1);

if strcmp(ScanType,'CEST') && ~flagDataDriven && lowmem == 2
    Phi = icollapse(permute(reshape(collapse(Phi, 2, sizes) * curvePhi_orig.',[],Nseg_orig,sizes_orig(5)),[1 3 2]), 2, sizes_orig);
elseif strcmp(ScanType,'CEST') && ~flagDataDriven
    Phi = icollapse(permute(reshape(collapse(Phi, 2, sizes) * pinv(curvePhi.') * curvePhi_orig.',[],Nseg_orig,sizes_orig(5)),[1 3 2]), 2, sizes_orig);
elseif strcmp(ScanType,'Cine')
    %Phi = Phi;
elseif lowmem == 2
    Phi = icollapse(collapse(Phi, 2, sizes) * curvePhi_orig.', 2, sizes_orig);
elseif ~isempty(curvePhi)
    Phi = icollapse(collapse(Phi, 2, sizes) * pinv(curvePhi.') * curvePhi_orig.', 2, sizes_orig);
elseif isempty(curvePhi) && lowmem == 1
    Phi = icollapse(interp1Segmented((collapse(Phi, 2, sizes)).',uni,linesPerShot,'cols','pchip').', 2, sizes_orig);
elseif isempty(curvePhi) && lowmem == 0
    sizes = size(Phi(:,uni,:,:,:));
    Phi = icollapse(interp1Segmented((collapse(Phi(:,uni,:,:,:), 2, sizes)).',uni,linesPerShot,'cols','pchip').', 2, sizes_orig);
end

% wallClock = repmat(vec(repmat(1:moduleLength, [Nseg_orig/moduleLength 1])), [ceil(Ntpoint/Nseg_orig) 1]);
% wallClock = wallClock(1:Ntpoint);  
% Phi = permute(reshape(Phi,[sizes(1) Nseg_orig/moduleLength moduleLength sizes(3:end)]),[1 2 4 5 3]);
sizes = size(Phi);
sizes(end+1:6) = 1;

% % interpolate full timings if necessary.
% if ~doBloch
%     Phi = ipermute(reshape(interp1(unique(mod(navIndices-1,Nseg)+1),reshape(permute(Phi,[2 1 3 4 5]),sizes(2),[]),1:sizes(2),'linear','extrap'),[Nseg L sizes(3:end)]),[2 1 3 4 5]);
%     Phi = reshape(Phi,[size(C,1) sizes(2:end)]);
% end


%% degate from Phi
tempNecho = size(Phi,6);
if ~isempty(find(Ridx==0,1))
    Ridx = interp1(navTiming(Ridx>0),Ridx(Ridx>0),navTiming,'linear','extrap').';
    Ridx = round(Ridx);
    Ridx(Ridx>rbins) = rbins;
    Ridx(Ridx<1) = 1;
end
if strcmp(ScanType,'Cine')
    Ridx_full = interp1(navTiming,Ridx(:),ACQTiming,'nearest','extrap').';
    %Hidx_full = interp1(navIndices,Hidx(:),1:Ntpoint,'nearest','extrap').';
    Hidx_full = round(cbins/(2*pi)*(pi+angle(interp1(navTiming,exp(1i*2*pi*(Hidx(:)/cbins-0.5)),ACQTiming,'linear','extrap')))).';
    Hidx_full(Hidx_full==0) = cbins;
elseif strcmp(ScanType,'CEST')
    Hidx_full = row(repmat(row(Hidx),[SGBlock 1]));
    Ridx_full = row(repmat(row(Ridx),[SGBlock 1]));
else
    Hidx_full = interp1Segmented(Hidx(:),navIndices,linesPerShot,'cols','nearest',navTiming,ACQTiming).';
    Ridx_full = interp1Segmented(Ridx(:),navIndices,linesPerShot,'cols','nearest',navTiming,ACQTiming).';
    
    if ~flagUseFirstEcho
        Hidx_full = row(repmat(Hidx_full,[tempNecho 1]));
        Ridx_full = row(repmat(Ridx_full,[tempNecho 1]));
    end
end

% switch ScanType
%     case 'SR'
%         wallClock_full = ceil((1:Ntpoint)/(size(Phi,2)*NDCEAve));
% %    case 'IR'
% %         wallClock_full = vec(repmat(1:Npreps, [Nseg Nrep])).';
% %         wallClock_full = wallClock_full(:).';
% %         wallClock_full = wallClock_full(1:Ntpoint);
%     case 'T2IR'
%         wallClock_full = row(repmat(1:moduleLength,[linesPerShot*Necho 1 ceil(Ntpoint/linesPerShot/moduleLength)]));
%         wallClock_full = wallClock_full(1:numel(ACQTiming)*Necho);
%     otherwise
%         wallClock_full = wallClock;
% end

wallClock_full = row(repmat(reshape(wallClock,tempNecho,[]),[SGBlock 1]));    
wallClock_full = wallClock_full(1:numel(ACQTiming)*tempNecho);

%XZ 11/14/2024
%Segidx_temp = [Segidx;Segidx+1];
%Segidx_full = row(repmat(Segidx_temp, [tempNecho 1]));

Segidx_full  = row(repmat(mod(0:Ntpoint-1,Nseg_orig) + 1,[tempNecho 1]));

echoIdx_full = repmat(1:tempNecho,[1 Ntpoint]);
Phi_rt_full  = degate(Phi(:,:).',sizes,Segidx_full,Hidx_full,Ridx_full,wallClock_full,echoIdx_full); %(C*UU').'
if strcmp(ScanType,'CEST')
    Phi_rt_full = Phi_rt_full.*row(curveMaskFull);
end
if flagUseFirstEcho
    Phi_rt_small = Phi_rt_full(:,navIndices);
    Phi_rt       = Phi_rt_full;
    Phi_rt(:,navIndices) = [];
else
    Phi_rt_small = Phi_rt_full(:,navIndices_full);
    Phi_rt       = Phi_rt_full;
    Phi_rt(:,navIndices_full) = [];
end


%%

if ~flagDataDriven
    if ~strcmp(ScanType,'CEST')
        curveMaskNav = ones(size(wallClock));
        curveMaskFull = ones(size(wallClock_full));
    end

    if flagUseFirstEcho
        navData = navData(:,:,:,:,1);
    else
        navData = reshape(permute(navData,[5 1 2 3 4]),size(navData,1)*size(navData,5),[]);
    end
    lrw = double(sqrt(sum(abs(navData(:,:).'-(navData(:,:).'*pinv(Phi_rt_small))*Phi_rt_small).^2)));
    h   = fminsearch(@(h)abs(norm(exp(-lrw.^2/h))^2/sqrt(numel(lrw))/norm(exp(-lrw.^2/h).^2)-sqrt(.95)),median(lrw).^2);
    
    x1 = median(lrw)./lrw;
    y1 = exp(-lrw.^2/h)/exp(-median(lrw).^2/h);
    S = @(a,b) a*x1 + b;
    [tempfit, res] = lsqnonlin(@(u)abs(S(u(1),u(2))-y1),[1,1],[0 0],[15 2]);  
    res_norm = res/norm(y1);
    
    try
        if strcmp(ScanType,'CEST')
            x2 = vec(realify(vec(Phi(1,:,1,1,:,1)./Phi(1,end,1,1,end,1))));
        else
            x2 = vec(realify(Phi(1,:,1,1,1,1)./Phi(1,end,1,1,1,1)));
        end
        y2 = vec(realify(curvePhi_orig(:,1)/curvePhi_orig(end,1)));
        res_phi = norm(x2-y2)/norm(y2);
        
        if flagCommandLine
            figure;
            subplot(1,2,1);plot(median(lrw)./lrw,exp(-lrw.^2/h)/exp(-median(lrw).^2/h),'.',x1,tempfit(1)*x1+tempfit(2),'.');axis([0 inf 0 inf]);title(['weight distribution (res = ' num2str(res_norm,'%.4f') ')']);drawnow;
            subplot(1,2,2);plot(x2);title('Phi1');axis([-inf inf -inf max(1,max(x2))]);drawnow;%title(['Phi1/curvePhi1 (res = ' num2str(res_phi,'%.4f') ')']);drawnow;
        end
    catch
    end
    
    if ~exist('NDCEAve','var')
        NDCEAve = 2;
    end
    
    lrw = exp(-lrw.^2/h);
    fprintf('\nmedian(lrw) = %f... ', median(lrw));
    
    lrw = lrw/median(lrw);
    if strcmp(ScanType,'Cine')
        lrw = interp1(navIndices,lrw,1:Ntpoint,'linear','extrap');
    elseif strcmp(ScanType,'CEST')
        lrw = row(repmat(row(lrw),[SGBlock 1]));
    else
        if flagUseFirstEcho
            lrw = interp1Segmented(lrw,navIndices,linesPerShot,'rows','linear');
        else
            lrwtemp = zeros(Necho,ceil(navIndices(end)/linesPerShot)*linesPerShot);
            for eIdx = 1:Necho
                lrwtemp(eIdx,:) = interp1Segmented(lrw(eIdx:Necho:end),navIndices,linesPerShot,'rows','linear');
            end
            lrw = row(lrwtemp);
        end
    end
    lrw(:) = 1;
    % XZ 11/14/2024
    lrw_full = row(lrw(1:floor(Ntpoint/Nseg_orig/NDCEAve)*Nseg_orig*NDCEAve)).*row(curveMaskFull(1:floor(Ntpoint/Nseg_orig/NDCEAve)*Nseg_orig*NDCEAve));
    %lrw_full = row(lrw(Segidx_full)).*row(curveMaskFull(Segidx_full));
    if flagUseFirstEcho
        lrw(navIndices) = [];
    else
        lrw(navIndices_full) = [];
    end
    
    wallClockcount = max(wallClock_full(1:numel(lrw_full)));      
    fitw_full = regate(lrw_full(:).^2,Segidx_full(1:numel(lrw_full)),Hidx_full(1:numel(lrw_full)),Ridx_full(1:numel(lrw_full)),row(wallClock_full(1:numel(lrw_full))).*row(curveMaskFull(1:numel(lrw_full))));
    fitw_full = sqrt(sum(sum(fitw_full,3),4));
%     fitw_full = sqrt(sum(reshape(lrw_full.^2,Nseg_orig,[],wallClockcount),2));
    fitw_full = reshape(fitw_full,1,Nseg,wallClockcount)/median(fitw_full(:));
    %fitw_full = reshape(fitw_full,1,Nseg_orig,wallClockcount)/median(fitw_full(:));
    
    [Phi_rt,Gr1,Gr2] = svd(bsxfun(@times,Phi_rt,lrw)','econ');
else
    [Phi_rt,Gr1,Gr2] = svd(Phi_rt','econ');
end
Phi_rt = Phi_rt';
Gr = Gr2*Gr1;

Phi_rt_full  = Gr\Phi_rt_full;
Phi_rt_small = Gr\Phi_rt_small;

if flagUseFirstEcho
    Phi_rt_new = zeros(size(Phi_rt,1)*Necho,size(Phi_rt,2),Necho);
    for n = 1:Necho
        Phi_rt_new(n:Necho:end,:,n) = Phi_rt;
    end
    Phi_rt = Phi_rt_new;
     
    Phi_new = zeros([sizes(1)*Necho,sizes(2:5),Necho]);
    for n = 1:Necho
        Phi_new(n:Necho:end,:,:,:,:,n) = Phi;
    end
    Phi = Phi_new;
end

fprintf('done. \n');

%%

dataArray.Hidx_full = Hidx_full;
dataArray.Ridx_full = Ridx_full;
dataArray.wallClock_full = wallClock_full;
if ~flagDataDriven
    dataArray.lrw       = lrw;
    dataArray.lrw_full  = lrw_full;
    dataArray.fitw_full = fitw_full;
end

temporalBasis.Phi          = Phi;
temporalBasis.Phi_rt       = Phi_rt;
temporalBasis.Phi_rt_full  = Phi_rt_full;
temporalBasis.Phi_rt_small = Phi_rt_small;
temporalBasis.Gr           = Gr;


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