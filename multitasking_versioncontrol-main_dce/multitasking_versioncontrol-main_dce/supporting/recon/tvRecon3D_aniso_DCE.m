function [reconOptions,spatialCoeff] = tvRecon3D_aniso_DCE(params,reconOptions,dataArray,temporalBasis,spatialCoeff,tempStructTensor)

extractVarFromStruct(params);
extractVarFromStruct(reconOptions);
%extractVarFromStruct(dataArray);
%extractVarFromStruct(temporalBasis);
%extractVarFromStruct(spatialCoeff);
    
vec = @(x) x(:);

Phi_rt = temporalBasis.Phi_rt;
L = size(Phi_rt,1);
linOrder = st.linOrder_shift;
parOrder = st.parOrder_shift;

if MBfactor == 1
    SEs = fftshift(dataArray.SEs,3);
else
    SEs = dataArray.SEs;
end
if isCartesian
    SEs = fftshift(SEs,1);
end

SEs(~isfinite(SEs))=0;
msdev = dataArray.msdev;

maxIter = reconOptions.wavelet.maxIter;

% try
% Progress bar
progress = waitbar(0,'','Name', sprintf('Progress'),...
        'CreateCancelBtn',...
        'setappdata(gcbf,''canceling'',1)');
setappdata(progress,'canceling',0)
waitbar(0, progress, sprintf('Wavelet recon iteration 1/%d...',maxIter));

if isfield(temporalBasis,'Phi2')
    Phi2 = temporalBasis.Phi2;
else
    fprintf('Calculating Phi^2 ... ')
    Phi2 = zeros(L, L, Ny, Nz);  
    for npy = 1:Ny
        for npz = 1:Nz
            t_ind = (linOrder==npy) & (parOrder==npz);
            Phi2(:,:,npy,npz) = Phi_rt(:,t_ind) * Phi_rt(:,t_ind)';
        end
    end
end

Ahb = spatialCoeff.Ahb;
U   = spatialCoeff.U_tensor;

clear dataArray temporalBasis;

% Whiten U
Nshift = floor(N/2) - floor(Norig/2);
Wt = reshape(U,Nx,Ny,Nz,[]);
Wt = reshape(Wt(Nshift+(1:Norig),Nshift+(1:Norig),:,:),[],L);
Wt=Wt'*Wt/size(Wt,1);
Wt = inv(sqrtm(Wt));
WtWh = Wt*Wt';
morozov=tempStructTensor.morozov;

mov = @(x)reshape(x,Nx,Ny,Nz,[]);
% if ~isCartesian && ~strcmp(ScanType,'T2prep')
%   cropTV=@(x)x(ovsRO/2+(1:Norig),ovsSL/2+(1:Norig),:,:,:);
%   uncropTV=@(x)padarray(x,[ovsRO/2, ovsSL/2, 0, 0, 0]);
% else
%   cropTV=@(x)x;
%   uncropTV=@(x)x;
% end
  cropTV=@(x)x;
  uncropTV=@(x)x;

zweight=2*params.dReadoutFOV_mm*Nz/(params.dThickness_mm*Nx);
mov = @(x)reshape(x,Ny,Nx,Nz,[]);
Wf=@(U)cropTV(cat(5,mov(U)-circshift(mov(U),[1 0 0 0]),mov(U)-circshift(mov(U),[0 1 0 0]),zweight*(mov(U)-circshift(mov(U),[0 0 1 0]))));
W=@(U)Wf(reshape(U,[],L)*Wt);
Whf=@(U)U(:,:,:,:,1)-circshift(U(:,:,:,:,1),[-1 0 0 0])+U(:,:,:,:,2)-circshift(U(:,:,:,:,2),[0 -1 0 0])+zweight*(U(:,:,:,:,3)-circshift(U(:,:,:,:,3),[0 0 -1 0]));
Wh=@(U)vec(reshape(Whf(uncropTV(U)),[],size(Wt,2))*Wt');
WhW=@(U)vec(Whf(uncropTV(Wf(reshape(U,[],L)*WtWh))));
group=@(x)cat(5,sqrt(sum(abs(x(:,:,:,:,1:2)).^2,5)),abs(x(:,:,:,:,3)));
ungroup=@(x)cat(5,x(:,:,:,:,1),x(:,:,:,:,1),x(:,:,:,:,2));
figmake=@(WU)sum(sum(group(WU(:,:,1,:,:)),4),5);

WU=W(U);
temp=group(WU);
if flagAutoLambdaWavelet
lambda=2*msdev^2/mean(abs(temp(:)-median(real(temp(:)))-1i*median(imag(temp(:)))))
else 
    lambda = reconOptions.wavelet.lambda;
end
Y=zeros(size(WU));
% alpha=max(vec(temp))
alpha=1;
figure,hist(vec(temp),1000);
im0 = log(figmake(WU));
figure,imshow(im0,[]);
fprintf('Adjust lambda and alpha if desired, then return.\n')

if flagUseGPU
  delete(gcp('nocreate')); %delete current parpool
  parpool('local',8); %open parpool with 8 workers
  gpuWorkerReset(); %reassign workers
end

rho=lambda/alpha;

tic;
fprintf('tv recon iteration ');
for it = 2:maxIter  
    fprintf('\n -- %d/%d ',it,maxIter);
    if getappdata(progress,'canceling')
        break
    else
        waitbar(it/maxIter, progress, sprintf('tv recon iteration %d/%d...',it,maxIter));
    end
    keyboard;
    figTitle = ['tv recon iteration ' num2str(it) '/' num2str(maxIter)];
  WU = W(U);
  figure(100),imshow([im0, log(figmake(WU))],[]),drawnow;
  figure(101),imshow(5*[exp(im0),figmake(WU)]/max(exp(im0(:)))),drawnow;
  Z=WU+Y/rho;
  Zg = group(Z);
  Zg = max(abs(Zg)-alpha,0)./Zg;
  Zg(isnan(Zg))=0;
  Z=ungroup(Zg).*Z;
  Y=Y+rho*(WU-Z);
  
  step = 1.25;
  alpha=alpha/step;
  rho=lambda/alpha;
  
  Uold=U;
  switch Trajectory
    case 'Radial'
%       U=pcg2(@(x)AhA_ps3D(x,st,SEs,Phi2)+rho/2*WhW(x),...
%         Ahb(:)+rho/2*Wh(Z-Y/rho),[],median([5 it 10]),[],[],U);
%           U=pcg1(@(x)AhA_ps3D(x,st,SEs,Phi2)+rho/2*WhW(x),...
%         Ahb(:)+rho/2*Wh(Z-Y/rho),[],1,[],[],U);
      U=pcg2(@(x)AhA_ps3D(x,st,SEs,Phi2),...
        Ahb(:),[],1,[],[],U);
      case 'Linogram'
      U=pcg2(@(x)AhA_ps_lin(x,st,SEs,Phi_rt,A,Ah,thetas)+rho/2*WhW(x),...
        Ahb(:)+rho/2*Wh(wavelet_applyop(@minus,Z,wavelet_applyop(@rdivide,Y,rho))),[],median([5 it 10]),[],[],U);
    case 'Cartesian'
      U=pcg2(@(x)AhA_ps_cart(x,st,SEs,Phi_rt,Ny,Nx,Nz)+rho/2*WhW(x),...
        Ahb(:)+rho/2*Wh(wavelet_applyop(@minus,Z,wavelet_applyop(@rdivide,Y,rho))),[],median([5 it 10]),[],[],U);
  end
  eps(it) = norm(U(:)-Uold(:))/norm(Uold(:))

  clear Uold Ax;
  if (eps(it) < 1e-3) && (eps(it) ~= 0)
    break;
  end
end
delete(progress);
fprintf('\n Done. ');
% catch
%     delete(progress);
%     fprintf('\n Wavelet recon failed. ');
% end

toc;

%% Store results

reconOptions.wavelet.alpha  = alpha;
reconOptions.wavelet.lambda = lambda;

spatialCoeff.U = U;
spatialCoeff.U_wavelet = U;


%% License check
function isValid = checkLicense(reconOptions)

% get machine ID
if ispc
    [~,info] = system('getmac');
elseif ismac
    [~,info] = system('ifconfig en0 | grep ether');
elseif isunix
    [~,info] = system('ip addr | grep ether');
else
    error('OS not recognized.');
end
expression = ('\w\w:\w\w:\w\w:\w\w:\w\w:\w\w|\w\w-\w\w-\w\w-\w\w-\w\w-\w\w');
strAdd = regexp(info,expression,'match');
for n = length(strAdd):-1:1
    if strcmp(strAdd{n},'ff:ff:ff:ff:ff:ff') || strcmp(strAdd{n},'ff-ff-ff-ff-ff-ff')
        strAdd(n) = [];
    end
end

md = java.security.MessageDigest.getInstance('MD5');    
try
    % load license file
    [licenseID,flag] = loadLicenseFile(reconOptions);
    licenseID = licenseID/sum(double('multitasking'));
    if flag == 1        % offline license check
        datestr = licenseID(end);
        licenseID(end) = [];
        licenseDays = datestr - datenum(date);
        if licenseDays > 0 
            for n = 1:length(strAdd)
                ID = double(strAdd{n}) * datestr;
                hash = dec2hex(uint8(double(md.digest(ID))+128));
                hash = hash(:).';
                isValid = contains(char(licenseID),hash);
                if isValid; break; end
            end
        else
            error('License expired.');
        end
    else                % online license check
        % get hashes from server
        % hashes  = webread(['https://agchristodoulou.github.io/MTcheck/' htmlID '.html']);
        URL = char(licenseID);
        hashes  = webread(URL);
        
        % get Hash
        for n = 1:length(strAdd)
            ID   = double(strAdd{n})*sum(double(date));
            hash = dec2hex(uint8(double(md.digest(ID))+128));
            hash = hash(:).';
            isValid = contains(hashes,hash);
            if isValid; break; end
        end
    end
catch
    isValid = false;
end
