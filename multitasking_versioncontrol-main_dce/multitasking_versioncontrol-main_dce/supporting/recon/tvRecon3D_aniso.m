function [reconOptions,spatialCoeff] = tvRecon3D_aniso(params,reconOptions,dataArray,temporalBasis,spatialCoeff,flagContinue)

% License check
if ~checkLicense(reconOptions)
    dlg = errordlg('Multitasking license check failed');
    waitfor(dlg);
    return;
end

% use spatialCoeff.U as initial U?
% flagContinue = 1: use U; 0: use U_tensor
if nargin < 6
    flagContinue = 0;
end

alpha = 1;

extractVarFromStruct(params);
extractVarFromStruct(reconOptions);
extractVarFromStruct(tv);

vec = @(x) x(:);
dispim = @(x) x(floor(Ny/2)-floor(Nydisp/2) + (1:Nydisp), floor(Nx/2)-floor(Nxdisp/2) + (1:Nxdisp), dispSlice, :);

if MBfactor == 1
    SEs = fftshift(dataArray.SEs,3);
else
    SEs = dataArray.SEs;
end
if isCartesian
    SEs = fftshift(SEs,1);
end

msdev = dataArray.msdev;

tic;

Phi_rt = temporalBasis.Phi_rt;
Phi2   = temporalBasis.Phi2;

if isfield(temporalBasis,'Gr')
    Gr     = temporalBasis.Gr;
    Phi    = temporalBasis.Phi;
    L      = size(Gr,1);
    tempNecho = size(Phi,1)/L;
    Phi_temp = zeros(L*tempNecho,min(cbins,2));
    if tempNecho == Necho
        for n = 1:Necho
            Phi_temp(n:Necho:end,:) = Gr\reshape(Phi(n:Necho:end,floor(min(size(Phi,2),linesPerShot)/2)+1,1:ceil(cbins/2):end,1,1,1),L,[]);
        end
    else
        Phi_temp = Gr\reshape(Phi(:,floor(min(size(Phi,2),linesPerShot)/2)+1,1:ceil(cbins/2):end,1,1,1),L,[]);
    end
else
    L = size(Phi,1);
    Phi_temp = zeros(L,min(cbins,2));
    Phi_temp(:,1) = mean(Phi_rt(:,dataArray.Ridx == 1 & dataArray.Hidx == 1),2);
    if size(Phi_temp,2) == 2
        Phi_temp(:,2) = mean(Phi_rt(:,dataArray.Ridx == 1 & dataArray.Hidx == floor(cbins/2)+1),2);
    end
    Phi_temp = reshape(Phi_temp,L,[]);
end
    
Ahb = spatialCoeff.Ahb;

if flagContinue
    U = spatialCoeff.U;
else
    try
        U = spatialCoeff.U_tensor;
    catch errormsg
        fprintf('%s\n ', errormsg.message);
        fprintf('\n Using "U" instead of "U_tensor". \n');
        U = spatialCoeff.U;
    end
end

clear dataArray temporalBasis;

% initial images
U_temp = reshape(U,Ny,Nx,Nz,[]);
if MBfactor == 1
    U_temp = fftshift(U_temp,3);
end
if isCartesian
    U_temp = fftshift(U_temp,1);
end
U_temp     = reshape(dispim(U_temp),[],L*tempNecho);
recon_init = reshape(abs(U_temp*Phi_temp),Nydisp,Nxdisp,[]);
recon_init = recon_init/prctile(recon_init(:),99);
recon_init = reshape(permute(recon_init,[2 1 3]),Nxdisp,[])';
recon_init = imageOrientLPS(recon_init,params);

if ~isCartesian
    Nyshift = floor(Ny/2) - floor(Nydisp/2);
    Nxshift = floor(Nx/2) - floor(Nxdisp/2);
    Wt = reshape(U,Ny,Nx,Nz,[]);
    Wt = reshape(Wt(Nyshift+(1:Nydisp),Nxshift+(1:Nxdisp),:,:),[],L);
else
    Wt = reshape(U,[],L);
end
Wt = Wt'*Wt/size(Wt,1);
Wt = inv(sqrtm(Wt));
WtWh = Wt*Wt';

%morozov = msdev^2 * numel(kspaceData);

if ~isCartesian && ~strcmp(ScanType,'T2prep')
    cropTV   = @(x) x(floor(Ny/2)-floor(Nydisp/2) + (1:Nydisp), floor(Nx/2)-floor(Nxdisp/2) + (1:Nxdisp), :, :, :);
    uncropTV = @(x) padarray(x,[floor(Ny/2)-floor(Nydisp/2), floor(Nx/2)-floor(Nxdisp/2), 0, 0, 0]);
else
    cropTV   = @(x) x;
    uncropTV = @(x) x;
end

mov = @(x) reshape(x,Ny,Nx,Nz,[]);
Wf  = @(U) cropTV(cat(5,(mov(U)-circshift(mov(U),[1 0 0 0]))/voxelSpacing(1),(mov(U)-circshift(mov(U),[0 1 0 0]))/voxelSpacing(2),(mov(U)-circshift(mov(U),[0 0 1 0]))/voxelSpacing(3)));
W   = @(U) Wf(reshape(U,[],L)*Wt);
Whf = @(U) (U(:,:,:,:,1)-circshift(U(:,:,:,:,1),[-1 0 0 0]))/voxelSpacing(1)+(U(:,:,:,:,2)-circshift(U(:,:,:,:,2),[0 -1 0 0]))/voxelSpacing(2)+(U(:,:,:,:,3)-circshift(U(:,:,:,:,3),[0 0 -1 0])/voxelSpacing(3));
Wh  = @(U) vec(reshape(Whf(uncropTV(U)),[],size(Wt,2))*Wt');
WhW = @(U) vec(Whf(uncropTV(Wf(reshape(U,[],L)*WtWh))));
group   = @(x) cat(5,sqrt(sum(abs(x(:,:,:,:,1:2)).^2,5)),abs(x(:,:,:,:,3)));
ungroup = @(x) cat(5,x(:,:,:,:,1),x(:,:,:,:,1),x(:,:,:,:,2));

if isCartesian
    figmake = @(WU) dispim(fftshift(fftshift(sum(sum(group(WU),4),5),1),3));
elseif MBfactor == 1
    figmake = @(WU) dispim(fftshift(sum(sum(group(WU),4),5),3));
else
    figmake = @(WU) dispim(sum(sum(group(WU),4),5));
end

WU = W(U);
Y = zeros(size(WU));

lambda = reconOptions.tv.lambda;
alpha  = reconOptions.tv.alpha;

temppath = pwd;

try
temp = group(WU);
lambda = 2*msdev^2/mean(abs(temp(:)-median(real(temp(:)))-1i*median(imag(temp(:)))))

im0 = log(figmake(WU));
figTitle = ['TV recon iteration 1/' num2str(maxIter)];
figure(100),imshow([im0, log(figmake(WU))],[]),title(figTitle);drawnow;
figure(101),imshow(1*[exp(im0),figmake(WU)]/max(abs(exp(im0(:))))),title(figTitle);drawnow;

% show intermediate images
figTitle = ['L:init, R: 1/' num2str(maxIter) ' (close this window to end TV recon)'];
figure(102),imshow([recon_init recon_init],'InitialMagnification',200),title(figTitle);drawnow;

% if flagUseGPU
%     delete(gcp('nocreate')); %delete current parpool
%     parpool('local',3); %open parpool with 8 workers
%     gpuWorkerReset(); %reassign workers
% end

rho = lambda / alpha;
tic;
fprintf('Begin iteration ');
for it = 2:maxIter
    fprintf('\n -- %d/%d ',it,maxIter);
    
    WU = W(U);
    Z  = WU + Y/rho;
    Zg = group(Z);
    Zg = max(abs(Zg)-alpha,0)./Zg;
    Zg(isnan(Zg))=0;
    Z  = ungroup(Zg).*Z;
    Y  = Y + rho*(WU-Z);
  
    step  = 1.25;
    alpha = alpha/step;
    rho   = lambda/alpha;

    Uold = U;
    switch Trajectory
        case 'Radial'
            if flagUseToeplitz
                [U,flag] = pcg2(@(x) AhA_ps_Toeplitz(x,st,SEs,Phi2,reconOptions.Om,flagUseGPU)+rho/2*WhW(x),...
                    Ahb(:)+rho/2*Wh(Z-Y/rho),[],median([5 it 10]),[],[],U);
            elseif flagUsefinufft && flagUseGPU
                cd([reconOptions.mainpath '/supporting/nufft/finufft/build']);
                [U,flag] = pcg2(@(x) AhA_ps3D_cufinufft(x, st, SEs, Phi2)+rho/2*WhW(x),...
                    Ahb(:)+rho/2*Wh(Z-Y/rho),[],median([5 it 10]),[],[],U);
            elseif flagUsefinufft
                [U,flag] = pcg2(@(x) AhA_ps3D_finufft(x, st, SEs, Phi2)+rho/2*WhW(x),...
                    Ahb(:)+rho/2*Wh(Z-Y/rho),[],median([5 it 10]),[],[],U); 
            elseif flagUseGPU
                [U,flag] = pcg2(@(x) AhA_ps3D(x,st,SEs,Phi2)+rho/2*WhW(x),...
                    Ahb(:)+rho/2*Wh(Z-Y/rho),[],median([5 it 10]),[],[],U); 
            else
                [U,flag] = pcg2(@(x) AhA_ps3D_irt(x, st, SEs, Phi2)+rho/2*WhW(x),...
                    Ahb(:)+rho/2*Wh(Z-Y/rho),[],median([5 it 10]),[],[],U); 
            end
        case 'Linogram'
            [U,flag] = pcg2(@(x) AhA_ps_lin(x,st,SEs,Phi_rt,A,Ah,thetas)+rho/2*WhW(x),...
                            Ahb(:)+rho/2*Wh(Z-Y/rho),[],median([5 it 10]),[],[],U);
        case 'Cartesian'
            [U,flag] = pcg2(@(x) AhA_ps_cart(x,SEs,Phi2)+rho/2*WhW(x),...
                            Ahb(:)+rho/2*Wh(Z-Y/rho),[],median([5 it 10]),[],[],U);
    end
    eps = norm(U(:)-Uold(:))/norm(Uold(:));
    clear Uold;
    
    fprintf(' eps = %f',eps)
       
    drawnow; if ~ishandle(102); fprintf ( ', user terminated TV recon after %i iterations\n', it ); break; end
    
    WU = W(U);
    
    figTitle = ['TV recon iteration ' num2str(it) '/' num2str(maxIter)];
    figure(100),imshow([im0, log(figmake(WU))],[]),title(figTitle);drawnow;
    figure(101),imshow(1*[exp(im0),figmake(WU)]/max(exp(im0(:)))),title(figTitle);drawnow;
    
    % show intermediate images
    U_temp = reshape(U,Ny,Nx,Nz,[]);
    if MBfactor == 1
        U_temp = fftshift(U_temp,3);
    end
    if isCartesian
        U_temp = fftshift(U_temp,1);
    end
    U_temp     = reshape(dispim(U_temp),[],L*tempNecho);
    recon_temp = reshape(abs(U_temp*Phi_temp),Nydisp,Nxdisp,[]);
    recon_temp = recon_temp/prctile(recon_temp(:),99);
    recon_temp = reshape(permute(recon_temp,[2 1 3]),Nxdisp,[])';
    recon_temp = imageOrientLPS(recon_temp,params);
    figTitle = ['L:init, R:' num2str(it) '/' num2str(maxIter) ' (close this window to end TV recon)'];
    figure(102),imshow([recon_init recon_temp],'InitialMagnification',200),title(figTitle);drawnow;

    if (eps < 1e-3)
        fprintf(', stop.');
        break;
    end
    
    drawnow; if ~ishandle(102); fprintf ( ', user terminated TV recon after %i iterations\n', it ); break; end

end
fprintf('\n TV recon finished. '); 
catch errormsg
    fprintf(2,'\n TV recon failed. \n');
    fprintf(2,'%s\n ', errormsg.message);
end
toc;

cd(temppath);

reconOptions.tv.alpha  = alpha;

spatialCoeff.U = U;
spatialCoeff.U_TV = U;

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
