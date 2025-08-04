function [reconOptions,spatialCoeff] = waveletRecon2DSlices(params,reconOptions,dataArray,temporalBasis,spatialCoeff,flagContinue)

% License check
% if ~checkLicense(reconOptions)
%     dlg = errordlg('Multitasking license check failed');
%     waitfor(dlg);
%     return;
% end
    
% use spatialCoeff.U as initial U?
% 1: use U; 0: use U_tensor
if nargin < 6
    flagContinue = 0;
end

alpha = 1;

extractVarFromStruct(params);
extractVarFromStruct(reconOptions);
extractVarFromStruct(wavelet);
 
linOrder = st.linOrder_shift;
parOrder = st.parOrder_shift;  

alpha_init = alpha;

vec = @(x) x(:);


if MBfactor == 1
    SEs_all = fftshift(dataArray.SEs,3);
else
    SEs_all = dataArray.SEs;
end
if isCartesian
    SEs_all = fftshift(SEs_all,1);
end

msdev = dataArray.msdev;

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

Ahb_all = reshape(spatialCoeff.Ahb,Ny,Nx,Nz,[]);

if flagContinue
    U_all = spatialCoeff.U;
else
    U_all = spatialCoeff.U_tensor;
end
U_all = reshape(U_all,Ny,Nx,Nz,[]);
U_wavelet = zeros(size(U_all));

temppath = pwd;

tic;

try
    temp = sum(sum(abs(U_all(:,:,:,1,1)),1),3);
    xstart = find(temp>0,1,'first');
    xend = find(temp>0,1,'last');

for xIdx = xstart:xend
    alpha = alpha_init;
    
    U   = vec(U_all(:,xIdx,:,:));
    SEs = SEs_all(:,xIdx,:,:);
    Ahb = reshape(Ahb_all(:,xIdx,:,:),[],L);
    
    % initial images
    U_temp = reshape(U,Ny,Nz,[]);
    if MBfactor == 1
        U_temp = fftshift(U_temp,2);
    end
    if isCartesian
        U_temp = fftshift(U_temp,1);
    end
    U_temp   = reshape(U_temp,[],L*tempNecho);
    recon_init = reshape(abs(U_temp*Phi_temp),Ny,Nz,[]);
    recon_init = recon_init/prctile(recon_init(:),99);
    recon_init = reshape(permute(recon_init,[2 1 3]),Nz,[])';
        
    Wt = reshape(U,[],L);
    Wt = Wt'*Wt/size(Wt,1);
    Wt = inv(sqrtm(Wt));
    WtWh = Wt*Wt';

    %morozov = msdev^2*numel(kspaceData);

    %Set up wavelets
    dwtmode('per'); %periodic
    wlev   = 4; %wavelet level
    wname  = 'sym4'; %wavelet name
    padamt = ceil([Ny Nz]/2^wlev)*2^wlev-[Ny Nz]  %each dimension needs to be divisible by 2^wlev (16 for wlev=4)
    mov = @(x)reshape(x,Ny,Nz,[]);
    mat = @(x)reshape(x,[],L);
    if sum(padamt)>0 %if needs padding
        prep     = @(x)padarray(mov(x),[padamt 0],0,'post');
        prep_adj = @(x)mat(x(1:Ny,1:Nz,:));
    else
        prep     = @(x)mov(x);
        prep_adj = @(x)mat(x);
    end
    [~,waveS] = wavelet2D(prep(U),wlev,wname);
    W   = @(x)wavelet2D(prep(mat(x)*Wt),wlev,wname);
    Wh  = @(x)vec(prep_adj(iwavelet2D(x,waveS,wname))*Wt');
    WhW = @(U)vec(prep_adj(prep(mat(U)*WtWh)));

    figmake = @(WU)waveletfigmake2D(WU,waveS,wname,wlev);
    WU = W(U);
    Y  = zeros(size(WU));

    % lambda = 2*msdev^2/mean(abs(WU(:)-median(real(WU(:)))-1i*median(imag(WU(:)))));
    lambda = reconOptions.wavelet.lambda;

    if flagAutoLambdaWavelet
        nonzeros = @(x)x(x~=0);
        Wtnew = eye(L);
        for l = 1:L
            [~,kmc] = kmeans(nonzeros(abs(WU((l-1)*1+(1:1),prod(waveS(1,:))+1:end)).'),2);
            %[~,kmc] = kmeans(nonzeros(abs(WU(l,prod(waveS(1,:))+1:end)).'),2);
            Wtnew(l,l) = (lambertw(-2*exp(-2)) + 2)/mean(kmc); % b is approx mean(kmc)*0.6275. calculating 1/b here
        end

        if flagCommandLine
            %figure,plot(diag(Wtnew),'.')
        end

        %%
        Wt   = Wt*Wtnew;
        WtWh = Wt*Wt';

        W   = @(x)wavelet2D(prep(mat(gather(x))*Wt),wlev,wname);
        Wh  = @(x)vec(prep_adj(iwavelet2D(x,waveS,wname))*Wt');
        WhW = @(U)vec(prep_adj(prep(mat(U)*WtWh)));

        WU  = W(U);

        lambda = 2*msdev^2;
        reconOptions.wavelet.lambda = lambda*mean(diag(Wtnew));
        fprintf('Auto-lambda: new lambda is %0.1e, roughly equivalent to %0.1e using the old scaling; alpha = %f.\n',lambda,lambda*mean(diag(Wtnew)),alpha);
    else
        fprintf('Wavelet lambda = %g, alpha = %f.\n', lambda, alpha);
    end
    % figure,hist(vec(abs(WU)),1000);
    im0 = log(figmake(WU));
    %figure,imshow(im0,[]);
    %fprintf('Adjust lambda and alpha if desired, then return.\n')

    %keyboard;

    %%

    rho = lambda/alpha;
    % dcp=real(U(:)'*AhA_ps(U,st,SEs,Phi_rt,sp)-2*U(:)'*Ahb(:)+norm(kspaceData(:))^2);
    % l1=lambda*sum(abs(WU(:)));

%     if flagUseGPU
%         delete(gcp('nocreate')); %delete current parpool
%         parpool('local',4); %open parpool with 8 workers
%         gpuWorkerReset(); %reassign workers
%     end
        
    fprintf('Begin iteration ');
    for it = 2:maxIter
        fprintf('\n -- %d/%d ',it,maxIter);

        WU = W(U);
        Z = WU + Y/rho;
        Z = sign(Z).*max(abs(Z)-alpha,0);
        Y = Y+rho*(WU-Z);

        step  = 1.1;
        alpha = alpha/step;
        rho   = lambda/alpha;

        Uold = U;

        if isCartesian
            [U,flag] = pcg2(@(x)AhA_ps_cart(x,SEs,Phi2)+rho/2*WhW(x), Ahb(:)+rho/2*Wh(Z-Y/rho),[],median([5 it 10]),[],[],U);
        else
            if flagUseToeplitz
                [U,flag] = pcg2(@(x) AhA_ps_Toeplitz(x,st,SEs,Phi2,reconOptions.Om,flagUseGPU)+rho/2*WhW(x),Ahb(:)+rho/2*Wh(Z-Y/rho),[],median([5 it 10]),[],[],U);
            elseif flagUsefinufft && flagUseGPU
                cd([reconOptions.mainpath '/supporting/nufft/finufft/build']);
                [U,flag] = pcg2(@(x) AhA_ps3D_cufinufft(x,st,SEs,Phi2)+rho/2*WhW(x), Ahb(:)+rho/2*Wh(Z-Y/rho),[],median([5 it 10]),[],[],U);
            elseif flagUsefinufft
                [U,flag] = pcg2(@(x) AhA_ps3D_finufft(x,st,SEs,Phi2)+rho/2*WhW(x), Ahb(:)+rho/2*Wh(Z-Y/rho),[],median([5 it 10]),[],[],U);
            elseif flagUseGPU
                [U,flag] = pcg2(@(x) AhA_ps3D(x,st,SEs,Phi2)+rho/2*WhW(x), Ahb(:)+rho/2*Wh(Z-Y/rho),[],median([5 it 10]),[],[],U);
            else
                [U,flag] = pcg2(@(x) AhA_ps3D_irt(x,st,SEs,Phi2)+rho/2*WhW(x), Ahb(:)+rho/2*Wh(Z-Y/rho),[],median([5 it 10]),[],[],U);
            end
        end

        eps = norm(U(:)-Uold(:))/norm(Uold(:));
        %   dcp(it)=real(U(:)'*AhA_ps(U,st,SEs,Phi_rt,sp)-2*U(:)'*Ahb(:)+norm(kspaceData(:))^2);
        %   l1(it)=lambda*sum(abs(WU(:)));
        clear Uold;

        fprintf(' eps = %f',eps);

        % show intermediate images
        if (it == maxIter) || (eps < 1e-3)
            drawnow; if ~ishandle(102); break; end

            figTitle = ['Wavelet recon ' num2str(xIdx) '/' num2str(Nx) ': iteration ' num2str(it) '/' num2str(maxIter)];
            figure(100),imshow([im0, log(figmake(WU))],[]),title(figTitle);drawnow;
            figure(101),imshow(5*[exp(im0),figmake(WU)]/max(exp(im0(:)))),title(figTitle);drawnow;
            U_temp = reshape(U,Ny,Nz,[]);
            if isCartesian
                U_temp = fftshift(fftshift(U_temp,1),2);
            end
            U_temp   = reshape(U_temp,[],L*tempNecho);
            recon_temp = reshape(abs(U_temp*Phi_temp),Ny,Nz,[]);
            recon_temp = recon_temp/prctile(recon_temp(:),99);
            recon_temp = reshape(permute(recon_temp,[2 1 3]),Nz,[])';
            figTitle = ['Wavelet recon ' num2str(xIdx) '/' num2str(Nx) '- L:init, R:' num2str(it) '/' num2str(maxIter)];
            figure(102),imshow([recon_init recon_temp],'InitialMagnification',200),title(figTitle);drawnow;
        end
        
        if (eps < 1e-3)
            fprintf(', stop.\n');
            break;
        end
        
    end
    U_wavelet(:,xIdx,:,:) = reshape(U,Ny,1,Nz,[]);
    drawnow; if ~ishandle(102); fprintf ( ', user terminated wavelet recon after %d slices; recon not finished.\n', xIdx ); break; end
end
if ishandle(102); fprintf('\n Wavelet recon finished. '); end
catch errormsg
    fprintf(2,'\n Wavelet recon failed. \n');
    fprintf(2,'%s\n', errormsg.message);
end
toc;

cd(temppath);

%% Store results

reconOptions.wavelet.alpha  = alpha;

spatialCoeff.U = U_wavelet;
spatialCoeff.U_wavelet = U_wavelet;



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