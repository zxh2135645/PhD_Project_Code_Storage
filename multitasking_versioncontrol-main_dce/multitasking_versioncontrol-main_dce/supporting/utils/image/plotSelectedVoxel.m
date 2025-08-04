function curve = plotSelectedVoxel(recon,ylim,cmap,slice,recon2,cIdx,tNorm)
% ========================================================================
%
% Plot time curves of selected voxel
%
%   Input:
%
%       recon  - 4D image array (Ny x Nx x Nz x Nt)
%
%       ylim   - voxel value display range
%                default is [-1 1]
%                if not specified, image will be normalized
%
%       cmap   - image display colormap
%                default is grayscale
%
%       slice  - display slice number
%                default is center slice
%
%       recon2 - second image array or fitResult
%                if fitResult, plot fitted curve
%                (currently only works for T2IR_1FA)
%
%       cIdx   - cardiac phase number, only used when recon2 is fitResult
%                default is 1
%
%   Output:
%
%       curve  - time curve of the current voxel
%
%
%  - 2022-03-20
%  - Created by Hsu-Lei Lee @ BIRI, Cedars-Sinai Medical Center
%
% ========================================================================

if nargin < 4 || isempty(slice)
    slice = floor(size(recon,3)/2) + 1;
end

[Ny,Nx,Nz,Nt,Np] = size(recon);
if nargin > 4 && ~isempty(recon2)
    if isstruct(recon2)
        [Nymap,Nxmap,Nzmap,Npmap] = size(recon2.T1map);
        Ntmap = recon2.fitParams.Nseg*recon2.fitParams.moduleLength;
    else
        [Nymap,Nxmap,Nzmap,Ntmap,Npmap] = size(recon2);
    end
    if (Ny==Nymap) && (Nx==Nxmap) && (Nt==Ntmap) && (Nz==Nzmap)
        recon  = squeeze(recon(:,:,slice,:,:));
        if ~isstruct(recon2)
            recon2 = squeeze(recon2(:,:,slice,:,:));
        end
    elseif (Ny==Nymap) && (Nx==Nxmap) && (Nt==Ntmap) && (Nz==1)
        recon  = squeeze(recon);
        if ~isstruct(recon2)
            recon2 = squeeze(recon2);
        end
    else
        disp('Arrays do not match.');
        recon2 = [];
    end
else
    recon  = squeeze(recon(:,:,slice,:,:));
    recon2 = [];
end

if nargin < 6 || (isempty(recon2) && cIdx > Np) || (isstruct(recon2) && cIdx > Npmap)
    cIdx = 1;
end
recon = recon(:,:,:,:,min(cIdx,Np));
if ~isstruct(recon2)
    recon2 = recon2(:,:,:,cIdx);
end
        
if nargin < 7
    tNorm = size(recon(:,:,:),3);
end

figure; af = gcf;
if isempty(recon2) || ~isfield(recon2,'T1map')
    if nargin < 3 || isempty(cmap)
        cmap = 'gray';
    end
    if nargin < 2
        recon = recon/max(abs(recon(:)));
        ylim = [-1 1];
        subplot(1,2,1);displayImage(mean(abs(recon),3),[],cmap);
    elseif isempty(ylim)
        recon = recon/max(abs(recon(:)));
        ylim = [-1 1];
        subplot(1,2,1);displayImage(mean(abs(recon),3),[],cmap);    
    else
        subplot(1,2,1);displayImage(mean(abs(recon),3),ylim,cmap);
    end
else
    if nargin < 3 || isempty(cmap)
        cmap = warmmetal;
    end
    recon = recon/max(abs(recon(:)));
    ylim = [-1 1];
    subplot(1,2,1);displayImage(wmedfilt2(recon2.T1map(:,:,slice,cIdx))*1000,[0 3000],cmap);
end
title('Click on a pixel to show signal curve');drawnow;
ax = gca;
h = impoint(ax);
position = round(h.getPosition);
curve = realify(squeeze(recon(position(2),position(1),:)));
curve = curve./sign(curve(end));
curve = curve(:);
strTitle = sprintf('voxel = (%d,%d)',position(1),position(2));
if ~isempty(recon2)
    if isstruct(recon2)
        switch recon2.fitParams.ScanType
            case {'IR','SR'}
                curve2 = simCurve_T1_1FA(recon2,[position(2),position(1),slice],cIdx,curve(:).');
            case {'IR_VFA','SR_VFA'}
                curve2 = simCurve_T1_VFA(recon2,[position(2),position(1),slice],cIdx,curve(:).');
            case {'Cine'}
                curve2 = zeros(size(curve));
            otherwise
                curve2 = simCurve_T1rhoT2IRVFA(recon2,[position(2),position(1),slice],cIdx,curve(:).');
                %curve2 = simCurve_T1T2_1FA(recon2,[position(2),position(1),slice],cIdx,curve(:).');
                %curve2 = simCurve_T1T2_2FA(recon2,[position(2),position(1),slice],cIdx,curve(:).');
        end
        strTitle = sprintf('%s, T1 = %d',strTitle,round(1000*recon2.T1map(position(2),position(1),slice,cIdx)));
        if isfield(recon2,'fitParams') && recon2.fitParams.numT2prep > 0
            strTitle = sprintf('%s, T2 = %d',strTitle,round(1000*recon2.T2map(position(2),position(1),slice,cIdx)));
        end
        if isfield(recon2,'fitParams') && recon2.fitParams.numT1rhoPrep > 0
            strTitle = sprintf('%s, T1rho = %d',strTitle,round(1000*recon2.T1rhomap(position(2),position(1),slice,cIdx)));
        end
        if isfield(recon2,'BIRmap')
            strTitle = sprintf('%s, B = %.2f',strTitle,recon2.BIRmap(position(2),position(1),slice,cIdx));
        end
        if isfield(recon2,'B1map')
            strTitle = sprintf('%s, beta = %.2f',strTitle,recon2.B1map(position(2),position(1),slice,cIdx));
        end
        strTitle2 = 'Blue: data, Red: fitted curve';
    else
        recon2 = squeeze(recon2(:,:,:))/max(abs(recon2(:)));
        curve2 = realify(squeeze(recon2(position(2),position(1),:)));
        curve2 = curve2./sign(curve2(tNorm));
        strTitle2 = 'Blue: data1, Red: data2';
    end
    if isstruct(recon2)
        curve(:,2) = curve2(:);
    else
        curve(:,2) = curve(tNorm,1)*curve2(:)/curve2(tNorm);
    end
else
    strTitle2 = 'data curve';
end
figure(af);subplot(1,2,2);plot(curve,'.-');title(strTitle2);
axis([0 size(squeeze(recon),3) ylim]);grid on;
title(ax,strTitle);

addNewPositionCallback(h,@(h) updateCurve(h,squeeze(recon),af,ax,ylim,slice,recon2,cIdx,tNorm));


function curve = updateCurve(position,recon,af,ax,ylim,slice,recon2,cIdx,tNorm)
position = round(position);
if sum(~isfinite(squeeze(recon(position(2),position(1),:)))) == 0
    curve = realify(squeeze(recon(position(2),position(1),:)));
    curve = curve./sign(curve(tNorm));
else
    curve = zeros(squeeze(recon(position(2),position(1),:)));
end
curve = curve(:);
strTitle = sprintf('voxel = (%d,%d)',position(1),position(2));
if ~isempty(recon2)
    if isstruct(recon2)
        switch recon2.fitParams.ScanType
            case {'IR','SR'}
                curve2 = simCurve_T1_1FA(recon2,[position(2),position(1),slice],cIdx,curve(:).');
            case {'IR_VFA','SR_VFA'}
                curve2 = simCurve_T1_VFA(recon2,[position(2),position(1),slice],cIdx,curve(:).');
            case {'Cine'}
                curve2 = zeros(size(curve));
            otherwise
                curve2 = simCurve_T1rhoT2IRVFA(recon2,[position(2),position(1),slice],cIdx,curve(:).');
                %curve2 = simCurve_T1T2_1FA(recon2,[position(2),position(1),slice],cIdx,curve(:).');
                %curve2 = simCurve_T1T2_2FA(recon2,[position(2),position(1),slice],cIdx,curve(:).');
        end
        strTitle = sprintf('%s, T1 = %d',strTitle,round(1000*recon2.T1map(position(2),position(1),slice,cIdx)));
        if isfield(recon2,'T2map')
            strTitle = sprintf('%s, T2 = %d',strTitle,round(1000*recon2.T2map(position(2),position(1),slice,cIdx)));
        end    
        if isfield(recon2,'BIRmap')
            strTitle = sprintf('%s, B = %.2f',strTitle,recon2.BIRmap(position(2),position(1),slice,cIdx));
        end
        if isfield(recon2,'B1map')
            strTitle = sprintf('%s, beta = %.2f',strTitle,recon2.B1map(position(2),position(1),slice,cIdx));
        end
        strTitle2 = 'Blue: data, Red: fitted curve';
    else
        curve2 = realify(squeeze(recon2(position(2),position(1),:)));
        curve2 = curve2./sign(curve2(tNorm));
        strTitle2 = 'Blue: data1, Red: data2';
    end
    if isstruct(recon2)
        curve(:,2) = curve2(:);
    else
        curve(:,2) = curve(tNorm,1)*curve2(:)/curve2(tNorm);
    end
else
    strTitle2 = 'data curve';
end

figure(af);subplot(1,2,2);
ax2 = gca;
%titleStr = ax2.Title.String;
plot(curve,'.-');title(strTitle2);
axis([0 size(recon,3) ylim]);grid on;
%title(ax2,titleStr);
title(ax,strTitle);
