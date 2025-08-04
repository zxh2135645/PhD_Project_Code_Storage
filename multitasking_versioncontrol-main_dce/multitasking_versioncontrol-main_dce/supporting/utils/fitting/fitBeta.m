function [Betamap,T1map_corr] = fitBeta(R1map,TR,FAs,FAidx)
%% Correct for beta

FAs = FAs*pi/180.0;
if nargin < 4
    FAidx = 1:numel(FAs);
end

[Ny,Nx,Nz,~] = size(R1map);
fitsResult = zeros(Ny*Nx,Nz,2);
for cphase = 1
    R1curve = @(R1,beta) R1 - (log(cos(beta*FAs(FAidx)))/TR);
    for sl = 1:Nz
        fprintf('slice: %d / %d\n', sl, Nz);
        im_mask = prod(R1map(:,:,sl,:),4)>0;
        R1temp = reshape(R1map(:,:,:,FAidx),[],length(FAs(FAidx)));
        R1temp = R1temp(im_mask(:),:);
        
        opts = [];
        opts.MaxFunEvals = 1000;
        opts.Display = 'off';
       
        x0  = double([1.0, 0.5]); % R1, beta
        xlb = double([0.1, 0.05]);
        xub = double([10, 1]);
 
        fitmat = [];
        parfor j = 1:size(R1temp,1)
            curve = double(R1temp(j,:));
            cost = @(x) double(abs(R1curve(x(1),x(2))-curve));  
            [tempfit, res] = lsqnonlin(cost, x0, xlb, xub, opts);
            fitmat(j,:) = tempfit;
        end
        fitsResult(im_mask(:),sl,1) = 1./fitmat(:,1);
        fitsResult(im_mask(:),sl,2) = fitmat(:,2);
        
        T1map_corr(:,:,sl) = reshape(fitsResult(:,sl,1),size(im_mask,1),size(im_mask,2));
        Betamap(:,:,sl)  = reshape(fitsResult(:,sl,2),size(im_mask,1),size(im_mask,2));
        
        %figure(41),imagesc(imrotate(round(1000*wmedfilt2(T1map_corr(:,:,sl,cphase))),0),[0 2000]),axis image,colormap('inferno');title('T1');colorbar;
        figure(41),imagesc(round(1000*(T1map_corr(:,:,sl,1))),[0 3000]),axis image equal tight off,colormap('inferno');title('T1 (beta-corrected)');colorbar;
        figure(43),imagesc(Betamap(:,:,sl,1),[0.0 1.0]),axis image equal tight off,colormap('jet');title('beta');colorbar;
    end
end

