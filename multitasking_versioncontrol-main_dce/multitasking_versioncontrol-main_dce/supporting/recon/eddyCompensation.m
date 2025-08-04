function [kspaceData, navData] = eddyCompensation(params,kspaceData,navData,thetas,parOrder)
% Eddy current compensation

phasestd = @(SNR) 2*pi/sqrt(12)./(1+0.6834*SNR.^1.4744);  % approx. phase stddev as function of SNR

if ~exist(strcat(params.fileString,'_eddy_compensation'),'dir')
    mkdir(strcat(params.fileString,'_eddy_compensation'));
end

centerPars = (parOrder == params.DC_kz);
if strcmp(params.ScanType,'CEST')
    centerPars(1:floor(2*numel(centerPars)/3)) = false;
end

bipolarmode = min(params.ReadoutMode,params.Necho);

if params.isCartesian
    [~,DC_kx] = max(sum(abs(navData(:,:,1,1)),1));
else
    DC_kx = params.DC;
end

navIndices = params.navShift:params.SGBlock:params.Ntpoint;
% Segidx_nav = mod(navIndices,params.linesPerShot);
% Segidx_nav = logical(Segidx_nav>params.linesPerShot/2);

kspaceIndices = 1:params.Ntpoint;
kspaceIndices(navIndices) = [];
Segidx_k = mod(kspaceIndices,params.linesPerShot);
Segidx_k = logical(Segidx_k>params.linesPerShot/2);
centerPars = logical(centerPars(:).*Segidx_k(:));

h = figure('Name','Eddy Current Compensation');
for coil = 1:params.Ncoils
    for echo = 1:bipolarmode
%         phasenav = angle(mean(navData(Segidx_nav,DC_kx,1,coil,1)));
        msdev = std2(kspaceData(:,[1 end],1,coil,echo)); % noise std dev
        temp  = double(kspaceData(centerPars,DC_kx,1,coil,echo));
        thetas_center = reshape(thetas(centerPars),[],1);
        figure(h);subplot(2,1,1);plot(thetas_center,angle(temp),'.');ylim([-pi, pi]),xlabel('Spoke angle');ylabel('DC phase (original)');title(sprintf('Coil #%02d, echo #%d', coil, echo));

        % phase weights
        pw = 1./phasestd(abs(temp)/msdev);
        pw = double(min(pw/median(pw),sqrt(10)));

        % coarse search
        x1N = 15;
        xN  = 21;
        x1s = linspace(-pi,pi,x1N);
        xs  = linspace(-2*pi,2*pi,xN);
        costs = zeros(xN,xN);
        for j = 1:xN
            for k = 1:xN
                for l = 1:x1N
                    costs(j,k,l) = norm(pw.*angle(temp.*exp(1i*(x1s(l)+xs(j)*sin(thetas_center*pi/180)+xs(k)*cos(thetas_center*pi/180)))));
                end
            end
        end
        [~,argmin] = min(costs(:));
        [j,k,l]    = ind2sub([xN xN x1N],argmin);

        options = optimoptions('lsqnonlin','Display','off');
        temp = lsqnonlin(@(x)pw.*angle(temp.*exp(1i*(x(1)+x(2)*sin(thetas_center*pi/180)+x(3)*cos(thetas_center*pi/180)))),[x1s(l),xs(j),xs(k)],[-pi, -2*pi, -2*pi],[pi, 2*pi, 2*pi],options);

        %kspaceData(:,:,:,coil,echo:bipolarmode:end) = bsxfun(@times,kspaceData(:,:,:,coil,echo:bipolarmode:end),exp(1i*(temp(1)+temp(2)*sin(thetas*pi/180)+temp(3)*cos(thetas*pi/180))));
        kspaceData(:,:,:,coil,echo:bipolarmode:end) = bsxfun(@times,kspaceData(:,:,:,coil,echo:bipolarmode:end),exp(1i*(temp(2)*sin(thetas*pi/180)+temp(3)*cos(thetas*pi/180)))); % don't use the constant term
        figure(h),subplot(2,1,2),plot(thetas_center,angle(kspaceData(centerPars,DC_kx,1,coil,echo)),'.'),ylim([-pi, pi]),xlabel('Spoke angle'),ylabel('DC phase (corrected)'),drawnow;
        saveas(h,fullfile(strcat(params.fileString,'_eddy_compensation'),sprintf('coil%02d_echo%d.png',coil,echo)));
        
        if echo <= size(navData,5)
            %navData(:,:,:,coil,echo:bipolarmode:end) = bsxfun(@times,navData(:,:,:,coil,echo:bipolarmode:end),exp(1i*(temp(1)+temp(3))));
            navData(:,:,:,coil,echo:bipolarmode:end) = bsxfun(@times,navData(:,:,:,coil,echo:bipolarmode:end),exp(1i*(temp(3))));   % don't use the constant term
        end
%         phasek = angle(mean(kspaceData(centerk,DC_kx,1,coil,1)));
%         fprintf('Coil #%d, echo #%d, phasenav1 = %.4f, phasenav2 = %.4f, phasek = %.4f, temp1 = %.4f, temp2 = %.4f, temp3 = %.4f.\n ',coil,echo,phasenav,angle(mean(navData(Segidx_nav,DC_kx,1,coil,1))),phasek,temp(1),temp(2),temp(3));
    end
end
close(h);
clear costs j k l pw coil x1N xN xs argmin temp h thetas_center centerPars phasestd msdev

