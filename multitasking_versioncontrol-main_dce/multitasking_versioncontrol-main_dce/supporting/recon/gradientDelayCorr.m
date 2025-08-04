function [kspaceData,navData] = gradientDelayCorr(kspaceData,navData,Ntrajs,Norig,linOrder,parOrder,flagGradientDelayCorr,isVE,DC_kz)
%% Gradient delay correction

N = size(kspaceData,2);
Ncoils = size(kspaceData,4);
thetas = (linOrder-1) * 180 / Ntrajs;

fprintf("total %d echoes. Echo ",size(kspaceData,5));
for echo = 1:size(kspaceData,5)
    fprintf("#%d. ", echo);
    kspaceTemp = kspaceData(:,:,:,:,echo);
    if flagGradientDelayCorr == 1
        cand_coords = unique([linOrder, parOrder],'rows');    % set of (theta,kz)-coordinates
        cand_coords = intersect(cand_coords(cand_coords(:,1)<=Ntrajs,:),bsxfun(@minus,cand_coords(cand_coords(:,1)>Ntrajs,:),[Ntrajs 0]),'rows'); % set of coordinates which also have opposed pairings
        ksig = zeros(size(cand_coords,1),size(kspaceTemp,2),size(kspaceTemp,4));
        nsig = ksig;
        for j = 1:size(cand_coords,1)
            t_ind = (linOrder==(cand_coords(j,1)+Ntrajs)) & (parOrder==cand_coords(j,2));
            ksig(j,:,:) = mean(ifft(kspaceTemp(t_ind,end:-1:1,:),[],2),1);  % rather than mean, rank-1 correction would be better
            t_ind = (linOrder==cand_coords(j,1)) & (parOrder==cand_coords(j,2));
            nsig(j,:,:) = mean(ifft(kspaceTemp(t_ind,:,:),[],2),1);
        end
        ksig    = abs(fft(ksig,[],2)); % actually only use k-space magnitude
        x       = [0:(ceil(N/2)-1) (-floor(N/2)):-1];
        cost    = @(shift) norm(reshape(ksig-abs(fft(nsig.*repmat(exp(1i*(2*shift-pi/Norig)*x),[size(ksig,1) 1 Ncoils]),[],2)),[],1));
        shifts  = linspace(-.1,.1,201);
        mincost = inf;
        for j = 1:numel(shifts)
            tempcost = cost(shifts(j));
            if tempcost < mincost
                mincost = tempcost;
                shift   = shifts(j);
            end
        end
        shift = fminsearch(cost,shift);

        kspaceData(:,:,:,:,echo) = fft(ifft(kspaceTemp,[],2).*repmat(exp(1i*shift*x), [size(kspaceTemp,1) 1 size(kspaceTemp,3) Ncoils]),[],2);
        if size(navData,5) >= echo
            navData(:,:,:,:,echo)    = fft(ifft(navData(:,:,:,:,echo),[],2).*repmat(exp(1i*shift*x), [size(navData,1) 1 size(navData,3) Ncoils]),[],2);
        end
    elseif flagGradientDelayCorr == 2
        % New setup for different versions
        tic;
        if ~isVE
            % calculate initial angle
            dk = kspaceData(:,:,:,:,1) - navData(:,:,:,:,1);
            mk = kspaceData(:,:,:,:,1) + navData(:,:,:,:,1);
            sim = 4*sum(abs(dk(:,:)).^2,2)./sum(abs(mk(:,:)).^2,2);
            clear dk mk
        %     if strcmp(ScanType,'SR')
        %         sim(:)=1; %Doesn't seem to help for perfusion scans
        %     end

            %find 0 degree line
            simlast = floor(numel(sim)/Ntrajs)*Ntrajs;
            klast   = floor(size(kspaceTemp,1)/Ntrajs)*Ntrajs;
            temp    = squeeze(mean(reshape(kspaceTemp(1:klast,:,1,:),Ntrajs,[],N,Ncoils),2));
            err     = zeros(Ntrajs,1);
            parfor j = 1:Ntrajs
                [~,sorter] = sort(mod(theta*mod((1:Ntrajs)-j,Ntrajs),360));
                err(j)     = norm(reshape(diff(temp(sorter,:,:),1),[],1));
            end
            [~,theta0] = min(err.*mean(reshape(sim(1:simlast),Ntrajs,[]),2));
            %         figure,plot(1./err./mean(reshape(sim(1:simlast),Ntrajs,[]),2))

            %setup gradient delay correction
            oneeighties = abs(thetas-180)==min(abs(thetas-180));
            ksig = ifft(kspaceTemp(oneeighties,end:-1:1,:),[],2);
            nsig = ifft(navData(oneeighties,:,:),[],2);
            %       nsig = ifftshift(padarray(ifft(ifftshift(fft(navData(oneeighties,:,:),[],2),2),[],2),[0 ovs/2 0]),2)/sqrt(N);
        else
            %setup gradient delay correction
            cand_coords = unique(linOrder);
            cand_coords(cand_coords>Ntrajs) = [];  % keep only thetas<180
            success  = false;
            attempts = 0;
            while ~success && (attempts<3)
                attempts = attempts+1;
                ksig     = zeros(size(cand_coords,2),size(kspaceTemp,2),size(kspaceTemp,4));
                nsig     = ksig;
                nocands  = [];
                for j = 1:size(cand_coords,1)
                    t_ind = (linOrder==cand_coords(j)+Ntrajs);   % pick exact bipolar readout
                    if sum(t_ind)==0
                        nocands(end+1) = j;
                    else
                        ksig(j,:,:) = mean(ifft(kspaceTemp(t_ind,end:-1:1,:),[],2),1); % rather than mean, rank-1 correction would be better
                        t_ind = (linOrder==cand_coords(j));
                        nsig(j,:,:) = mean(ifft(kspaceTemp(t_ind,:,:),[],2),1);
                    end
                end          
                used_thetas = (cand_coords(:)-1)*180.0/Ntrajs;
                used_thetas(nocands) = [];
                ksig(nocands,:,:)    = [];
                nsig(nocands,:,:)    = [];
                success = ~isempty(ksig);
            end
        end

        % Gradient delay correction
        [~,kmax] = max(sqrt(sum(sum(abs(fft(ksig,[],2)).^2,1),3)));
        [~,nmax] = max(sqrt(sum(sum(abs(fft(nsig,[],2)).^2,1),3)));
        shift1 = (kmax-nmax+1)/2*pi/Norig;
        ksig = abs(fft(ksig,[],2)); %actually only use k-space magnitude

        x = [0:(Norig-1) (-Norig):-1];
        if isVE
            shiftfun=@(shift,thetas)shift(1)./sqrt(1-(shift(2)*cos(thetas(:)*pi/180+shift(3))).^2);
            cost=@(shift)norm(reshape(ksig-abs(fft(bsxfun(@times,nsig,exp(1i*(2*shiftfun(shift,used_thetas)-pi/Norig)*x)),[],2)),[],1));
            shifts1 = shift1 + pi/Norig*(-1:.5:1);
            shifts2 = linspace(0,3/4,11); shifts2(end)=[];
            shifts3 = linspace(-pi/2,pi/2,11); shifts3(end)=[];
            mincost = inf;
            for i = 1:numel(shifts1)
                for j = 1:numel(shifts2)
                    for k = 1:numel(shifts3)
                        tempcost = cost([shifts1(i),shifts2(j),shifts3(k)]);
                        if tempcost < mincost
                            mincost = tempcost;
                            shift   = [shifts1(i),shifts2(j),shifts3(k)];
                        end
                    end
                end
            end

            shiftfun=@(shift,thetas)shift(1)./sqrt(1-(shift(2)*cos(thetas(:)*pi/180+shift(3))).^2);
            cost=@(shift)abs(reshape(ksig-abs(fft(bsxfun(@times,nsig,exp(1i*(2*shiftfun(shift,used_thetas)-pi/Norig)*x)),[],2)),[],1));
            shift    = lsqnonlin(@(x)cost(x)/max(cost([0 0 0])),shift,[shift1-2*pi/Norig 0 -pi/2],[shift1+2*pi/Norig 3/4 pi/2]);

            rad=shiftfun(shift,(-180:179).');
            figure,plot(rad.*cosd(-180:179).',rad.*sind(-180:179).'),axis equal
            figure,imagesc(abs(ksig(:,:,1)));
            figure,imagesc(abs(fft(nsig(:,:,1),[],2)));
            figure,imagesc(abs(fft(bsxfun(@times,nsig(:,:,1),exp(1i*(2*shiftfun(shift,used_thetas)-pi/Norig)*x)),[],2)));

            kspaceData(:,:,:,:,echo) = fft(bsxfun(@times,ifft(kspaceTemp,[],2),exp(1i*shiftfun(shift,thetas)*x)),[],2);
            if size(navData,5) >= echo
                navData(:,:,:,:,echo)    = fft(bsxfun(@times,ifft(navData(:,:,:,:,echo),[],2),exp(1i*shiftfun(shift,0)*x)),[],2); %doesn't really matter
            end
        else
            cost  = @(shift)norm(reshape(ksig-bsxfun(@times,nsig,exp(1i*(2*shift-pi/Norig)*x)),[],1));
            shift = fminsearch(cost,shift1);
            kspaceData(:,:,:,:,echo) = fft(bsxfun(@times,ifft(kspaceTemp,[],2),exp(1i*shift*x)),[],2);
            if (size(navData,5) > 1) || (echo == 1)
                navData(:,:,:,:,echo)    = fft(bsxfun(@times,ifft(navData(:,:,:,:,echo),[],2),exp(1i*shift*x)),[],2); %doesn't really matter
            end
        end
        toc;
    elseif flagGradientDelayCorr == 3  % this section is for code testing 
        tic;
        %setup gradient delay correction
        cand_coords = unique(linOrder);
        cand_coords(cand_coords>Ntrajs) = [];  % keep only thetas<180
        success  = false;
        attempts = 0;
        while ~success && (attempts<3)
            attempts = attempts+1;
            ksig     = zeros(size(cand_coords,2),size(kspaceTemp,2),size(kspaceTemp,4));
            nsig     = ksig;
            nocands  = [];
            for j = 1:size(cand_coords,1)
                t_ind = (linOrder==cand_coords(j)+Ntrajs);   % pick exact bipolar readout
                if sum(t_ind)==0
                    nocands(end+1) = j;
                else
                    ksig(j,:,:) = mean(ifft(kspaceTemp(t_ind,end:-1:1,:),[],2),1); % rather than mean, rank-1 correction would be better
                    t_ind = (linOrder==cand_coords(j));
                    nsig(j,:,:) = mean(ifft(kspaceTemp(t_ind,:,:),[],2),1);
                end
            end          
            used_thetas = (cand_coords(:)-1)*180.0/Ntrajs;
            used_thetas(nocands) = [];
            ksig(nocands,:,:)    = [];
            nsig(nocands,:,:)    = [];
            success = ~isempty(ksig);
        end
            
        % Gradient delay correction
        [~,kmax] = max(sqrt(sum(sum(abs(fft(ksig,[],2)).^2,1),3)));
        [~,nmax] = max(sqrt(sum(sum(abs(fft(nsig,[],2)).^2,1),3)));
        shift1 = (kmax-nmax+1)/2*pi/Norig;
        
        ksig = (fft(ksig,[],2));
        x = [0:(Norig-1) (-Norig):-1];

        ksigtemp = ksig; nsigtemp = nsig;
        
        for n = 1:size(ksig,3)
            n
            ksig = ksigtemp(:,:,n);
            nsig = nsigtemp(:,:,n);
            anglefun = @(a) a - (a>2*pi)*2*pi + (a<-2*pi)*2*pi; 
            shiftfun = @(shift,thetas) shift(1)./sqrt(1-(shift(2)*cos(thetas(:)*pi/180+shift(3))).^2);
            cost     = @(shift) norm(anglefun(reshape(angle(ksig)-angle(fft(bsxfun(@times,nsig,exp(1i*(2*shiftfun(shift,used_thetas)-pi/Norig)*x)),[],2)),[],1)));
            shifts1 = shift1 + pi/Norig*(-1:.5:1);
            shifts2 = linspace(0,3/4,11); shifts2(end)=[];
            shifts3 = linspace(-pi/2,pi/2,11); shifts3(end)=[];
            mincost = inf;
            for i = 1:numel(shifts1)
                for j = 1:numel(shifts2)
                    for k = 1:numel(shifts3)
                        tempcost = cost([shifts1(i),shifts2(j),shifts3(k)]);
                        if tempcost < mincost
                            mincost = tempcost;
                            shift   = [shifts1(i),shifts2(j),shifts3(k)];
                        end
                    end
                end
            end

            cost  = @(shift) abs(anglefun(reshape(angle(ksig)-angle(fft(bsxfun(@times,nsig,exp(1i*(2*shiftfun(shift,used_thetas)-pi/Norig)*x)),[],2)),[],1)));
            shift = lsqnonlin(@(x)cost(x)/max(cost([0 0 0])),shift,[shift1-2*pi/Norig 0 -pi/2],[shift1+2*pi/Norig 3/4 pi/2]);

            rad = shiftfun(shift,(-180:179).');
            figure(1),plot(rad.*cosd(-180:179).',rad.*sind(-180:179).'),axis equal;
            figure(2),imagesc(abs(ksig));
            figure(3),imagesc(abs(fft(nsig,[],2)));
            figure(4),imagesc(abs(fft(bsxfun(@times,nsig,exp(1i*(2*shiftfun(shift,used_thetas))*x)),[],2)));

            kspaceData(:,:,:,n,echo) = fft(bsxfun(@times,ifft(kspaceTemp(:,:,:,n),[],2),exp(1i*shiftfun(shift,thetas)*x)),[],2);
            if size(navData,5) >= echo
                navData(:,:,:,n,echo) = fft(bsxfun(@times,ifft(navData(:,:,:,n,echo),[],2),exp(1i*shiftfun(shift,0)*x)),[],2); %doesn't really matter
            end
        end
        toc;
    elseif flagGradientDelayCorr == 9   % this section is for code testing 
        ksig = zeros(Ntrajs,size(kspaceTemp,2),size(kspaceTemp,4));
        nsig = ksig;
        for j = 1:Ntrajs
            t_ind = (linOrder==(j+Ntrajs)) & (parOrder==DC_kz);
            if sum(t_ind) > 0
                ksig(j,:,:) = mean(ifft(kspaceTemp(t_ind,end:-1:1,:),[],2),1);  % rather than mean, rank-1 correction would be better
                t_ind = (linOrder==j) & (parOrder==DC_kz);
                if sum(t_ind) > 0
                    nsig(j,:,:) = mean(ifft(kspaceTemp(t_ind,:,:),[],2),1);
                else
                    ksig(j,:,:) = 0;
                end
            end
        end
        ksig    = abs(fft(ksig,[],2)); % actually only use k-space magnitude
        x       = [0:(ceil(N/2)-1) (-floor(N/2)):-1];
        cost    = @(shift) norm(reshape(ksig-abs(fft(nsig.*repmat(exp(1i*(2*shift-pi/Norig)*x),[size(ksig,1) 1 Ncoils]),[],2)),[],1));
        shifts  = linspace(-.1,.1,201);
        mincost = inf;
        for j = 1:numel(shifts)
            tempcost = cost(shifts(j));
            if tempcost < mincost
                mincost = tempcost;
                shift   = shifts(j);
            end
        end
        shift = fminsearch(cost,shift);

        kspaceData(:,:,:,:,echo) = fft(ifft(kspaceTemp,[],2).*repmat(exp(1i*shift*x), [size(kspaceTemp,1) 1 size(kspaceTemp,3) Ncoils]),[],2);
        if size(navData,5) >= echo
            navData(:,:,:,:,echo)    = fft(ifft(navData(:,:,:,:,echo),[],2).*repmat(exp(1i*shift*x), [size(navData,1) 1 size(navData,3) Ncoils]),[],2);
        end
    end
end

