switch ScanType
  case {'IR','Cine','T2prep'}
    WallClock=false;
    wall_clock(:)=1;
  case {'T2IR','SR','T2star'}
    WallClock=true;
end

Nseg = 192;
[Navdata_tensor,mask]=regate(navdata_temp,Segidx,Hidx,Ridx,wall_clock);
if size(Navdata_tensor,2)<Nseg
  Navdata_tensor(:,Nseg,:,:)=0;
  mask(:,Nseg,:,:)=0;
end
sizes=size(Navdata_tensor);
if ~WallClock
  sizes((end+1):5)=1;
end
morozov=msdev^2*norm(double(mask(:)))^2

% doBloch=~isempty(curvePhi)
% doSpline=sum(abs(sms))~=0
doBloch = 1;
doSpline = 0;
% isUndersampled=~exp(sum(log(mask(:)))); %~prod(mask(:));

% initial guess
if ~doSpline && ~doBloch%if ~doSpline && isUndersampled && ~doBloch
  tempdim=3; %find(sizes~=1,1,'last');
  tempsizes=ones(size(sizes));
  tempsizes(tempdim)=sizes(tempdim); % TODO
  Navdata_sm=repmat(sum(Navdata_tensor.*mask,tempdim)./sum(mask,tempdim),tempsizes);
  Navdata_sm(logical(mask))=Navdata_tensor(logical(mask));
  Navdata_sm(isnan(Navdata_sm))=0;
else
  Navdata_sm=Navdata_tensor;
end
if doBloch
  Navdata_bloch=pcg(@(x)vec((reshape(permute(mask.^2,[1 3 4 5 2]),[],Nseg).*(reshape(x,[],cL)*curvePhi.'))*curvePhi),...
    vec(reshape(permute(mask.^2.*Navdata_tensor,[1 3 4 5 2]),[],Nseg)*curvePhi),[],200);
  Navdata_bloch=ipermute(reshape(reshape(Navdata_bloch,[],cL)*curvePhi.',sizes([1 3 4 5 2])),[1 3 4 5 2]);
  %     morozov_new=morozov-norm((Navdata_bloch(:)-Navdata_tensor(:)).*mask(:))^2
end

if true
     Navdata_sm=Navdata_bloch;
else
    if doBloch && doSpline
        if strcmp(ScanType,'T2prep')
            Navdata_sm=repmat(sum(Navdata_bloch.*mask,5)./sum(mask,5),[1 1 1 1 sizes(5)]);
            Navdata_sm(logical(mask))=Navdata_tensor(logical(mask));
        else
            Navdata_sm=LRTCp_splines(Navdata_bloch.*logical(mask),0,sms,mask,5,Navdata_bloch.*logical(mask),[],smorders,circs,ts_proj);
        end
        Navdata_sm=LRTCp_splines(Navdata_tensor,lr,sms,mask,20,Navdata_sm,curvePhi,smorders,circs,ts_proj);
        %     morozov_new
    elseif doSpline
        Navdata_sm=LRTCp_splines(Navdata_tensor,0,sms,mask,5,Navdata_tensor.*logical(mask),[],smorders,circs,ts_proj);
        Navdata_sm=LRTCp_splines(Navdata_tensor,lr,sms,mask,20,Navdata_sm,[],smorders,circs,ts_proj);
    elseif doBloch && isUndersampled
        Navdata_sm=LRTCp(Navdata_tensor,lr,mask,20,Navdata_bloch.*logical(mask),curvePhi);
    elseif doBloch && ~isUndersampled
        Navdata_sm=Navdata_bloch;
    elseif isUndersampled
        Navdata_sm=LRTCp(Navdata_tensor,lr,mask,20,Navdata_sm);
    else
        Navdata_sm=Navdata_tensor;
    end
end
% morozov
% TODO
% cL shouldn't be 1
% cL = 96;
ranks=[16, cL,cbins,rbins, params.NEco]; % Default ranks: # of basis ims, # of T1 recovery basis functions, # of T2* decay basis functions, # of "wall clock" basis functions
%ranks(ranks==1)=[];
% [C,UU,ranks]=choose_C(squeeze(Navdata_sm),ranks,squeeze(Navdata_tensor),squeeze(mask));
[C,UU,ranks]=choose_C(squeeze(Navdata_sm),ranks);
ranks

if doBloch %reinforce dictionary subspace
    UU=reshape(curvePhi*(pinv(curvePhi)*reshape(UU,Nseg,[])),size(UU));
end

Phi=C*UU';
clear UU;
L=size(C,1);
Phi=reshape(Phi,[L sizes(2:end)]);


% figure,imagesc(squeeze(real(Phi(1,:,:))))

U_sm=pcg(@(x)vec(reshape(vec(reshape(x,[],L)*Phi(:,:)).*mask(:).^2,[],prod(sizes(2:end)))*Phi(:,:)'),vec(reshape(Navdata_tensor(:).*mask(:).^2,[],prod(sizes(2:end)))*Phi(:,:)'),[],20);
temp=reshape(U_sm,[],L)*Phi(:,:);
if ~isempty(ts_proj)
  figure,imagesc(abs(conj(ts_proj)*temp))
else
  figure,imagesc(abs(temp))
end
[morozov norm((temp-Navdata_tensor(:,:)).*mask(:,:),'fro')^2]

% interpolate full timings if necessary.
Phi=reshape(Phi,[L sizes(2:end)]);

Hidx_full=interp1(nav_indices,Hidx(:),1:Nread,'nearest','extrap');

Ridx_full=interp1(nav_indices,Ridx(:),1:Nread,'nearest','extrap');

%xg
Hidx_full = circshift(Hidx_full,1);
Hidx_full(1,1) = Hidx_full(1,2);

Ridx_full = circshift(Ridx_full,1);
Ridx_full(1,1) = Ridx_full(1,2);


%xg

% Segidx_full=repmat(reshape(1:Nseg,SGblock,[]),[params.NEco, 1]);
% Segidx_full=repmat(Segidx_full(:),[numel(Hidx_full)/numel(Segidx_full), 1]);

% TODO navline/SGblock
% Segidx_full = Hidx_full;
% Segidx_full(:) = 1;
Segidx_full = interp1(nav_indices, Segidx(:), 1:Nread, 'nearest', 'extrap');
Segidx_full = circshift(Segidx_full, 1);
Segidx_full(1,1) = Segidx_full(1,2);
% interp1(nav_indices,Ridx(:),1:Nread,'nearest','extrap');
wall_clock_full=interp1(nav_indices,wall_clock(:),1:Nread,'nearest','extrap');
wall_clock_full = circshift(wall_clock_full,1);
wall_clock_full(1,1) = wall_clock_full(1,2);
% wall_clock_full = repmat(wall_clock',1,2);
% Ridx_full=repmat((1:params.NEco),[SGblock,1,numel(Segidx_full)/SGblock/params.NEco]);
% Ridx_full=Ridx_full(:);

Phi_rt_full = degate(Phi(:,:).',[size(C,1) sizes(2:end)],Segidx_full,Hidx_full,Ridx_full,wall_clock_full); %(C*UU').'
Phi_rt_small=Phi_rt_full(:,nav_indices);
Phi_rt=Phi_rt_full;

if ~iscartesian
  Phi_rt(:,nav_indices)=[];
else
  Phi_rt(:,nav_indices)=0;
end

% [Phi_rt,Gr1,Gr2] = svd(Phi_rt','econ');
% Phi_rt = Phi_rt';
% Gr = Gr2*Gr1;



lrw=sqrt(sum(abs(nav_data(:,:).'-(nav_data(:,:).'*pinv(Phi_rt_small))*Phi_rt_small).^2));
lrw=interp1(nav_indices,lrw,1:Nread,'linear','extrap');
h=fminsearch(@(h)abs(norm(exp(-lrw.^2/h))^2/sqrt(numel(lrw))/norm(exp(-lrw.^2/h).^2)-sqrt(.95)),median(lrw).^2)
figure,plot(median(lrw)./lrw,exp(-lrw.^2/h)/exp(-median(lrw).^2/h),'.');
lrw=exp(-lrw.^2/h);
median(lrw)
lrw=lrw/median(lrw);

if ~iscartesian
    lrw(nav_indices)=[];
else
    lrw(nav_indices)=0;
end

[Phi_rt,Gr1,Gr2] = svd(bsxfun(@times,Phi_rt,lrw)','econ');
Phi_rt = Phi_rt';
Gr = Gr2*Gr1;

Phi_rt_full=Gr\Phi_rt_full;
Phi_rt_small=Gr\Phi_rt_small;

kspace_data=bsxfun(@times,kspace_data,lrw(:));
