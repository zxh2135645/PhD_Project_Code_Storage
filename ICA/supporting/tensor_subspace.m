switch ScanType
  case {'IR','Cine'}
    WallClock=false;
    wall_clock(:)=1;
  case {'T2IR','SR','T2prep'}
    WallClock=true;
end

[Navdata_tensor,mask]=regate(navdata_temp,Segidx,Hidx,Ridx,wall_clock);
sizes=size(Navdata_tensor);
if ~WallClock
    sizes((end+1):5)=1;
end
morozov=msdev^2*norm(double(mask(:)))^2

doBloch=~isempty(curvePhi)
doSpline=sum(abs(sms))~=0

% initial guess
if ~doSpline
    tempdim=3; %find(sizes~=1,1,'last');
    tempsizes=ones(size(sizes));
    tempsizes(tempdim)=sizes(tempdim);
    Navdata_sm=repmat(sum(Navdata_tensor.*mask,tempdim)./sum(mask,tempdim),tempsizes);
    Navdata_sm(logical(mask))=Navdata_tensor(logical(mask));
    Navdata_sm(isnan(Navdata_sm))=0;
else
    Navdata_sm=Navdata_tensor;
end
if doBloch
    Navdata_bloch=pcg(@(x)vec((reshape(permute(mask.^2,[1 3 4 5 2]),[],Nseg/2).*(reshape(x,[],cL)*curvePhi.'))*curvePhi),vec(reshape(permute(mask.^2.*Navdata_tensor,[1 3 4 5 2]),[],Nseg/2)*curvePhi),[],200);
    Navdata_bloch=ipermute(reshape(reshape(Navdata_bloch,[],cL)*curvePhi.',sizes([1 3 4 5 2])),[1 3 4 5 2]);
%     morozov_new=morozov-norm((Navdata_bloch(:)-Navdata_tensor(:)).*mask(:))^2
end
%Low Rank Tensor Compeletion
if doBloch && doSpline
  if strcmp(ScanType,'T2prep')
    Navdata_sm=repmat(sum(Navdata_bloch.*mask,5)./sum(mask,5),[1 1 1 1 sizes(5)]);
    Navdata_sm(logical(mask))=Navdata_tensor(logical(mask));
  else
    Navdata_sm=LRTCp_splines(Navdata_bloch.*logical(mask),0,sms,mask,5,Navdata_bloch.*logical(mask),[],smorders,circs);
  end
    Navdata_sm=LRTCp_splines(Navdata_tensor,lr,sms,mask,20,Navdata_sm,curvePhi,smorders,circs);
%     morozov_new
elseif doSpline
    Navdata_sm=LRTCp_splines(Navdata_tensor,0,sms,mask,5,Navdata_tensor.*logical(mask),[],smorders,circs);
    Navdata_sm=LRTCp_splines(Navdata_tensor,lr,sms,mask,20,Navdata_sm,[],smorders,circs);
elseif doBloch
    Navdata_sm=LRTCp(Navdata_tensor,lr,mask,20,Navdata_bloch.*logical(mask),curvePhi);
else
    Navdata_sm=LRTCp(Navdata_tensor,lr,mask,20,Navdata_sm);
end
morozov
%watch_notification(mainpath,'Choose model orders');
switch ScanType
  case {'Cine','IR'}
    ranks=[32, cL, sizes(3:end)];
  case 'T2IR'
    ranks=[42, cL, sizes(3:end)];
  case 'T2prep'
    ranks=[36, cL, sizes(3:end-1), 32];
  case 'SR'
    ranks=[36, cL, sizes(3:end)];
end
ranks(ranks==1)=[];
%[C,UU,ranks]=choose_C(squeeze(Navdata_sm),ranks,squeeze(Navdata_tensor),squeeze(mask));
[C,UU,ranks]=choose_C(squeeze(Navdata_sm),ranks);
ranks

Phi=C*UU';
L=size(C,1);
Phi=reshape(Phi,[L sizes(2:end)]);
% figure,imagesc(squeeze(real(Phi(1,:,:))))

U_sm=pcg(@(x)vec(reshape(vec(reshape(x,[],L)*Phi(:,:)).*mask(:).^2,[],prod(sizes(2:end)))*Phi(:,:)'),vec(reshape(Navdata_tensor(:).*mask(:).^2,[],prod(sizes(2:end)))*Phi(:,:)'),[],20);
temp=reshape(U_sm,[],L)*Phi(:,:);
figure,imagesc(abs(temp))

if doBloch
  % use golden-angle timings now, not navigator timings.
  if ~isVE
    UU=reshape(curveU(1:2:end,1:cL)*pinv(curvePhi)*reshape(UU,Nseg/2,[]),size(UU));
  else
    UU=reshape(curveU(2:2:end,1:cL)*pinv(curvePhi)*reshape(UU,Nseg/2,[]),size(UU));
  end
end

Phi=C*UU';
Phi=reshape(Phi,[L sizes(2:end)]);

Phi_rt = degate(Phi(:,:).',[size(C,1) sizes(2:end)],Segidx,Hidx,Ridx,wall_clock); %(C*UU').'=conj(UU)*C.';

clear UU;

[Phi_rt,Gr1,Gr2] = svd(Phi_rt','econ');
Phi_rt = Phi_rt';
Gr = Gr2*Gr1';
