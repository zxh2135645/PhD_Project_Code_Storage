%Orthonormalize U (introduces image weighting into Phi_rt, so that binning
%it is identical to directly binning the real-time images)
if exist('U','var')
    
  %xg 09012020, to apply a mask on U for better binning results  
%   roi_mask = [];
%   slice_mask = [];
%   for i = 1:Nz
%       sample_image = ifftshift(reshape(U_init,Ny,Nx,Nz,[]),1);
%       maxV = max(sample_image(:));
%       imshow(abs(sample_image(:,:,i,1))/abs(maxV));
%       roi = impoly;
%       keyboard;
%       slice_mask = createMask(roi);
%       slice_mask = ifftshift(roi_mask,1);
%       roi_mask(:,:,i) = slice_mask;
% %       roi_mask = repmat(roi_mask, 1,1, Nz);
%       roi_mask = vec(roi_mask);
%   end
%   %xg 09012020
  Wti=reshape(U_init,[],L);
% 
%   Wti=Wti(roi_mask,:);
  Wti=Wti'*Wti;
  Wti = sqrtm(Wti); %actually inverse of Wt
  U_init=vec(reshape(U_init,[],L)/Wti);
  U=U_init;
else
  Wti = inv(S_rt(1:L,1:L)); %untested
  U=reshape(nav_data(:,:).'*pinv(Phi_rt_small)/Wti,[],Ncoils,L);
  U=padarray(complex(imresize(real(U),[Norig Norig],'nearest'),imresize(imag(U),[Norig Norig],'nearest')),[ovs/2 ovs/2 0]);
  U=vec(repmat(reshape(U,Ny,Nx,1,L),[1 1 Nz 1]));
  U_init=U;
end
Phi_rt_init=Wti*Phi_rt_init;
Phi_rt=Phi_rt_init;
Phi_rt_full_init=Wti*Phi_rt_full_init;
Phi_rt_full=Phi_rt_full_init;
Phi_rt_small_init=Wti*Phi_rt_small_init;
Phi_rt_small=Phi_rt_small_init;

if strcmp(ScanType,'T2IR') %temporarily regenerate Bloch subspace
  Nseg=Nseg*5; %5 different preps
  params.lSegments=params.lSegments*5; %5 different preps
  [curveU,curveS]=gen_curve_subspaceT1T2(Nseg/5,params.lEchoSpacing,alpha_deg,180,(params.alTR_seconds-params.lEchoSpacing*Nseg)/2);
  curvePhi = curveU(:,1:cL);
end

if rbins > 1 && resp == 1
  respml;
else
  Ridx = ones(size(Phi_rt_small,2),1);
end

if cbins > 1 && card == 1
  if strcmp(ScanType,'SR')
    cardiac_binning;
  else
    cardml;
  end
else
  Hidx = ones(size(Ridx));
end

switch ScanType
  case 'T2IR' %replace T1T2 subspace with T1 subspace (we will learn the T2 subspace)
    Nseg=Nseg/5;
    params.lSegments=params.lSegments/5;
    [curveU,curveS]=gen_curve_subspace(Nseg,params.lEchoSpacing,alpha_deg,alpha0_deg,(params.alTR_seconds-params.lEchoSpacing*Nseg)/2);
    curvePhi = curveU(:,1:cL);
    wall_clock=repmat(vec(repmat(1:5,[Nseg/2 1])),[numel(Hidx)/(Nseg/2*5) 1]);
    wall_clock=wall_clock(1:numel(Hidx));
  case 'Cine'
    wall_clock=1+cumsum(diff(Hidx)==(1-cbins));
    wall_clock(end+1)=wall_clock(end);
    wall_clock(wall_clock==wall_clock(end))=wall_clock(end)-1;
    curvePhi=[];
    Segidx(:)=1;
  case {'IR','SR'}
    wall_clock=1+cumsum(diff(Hidx)==(1-cbins));
    wall_clock(end+1)=wall_clock(end);
    wall_clock(wall_clock==wall_clock(end))=wall_clock(end)-1;
  case 'T2prep'
    wall_clock=1+cumsum(diff(Hidx)==(1-cbins));
    wall_clock(end+1)=wall_clock(end);
    wall_clock=ceil(wall_clock/5);
    wall_clock(wall_clock==wall_clock(end))=wall_clock(end)-1;
  case 'T2star'
%     Ridx=repmat((1:params.NEco).',[1,numel(Ridx)/params.NEco]);
%     Ridx=Ridx(:);
%     Segidx=vec(repmat(nav_indices(1:NnavsPerBlock),[params.NEco 1]));
%     Segidx=repmat(Segidx,[numel(Ridx)/numel(Segidx), 1]);
% %     wall_clock=ceil((1:numel(Ridx))/(2*params.NEco*NnavsPerBlock)).'; %every 2 SR periods (~1 sec)
% %     wall_clock=ceil((1:numel(Ridx))/(params.NEco*NnavsPerBlock)).'; %every SR period (~0.5 sec)
    
    % Are this 8 hard-coded? Number of Echoes?
    Ridx = reshape(repmat(Ridx,1,8).',[],1);
    
    % TODO
    Segidx = Ridx; % 
    Segidx(:) = 1;
    Hidx = reshape(repmat(Hidx,1,8).',[],1);
    wall_clock = vec(repmat((1:params.NEco).',[1,numel(Ridx)/params.NEco]));
end

figure,subplot(2,1,1),plot(Hidx)
subplot(2,1,2),plot(Ridx)

% clear temp;
% for j=1:rbins
%   temp(:,:,j)=abs(reshape(reshape(dispim(reshape(U_init,Ny,Nx,Nz,[])),[],L_init)...
%     *mean(Phi_rt_small_init(:,Ridx==j),2),Norig,Norig,[]));
% end
% implay(abs(temp)/cw)
% % implay(abs(dispim(reshape(reshape(U_init,[],L_init)*Phi_rt_init(:,Ridx==1),Ny,Nx,[])))/cw)
% 
% 
% clear temp;
% for j=1:cbins
%   temp(:,:,j)=abs(reshape(reshape(dispim(reshape(U_init,Ny,Nx,Nz,[])),[],L_init)...
%     *mean(Phi_rt_small_init(:,Hidx==j & Ridx==1),2),Norig,Norig,[]));
% end
% implay(abs(temp)/cw)
% implay(abs(dispim(reshape(reshape(U_init,[],L_init)*Phi_rt_init(:,(Ridx(:)==1) & (Hidx(:)==1)),Ny,Nx,[])))/cw)
