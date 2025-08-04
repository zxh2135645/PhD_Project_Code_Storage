function [realtime_sig, temp_LV] = realtime_mean_sig(Params, SpatialCoeff, TemporalBasis, ReconstructedImages)

Ny=Params.Ny;
Nx=Params.Nx;
Nydisp=Params.Nydisp;
Nxdisp=Params.Nxdisp;
numFA  = 1;
Nz=Params.Nz;
slices=1:Nz;
Nseg = Params.linesPerShot;
moduleLength = Params.moduleLength;
MBfactor = Params.MBfactor;
dispSlice=1;
dispim = @(x) x(floor(Ny/2)-floor(Nydisp/2) + (1:Nydisp), floor(Nx/2)-floor(Nxdisp/2) + (1:Nxdisp),:);
%dispim = @(x) x(floor(Ny/2)-floor(Nydisp/2) + (1:Nydisp), floor(Nx/2)-floor(Nxdisp/2) + (1:Nxdisp), dispSlice, :);
U_init=SpatialCoeff.U_init;
Phi_rt_init=TemporalBasis.Phi_rt_full_init;
reconreal_cell=cell(ceil(size(Phi_rt_init,2)/(Nseg*moduleLength*numFA)),1);
%reconreal_cell=cell(1,1);
Utemp = reshape(U_init,Ny,Nx,Nz,[]);
if MBfactor == 1
    Utemp = fftshift(Utemp,1);
end
%Utemp = dispim(Utemp);
%Nd = [Nydisp Nxdisp numel(slices)];
for ii = 1:size(reconreal_cell,1)
    reconreal = reshape(Utemp,[],size(Phi_rt_init,1));
    reconreal = reconreal(:,1:size(Phi_rt_init,1));
    reconreal = (reshape(reconreal*Phi_rt_init(1:size(Phi_rt_init,1),(ii-1)*(Nseg*moduleLength*numFA)+1:min(ii*(Nseg*moduleLength*numFA),size(Phi_rt_init,2))),Ny,Nx,Nz,[]));
    %reconreal = (reshape(reconreal*Phi_rt_init(:,:),Ny,Nx,Nz,[]));
    reconreal = reconreal(:,:,1,:);
    reconreal = realify(dispim(reconreal));
    reconreal_cell{ii} = reconreal;
end

cw = prctile(abs(ReconstructedImages.reconRealtime(:)),99);
figure;
imshow(abs(ReconstructedImages.reconRealtime(:,:,1,340/2))./cw);
title('LV Contour');
set(gcf, 'units', 'normalized', 'Position', [0.045, 0.3, 0.4, 0.5])
%truesize([3000 3000])
[temp_LV,tmp_xi,tmp_yi]=roipoly;

%reconreal = squeeze(ReconstructedImages.reconRealtime.*cw);
for ii = 1:size(reconreal_cell,1)
    reconreal_cell{ii} = reshape(reconreal_cell{ii},size(reconreal_cell{ii},1)*size(reconreal_cell{ii},2),[]);
    reconreal_cell{ii} = reconreal_cell{ii}(temp_LV(:),:);
end

realtime_sig=cellfun(@mean, arrayfun(@(i) cell2mat(reconreal_cell(i:min(i+32, end))), 1:22:length(reconreal_cell), 'UniformOutput', false), 'UniformOutput', false);
realtime_sig=cell2mat(realtime_sig);
% realtime_sig = mean(cat(2,reconreal_cell{1:end-1}),1);
% realtime_sig = reshape(realtime_sig,(Nseg*moduleLength*numFA),[]);
% realtime_sig = mean(realtime_sig',1);
% realtime_sig = realtime_sig./realtime_sig(end);

end