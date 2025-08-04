L=32;sizes=size(TemporalBasis.Phi);
Phi=reshape(TemporalBasis.Phi,[L sizes(2:end)]);

%%

temp=TemporalBasis.Gr\reshape(Phi(:,:,7,1,1),L,[]);
% temp = reshape(reshape(dispim(reshape(U,Nx,Ny,Nz,[])),[],L)*temp,Nx, Ny,[]);
dispim = @(x,st)fftshift(x(:,:,1,:),1);
temp=(reshape(reshape(dispim(reshape(SpatialCoeff.U,Params.Ny,Params.Nx,Params.Nzdisp,[])),[],L)*temp,Params.Ny,Params.Nx,[]));


% if ~exist('cw','var')
    cw=max((abs(temp(:))));
% end
%
h = implay([abs(temp)/cw angle(temp)/2/pi+0.5]);
set(h.Parent,'Name','old_basal_echo1');
 save_path='/home/biri/';

% Han_CreateVideo(abs(temp(:,:,:))/cw,10,fullfile(save_path,'result_reconWT.avi'),'avi')

% clear temp;
