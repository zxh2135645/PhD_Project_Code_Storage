Phi=reshape(Phi,[L sizes(2:end)]);
switch ScanType
  case 'Cine'
    temp=Gr\reshape(Phi(:,1,:,1,ceil(end/2)),L,[]);
  case 'IR'
    temp=Gr\reshape(Phi(:,:,1,1,ceil(end/2)),L,[]);
  case 'T2prep'
    temp=Gr\reshape(Phi(:,6,1,:),L,[]);
  case 'T2IR'
    temp=Gr\reshape(Phi(:,:,1,1,:),L,[]);
  case 'SR'
    temp=Gr\reshape(Phi(:,ceil(end/2),1,1,:),L,[]);
  case 'T2star'
    temp=Gr\reshape(Phi(:,end,1,1,:),L,[]); %first to last segment; all time points
end

% % ZH
% % add filter to eliminate ringing artifacts
% using_win=input('Using Filter Window? Yes [1] or No [0]: ')
% if using_win==1
%     tempwin_U=zeros(Ny,Nx,Nz,L);
%     for sli=1:Nz
%         temp_U=reshape(U,Ny,Nx,Nz,L,[]);
%         temp_U=fftshift(temp_U(:,:,sli,:),1);
%         tempwin=fft2(temp_U);
%         % Gaussian window
%         tempwin=ifft2(bsxfun(@times,tempwin,ifftshift(gausswin(Ny,1.2)*gausswin(Nx,1.2).')));
%         % Hanning window
% %         tempwin=ifft2(bsxfun(@times,tempwin,ifftshift(hann(Ny,'periodic')*hann(Nx,'periodic').')));
%         tempwin_U(:,:,sli,:)=tempwin;
%     end
%     figure,imshow(abs([temp_U(:,:,1,1) tempwin(:,:,1,1)]),[]);
%     U_temp=reshape(tempwin_U,[],1);
%     U=U_temp;
%     
%     temp=fftshift(ifftshift(reshape(reshape(dispim(reshape(U,Ny,Nx,Nz,[])),[],L)*temp,Ny,Nx,[]),1),1);
% else
%     temp=ifftshift(reshape(reshape(dispim(reshape(U,Ny,Nx,Nz,[])),[],L)*temp,Ny,Nx,[]),1);
% end
% % ZH
% 


% for j = 1:size(temp,2)
%     rbin = mode(Ridx(wall_clock==j));
%     temp(:,j) = Gr\Phi(:,ceil(end/2),ceil(end/2),rbin,j);
% end
temp=Gr\reshape(Phi(:,:,1,end),L,[]);
% temp = reshape(reshape(dispim(reshape(U,Nx,Ny,Nz,[])),[],L)*temp,Nx, Ny,[]);
dispim = @(x,st)fftshift(x(:,:,4,:),1);
temp=(reshape(reshape(dispim(reshape(U,Ny,Nx,Nz,[])),[],L)*temp,Ny,Nx,[]));


% if ~exist('cw','var')
    cw=0.5*max(vec(abs(temp)));
% end
%%
h = implay(abs(temp)/cw);
set(h.Parent,'Name','old_basal_echo1');
% clear temp;
%%
% figure('NumberTitle','off','Name','SliceThickness = 1.3mm '),
% for i = 1:8
%     subplot(1,8,i), imshow(abs(temp(:,:,i))/cw);
% end
% % 
%%
% % figure,
% % for i = 1:8
% %     subplot(2,4,i),imshow(images(:,:,i)/max(vec(images)));
% % end
% 
%%
% clear all 
% load('20P48_4wk_3mm_5meas_binningResults_0214.mat', 'Gr','L','Nx', 'Ny', 'Nz','Phi','U','dispim', 'vec','sizes','ScanType');
%     
%     BinningModified
%     
    
    
    
    
    
    
    


