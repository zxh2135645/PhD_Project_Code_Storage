Phi=reshape(Phi,[L sizes(2:end)]);
switch ScanType
  case 'Cine'
    temp=Gr\reshape(Phi(:,1,:,1,ceil(end/2)),L,[]);
  case {'IR','T2prep'}
    %temp=Gr\reshape(Phi(:,:,1,1,ceil(end/2)),L,[]);
    temp=Gr\reshape(Phi(:,3,1,4,:),L,[]);
    
  case 'T2IR'
    temp=Gr\reshape(Phi(:,:,2,1,:),L,[]);
  case 'SR'
    temp=Gr\reshape(Phi(:,ceil(end/2),1,1,:),L,[]);
    
end
temp=reshape(reshape(U,Ny*Nx,[])*temp,Ny,Nx,[]);
implay(abs(dispim(temp))/max(abs(vec(dispim(temp)))));
%clear temp;