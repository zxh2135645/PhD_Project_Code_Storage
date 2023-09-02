Phi=reshape(Phi,[L sizes(2:end)]);
switch ScanType
  case 'Cine'
    temp=Gr\reshape(Phi(:,1,:,1,ceil(end/2)),L,[]);
  case {'IR','T2prep'}
    %temp=Gr\reshape(Phi(:,:,1,1,ceil(end/2)),L,[]);
    %temp=Gr\reshape(Phi(:,3,6,5,:),L,[]);
    temp1=Gr\reshape(Phi(:,3,5,3,:),L,[]);
    temp2=Gr\reshape(Phi(:,end-3,5,3,:),L,[]);
  case 'T2IR'
    temp=Gr\reshape(Phi(:,:,2,1,:),L,[]);
  case 'SR'
    temp=Gr\reshape(Phi(:,ceil(end/2),1,1,:),L,[]);
end
temp1=reshape(reshape(U,Ny*Nx,[])*temp1,Ny,Nx,[]);
temp2=reshape(reshape(U,Ny*Nx,[])*temp2,Ny,Nx,[]);
%implay(abs(dispim(temp))/max(abs(vec(dispim(temp)))));
T2=(abs(30/log(abs(temp2)./abs(temp1))));
%clear temp;