dispim = @(x)fftshift(x(:,:,1,:),1);

temp2 = cell(12,1);
for j = 7
    temp2{j,1} = abs(reshape(reshape(dispim(reshape(SpatialCoeff.U,Params.Ny,Params.Nx,Params.Nz,[])),[],32)*TemporalBasis.Phi_rt_small(:,DataArray.Ridx==1&DataArray.Hidx==j),Params.Ny,Params.Nx,[],1));
    temp2{j,1} = temp2{j,1}./prctile(temp2{j,1}(:),99.9);
end
dataArray.binsCardone =temp2;
implayZoom(dataArray.binsCardone{j});