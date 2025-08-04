function dataArray = genRespBinsImage(params,reconOptions,dataArray,temporalBasis,spatialCoeff)
%v0.2

extractVarFromStruct(params);
extractVarFromStruct(reconOptions);
extractVarFromStruct(dataArray);
extractVarFromStruct(temporalBasis);
extractVarFromStruct(spatialCoeff);


%% Images for individual bins

dispim = @(x) x(floor(Ny/2)-floor(Nydisp/2) + (1:Nydisp), floor(Nx/2)-floor(Nxdisp/2) + (1:Nxdisp), :, :);

temp = zeros(Nydisp,Nxdisp,Nz,rbins);
for j = 1:rbins
    temp(:,:,:,j) = abs(reshape(reshape(dispim(reshape(U_init,Ny,Nx,Nz,[])),[],L_init)...
                    *mean(Phi_rt_small_init(:,Ridx==j),2),Nydisp,Nxdisp,Nz,[]));
end
cw = prctile(temp(:),99);
dataArray.binsResp = abs(temp)/cw;

if flagCommandLine
    figure(11); subplot(2,1,2),plot(Ridx);axis([-inf inf 0 rbins+1]);
    implay(dataArray.binsResp,2);
end


