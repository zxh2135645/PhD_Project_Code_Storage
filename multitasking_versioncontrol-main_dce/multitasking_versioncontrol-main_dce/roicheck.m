dispim = @(x)fftshift(x(:,:,1,:),1);
vec = @(x) x(:);
Phi = TemporalBasis.Phi;
Ny = Params.Ny;
Nx = Params.Nx;
Nz = Params.Nz;
L = ReconOptions.L;
U = SpatialCoeff.U;
Gr = TemporalBasis.Gr;

for i=1:size(Phi,5)
    temp = Gr\reshape(Phi(:,:,10,1,i), L, []);
    temp = reshape(reshape(dispim(reshape(U,Ny,Nx,Nz,[])),[],L)*temp, Ny, Nx, [], Params.Necho);

    figure;imshow(temp(:,:,100),[])
    roi=roipoly;
    for j=1:size(temp,3)
        dat=temp(:,:,j);
        mean_val(i,j)=mean(dat(roi));
        std_val(i,j)=std(dat(roi));
    end
end

figure;
for i=1:size(Phi,5)
    plot(1:680,mean_val(i,:))
    hold on;
    lege{i}=['temept',num2str(i)];

end
legend(lege)