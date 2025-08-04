function im = cufinufft_SoS_adj(x, st)

if ~isfield(st, 'cufinufft_adj_scalar')
    st.cufinufft_adj_scalar = 1;
end

Ny = st.Nd(1);
Nx = st.Nd(2); 
Nz = st.Nz;
NM = st.M;

x = reshape(x, NM, []);
Nim = size(x,2);
batch = 50;
Nseg = floor(Nim/batch);

opts.upsampfac = st.osf;

im = zeros(Ny,Nx,Nim);
if Nim == 1
    im = cufinufftf2d1(single(st.om(:,1)),single(st.om(:,2)),single([x x]),1,1e-6,Ny,Nx,NM,2);
    im = im(:,:,1);
else
    for n = 1:Nseg
        im(:,:,(1:batch)+batch*(n-1)) = cufinufftf2d1(single(st.om(:,1)),single(st.om(:,2)),single(x(:,(1:batch)+batch*(n-1))),+1,1e-6,Ny,Nx,NM,size(x(:,(1:batch)+batch*(n-1)),2));
    end
    if Nim > batch*Nseg
        im(:,:,batch*Nseg+1:end) = cufinufftf2d1(single(st.om(:,1)),single(st.om(:,2)),single(x(:,batch*Nseg+1:end)),+1,1e-6,Ny,Nx,NM,size(x(:,batch*Nseg+1:end),2));
    end
    %im = finufft2d1(single(st.om(:,1)),single(st.om(:,2)),single(x),+1,1e-6,Ny,Nx,opts);
end
im = reshape(im, Ny, Nx, Nz, []) * st.cufinufft_adj_scalar;
if Nz > 1
    im = ifft(im,[],3) * sqrt(Nz);
end

