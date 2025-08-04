function im = finufft_SoS(x, st)

if ~isfield(st, 'finufft_fwd_scalar')
    st.finufft_fwd_scalar = 1;
end
    
Ny = st.Nd(1);
Nx = st.Nd(2); 
Nz = st.Nz;
NM = st.M;

if Nz > 1
    x = fft(x,[],3) / sqrt(Nz);
end
x = reshape(x, Ny, Nx, []);
Nim = size(x,3);
batch = 50;
Nseg = floor(Nim/batch);

opts.upsampfac = st.osf;

im = zeros(NM,Nim);
for n = 1:Nseg
    im(:,(1:batch)+batch*(n-1)) = finufft2d2(single(st.om(:,1)),single(st.om(:,2)),-1,1e-6,single(x(:,:,(1:batch)+batch*(n-1))),opts);
end
if Nim > batch*Nseg
    im(:,batch*Nseg+1:end) = finufft2d2(single(st.om(:,1)),single(st.om(:,2)),-1,1e-6,single(x(:,:,batch*Nseg+1:end)),opts);
end
%im = finufft2d2(single(st.om(:,1)),single(st.om(:,2)),-1,1e-6,single(x),opts);
im = reshape(im, NM, Nz, []) * st.finufft_fwd_scalar;

