function k = irt_nufft_SoS(x, st)

if ~isfield(st, 'irtnufft_fwd_scalar')
    st.irtnufft_fwd_scalar = 1;
end

Ny = st.Nd(1);
Nx = st.Nd(2); 
Nz = st.Nz;
NM = st.M;

if Nz > 1
    x = fft(x,[],3) / sqrt(Nz);
end
x = reshape(x, Ny, Nx, []);
% k = nufft(x,st)/sqrt(NM);
k = nufft(x,st) * st.irtnufft_fwd_scalar;
k = reshape(k, NM, Nz, []);

