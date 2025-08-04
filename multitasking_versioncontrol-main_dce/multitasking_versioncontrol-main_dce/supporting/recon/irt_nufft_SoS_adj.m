function im = irt_nufft_SoS_adj(x, st)

if ~isfield(st, 'irtnufft_adj_scalar')
    st.irtnufft_adj_scalar = 1;
end

Ny = st.Nd(1);
Nx = st.Nd(2); 
Nz = st.Nz;
NM = st.M;

x = reshape(x, NM, []);

% im = nufft_adj(x,st)/sqrt(NM);
im = nufft_adj(x,st) * st.irtnufft_adj_scalar;
im = reshape(im, Ny, Nx, Nz, []);
if Nz > 1
    im = ifft(im,[],3) * sqrt(Nz);
end

