
vec=@(x)x(:);
prep=@(x,st)reshape(x,st.Nd(1),st.Nd(2),[]);
prep_adj=@(x,st)reshape(x,st.M,[]);
dcff=@(k)(k+(k==0)*min(k(k~=0))/8);
dcf=@(st,c)c*dcff(sqrt(sum(abs(st.om).^2,2)));
crop=@(x,st)x(1:st.Nd(1),1:st.Nd(2),:,:);
dispim=@(x)x(ovs/2 + (1:Norig),ovs/2 + (1:Norig),:);
if useGPU
  Ainv=@(x,st,c)st.F'*(diag_sp(dcf(st,c))*prep_adj(x,st));
  A=@(x,st)st.F*reshape(x,prod(st.Nd),[]);
  Ah=@(x,st)st.F'*prep_adj(x,st);
  AhA=@(x,st,GhG)vec(st.F'*(st.F*reshape(x,prod(st.Nd),[])));
  A_sense=@(x,st,SEs,Nt)st.FS*reshape(x,prod(st.Nd),[]);
  Ah_sense=@(x,st,SEs,Nt)vec(st.FS'*prep_adj(x,st));
  AhA_sense=@(x,st,SEs,GhG,Nt)vec(st.FS'*(st.FS*reshape(x,prod(st.Nd),[])));
else
  Ainv=@(x,st,c)nufft_adj(diag_sp(dcf(st,c))*prep_adj(x,st),st);
  A=@(x,st)nufft(prep(x,st),st);
  Ah=@(x,st)vec(nufft_adj(prep_adj(x,st),st));
  AhA=@(x,st,GhG)vec(conj(repmat(st.sn,[1 1 numel(x)/prod(st.Nd)])).*...
    crop(ifft2(reshape(...
    GhG*reshape(fft2(padarray(repmat(st.sn,[1 1 numel(x)/prod(st.Nd)]).*prep(x,st),...
    st.Kd-st.Nd,0,'post')),prod(st.Kd),[]),st.Kd(1),st.Kd(2),[])),st)*prod(st.Kd));
  A_sense=@(x,st,SEs,Nt)nufft(repmat(prep(x,st),[1 1 1 size(SEs,3)]).*repmat(reshape(SEs,size(SEs,1),size(SEs,2),1,size(SEs,3)),[1 1 Nt 1]),st);
  Ah_sense=@(x,st,SEs,Nt)vec(sum(repmat(reshape(conj(SEs),[],1,size(SEs,3)),[1 Nt 1]).*reshape(nufft_adj(prep_adj(x,st),st),[],Nt,size(SEs,3)),3));
  AhA_sense=@(x,st,SEs,GhG,Nt)vec(repmat(st.sn,[1 1 Nt]).*sum(repmat(reshape(conj(SEs),size(SEs,1),size(SEs,2),1,size(SEs,3)),[1 1 Nt]).*crop(ifft2(reshape(...
    GhG*reshape(fft2(padarray(repmat(reshape(SEs,size(SEs,1),size(SEs,2),1,size(SEs,3)),[1 1 Nt 1]).*repmat(repmat(st.sn,[1 1 Nt]).*prep(x,st),[1 1 1 Ncoils]),...
    st.Kd-st.Nd,0,'post')),prod(st.Kd),[]),st.Kd(1),st.Kd(2),Nt,[])),st),4)*prod(st.Kd));
end

% set up preconditioner
[kx, ky]=ndgrid(-(N/2):(N/2-1),-(N/2):(N/2-1));
k = sqrt(kx.^2+ky.^2);
window = 1./fftshift(abs(fft2(Ah(ones(st.M,1),st))));%Randy modification
%original
%window=ones(N,N);%Randy modification
window=reshape(window,[N N]);
window=window*(floor(3*N/4-1)-N/2-1)/window(N/2+1,floor(3*N/4-1));
window(k~=0) = k(k~=0);
window=ifftshift(window);
window(window>(N/2-1))=min(window(:));
window=window/N;

% window=1./abs(fft2(prep(Ah(ones(st.M,1),st),st)));
% window(ifftshift(k)>(N/2-1))=min(window(:));

Mf=@(x)ifft2(fft2(x).*repmat(window,[1 1 size(x,3)]));
M=@(x)vec(Mf(prep(x,st)));
