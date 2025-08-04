load AhAx_test.mat;
tic;
AhAx = AhA_ps_cufinufft(filtAhb, st, SEs, Phi2);
AhAx = reshape(AhAx, st.Nd(1), st.Nd(2), st.Nz, []);
toc;

%%
U0 = AhA_ps_cufinufft(filtAhb, st, SEs, Phi2);
c2 = mean(abs(Ahb(:)),'all')./mean(abs(U0),'all');%real(pinv(U0(:))*Ahb(:));
U0 = filtAhb*c2;

tic;
U = pcg(@(x) AhA_ps_cufinufft(x, st, SEs, Phi2), Ahb(:), [], 10, [], [], U0(:));
U = reshape(U, st.Nd(1), st.Nd(2), st.Nz, []);
toc;

%%
% load('FSU_cu_test_rank2.mat');
% tic;
% SU_rank2= cufinufftf2d1(single(st.om(:,1)),single(st.om(:,2)), filtFSU_rank2, 1, 1e-6, 320,320, 280*320, 56*12 * 2); 
% SU_rank2 = ifft(reshape(SU_rank2, 320,320,56,12,2),[],3);
% toc;
% 
% U_rank2 = squeeze(sum(SU_rank2 .*conj(SEs),4));
% 
% L = size(U_rank2,4);
% 
% Phi2 = zeros(L,L,280,56,'like',single(1));
% st.om = single(st.om);
% 
% U_SE_rank2 = AhA_ps_cufinufft(U_rank2, st, SEs, Phi2);
