fprintf('Estimating real-time subspace ... ');
navdata_temp=reshape(reshape(nav_data,[],Ncoils)*mixer(:,1:(Ncoils/3)),size(nav_data,1),[]);
[Phi_rt,S_rt,~]=svde(navdata_temp);
L = 32; %min(find(diag(S_rt)/S_rt(1)>0.01,1,'last'),32);
%figure,plot(diff(log(diag(S_rt)),2))
%figure,plot(20*log10(diag(S_rt)/S_rt(1)))
Phi_rt=Phi_rt(:,1:L).';
fprintf('done\n');