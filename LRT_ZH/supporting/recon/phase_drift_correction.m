fprintf('Phase drift correction...');

% for i = 1:params.NEco
% 
%     ph_corrs=linspace(-pi,pi,101)/(Nread/params.NEco);
%     cost=@(x)std(exp(-1i*x*nav_indices(i:params.NEco:end)).*sign(Phi_rt_small(1,i:params.NEco:end)));
%     costs=zeros(numel(ph_corrs),1);
%     for j=1:numel(ph_corrs)
%       costs(j)=cost(ph_corrs(j));
%     end
%     [~,j]=min(costs);
%     ph_corr=fminsearch(cost,ph_corrs(j));
%     ph_corr=exp(-1i*ph_corr*(1:(Nread)));
% 
% %     ph_corr_k=ph_corr;
% %     ph_corr_k(nav_indices(i:params.NEco:end))=[];
% 
%     nav_data(i:params.NEco:end,:,:,:)=bsxfun(@times,nav_data(i:params.NEco:end,:,:,:),vec(ph_corr(nav_indices(i:params.NEco:end))));
%     if iscartesian
%         kspace_data(ss_index+i-1,:,:,:)=bsxfun(@times,kspace_data(ss_index+i-1,:,:,:),ph_corr(i:params.NEco:end)');
%     else
%         kspace_data(ss_index+i-1,:,:,:)=bsxfun(@times,kspace_data(ss_index+i-1,:,:,:),ph_corr_k(:));
%     end
% end


%%


ph_corrs=linspace(-pi,pi,101)/Nread;
cost=@(x)std(exp(-1i*x*nav_indices).*sign(Phi_rt_small(1,:)));
costs=zeros(numel(ph_corrs),1);
for j=1:numel(ph_corrs)
  costs(j)=cost(ph_corrs(j));
end
[~,j]=min(costs);
ph_corr=fminsearch(cost,ph_corrs(j));
ph_corr=exp(-1i*ph_corr*(1:Nread));

ph_corr_k=ph_corr;
ph_corr_k(nav_indices)=[];

nav_data=bsxfun(@times,nav_data,vec(ph_corr(nav_indices)));
if iscartesian
    kspace_data=bsxfun(@times,kspace_data,ph_corr(:));
else
    kspace_data=bsxfun(@times,kspace_data,ph_corr_k(:));
end
