fprintf('Estimating real-time subspace ... ');


[Phi_rt,S_rt,~]=svde(nav_data(:,:));%Phi_rt based on all echoes
%XG, Phi_rt based on first echo;
% [Phi_rt,S_rt,~]=svde(nav_data(1:params.NEco:end,:));
%XG end

L = min(find(diag(S_rt)/S_rt(1)>0.001,1,'last'),32)
% fprintf('What L do you really want? \n');
% keyboard;

figure,plot(diff(log(diag(S_rt)),2))
figure,plot(20*log10(diag(S_rt)/S_rt(1)))
Phi_rt_small=Phi_rt(:,1:L).';

NnavsPerBlock=params.lSegments/SGblock;


temp=permute(reshape(Phi_rt_small,L,params.NEco,NnavsPerBlock,[]),[3 1 2 4]);
% temp = Phi_rt_small';


%XG Phi_rt based on first echo
% temp=permute(reshape(Phi_rt_small,L,NnavsPerBlock,[]),[2 1 3]);
%XG end

clear Phi_rt;
%XG Phi_rt extrap instead of linear.
Phi_rt=interp1(nav_indices(1:NnavsPerBlock),temp,1:params.lSegments,'pchip','extrap');
% Phi_rt=interp1(nav_indices,temp,1:Nread,'linear');
% Phi_rt=interp1(nav_indices,temp,1:Nread,'previous','extrap');

temp_hold=interp1(nav_indices(1:NnavsPerBlock),temp,1:params.lSegments,'previous','extrap');
% temp_hold=interp1(nav_indices,temp,1:Nread,'previous','extrap');

Phi_rt(isnan(Phi_rt))=temp_hold(isnan(Phi_rt));
clear temp temp_hold;

Phi_rt=ipermute(Phi_rt,[3 1 2 4]); %restore dimension ordering
Phi_rt=reshape(permute(reshape(Phi_rt,L,params.NEco,SGblock,[]),[1 3 2 4]),L,[]); %swap echos and SGblocks

% Phi_rt=reshape(repmat(reshape(Phi_rt_small,L,1,params.NEco,[]),[1 SGblock 1 1]),L,[]);

% Phi_rt = Phi_rt';
Phi_rt_full=Phi_rt;

% Phi_rt=Phi_rt_full;

if ~iscartesian
  Phi_rt(:,nav_indices)=[];
end

fprintf('done\n');


%% iterative single echo estimate

% % Phi_rt = nav_data(:,:);
% 
% [Phi_rt,S_rt,~]=svde(nav_data(1:params.NEco:end,:));
% 
% L = 7;
%             
% Phi_rt_small=Phi_rt(:,1:L).';
% temp = Phi_rt_small';
% 
% Phi_rt=interp1(nav_indices(1:params.NEco:end),temp,1:Nread,'previous','extrap');
% % temp_hold=interp1(nav_indices,temp,1:Nread,'previous','extrap');
% % Phi_rt(isnan(Phi_rt))=temp_hold(isnan(Phi_rt));
% % clear temp temp_hold;
%         
% Phi_rt = Phi_rt';
% Phi_rt_full=Phi_rt;
% 
% Phi_rt=Phi_rt_full;
% if ~iscartesian
%   Phi_rt(:,nav_indices)=[];
% end
% 
% fprintf('done\n');        
% %         
% %         
% %         
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        


