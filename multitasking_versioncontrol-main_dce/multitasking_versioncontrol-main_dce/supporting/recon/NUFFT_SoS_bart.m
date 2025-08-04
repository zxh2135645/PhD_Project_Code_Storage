function out = NUFFT_SoS_bart(x,st)
% x    - image to be transform back to k-space with stack-of spirals or 
%        stack-of-stars trajectory, size(x) = [Ny Nx Nz Nims]
%
% st.traj_bart - k-space trajectory, size(traj) = [3 Ntrajs Nkx]
% st.Nd        - target image, Nd = [Ny Nx]
% st.Nz        - taret image Nz
%
% out  - k-space data in stack-of spirals or stack-of-stars trajectory,
%        size(out) = [Ntrajs Nkx Nz Nims]
%

if ~isfield(st, 'bartnufft_fwd_scalar')
    st.bartnufft_fwd_scalar = 1;
end
    
x = x(:,:,:,:);

try
    [Ny,Nx,Nz,Nims]  = size(x);
    [~,Ntrajs,Nkx] = size(st.traj_bart);
    
    %fprintf('nufft...\n');
    x   = fft(x,[],3);
    
    bart_str = 'nufft';
    try
        x = reshape(x,Ny,Nx,1,[]);
        out = bart(bart_str,st.traj_bart,x);
        out = reshape(out,Ntrajs,Nkx,Nz,[]);
    catch
        x   = permute(x,[1 2 4 3]); % swap Nz and Nims dimension
        out = complex(zeros(1,Ntrajs,Nkx,Nims,Nz),0);
        for n = 1:Nims
            %fprintf('volume %d/%d...',n,Nims);
            %out(1,:,:,n,:) = bart(bart_str,st.traj_bart,x(:,:,n,:));
            evalc('out(1,:,:,n,:) = bart(bart_str,st.traj_bart,x(:,:,n,:));');
        end
        out = permute(out,[2 3 5 4 1]); % [Nkx Ntrajs Nz Nims]
    end
    out = out * st.bartnufft_fwd_scalar;
catch errormsg
    fprintf(2,'\n NUFFT_SoS_bart failed. \n');
    fprintf(2,'%s\n ', errormsg.message);
end

