function out = NUFFT_SoS_adj_bart(x,st,useDCF)
% x      - k-space data collected using stack-of spirals or stack-of-stars
%          trajectory, size(x) = [Ntrajs Nkx Nz Nims]
%
% st.traj_bart - k-space trajectory, size(traj) = [3 Ntrajs Nkx]
% st.Nd        - target image, Nd = [Ny Nx]
% st.Nz        - taret image Nz
% 
% useDCF - 1: inverse nufft, estimate density compensation function
%          0: adjoint nufft, no density compensation
%
% out    - 3D image series, size(out) = [Ny Nx Nz Nims]
%

if nargin < 3
    useDCF = 0;
end

if ~isfield(st, 'bartnufft_adj_scalar')
    st.bartnufft_adj_scalar = 1;
end

Nd = st.Nd;
[Ntrajs, Nkx, Nz, Nims] = size(x);

try
    if st.flagUseGPU
        if useDCF
            %fprintf('inverse nufft...\n');
            bart_str = ['nufft -i -g -d ' num2str(Nd(1)) ':' num2str(Nd(2)) ':1 -t'];
        else
            %fprintf('adjoint nufft...\n');
            bart_str = ['nufft -a -g -d ' num2str(Nd(1)) ':' num2str(Nd(2)) ':1 -t'];
        end
    else
        if useDCF
            %fprintf('inverse nufft...\n');
            bart_str = ['nufft -i -d ' num2str(Nd(1)) ':' num2str(Nd(2)) ':1 -t'];
        else
            %fprintf('adjoint nufft...\n');
            bart_str = ['nufft -a -d ' num2str(Nd(1)) ':' num2str(Nd(2)) ':1 -t'];
        end
    end
    
    try
        x = reshape(x,1,Ntrajs,Nkx,1,[]);
%         evalc('out = bart(bart_str,st.traj_bart,x);');
        out = bart(bart_str,st.traj_bart,x);
        out = reshape(out,Nd(1),Nd(2),Nz,[]);
    catch
        x   = permute(x,[5 1 2 4 3]);   % swap Nz and Nims dimension
        out = complex(zeros(Nd(1),Nd(2),Nims,Nz),0);
        for n = 1:Nims
            out(:,:,n,:) = bart(bart_str,st.traj_bart,x(1,:,:,n,:));
%             evalc('out(:,:,n,:) = bart(bart_str,st.traj_bart,x(1,:,:,n,:));');
        end
        out = permute(out,[1 2 4 3]);   % swap Nz and Nims dimension back
    end
    out = out*st.bartnufft_adj_scalar;
    out = ifft(out,[],3);

catch errormsg
    fprintf(2,'\n NUFFT_SoS_adj_bart failed. \n');
    fprintf(2,'%s\n ', errormsg.message);
end

