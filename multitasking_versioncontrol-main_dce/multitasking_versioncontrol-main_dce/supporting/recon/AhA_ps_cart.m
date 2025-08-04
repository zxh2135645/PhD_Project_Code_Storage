function [AhAx, lowmem] = AhA_ps_cart(x, SEs, Phi2, lowmem)

if nargin < 4
    lowmem = 0;
end

L  = size(Phi2,1);
Ny = size(SEs,1);
Nz = size(SEs,3);
Ncoils = size(SEs,4);

% try
FU = fft(fft(bsxfun(@times,reshape(x,Ny,[],Nz,1,L),SEs),[],1),[],3)/sqrt(Ny)/sqrt(Nz);

try
    if lowmem == 0
        tempPhi2 = permute(Phi2,[3 5 4 6 1 2]);
        FU   = reshape(sum(FU.*tempPhi2,5),Ny,[],Nz,Ncoils,L);
        AhAx = reshape(sum(bsxfun(@times, ifft(ifft(FU,[],1),[],3)*sqrt(Ny)*sqrt(Nz),conj(SEs)),4),[],1);  
    end
catch errormsg
    fprintf('%s\n ', errormsg.message);
    fprintf('\n Try looping through trajectories... \n');
    lowmem = 1;
end

try
    if lowmem == 1
        for npy = 1:Ny
            for npz = 1:Nz
                FU(npy,:,npz,:,:) = reshape(reshape(FU(npy,:,npz,:,:),[],L) * Phi2(:,:,npy,npz),[],Ncoils,L);
            end
        end
        AhAx = reshape(sum(bsxfun(@times, ifft(ifft(FU,[],1),[],3)*sqrt(Ny)*sqrt(Nz),conj(SEs)),4),[],1);    
    end
catch errormsg
    fprintf('%s\n ', errormsg.message);
    fprintf('\n Try looping through coils and trajectories... \n');
    lowmem = 2;
end

try
    if lowmem == 2
        for ncoil = 1:Ncoils
            for npy = 1:Ny
                for npz = 1:Nz
                    FU(npy,:,npz,ncoil,:) = reshape(reshape(FU(npy,:,npz,ncoil,:),[],L) * Phi2(:,:,npy,npz),[],L);
                end
            end
        end
        AhAx = reshape(sum(bsxfun(@times, ifft(ifft(FU,[],1),[],3)*sqrt(Ny)*sqrt(Nz),conj(SEs)),4),[],1); 
    end
catch errormsg
    fprintf(2,'\n AhA_ps_cart failed. \n');
    fprintf(2,'%s\n ', errormsg.message);
end


% function AhAx = AhA_ps_cart(x,SEs,Phi2)
% 
% L  = size(Phi2,1);
% Ny = size(Phi2,3);
% Nz = size(Phi2,4);
% Ncoils = size(SEs,4);
% 
% % try
% FU = fft(fft(bsxfun(@times,reshape(x,Ny,[],Nz,1,L),SEs),[],1),[],3)/sqrt(Ny)/sqrt(Nz);
% 
% for npy = 1:Ny
%     for npz = 1:Nz
%         FU(npy,:,npz,:,:) = reshape(reshape(FU(npy,:,npz,:,:),[],L) * Phi2(:,:,npy,npz),[],Ncoils,L);
%     end
% end
% 
% AhAx = reshape(sum(bsxfun(@times, ifft(ifft(FU,[],1),[],3)*sqrt(Ny)*sqrt(Nz),conj(SEs)),4),[],1);
% 
% end