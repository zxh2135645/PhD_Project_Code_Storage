function im = NUFFT_SoS_adj(x,st)

try
  am=getGPUmem(); %get available memory (am)
  if min(am) < 4e3
    gpuWorkerReset(); %reset if < 4e3 MB
  end
  sz=size(x);
  x=x(:,:,:);
  %im = zeros(st.Nd);
  % size(x)
  Nims=size(x,3);
  
  % delete(gcp('nocreate')); %delete any open pool
  % parpool('local',6); %start a new poool
  parfor np=1:st.Nd(3)
    F=gpuNUFFT(st.om.'/(2*pi),ones(1,st.M),2,6,8,st.Nd(1:2),[],true);
    for nc=1:Nims %Ncoils
      im(:,:,np,nc)=F'*x(:,np,nc);
    end
  end
  
  im=ifft(im,[],3)*sqrt(st.Nd(3));
  im=reshape(im,[st.Nd,sz(3:end)]);
  
catch errormsg
  errormsg
  delete(gcp('nocreate'));
  parpool('local',8)
  im=NUFFT_SoS_adj(x,st);
end

end
