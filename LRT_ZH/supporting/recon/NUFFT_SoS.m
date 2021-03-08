function out = NUFFT_SoS(x,st)
%10/16/2016 added try/catch for worker pool issue
try
am=getGPUmem(); %get available memory (am)
if min(am) < 4e3
  gpuWorkerReset(); %reset if < 4e3 MB
end
sz=size(x);
x=x(:,:,:,:);
ks=fft(x,[],3)/sqrt(st.Nd(3));

out = zeros(st.M,st.Nd(3),size(ks,4));
Nims=size(ks,4);


parfor np=1:st.Nd(3)
  F=gpuNUFFT(st.om.'/(2*pi),ones(1,st.M),2,6,8,st.Nd(1:2),[],true);
  for nc=1:Nims
    out(:,np,nc)=F*ks(:,:,np,nc);
  end
  
end

out=reshape(out,[st.M,st.Nd(3),sz(4:end)]);
catch errormsg
    errormsg
    delete(gcp('nocreate'));
    parpool('local',8)
    out=NUFFT_SoS(x,st);
end

end
