if ~isVE
  %calculate initial angle
  dk=kspace_data-nav_data;
  mk=kspace_data+nav_data;
  sim=4*sum(abs(dk(:,:)).^2,2)./sum(abs(mk(:,:)).^2,2);
  clear dk mk
  if strcmp(ScanType,'SR')
    sim(:)=1; %Doesn't seem to help for perfusion scans
  end
end

%%
fib=[1 1];
while fib(2) < params.lBaseResolution * pi / 2
  fib = [fib(2), sum(fib)];
end
trajs = fib(2);

% phi=(1+sqrt(5))/2; %exact
phi = fib(2)/fib(1);
theta=180/phi;

%% setup for different versions
if ~isVE
  %find 0 degree line
  simlast=floor(numel(sim)/trajs)*trajs;
  klast = floor(size(kspace_data,1)/trajs)*trajs;
  temp=squeeze(mean(reshape(kspace_data(1:klast,:,1,:),trajs,[],N,Ncoils),2));
  err=zeros(trajs,1);
  parfor j=1:trajs
    [~,sorter]=sort(mod(theta*mod((1:trajs)-j,trajs),360));
    err(j)=norm(reshape(diff(temp(sorter,:,:),1),[],1));
  end
  [~,theta0]=min(err.*mean(reshape(sim(1:simlast),trajs,[]),2));
  %         figure,plot(1./err./mean(reshape(sim(1:simlast),trajs,[]),2))
  thetas=mod(theta*mod((1:size(kspace_data,1))-theta0,trajs),360);
  
  %setup gradient delay correction
  oneeighties=abs(thetas-180)==min(abs(thetas-180));
  ksig = ifft(kspace_data(oneeighties,end:-1:1,:),[],2);
  nsig = ifft(nav_data(oneeighties,:,:),[],2);
  %       nsig = ifftshift(padarray(ifft(ifftshift(fft(nav_data(oneeighties,:,:),[],2),2),[],2),[0 ovs/2 0]),2)/sqrt(N);
else
  thetas=mod(theta*mod((1:size(kspace_data,1))-1,2*trajs),360);
  
  %setup gradient delay correction
  cand_coords=unique(thetas);
  cand_coords(cand_coords>=180)=[];  %keep only thetas<180
  ksig=zeros(size(cand_coords,1),size(kspace_data,2),size(kspace_data,4));
  nsig=ksig;
  for j=1:size(cand_coords,1)
    t_ind=(thetas==cand_coords(j)+180);
    ksig(j,:,:)=mean(ifft(kspace_data(t_ind,end:-1:1,:),[],2),1); %rather than mean, rank-1 correction would be better
    t_ind=(thetas==cand_coords(j));
    nsig(j,:,:)=mean(ifft(kspace_data(t_ind,:,:),[],2),1);
  end
end

%% Gradient delay correction
x=[0:(Norig-1) (-Norig):-1];
cost=@(shift)norm(reshape(ksig-nsig.*repmat(exp(1i*(2*shift-pi/Norig)*x),[size(ksig,1) 1 Ncoils]),[],1));
shifts=linspace(-.1,.1,201);
mincost=inf;
for j=1:numel(shifts)
  tempcost=cost(shifts(j));
  if tempcost < mincost
    mincost=tempcost;
    shift=shifts(j);
  end
end
shift=fminsearch(cost,shift);
kspace_data=fft(ifft(kspace_data,[],2).*repmat(exp(1i*shift*x),...
  [size(kspace_data,1) 1 size(kspace_data,3) Ncoils]),[],2);
nav_data=fft(ifft(nav_data,[],2).*repmat(exp(1i*shift*x),...
  [size(nav_data,1) 1 size(nav_data,3) Ncoils]),[],2);

if isVE %force all angles into [0, 180)-degree range
  kspace_data(thetas>=180,:,:,:)=kspace_data(thetas>=180,[1 end:-1:2],:,:);
  thetas=mod(theta*mod((1:size(kspace_data,1))-1,trajs),180);
end

r=linspace(-pi,pi,2*Norig+1); r(end)=[];
om1=cos(thetas(1:trajs)*pi/180).'*r;
om2=sin(thetas(1:trajs)*pi/180).'*r;
om=[om1(:), om2(:)];

fprintf('Setting up NUFFT...')
if useGPU
  st.Nd=[N N];
  st.M=size(om,1);
  st.om=om;
  st.F=gpuNUFFT(om.'/(2*pi),ones(1,st.M),2,6,8,st.Nd,[],true);
  GhG=[];
else
  st=nufft_init(om,[N N],[6 6],2*[N N],[N N]/2);
  GhG = st.p.G'*st.p.G;
end
clear om1 om2 om
fprintf('done.\n')
