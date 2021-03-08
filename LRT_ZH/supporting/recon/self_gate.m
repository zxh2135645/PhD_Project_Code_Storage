function [Navdata_tensor,mask,Segidx,Hidx,Ridx,wall_clock]=self_gate(Navdata,Hfreq_range,Hframes,Rfreq_range,Rframes,frame_rate,segments)

if nargin < 7
  segments=1;
end

isMATLAB=logical(numel(ver('matlab')));
between=@(x,lb,ub)(x>lb)&(x<ub);
iseven=@(x)mod(x,2)==0;
vec=@(x)x(:);

[U,S,~]=svd(Navdata,'econ');
N=size(U,1);
L=min(3,size(U,2));

USf=fftshift(fft(U(:,1:L)),1)*S(1:L,1:L);
fs=linspace(-frame_rate/2,frame_rate/2,floor(N/2)*2+1);
if iseven(N)
  fs=fs(1:(end-1));
end

%% Cardiac gate
Hfreq_range=between(fs,Hfreq_range(1),Hfreq_range(2)) ...
  | between(-fs,Hfreq_range(1),Hfreq_range(2));
[~,Hargmax]=max(vec(abs(USf.*repmat(Hfreq_range(:),[1 L]))));
[Hfreq, Hcomp]=ind2sub(size(U),Hargmax);
Hfreq=abs(fs(Hfreq))

%NOT ROBUST
features=[abs(USf(Hfreq_range,Hcomp)), abs(abs(fs(Hfreq_range).')-Hfreq)];
for j=1:size(features,2)
  features(:,j)=features(:,j)/std(features(:,j));
end
features(:,2)=2*features(:,2);
if isMATLAB
  try
    [idx,c]=kmeans(features,3);
  catch
    keyboard
  end
else
  [c,~,idx]=newKmeans(features,2);
end
figure, hold all
for j=1:size(c,1)
  plot(features(idx==j,1),features(idx==j,2),'.');
end
hold off
drawnow
[~,Hgroup]=min(sum([-c(:,1) c(:,2)],2));

Hsig=zeros(N,1);
Hsig(Hfreq_range)=USf(Hfreq_range,Hcomp).*(idx==Hgroup);
Hsig=ifft(ifftshift(Hsig,1));
[Hsig,~,~]=svd([real(Hsig) imag(Hsig)],'econ');
Hsig=Hsig(:,1);
Hsig=angle(hilbert(Hsig));
[~,Hidx]=histc(Hsig(:),linspace(-pi,pi,Hframes+1));

%% Respiratory gate
if prod(Rfreq_range==0)
  Ridx=mod(0:(numel(Hidx)-1),Rframes).'+1; %allows Dchen data override
else
  Rfreq_range=between(fs,Rfreq_range(1),Rfreq_range(2)) ...
    | between(-fs,Rfreq_range(1),Rfreq_range(2));
  [~,Rargmax]=max(vec(abs(USf.*repmat(Rfreq_range(:),[1 L]))));
  [Rfreq, Rcomp]=ind2sub(size(U),Rargmax);
  Rfreq=abs(fs(Rfreq))
  
  %NOT ROBUST
  features=[abs(USf(Rfreq_range,Rcomp)), abs(abs(fs(Rfreq_range).')-Rfreq)];
  for j=1:size(features,2)
    features(:,j)=features(:,j)/std(features(:,j));
  end
  if isMATLAB
    [idx,c]=kmeans(features,2);
  else
    [c,~,idx]=newKmeans(features,2);
  end
  figure, hold all
  for j=1:size(c,1)
    plot(features(idx==j,1),features(idx==j,2),'.');
  end
  hold off
  drawnow
  [~,Rgroup]=min(sum([-c(:,1) c(:,2)],2));
  
  Rsig=zeros(N,1);
  Rsig(Rfreq_range)=USf(Rfreq_range,Rcomp).*(idx==Rgroup);
  Rsig=ifft(ifftshift(Rsig,1));
  [Rsig,~,~]=svd([real(Rsig) imag(Rsig)],'econ');
  Rsig=Rsig(:,1);
  Rsig=angle(hilbert(Rsig));
  [~,Ridx]=histc(Rsig(:),linspace(-pi,pi,Rframes+1));
end
%% Form tensor
wall_clock=ceil((1:N)*Hfreq/frame_rate);
% wall_clock(wall_clock==1)=2; %group first two HBs so no extrapolation
% wall_clock = wall_clock-1; %adjust
wall_clock(wall_clock==max(wall_clock))=max(wall_clock)-1; %group last two HBs so no extrapolation

Segidx=mod(0:(N-1),segments)+1;

%below is included only for memory calculations
try
  Navdata_tensor=complex(zeros(size(Navdata,2),segments,Hframes,Rframes,max(wall_clock)));
catch
  wall_clock(:)=1;
  disp('Memory warning: Collapsing along heart beat dimension!')
  Navdata_tensor=complex(zeros(size(Navdata,2),segments,Hframes,Rframes,max(wall_clock)));
end

try
  mask=zeros(size(Navdata,2),segments,Hframes,Rframes,max(wall_clock));
catch
  wall_clock(:)=1;
  disp('Memory warning: Collapsing along heart beat dimension!')
  Navdata_tensor=complex(zeros(size(Navdata,2),segments,Hframes,Rframes,max(wall_clock)));
  mask=zeros(size(Navdata,2),segments,Hframes,Rframes,max(wall_clock));
end

return
