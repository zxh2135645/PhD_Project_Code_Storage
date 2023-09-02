if ~exist('card_roi','var')
    h=figure;
    imshow(dispim(abs(composite_fbp)),[]);
    title('Select Cardiac ROI')
    card_roi=imrect;
    card_roi=padarray(createMask(card_roi),[ovs/2 ovs/2]);
    close(h)
    drawnow
end

missing=round((params.alTR_seconds+[.012 .02 .03 .04 .05]-params.lSegments*params.lEchoSpacing)/params.lEchoSpacing/2)
Phi_rt=reshape(Phi_rt_init,size(Phi_rt_init,1),Nseg/2,5,[]);
Phi_rt(:,end+(1:max(missing)),:,:)=repmat(Phi_rt(:,end,:,:),[1 max(missing) 1 1]);
for j=1:5
    Phi_rt(:,end-max(missing)+missing(j)+1:end,j,:)=0;
end
Phi_rt=Phi_rt(:,:);
Phi_rt=Phi_rt(:,Phi_rt(1,:)~=0);

Nseg2 = 5*params.lSegments+sum(missing)*2;
temp=reshape(U_init,[],L_init);
temp=Phi_rt.'*temp(card_roi,:).';
%   temp=navdata_temp;
frame_rate = 1/(params.lEchoSpacing*2);
fs=linspace(-frame_rate/2,frame_rate/2,floor(size(temp,1)/2)*2+1);
if mod(size(temp,1),2)==0 %iseven
    fs=fs(1:(end-1));
end
if false %alternative option
  temp=reshape(bsxfun(@minus,reshape(temp,[],Nseg2/2,size(temp,2)),mean(reshape(temp,[],Nseg2/2,size(temp,2)),2)),size(temp));
else
  notchfilt = mod(fs,1/(Nseg2*params.lEchoSpacing));
  notchfilt = notchfilt < (fs(2)-fs(1)) | notchfilt > (1/(Nseg2*params.lEchoSpacing)-(fs(2)-fs(1)));
  notchfilt = ifftshift(~notchfilt);
  temp=ifft(fft(temp).*repmat(notchfilt(:),[1 size(temp,2)]));
end
%plot(fs*60,sqrt(sum(abs(fftshift(fft(temp))).^2,2))),xlim([0 120])
[~,~,~,Hidx,~,~]=self_gate(temp,[55 109]/60,cbins,[10 49]/60,1,frame_rate,Nseg2/2);
takeme=[];
for j=1:5
    takeme=[takeme; ones(Nseg/2,1); zeros(missing(j),1)];
end
takeme=repmat(takeme,[ceil(numel(Hidx)/numel(takeme)) 1]);
takeme=logical(takeme(1:numel(Hidx)));
Hidx=Hidx(takeme);
wall_clock=repmat(vec(repmat(1:5,[Nseg/2 1])),[numel(Hidx)/(Nseg/2*5) 1]);
wall_clock=wall_clock(1:numel(Hidx));
Segidx = repmat(1:(Nseg/2),[1 numel(Hidx)/Nseg*2]).';