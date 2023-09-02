function [sf,I,maxent]=imdisp(I)
% [sf,I,maxent]=imdisp(I)
% Maximum entropy scaling of input I

I=double(I);

entropy=@(x)-sum(x(:).*log2(x(:)),'omitnan'); %entropy function
pixent=@(x)entropy(hist(uint8(255*I(:)/x),0:255)/numel(I)); %pixelwise entropy for 8-bit integer image

%Coarse search
sfs=linspace(median(I(:)),max(I(:)),21);
ent=zeros(size(sfs));
for j=1:numel(sfs)
  ent(j)=pixent(sfs(j));
end
[~,argmax]=max(ent);
sf=sfs(argmax);

%Fine-tuning
sf=fminsearch(@(x)-pixent(x),sf);
maxent=pixent(sf);

figure,plot(sfs,ent),xlabel('Max. value'),ylabel('Bits/pixel'),ylim([0 8])
hold all
plot(sf,maxent,'+');
hold off

I=uint8(255*I/sf);

if ismatrix(squeeze(I))
  figure,imshow(I)
else
  try
    implay(I)
  catch
  end
end

figure,hist(I(:),0:255)

end
