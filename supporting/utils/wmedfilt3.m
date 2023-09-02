function y=wmedfilt3(x,h)
% y=wmedfilt3(x)
% y=wmedfilt3(x,h)
%
% x: 3D volume
% h: smoothness parameter (higher is smoother), default = 0.5

if nargin<2
    h=0.5;
end

y=zeros(size(x));

w=sqrt(bsxfun(@plus,bsxfun(@plus,(-1:1).'.^2,(-1:1).^2),reshape((-1:1).^2,1,1,[])));
w=normpdf(w,0,max(abs(w(:)))*h); %Gaussian weighting function
w=w/sum(w(:));

for frame=1:size(x,4)
  
  [sorted,wind]=sort(cat(4,...
    circshift(x(:,:,:,frame),[-1 -1 -1]),circshift(x(:,:,:,frame),[0 -1 -1]),circshift(x(:,:,:,frame),[1 -1 -1]),...
    circshift(x(:,:,:,frame),[-1  0 -1]),circshift(x(:,:,:,frame),[0  0 -1]),circshift(x(:,:,:,frame),[1  0 -1]),...
    circshift(x(:,:,:,frame),[-1  1 -1]),circshift(x(:,:,:,frame),[0  1 -1]),circshift(x(:,:,:,frame),[1  1 -1]),...
    circshift(x(:,:,:,frame),[-1 -1 0]),circshift(x(:,:,:,frame),[0 -1 0]),circshift(x(:,:,:,frame),[1 -1 0]),...
    circshift(x(:,:,:,frame),[-1  0 0]),circshift(x(:,:,:,frame),[0  0 0]),circshift(x(:,:,:,frame),[1  0 0]),...
    circshift(x(:,:,:,frame),[-1  1 0]),circshift(x(:,:,:,frame),[0  1 0]),circshift(x(:,:,:,frame),[1  1 0]),...
    circshift(x(:,:,:,frame),[-1 -1 1]),circshift(x(:,:,:,frame),[0 -1 1]),circshift(x(:,:,:,frame),[1 -1 1]),...
    circshift(x(:,:,:,frame),[-1  0 1]),circshift(x(:,:,:,frame),[0  0 1]),circshift(x(:,:,:,frame),[1  0 1]),...
    circshift(x(:,:,:,frame),[-1  1 1]),circshift(x(:,:,:,frame),[0  1 1]),circshift(x(:,:,:,frame),[1  1 1])),4);
  
  wmed=cumsum(w(wind),4);
  
  for j=1:size(x,1)
    for k=1:size(x,2)
        for l=1:size(x,3)
      temp=find(wmed(j,k,l,:)==.5,1);
      if isempty(temp)
        y(j,k,l,frame)=sorted(j,k,l,find(wmed(j,k,l,:)>.5,1));
      else
        y(j,k,l,frame)=sorted(j,k,l,temp)+sorted(j,k,l,temp+1);
      end
    end
  end
  
end

end