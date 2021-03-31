function y=wmedfilt2(x)

y=zeros(size(x));

w=sqrt(bsxfun(@plus,(-1:1).'.^2,(-1:1).^2));
w=normpdf(w,0,max(abs(w(:)))/2); %Gaussian weighting function
w=w/sum(w(:));

for frame=1:size(x,3)
  
  [sorted,wind]=sort(cat(3,...
    circshift(x(:,:,frame),[-1 -1]),circshift(x(:,:,frame),[0 -1]),circshift(x(:,:,frame),[1 -1]),...
    circshift(x(:,:,frame),[-1  0]),circshift(x(:,:,frame),[0  0]),circshift(x(:,:,frame),[1  0]),...
    circshift(x(:,:,frame),[-1  1]),circshift(x(:,:,frame),[0  1]),circshift(x(:,:,frame),[1  1])),3);
  
  wmed=cumsum(w(wind),3);
  
  for j=1:size(x,1)
    for k=1:size(x,2)
      temp=find(wmed(j,k,:)==.5,1);
      if isempty(temp)
        y(j,k,frame)=sorted(j,k,find(wmed(j,k,:)>.5,1));
      else
        y(j,k,frame)=sorted(j,k,temp)+sorted(j,k,temp+1);
      end
    end
  end
  
end

end