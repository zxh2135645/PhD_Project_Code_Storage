function [cc mask]=fuzzybiassmooth(image,k) 
global biasout
x=image;
original=x;
sigmas=1;  %%% Spatial (Geometric) Sigma.
ksize=3;   %%% Kernal Size.
gk=gauss_ker2D(sigmas,ksize);
image=x;
[m n]=size(image);
image=image(image>5);
inidi=isempty(image);
if inidi==1
    cc=zeros(1,k);
    s=size(image,1);
    bias=zeros(s,1);
    mask=zeros(size(original));
       return
end
image=double(image);
s=size(image,1);
mask=zeros(s,1);
distance=zeros(s,k);
dist=zeros(s,k);
newdist=zeros(s,k);
u=zeros(s,k);
v=zeros(s,k);
bias=zeros(s,1);
norm=0.08;
bias1=zeros(s,k);
num=zeros(s,k);
den=zeros(s,k);
bias2=zeros(s,k);
gk1=[1/9 1/9 1/9;1/9 1/9 1/9;1/9 1/9 1/9];
%gk1=[1/9 1/9 1/9 1/9 1/9 1/9;1/9 1/9 1/9 1/9 1/9 1/9;1/9 1/9 1/9 1/9 1/9 1/9];
[muk,outk]=kmeans(uint8(image),k);
         spatial=image;
    init=zeros(1,k);
    for i=1:s
        bias(i,1)=0.00000000000000000000000000001;
    end
     cc=zeros(1,k);
      ccc=zeros(1,k);
      tmpMatrix=zeros(1,k);
    for i=1:k
        init(i)=muk(i);
    end
   ttFcm=0;
while(ttFcm<500)
    ttFcm=ttFcm+1;
    if ttFcm==1 
       for i=1:k
        cc(i)=init(i);
       end
    end
    for i=1:k
        distance(:,i)=(spatial)-bias-cc(i);
        distance(:,i)=distance(:,i).*distance(:,i)+0.000000000000000000001;
        dist(:,i)=1./distance(:,i);
    end
        distsum=sum(dist,2);
        for z=1:k
            newdist(:,z)=distance(:,z).*distsum;
            u(:,z)=1./newdist(:,z);
            smooth=bias;
            t=1;
        for j=1:n
          for i=1:m
            if original(i,j)>5
                biassmooth(i,j)=smooth(t,1);
                t=t+1;
            end 
            if original(i,j)<=5
                biassmooth(i,j)=0;
            end
           end
        end
     % biassmooth=conv2(biassmooth,gk1,'same');
    %   biassmooth= FAST_NLM_II(biassmooth,2,5,10);
        for i=1:m
            for j=1:n
                 if original(i,j)<=5
                    biassmooth(i,j)=0;
                 end
            end
        end
        biasout=biassmooth;
        bias=biassmooth(biassmooth~=0);
        smooth=u(:,z);
            t=1;
         for j=1:n
          for i=1:m
            if original(i,j)>5
                usmooth(i,j)=smooth(t,1);
                t=t+1;
            end 
            if original(i,j)<=5
                usmooth(i,j)=0;
            end
           end
        end
   %   usmooth=conv2(usmooth,gk1,'same');
     % usmooth= FAST_NLM_II(usmooth,2,5,0.38);
        for i=1:m
            for j=1:n
                 if original(i,j)<=5
                    usmooth(i,j)=0;
                 end
            end
        end
        u(:,z)=usmooth(usmooth~=0);
      

            % u(:,i) = FAST_NLM_II(u(:,i),2,5,0.36357);
           % u(:,z)=conv2(u(:,z),gk1,'same');
          % num(:,z)=(u(:,z).*u(:,z).*(spatial));
           %den(:,z)=(u(:,z).*u(:,z));
         %  v(:,z)=num(:,z)./den(:,z);
            ccc(z)=sum(sum(u(:,z).*u(:,z).*(spatial-bias)))/sum(sum(u(:,z).*u(:,z)));
            cccmat(:,z)=repmat(ccc(z),[s 1]);
                bias1(:,z)=u(:,z).*u(:,z).*(spatial-cccmat(:,z));
           bias2(:,z)=u(:,z).*u(:,z);
        end
        bias3=sum(bias1,2)./sum(bias2,2);
                 % bias=abs(bias);
             tmpMatrix=abs(cc-ccc).*abs(cc-ccc);
                
    %%%%%%%%%%%%%%%
  if max(tmpMatrix)<0.001
       break;
  else
     cc=ccc;
     bias=bias3;
  end
end
for z=1:k
  smooth=u(:,z);
            t=1;
         for j=1:n
          for i=1:m
            if original(i,j)>5
                usmooth(i,j)=smooth(t,1);
                t=t+1;
            end 
            if original(i,j)<=5
                usmooth(i,j)=0;
            end
           end
        end
      usmooth=conv2(usmooth,gk1,'same');
     
      

   % usmooth= FAST_NLM_II(usmooth,2,5,0.6);
        for i=1:m
            for j=1:n
                 if original(i,j)<=5
                    usmooth(i,j)=0;
                 end
            end
        end
        u(:,z)=usmooth(usmooth~=0);
end
%%%%%%%%%%%%%%%%%%
for j=1:s
    c=u(j,:);
    a=find(c==max(c));  
    mask(j,1)=a(1);
end
copy=original;
w=1;
for j=1:n
    for i=1:m
            if original(i,j)>5
                copy(i,j)=mask(w,1);
                w=w+1;
            end
            if original(i,j)<=5
                copy(i,j)=0;
             end
     end
end
mask=copy;
old=cc;
new=sort(cc);
h=1;
for i=1:k
    for j=1:k
         t=old(1,i)-new(1,j);
        if t==0
            label(h)=j;
            h=h+1;
        end
    end
end
for p=1:k
for i=1:m
    for j=1:n
        if mask(i,j)==p
            mask(i,j)=label(p);
        end
       
    end
end
end


