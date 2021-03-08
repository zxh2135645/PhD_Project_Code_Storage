function Y=wavelet_figmake(WU)

j=size(WU,1);
wsl=floor(size(WU(j,1).dec{1},3)/2+1);
for k=1:size(WU,2)
  Y(:,:,k)=[WU(j,k).dec{1}(:,:,wsl), WU(j,k).dec{1,2,1}(:,:,wsl); WU(j,k).dec{2,1,1}(:,:,wsl), WU(j,k).dec{1,1,2}(:,:,wsl)];
end

for j=size(WU,1)-1:-1:1
  wsl=floor(size(WU(j,1).dec{1},3)/2+1);
  Yold=Y;
  clear Y;
  for k=1:size(WU,2)
    Y(:,:,k)=[Yold(:,:,k), WU(j,k).dec{1,2,1}(:,:,wsl); WU(j,k).dec{2,1,1}(:,:,wsl), WU(j,k).dec{1,1,2}(:,:,wsl)];
  end
end

Y=sum(abs(Y),3);

return;