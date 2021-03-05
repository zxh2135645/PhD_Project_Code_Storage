function y=wave3d(x,level,wname)

if nargin < 2
  level=3;
elseif isempty(level)
  level=3;
end

checkLevel=true;
while checkLevel
  if sum(mod(size(x(:,:,:,1)),2^level))>0
    level=level-1;
  else
    checkLevel=false;
  end
end

if nargin < 3
  wname='db4';
elseif isempty(wname)
  wname='db4';
end

for j=1:size(x,4)
  y(1,j)=dwt3(x(:,:,:,j),wname,'mode','per');
  for l=2:level
    y(l,j)=dwt3(y(l-1,j).dec{1},wname,'mode','per');
  end
end

return