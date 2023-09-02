function realfun=realify(x,type)
% realify(x)
% realify(x,'cols')
% realify(x,'rows')

if nargin < 2
  type='standard';
elseif isempty(type)
  type='standard';
end
type=lower(type);

switch type
  case 'standard'
    [~,~,v]=svd([real(x(:)), imag(x(:))],'econ');
    theta=angle(complex(v(1,1),v(2,1)));
    realfun=real(x*exp(-1i*theta));
    realfun=realfun/sign(mean(realfun(:))); %mostly positive
  case 'cols'
    realfun=zeros(size(x));
    for j=1:size(x,2)
      realfun(:,j)=realify(x(:,j));
    end
  case 'rows'
    realfun=zeros(size(x));
    for j=1:size(x,1)
      realfun(j,:)=realify(x(j,:));
    end
  otherwise %unrecognized
    realfun=realify(x);
end

end