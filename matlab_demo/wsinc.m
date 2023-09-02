function h = wsinc(tbw, ns)

% rf = wsinc(tbw, ns)
%   
%   tbw  --  time bandwidth product
%   ns    --  number of samples
%   h   -- windowed sinc function, normalized so that sum(h) = 1

xm = (ns-1)/2;
x = [-xm:xm]/xm;
h = sinc(x*tbw/2).*(0.54+0.46*cos(pi*x));

h = h/sum(h);
