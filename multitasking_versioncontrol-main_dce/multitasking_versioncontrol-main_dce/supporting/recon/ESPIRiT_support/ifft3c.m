function res = ifft3c(x)

S = size(x);
fctr = S(1)*S(2)*S(3);

x = reshape(x,S(1),S(2),S(3),prod(S(4:end)));

res = zeros(size(x));
for n=1:size(x,4)
    res(:,:,:,n) = sqrt(fctr)*fftshift(fftshift(fftshift(ifft(ifft(ifft(ifftshift(ifftshift(ifftshift(x(:,:,:,n),1),2),3),[],1),[],2),[],3),1),2),3);
end

res = reshape(res,S);

