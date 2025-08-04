function plotPhi(Phi,dt,idx,titleText)

if nargin < 3
    idx = 1:min(size(Phi,1),4);
end
if nargin < 2
    dt = 0.0035;
end

Nrows   = min(length(idx),size(Phi,1));
Npoints = min(3600,size(Phi,2));
df      = 1/(dt*size(Phi,2));

figure;

for n = 1:Nrows
    subplot(Nrows,2,2*n-1); plot(dt:dt:dt*Npoints,realify(Phi(idx(n),1:Npoints))); axis tight;
    subplot(Nrows,2,2*n);   plot(-df*size(Phi,2)/2:df:df*size(Phi,2)/2-df,abs(fftshift(fft(Phi(idx(n),:)))));axis([-3 3 0 Inf]);
end

if nargin > 3
    subplot(Nrows,2,1);title(['Temporal Basis (Sec): ' titleText]);
    subplot(Nrows,2,2);title(['Spectrum (Hz): ' titleText]);
else
    subplot(Nrows,2,1);title('Temporal Basis (Sec)');
    subplot(Nrows,2,2);title(['Spectrum (Hz)']);
end