%% Gaussian, interpolated as in Bieri and Scheffler, 2006
% return LUT interpolated from -1.5kHz to +1.5KHz
% Units of G are us (microseconds)

function [ff,G] = GaussianLineShape(T2b)
% T2b is us

%%% define frequency range
n=512;
ff = linspace(-30e3,30e3,n);

%%% compute G for this range
G = zeros([n 1]);

for ii=1:n
    G(ii) = Gaussian(ff(ii));
end

%%% interpolate
po = find(abs(ff)<1.5e3); % points to interpolate
% Because it's close to 0...
pu = find((abs(ff)>1.5e3)&(abs(ff)<2e3)); % points to use

Gi = spline(ff(pu),G(pu),ff(po));
G(po) = Gi;

G = G*1e6; % us

    function gg = Gaussian(f) % Gaussian Line
        g = exp(-(2*pi*f*T2b).^2 / 2);
        gg = T2b / sqrt(2*pi) .* g;
    end

end