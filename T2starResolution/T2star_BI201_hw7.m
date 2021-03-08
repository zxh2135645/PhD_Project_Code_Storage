% BI201 - HW7 problem 1
% 1. Estimating T2*
%    T2* is the result of a range of off-resonance frequencies being present within a voxel, and can be
%    modeled by averaging Net Magnetization vectors which have a range of off-resonances. For this
%    problem, you will need to design a T2* simulator by creating a set of net magnetizations that have
%    the same T2 value but a range of off-resonances, and average together to model the voxel signal.
%    Include the time evolution of the spins, so you can plot the signal magnitude as a function of time.
%    Plot the voxel signal and estimate the T2* for the following conditions at 3T:
%    (a) T2 = 80 ms, off-resonance range of ±0.02 ppm
%    (b) T2 = 80 ms, off-resonance range of ±0.05 ppm
%    (c) T2 = 80 ms, off-resonance range of ±0.1 ppm
%    (d) Bonus question: A more realistic simulation would be to use a Gaussian distribution of offresonance
%    frequencies within a voxel to average over. Redesign your simulator and repeat the
%    above plots/estimates, but with the ranges representing the standard deviation of the Gaussian
%    distribution.

clear,clc,close all;
t = linspace(0, 0.2, 401);
B0 = 3e4;%G
T2= [0.08 0.08 0.08];
dfmax = [0.02e-6 0.05e-6 0.1e-6];
%uncomment for 1d
% SD = [0.02 0.05 0.1];
% fit T2* values: T2s = zeros(1,3);%
for n=1:3
maxdf = dfmax(n) * B0 * 4257;
df = linspace(-maxdf, maxdf);
%uncomment for 1d
%dx = linspace(-SD(n),SD(n));
%gaussw = 1/(SD(n)*sqrt(2*pi))*exp(-dx.^2/(2*SD(n)^2));
dfmat = repmat(df.', 1, length(t));
tmat = repmat(t, length(df), 1);
Mxy = exp(-tmat/T2(n)) .* exp(1i*2*pi*dfmat.*tmat);
%uncomment for 1d
%Mxy = Mxy .* repmat(gaussw.', 1, length(t));

objfcn = @(v)v(1)+v(2)*exp(-1/v(3)*t)-real(sum(Mxy,1))/max(real(sum(Mxy,1)));
x0 = 1e-3*[1;1;1];
[v, resnorm] = lsqnonlin(objfcn,x0,0,0.1);
T2s(n) = v(3);
subplot(3,1,n), plot(t, real(sum(Mxy,1))/max(real(sum(Mxy,1))), t, exp(-t/T2s(n)), '--')
legend('signal', 'T2* fit')
title(['Estimated T_2* = ' num2str(T2s(n)*1e3) ' ms'])
end


%% 1d
t = linspace(0, 0.2, 401);
B0 = 3e4;%G
T2= [0.08 0.08 0.08];
dfmax = [0.02e-6 0.05e-6 0.1e-6];
SD = [0.02 0.05 0.1];
% fit T2* values: T2s = zeros(1,3);%
for n=1:3
maxdf = dfmax(n) * B0 * 4257;
df = linspace(-maxdf, maxdf);
%uncomment for 1d
dx = linspace(-SD(n),SD(n));
gaussw = 1/(SD(n)*sqrt(2*pi))*exp(-dx.^2/(2*SD(n)^2));
dfmat = repmat(df.', 1, length(t));
tmat = repmat(t, length(df), 1);
Mxy = exp(-tmat/T2(n)) .* exp(1i*2*pi*dfmat.*tmat);
%uncomment for 1d
Mxy = Mxy .* repmat(gaussw.', 1, length(t));

% Mxy is then normalized
objfcn = @(v)v(1)+v(2)*exp(-1/v(3)*t)-real(sum(Mxy,1))/max(real(sum(Mxy,1)));
x0 = 1e-3*[1;1;1];
[v, resnorm] = lsqnonlin(objfcn,x0,0,0.1);
T2s(n) = v(3);
subplot(3,1,n), plot(t, real(sum(Mxy,1))/max(real(sum(Mxy,1))), t, exp(-t/T2s(n)), '--');
legend('signal', 'T2* fit');
title(['Estimated T_2* = ' num2str(T2s(n)*1e3) ' ms']);
end

%% What if lorenzian
t = linspace(0, 0.2, 401);
B0 = 3e4;%G
T2= [0.08 0.08 0.08];
dfmax = [0.02e-6 0.05e-6 0.1e-6];
gamma = [0.02 0.05 0.1]/2;
% fit T2* values: T2s = zeros(1,3);%
for n = 1:3
maxdf = dfmax(n) * B0 * 4257;
df = linspace(-maxdf, maxdf);
%uncomment for 1d
dx = linspace(-gamma(n),gamma(n));
% gaussw = 1/(gamma(n)*sqrt(2*pi))*exp(-dx.^2/(2*gamma(n)^2));
lorenz = 1 ./ (pi * gamma(n) .* (1 + (dx.^2 ./ gamma(n).^2)));
dfmat = repmat(df.', 1, length(t));
tmat = repmat(t, length(df), 1);
Mxy = exp(-tmat/T2(n)) .* exp(1i*2*pi*dfmat.*tmat);

%uncomment for 1d
Mxy = Mxy .* repmat(lorenz.', 1, length(t));

% Mxy is then normalized
objfcn = @(v)v(1)+v(2)*exp(-1/v(3)*t)-real(sum(Mxy,1))/max(real(sum(Mxy,1)));
x0 = 1e-3*[1;1;1];
[v, resnorm] = lsqnonlin(objfcn,x0,0,0.1);
T2s(n) = v(3);
subplot(3,1,n), plot(t, real(sum(Mxy,1))/max(real(sum(Mxy,1))), t, exp(-t/T2s(n)), '--');
legend('signal', 'T2* fit');
title(['Estimated T_2* = ' num2str(T2s(n)*1e3) ' ms']);
end