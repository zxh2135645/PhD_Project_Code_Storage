%% Bloch simulation
TR = 0.0035;
TEs = [0 30 40 50 60 0 30 40 50 60]*1e-3;
TEs = 0;

beta = 0.5;
alphaArray = beta*[3 10 3 10 3 10 3 10 3 10]*pi/180;
alphaArray = beta*[5 5 5 5 5 5 5 5 5 5]*pi/180;

Nseg = 710;
moduleLength = 1;

T1 = 1;
T2 = 0.05;
BIR = 1;
BT2 = 1;
e1 = @(T1) exp(-TR/T1);

repetitions = 100;

Msig = zeros(1,Nseg*moduleLength*repetitions);

M = 1;
count = 0;
for rep = 1: repetitions
    for shot = 1:moduleLength
        FA = alphaArray(mod(shot-1,moduleLength)+1);
        M  = M * cos(BIR*pi)*((cos(BT2*pi/2)^2+sin(BT2*pi/2)^2*exp(-TEs(mod(shot-1,moduleLength)+1)/T2)));
        for n = 1:Nseg
            count = count + 1;
            Msig(count) = M*sin(FA);
            M = M * cos(FA)*e1(T1) + (1-e1(T1));
        end
    end
end

figure;plot(Msig(end-Nseg*moduleLength+1:end));