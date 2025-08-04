% this is to test if their setup is:
% S = (w+F*e^(f*t))*e^(j*f_B*t)
row = 8;
col = 98;
S_ref = squeeze(I(row,col,:));
w = water(row, col);
f = fat(row, col);
r = R2s(row, col);
deltaf_B = delta_B0(row, col);
Deltaf = algoParams.species(2).frequency;
Deltaf = 42.58*3*Deltaf;
perc = algoParams.species(2).relAmps;
Expf = exp(1i*2*pi*Deltaf.'*TE');
Fat_spect = perc*Expf;
S = (w + f*Fat_spect).*exp(-r*TE').*exp(1i*2*pi*deltaf_B*TE');
figure, plot(TE,abs(S),TE,abs(S_ref))
figure, plot(TE, atan2(imag(S),real(S)),TE,atan2(imag(S_ref),real(S_ref)))
% OK, the model is
% S = (w+F*e^(f*t))*e^(j*f_B*t)
% i.e. water resonante faster, and fat is slower (-420Hz)