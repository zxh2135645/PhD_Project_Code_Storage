function I = waveletfigmake2D(WU,waveS,wname,wlev)

WU = sum(abs(WU),1);

I = appcoef2(WU,waveS,wname);
for j = wlev:-1:1
    I = [I, detcoef2('h',WU,waveS,j); detcoef2('v',WU,waveS,j), detcoef2('d',WU,waveS,j)];
end
