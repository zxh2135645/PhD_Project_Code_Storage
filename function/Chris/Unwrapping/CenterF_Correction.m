function [unwph_uf_CenF] = CenterF_Correction(unwph_uf,iFreq_raw,Mask)
[x,y,z] = size(unwph_uf);
[counts_unwph,centers_unwph] = hist(unwph_uf(Mask(:)),20);
[~,count_unwph] = max(counts_unwph);
center_unwph = centers_unwph(count_unwph);
[counts_iFreq,centers_iFreq] = hist(iFreq_raw(Mask(:)),20);
[~,count_iFreq] = max(counts_iFreq);
center_iFreq = centers_iFreq(count_iFreq);
Dev = center_unwph - center_iFreq;
jump = round(Dev/2/pi,0);
unwph_uf_CenF = -2*pi*jump.*ones(x,y,z)+unwph_uf;
disp(['Center Frequence Correction ',num2str(jump),' * 2 pi Jump']);
end