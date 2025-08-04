function [pSWI,nSWI,tSWI,phase_hp] = SWI(reconMag,reconPhase,filterSize,threshold,power_weighting,QSM,chi_1,chi_2)

if nargin > 5 && ~isempty(QSM)
    flagDotSWI = true;
    if nargin < 7
        chi_1 = 0;
    end
    if nargin < 8
        chi_2 = 550;
    end
else
    flagDotSWI = false;
end

if nargin < 5
    power_weighting = 4;
end
if nargin < 4 || threshold < 0
    threshold = pi;
end
if nargin < 3 || filterSize < 0
    filterSize = 65;
end

% convert input to double
reconMag   = double(reconMag);
reconPhase = double(reconPhase);
recon = reconMag .* exp(1j*reconPhase);

% phase filtering
disp('Calculating SWI... ');
reconTemp_lowpass = zeros(size(reconPhase));
lowpass_filter = zeros(size(reconPhase(:,:,1,1)));
lowpass_filter(1:filterSize,1:filterSize) = hamming(filterSize)*hamming(filterSize)'; 
lowpass_filter = circshift(lowpass_filter,[-floor(filterSize/2) -floor(filterSize/2)]);
for n = 1:size(reconPhase,4)
    for slice = 1:size(reconPhase,3)
        reconTemp = exp(1j.*reconPhase(:,:,slice,n));
%         reconTemp_lowpass(:,:,slice,n) = filter2(lowpass_filter,reconTemp);   
        reconTemp_lowpass(:,:,slice,n) = fftshift(ifft2(fft2(ifftshift(reconTemp)).*lowpass_filter));   
    end
end
phase_hp = exp(1j.*reconPhase) ./ reconTemp_lowpass;
phase = angle(phase_hp);

% positive phase mask
pSmask          = (threshold-phase)/threshold;
pSmask(phase>threshold)	= 0; 
pSmask(phase<0) = 1;
% negative phase mask
dSmask          = (threshold+phase)/threshold;
dSmask(phase>0) = 1; 
dSmask(phase<-threshold)= 0;

% generating SWI
pSWI = bsxfun(@times,reconMag,pSmask.^power_weighting);
nSWI = bsxfun(@times,reconMag,dSmask.^power_weighting);
SWI = bsxfun(@times,reconMag,phase.^power_weighting);

% generating tSWI
if flagDotSWI
    disp('Calculating tSWI... ');
    % process QSM_mask
    QSM_mask = 1 - (QSM-chi_1)/(chi_2-chi_1);
    QSM_mask(QSM_mask>1) = 1;
    QSM_mask(QSM_mask<0) = 0;
    tSWI = abs(reconMag).*(QSM_mask.^2);
else
    disp('No QSM result. tSWI not available.')
    tSWI = [];
end


% % convert phase to complex value
% phase_cplx = exp(1i*phase);
% 
% dim = size(phase_cplx);
% 
% % for single echo
% if length(dim) < 4
%     dim(4) = 1;
% end
% 
% % create a 2D low-pass filtre
% lowpass_filter = hamming(filterSize)*hamming(filterSize)'; 
% 
% swi_phase = zeros(dim);
% % loops all echoes
% for kt = 1:dim(4)
%     disp(['Processing echo ' num2str(kt) ' ...']);
%     % loops all slices, the thresholding is actually performed on
%     % slice-by-slice basis
%     for kz = 1:dim(3)
%         % get slice with complex-valued data
%         c = phase_cplx(:,:,kz, kt);
%         
%         % filtre the complex-valued slice with the above low-pass filtre
%         c_lowpass = filter2(lowpass_filter, c);   
%         
%         % complex-valued division between the original data and low-passed
%         % data is equivalent to high-pass filtred the original data
%         swi_highpass = c ./ c_lowpass;
%         % KC: compute the high-pass filtred phase from the complex-valued slice
%         swi_phase(:, :, kz, kt) = angle(swi_highpass);
%     end 
% end
% 
% % switch method
% %     case 'default'
% %         
% %     case 'multiecho (testing)'
% %         swi_phase = cumsum(swi_phase,4);
% % end
% 
% % positive phase mask
% pSmask                  = (thres-swi_phase)/thres;
% pSmask(swi_phase>thres)	= 0; 
% pSmask(swi_phase<0)   	= 1;
% % negative phase mask
% dSmask                 	= (thres+swi_phase)/thres;
% dSmask(swi_phase>0)   	= 1; 
% dSmask(swi_phase<-thres)= 0;
% 
% pSWI = bsxfun(@times,magn,pSmask.^m);
% nSWI = bsxfun(@times,magn,dSmask.^m);
% 
% end