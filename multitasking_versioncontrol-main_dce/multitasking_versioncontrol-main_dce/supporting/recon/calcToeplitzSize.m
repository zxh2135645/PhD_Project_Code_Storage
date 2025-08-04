function [reconOptions,dataArray] = calcToeplitzSize(params,reconOptions,dataArray)

extractVarFromStruct(params);
extractVarFromStruct(reconOptions);
extractVarFromStruct(dataArray);

vec = @(x) x(:);
energyk = fft2(abs(fbpComposite(:,:,1)).^2)/norm(fbpComposite(:,:,1))^2;
list = [];
minv = log2(Norig);
maxv = log2(N);
vals = round(2.^unique([minv ceil(minv):maxv maxv]));
for height = vals
    for width = vals
        box = padarray(ones(height,width),[N-height, N-width]/2);
        energymap = ifft2(fft2(box).*energyk);
        try
            if flagUseGPUToeplitz
                tic;
                temp = fft2(zeros(2*height,2*width,L*Ncoils,'single','gpuArray'));
                wait(gpudev);
                calctime = 1000*toc;
            else
                tic;
                temp = fft2(zeros(2*height,2*width,L*Ncoils,'single'));
                calctime = 1000*toc;
            end
        catch
            calctime = inf;
        end
        clear temp;
        [maxe,argmax]=max(energymap(:));
        list = [list; 2*height, 2*width, 100*sqrt(maxe), calctime, argmax];
    end
end

cands = find(list(:,3)>99);
[~,argmin] = min(list(cands,4));
Om_size = double(list(cands(argmin),1:2));
fprintf('Toeplitz matrix size = [%dx%d]\n',Om_size(1),Om_size(2));
if flagCommandLine
    figure,plot(list(:,4),list(:,3),'o'),xlabel('Time (ms)'),ylabel('Accuracy (%)');
    hold all;
    plot(list(cands(argmin),4),list(cands(argmin),3),'or'),title(sprintf('Toeplitz size: %d x %d',Om_size(1),Om_size(2)));
end

height = Om_size(1)/2;
width  = Om_size(2)/2;
SE_mask = padarray(ones(height,width),[N-height, N-width]/2);
[vshift,hshift] = ind2sub([N N],list(cands(argmin),5));
SE_mask = circshift(SE_mask,[vshift, hshift]-1);

reconOptions.Om_size = Om_size;
dataArray.SEs_unmasked = SEs;
dataArray.SEs = bsxfun(@times,SEs,SE_mask);

