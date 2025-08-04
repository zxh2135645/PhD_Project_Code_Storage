%% DistorCor_demo
%
%  Sample code for image-based 2D distortion correction using external
%  application DistorCor and Siemens MrDistorCor libraries.
%
%  To run this demo two variables are needed in the Matlab Workspace:
%  (DemoData/DistorCorDemo_vars.mat contains a pair of twix_obj/img).
%
%  The path set up section at the beginning of multitasking_reon.m needs to
%  be run first in order for twix_obj to be read correctly.
%
%  twix_obj - rawdata header loaded with mapVBVD()
%
%  img      - 2D/3D/4D image matrix [Npe x Nro x Nz x NR]
%             If matrix is bigger than 1024x1024x64, split it into multiple
%             matrices.
%
%   To have the correct image orientation, the data should be reconstructed 
% 	as described below:
%
%   1)  Apply complex conjugation to kspaceData before recon if patient 
% 		position is head-first. Do kspaceData  = conj(kspaceData) right 
%       after loading data from twix object.
%       
%   2)  For radial data, the trajectory om for NUFFT should be set as
%
%           r   = linspace(-pi, pi, Nkx+1); r(end)=[];
%           om1 = sin(thetas(1:Ntrajs)*pi/180)*r;
%           om2 = cos(thetas(1:Ntrajs)*pi/180)*r;
%           om  = [om1(:), om2(:)];
%
%       where thetas is the dAzimuthalAngle during data acquisition
%
%   -- Hsu-Lei Lee, last edited on 11/06/2020
%


%% absolute path + filename prefix for the temporary files 

prefix_path = '/tmp/DistorCorDemoTemp';


%% add DistorCor folder to library path 

% strLibPath = absolute path to the folder that contains the executable and .so libraries
path_temp  = what(mainpath);                                    	    % path to the multitasking recon tool, set at the beginning of multitasking_reon.m 
strLibPath = [path_temp.path '/supporting/utils/DistorCor/lib'];    	% I put the DistorCor folder under multitasking/supporting/utils/ hence 'mainpath' is used here

% add strLibPath to the Matlab library path
MatlabPath = getenv('LD_LIBRARY_PATH'); 
setenv('LD_LIBRARY_PATH',[strLibPath ';' MatlabPath]);

% if you are running this on recon2 and it complains about cannot load libcrypto.so.1.1, 
% run the following two lines once to create a symbolic link to libcryptopp.so
% strSymLink = ['ln -s /opt/anaconda3/envs/python37/lib/libcrypto.so.1.1 ' strLibPath '/libcrypto.so.1.1'];
% system(strSymLink);

%% export image/header and run distortion correction

% create binary for the real image (prefix_path_real_Img) and save header file (prefix_path_real_DistorCorHdr)
prefix_real = [prefix_path '_real'];
prepDistorCor(prefix_real,single(real(img)),twix_obj);

% run correction. there will be error messages regarding UTrace, please ignore
strRunDistorCor = [strLibPath '/DistorCor ' prefix_real];  % Linux command >> DistorCor [file path+name] 
[~] = system(strRunDistorCor);

% if image is complex, repeat for the imaginary part
if ~isreal(img)
    prefix_imag = [prefix_path '_imag'];
    prepDistorCor(prefix_imag,single(imag(img)),twix_obj);
    strRunDistorCor = [strLibPath '/DistorCor ' prefix_imag];
    [~] = system(strRunDistorCor);
end

clear strRunDistorCor;

%% read corrected image binary back into Matlab

[tempNx,tempNy,tempNz,tempNR] = size(img);
fnameImgCor   = [prefix_real '_Img_DistCor'];
img_distorcor = readImgBinary(fnameImgCor,tempNx*tempNy*tempNz*tempNR); 
img_distorcor = reshape(img_distorcor,tempNx,tempNy,tempNz,tempNR);

if ~isreal(img)
    fnameImgCor        = [prefix_imag '_Img_DistCor'];
    img_distorcor_imag = readImgBinary(fnameImgCor,tempNx*tempNy*tempNz*tempNR); 
    img_distorcor_imag = reshape(img_distorcor_imag,tempNx,tempNy,tempNz,tempNR);
    img_distorcor      = img_distorcor + 1i*img_distorcor_imag;
end

clear fnameImgCor img_distorcor_imag tempNx tempNy tempNz tempNR

%% reset Matlab library path 

setenv('LD_LIBRARY_PATH', MatlabPath );

clear path_temp strLibPath MatlabPath

%% display corrected images

slice = ceil(size(img,3)/2);

figure;
subplot(1,3,1);imagesc(abs(img(:,:,slice)));colormap('gray');axis equal tight off;title('No Correction');
subplot(1,3,2);imagesc(abs(img_distorcor(:,:,slice)));colormap('gray');axis equal tight off;title('DistorCor2D');
subplot(1,3,3);imagesc(abs(img_distorcor(:,:,slice))-abs(img(:,:,slice)));colormap('gray');axis equal tight off;title('Difference');

clear slice

%% delete temporary files

temp1 = [prefix_real '_Img'];
temp2 = [prefix_real '_DistorCorHdr'];
temp3 = [prefix_real '_Img_DistCor'];
delete(temp1,temp2,temp3);

if ~isreal(img)
    temp1 = [prefix_imag '_Img'];
    temp2 = [prefix_imag '_DistorCorHdr'];
    temp3 = [prefix_imag '_Img_DistCor'];
    delete(temp1,temp2,temp3);
end

clear fnameImgCor prefix_path prefix_real prefix_imag img_distorcor_imag temp1 temp2 temp3
