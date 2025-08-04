function img_distorcor = runDistorCor(img, mainpath, twix_obj, hdr)

if nargin < 4
    hdr = [];
end

%% absolute path + filename prefix for the temporary files 

prefix_path = '/tmp/DistorCorTemp';


%% add DistorCor folder to library path 

% strLibPath = absolute path to the folder that contains the executable and .so libraries
path_temp  = what(mainpath);                                    % path to the multitasking recon tool, set at the beginning of multitasking_reon.m 
strLibPath = [path_temp.path '/supporting/utils/DistorCor/lib'];    % I put the DistorCor folder under multitasking/supporting/utils/ hence 'mainpath' is used here

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
prepDistorCor(prefix_real,single(real(img)),twix_obj, hdr);

% run correction. there will be error messages regarding UTrace, please ignore
strRunDistorCor = [strLibPath '/DistorCor ' prefix_real];  % Linux command >> DistorCor [file path+name] 
[~] = system(strRunDistorCor);

% if image is complex, repeat for the imaginary part
if ~isreal(img)
    prefix_imag = [prefix_path '_imag'];
    prepDistorCor(prefix_imag,single(imag(img)),twix_obj, hdr);
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
