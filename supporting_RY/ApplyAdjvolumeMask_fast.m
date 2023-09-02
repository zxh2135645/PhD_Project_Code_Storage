% get the shimming box from the struct sAdjData and then generate the mask

% 08/16/2019 Mona,Tim
function Mask_matrix = ApplyAdjvolumeMask(twix)
%% extract the info from the data file  and create the new blank matrix, saving useful information in the struct 'info'

info = get_info(twix);
x = twix.image.NCol;
y = twix.image.NLin;
z = twix.image.NPar*twix.image.NSli;

%% get affine matrix 
% get the left up corner cordf of the img
voxel = [info.img.voxel_size(1) info.img.voxel_size(2) info.img.voxel_size(3)];
voxel = repmat(voxel,3,1);
rotation_img = rotation_matrix(info.img.dSag,info.img.dCor,info.img.dTra,info.img.dRot);
center_img = [info.img.dSag_center info.img.dCor_center info.img.dTra_center];
temp_affine = [[rotation_img.*voxel center_img'];[0 0 0 1]];
% corner = temp_affine*[-info.img.arraysize(1)/2 info.img.arraysize(2)/2 -info.img.arraysize(3)/2 1]';
corner = temp_affine*[-x/2 y/2 z/2 1]';

% get the img affine matrix
corner_img = [corner(1) corner(2) corner(3)];
img_affine = [[rotation_img.*voxel corner_img'];[0 0 0 1]];

%% patient position

% get the shimmingbox affine matrix 
rotation_shimbox = rotation_matrix(info.shimbox.dSag,info.shimbox.dCor,info.shimbox.dTra,info.img.dRot);
center_shimbox = [info.shimbox.dSag_center info.shimbox.dCor_center info.shimbox.dTra_center];
shim_affine = [[rotation_shimbox center_shimbox'];[0 0 0 1]];

% calculate the vertex of the shimmingbox
vertex(1).x = 0 + info.shimbox.dReadoutFOV/2;
vertex(1).y = 0 - info.shimbox.dPhaseFOV/2;
vertex(1).z = 0 + info.shimbox.dThickness/2;

vertex(2).x = 0 + info.shimbox.dReadoutFOV/2;
vertex(2).y = 0 + info.shimbox.dPhaseFOV/2;
vertex(2).z = 0 + info.shimbox.dThickness/2;

vertex(3).x = 0 + info.shimbox.dReadoutFOV/2;
vertex(3).y = 0 - info.shimbox.dPhaseFOV/2;
vertex(3).z = 0 - info.shimbox.dThickness/2;

vertex(4).x = 0 - info.shimbox.dReadoutFOV/2;
vertex(4).y = 0 - info.shimbox.dPhaseFOV/2;
vertex(4).z = 0 + info.shimbox.dThickness/2;   
for i = 1:4
    vertex(i).cordf = shim_affine*[vertex(i).x vertex(i).y vertex(i).z 1]';
    vertex(i).cordf = vertex(i).cordf(1:3);
end

box_center = [info.shimbox.dSag_center info.shimbox.dCor_center info.shimbox.dTra_center]';

%% decide and generate
% for each voxel in the img matrix, decide whether it's in the shimming box
tic
word_cordf = zeros([4 x*y*z]);
a = 1;
for i = 1:z
    for j = 1:y
        % slicer and phase direction inversed
        temp = [[0:1:x-1]
            [repmat(-j+1,1,x)]
            [repmat(-i+1,1,x)]
            [repmat(1,1,x)]];
        word_cordf(:,(a-1)*x+1:a*x) = temp;
        a = a+1;
    end
end

word_cordf = img_affine*word_cordf;
result = in_shimbox(word_cordf(1:3,:),box_center,info.shimbox.dReadoutFOV/2,info.shimbox.dPhaseFOV/2,info.shimbox.dThickness/2,vertex(1).cordf,vertex(2).cordf,vertex(3).cordf,vertex(4).cordf);
Mask_matrix = reshape(result,[x,y,z]);
toc

end

%% additional function

% calculate the rotation_matrix 
function rotationMatrix = rotation_matrix(dsag,dcor,dtra,dRot)
    B = [dsag,dcor,dtra]';
    A = [0 0 1]';
    C = cross(A, B);
    if(norm(C) == 0)
        rotationtemp1 = [[1 0 0]
            [0 1 0]
            [0 0 1]];
    else
        theta = acos((A'*B) / ( norm(A)*norm(B) ));
        r = C / norm(C) * theta;
       % rotationtemp1 = rotationVectorToMatrix(r);
        rotationtemp1 = Vector2RotationMatrix(r)';
    end
    
    % do the drot
    if dRot ~= 0
        r = B / norm(B) * dRot;
%         rotationtemp2 = rotationVectorToMatrix(r);
        rotationtemp2 = Vector2RotationMatrix(r)';
        rotationMatrix = rotationtemp2*rotationtemp1;
    else
        rotationMatrix = rotationtemp1;
    end
    
end


% decide whether the point is in a volume
function result = in_shimbox(x,y,x0,y0,z0,vertex0,vertex1,vertex2,vertex3)
    % x is point to judge
    % y is center point
    % x0 is dis from center along x
    % y0 is dis from center along y
    % z0 is dis from center along z
    
    normX = cross(vertex1-vertex0,vertex2-vertex0);
    normZ = cross(vertex1-vertex0,vertex3-vertex0);
    normY = cross(vertex2-vertex0,vertex3-vertex0);
    normX1 = norm(normX);
    normY1 = norm(normY);
    normZ1 = norm(normZ);
    % for fast calculation
    normX = repmat(normX,1,size(x,2));
    normY = repmat(normY,1,size(x,2));
    normZ = repmat(normZ,1,size(x,2));
    y = repmat(y,1,size(x,2));
    
    lenX = abs(dot(x-y,normX)/norm(normX1));
    lenY = abs(dot(x-y,normY)/norm(normY1));
    lenZ = abs(dot(x-y,normZ)/norm(normZ1));
    
    x0 = repmat(x0,1,size(x,2));
    y0 = repmat(y0,1,size(x,2));
    z0 = repmat(z0,1,size(x,2));
    
    resultx = (lenX <= x0);
    resulty = (lenY <= y0);
    resultz = (lenZ <= z0);
    result = resultx & resulty & resultz;
end

% get the necessary info from the twix_obj
function info = get_info(twix_in)
    prot = twix_in.hdr.MeasYaps;
    
    % information about the image
    img_info = prot.sSliceArray.asSlice{1,1};
    if isfield(img_info,'dThickness')
        info.img.dThickness = img_info.dThickness;
    else
        info.img.dThickness = 0;
    end
    
    if isfield(img_info,'dPhaseFOV')
        info.img.dPhaseFOV = img_info.dPhaseFOV;
    else
        info.img.dPhaseFOV = 0;
    end
    
    if isfield(img_info,'dReadoutFOV')
        info.img.dReadoutFOV = img_info.dReadoutFOV;
    else
        info.img.dReadoutFOV = 0;
    end
    
    if isfield(img_info,'sNormal')
        if isfield(img_info.sNormal,'dTra')
            info.img.dTra = img_info.sNormal.dTra;
        else
            info.img.dTra = 0;
        end

        if isfield(img_info.sNormal,'dSag')
            info.img.dSag = img_info.sNormal.dSag;
        else
            info.img.dSag = 0;
        end

        if isfield(img_info.sNormal,'dCor')
            info.img.dCor = img_info.sNormal.dCor;
        else
            info.img.dCor = 0;
        end
    else
        info.img.dTra = 0;
        info.img.dSag = 0;
        info.img.dCor = 0;
    end
    
    if isfield(img_info,'sPosition')
        if isfield(img_info.sNormal,'dTra')
            info.img.dTra_center = img_info.sPosition.dTra;
        else
            info.img.dTra_center = 0;
        end

        if isfield(img_info.sNormal,'dSag')
            info.img.dSag_center = img_info.sPosition.dSag;
        else
            info.img.dSag_center = 0;
        end

        if isfield(img_info.sNormal,'dCor')
            info.img.dCor_center = img_info.sPosition.dCor;
        else
            info.img.dCor_center = 0;
        end
    else
        info.img.dTra_center = 0;
        info.img.dSag_center = 0;
        info.img.dCor_center = 0;
    end  
    
    if isfield(img_info,'dInPlaneRot')
        info.img.dRot = img_info.dInPlaneRot;
    else
        info.img.dRot = 0;
    end
    
    % !!!!!!!!!!!!!!!!!!check the calculation of voxel size!!!!!!!!!!!!!!!!
    if isfield(twix_in.hdr.Meas,'NImageCols')
        info.img.arraysize(1) = twix_in.hdr.Meas.NImageCols;
    else
        info.img.arraysize(1) = twix_in.image.NCol/2;
    end
    if isfield(twix_in.hdr.Meas,'NImageLins')
        info.img.arraysize(2) = twix_in.hdr.Meas.NImageLins;
    else
        info.img.arraysize(2) = twix_in.image.NLin;
    end
    if isfield(twix_in.hdr.Meas,'NImagePar')
        info.img.arraysize(3) = twix_in.hdr.Meas.NImagePar;    
    elseif isfield(twix_in.hdr.Meas,'NProtPar')
        info.img.arraysize(3) = twix_in.hdr.Meas.NProtPar*twix_in.hdr.Meas.NProtSlc;
    else
        info.img.arraysize(3) = twix_in.image.NPar*twix_in.image.NSli;
    end
    info.img.voxel_size(1) = info.img.dReadoutFOV/info.img.arraysize(1);
    info.img.voxel_size(2) = info.img.dPhaseFOV/info.img.arraysize(2);
    info.img.voxel_size(3) = info.img.dThickness/info.img.arraysize(3);
    
    % information about the shimming box
    shim_info = prot.sAdjData.sAdjVolume;
    if isfield(shim_info,'dThickness')
        info.shimbox.dThickness = shim_info.dThickness;
    else
        info.shimbox.dThickness = 0;
    end
    
    if isfield(shim_info,'dPhaseFOV')
        info.shimbox.dPhaseFOV = shim_info.dPhaseFOV;
    else
        info.shimbox.dPhaseFOV = 0;
    end
    
    if isfield(shim_info,'dReadoutFOV')
        info.shimbox.dReadoutFOV = shim_info.dReadoutFOV;
    else
        info.shimbox.dReadoutFOV = 0;
    end
    
    if isfield(shim_info,'sNormal')
        if isfield(shim_info.sNormal,'dTra')
            info.shimbox.dTra = shim_info.sNormal.dTra;
        else
            info.shimbox.dTra = 0;
        end

        if isfield(shim_info.sNormal,'dSag')
            info.shimbox.dSag = shim_info.sNormal.dSag;
        else
            info.shimbox.dSag = 0;
        end

        if isfield(shim_info.sNormal,'dCor')
            info.shimbox.dCor = shim_info.sNormal.dCor;
        else
            info.shimbox.dCor = 0;
        end
    else
        info.shimbox.dTra = 0;
        info.shimbox.dSag = 0;
        info.shimbox.dCor = 0;
    end
    
    if isfield(shim_info,'sPosition')
        if isfield(shim_info.sPosition,'dTra')
            info.shimbox.dTra_center = shim_info.sPosition.dTra;
        else
            info.shimbox.dTra_center = 0;
        end

        if isfield(shim_info.sPosition,'dSag')
            info.shimbox.dSag_center = shim_info.sPosition.dSag;
        else
            info.shimbox.dSag_center = 0;
        end

        if isfield(shim_info.sPosition,'dCor')
            info.shimbox.dCor_center = shim_info.sPosition.dCor;
        else
            info.shimbox.dCor_center = 0;
        end
    else
        info.shimbox.dTra_center = 0;
        info.shimbox.dSag_center = 0;
        info.shimbox.dCor_center = 0;
    end  
        
    if isfield(shim_info,'dInPlaneRot')
        info.shimbox.dRot = shim_info.dInPlaneRot;
    else
        info.shimbox.dRot = 0;
    end
end
