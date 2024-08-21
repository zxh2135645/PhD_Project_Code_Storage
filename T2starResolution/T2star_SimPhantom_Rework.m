clear all;
close all;

addpath('../function/create-a-heart-phantom-master/')
addpath('../function/')
%%
%=====================================
%Create phantom for MRiLab simulation
%version 1.0
%This script is adapted from MRiLab
%Author: Hongzhang Chen
%=====================================


% Xinheng Zhang
% Used T2 as a label for color rendering
transmural_array = [0.025, 0.05, 0.10, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8];

theta = pi/2;
Groove = 180;

r_endo = 0.3; % endo -> 30 mm
r_epi =  0.5; % epi  -> 50 mm

Phantom_shape_cell = cell(1,length(transmural_array));

for j = 1:length(transmural_array) %loop the code 4 times to produce 4 different phantoms
    %parameters of the virtual object are stored in the structure array, VObj

    VObj = struct;
    transmural = transmural_array(j);

    VObj.Gyro=[267538030.3797];
    AttributeOpt={'Normal','MT','ME','GM'};
    VObj.Model=AttributeOpt{1};
    VObj.Name='VObj_Custom';
    VObj.Notes=' 5 sphere phantom';
    VObj.Type='Muscle';
    VObj.TypeNum=[1];
    VObj.XDim=[512];
    VObj.XDimRes=[0.002]; % decimeter
    VObj.YDim=[512];
    VObj.YDimRes=[0.002]; % decimeter
    VObj.ZDim=[24];
    VObj.ZDimRes=[0.1]; % centimeter
    [x,y,z]=meshgrid(((-(VObj.XDim-1)/2):((VObj.XDim-1)/2))*VObj.XDimRes,((-(VObj.YDim-1)/2):((VObj.YDim-1)/2))*VObj.YDimRes,((-(VObj.ZDim-1)/2):((VObj.ZDim-1)/2))*VObj.ZDimRes);
    
    %create matrices to store properties (e.g relaxation time t1 and t2)
    rho=zeros([size(x) VObj.TypeNum]);
    t1=zeros([size(x) VObj.TypeNum]);
    t2=zeros([size(x) VObj.TypeNum]);
    t2star=zeros([size(x) VObj.TypeNum]);
    chemshift=zeros(1,VObj.TypeNum);
    typeflag=zeros(1,VObj.TypeNum);
    lineshapeflag=zeros(1,VObj.TypeNum);
    k=zeros([size(x) VObj.TypeNum^2]);
    econ=zeros([size(x) 3]);  
    massden=zeros(size(x));
     
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Cylinder 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %This part has the parameters of the largest sphere 
    %important parameters are set by using 'sturct master'
    Alpha=[1];
    Color='r';
    Notes='Prolate Ellipsoid (Muscle)';
    p.CenterX=[0];
    p.CenterY=[0];
    p.CenterZ=[0];
    p.FaceNum=[100];
    p.Length = [4.8]; % XZ
    p.Radius = r_epi;
    ChemShift=[0];
    K=[10];
    LineShapeFlag=[0];
    master.ECon=[0.29 0.29 0.29];
    master.MassDen=[1090];
    master.Rho=[0.7];
    master.T1=[1.1];
    master.T2=[1];
    master.T2Star=[40];
    TypeFlag=[0];
    TypeIdx=[1];



    %create 3D Cylindral virtual object
    [mask, fvc]=VObjCylinder(p,x,y,z,flag,theta,Groove);
    Obj{1,1}=fvc;
    Obj{2,1}=Color;
    Obj{3,1}=Alpha;
    
    %tissue properties assignment---------------------------------
    %element-wise multiplication for matrices
    rho(:,:,:,TypeIdx)=rho(:,:,:,TypeIdx)+mask.*master.Rho; 
    t1(:,:,:,TypeIdx)=t1(:,:,:,TypeIdx)+mask.*master.T1;
    t2(:,:,:,TypeIdx)=t2(:,:,:,TypeIdx)+mask.*master.T2;
    t2star(:,:,:,TypeIdx)=t2star(:,:,:,TypeIdx)+mask.*master.T2Star;
    massden=massden+mask.*master.MassDen;
    econ(:,:,:,1)=econ(:,:,:,1)+mask.*master.ECon(1);
    econ(:,:,:,2)=econ(:,:,:,2)+mask.*master.ECon(2);
    econ(:,:,:,3)=econ(:,:,:,3)+mask.*master.ECon(3);
    %direct assignment for single values
    chemshift(TypeIdx)=ChemShift;
    typeflag(TypeIdx)=TypeFlag;
    lineshapeflag(TypeIdx)=LineShapeFlag;
    %-------------------------------------------------------------
    
    if ~strcmp(VObj.Model,'Normal')
        for i=1:VObj.TypeNum
            k(:,:,:,(TypeIdx-1)*VObj.TypeNum+i)=k(:,:,:,(TypeIdx-1)*VObj.TypeNum+i)+mask.*K(i);
        end
    end
    p=[];


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Cylinder 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Sphere 2 is radius-variable 
    Alpha=[1];
    Color='b';
    Notes='Sphere 2';
    p.CenterX=[0];
    p.CenterY=[0];
    p.CenterZ=[0];
    p.FaceNum=[100];
    p.Length = [4.8]; % XZ
    p.Radius = r_endo;% increment of radius between 0.01 to 0.018 at a step of 0.002
    ChemShift=[0];
    ECon=[0.5 0.5 0.5];
    K=[10];
    LineShapeFlag=[0];
    MassDen=[1000];
    Rho=[1];
    T1=[1];
    T2=[2];
    T2Star=[0.0175];
    TypeFlag=[0];
    TypeIdx=[1];

    %create 3D sphere virtual object
    [mask, fvc]=VObjCylinder(p,x,y,z,flag,theta,Groove);
    Obj{1,2}=fvc;
    Obj{2,2}=Color;
    Obj{3,2}=Alpha;

    %tissue properties assignment---------------------------------
    %element-wise multiplication for matrices
    rho(:,:,:,TypeIdx)=rho(:,:,:,TypeIdx)+mask.*Rho-mask.*master.Rho; 
    t1(:,:,:,TypeIdx)=t1(:,:,:,TypeIdx)+mask.*T1-mask.*master.T1;
    t2(:,:,:,TypeIdx)=t2(:,:,:,TypeIdx)+mask.*T2-mask.*master.T2;
    t2star(:,:,:,TypeIdx)=t2star(:,:,:,TypeIdx)+mask.*T2Star-mask.*master.T2Star;
    massden=massden+mask.*MassDen-mask.*master.MassDen;
    econ(:,:,:,1)=econ(:,:,:,1)+mask.*ECon(1)-mask.*master.ECon(1);
    econ(:,:,:,2)=econ(:,:,:,2)+mask.*ECon(2)-mask.*master.ECon(2);
    econ(:,:,:,3)=econ(:,:,:,3)+mask.*ECon(3)-mask.*master.ECon(3);
    %direct assignment for single values
    chemshift(TypeIdx)=ChemShift;
    typeflag(TypeIdx)=TypeFlag;
    lineshapeflag(TypeIdx)=LineShapeFlag;
    %--------------------------------------------------------------
    if ~strcmp(VObj.Model,'Normal')
        for i=1:VObj.TypeNum
            k(:,:,:,(TypeIdx-1)*VObj.TypeNum+i)=k(:,:,:,(TypeIdx-1)*VObj.TypeNum+i)+mask.*K(i);
        end
    end
    p=[];



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Cylinder Hemo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Alpha=[1];
    Color='b';
    Notes='Sphere 3';
    p.CenterX=[0.00];
    p.CenterY=[0.00];
    p.CenterZ=[0.00];
    p.FaceNum=[100];
    p.Length = [4.8]; % XZ
    p.Radius1 = r_endo + (r_epi - r_endo) * 0.1; % decimeter
    p.Radius2 = r_endo + (r_epi - r_endo) * (0.1 + transmural) ; % decimeter
    ChemShift=[0];
    ECon=[0.5 0.5 0.5];
    K=[10];
    LineShapeFlag=[0];
    MassDen=[1000];
    Rho=[1];
    T1=[1];
    T2=[3];
    T2Star=[0.005];
    TypeFlag=[0];
    TypeIdx=[1]; 

    [mask, fvc1, fvc2]=VObjAnnulusCylinder(p,x,y,z,theta,Groove,flag);
    Obj{1,2}=fvc1;
    Obj{2,2}=Color;
    Obj{3,2}=Alpha;

    rho(:,:,:,TypeIdx)=rho(:,:,:,TypeIdx)+mask.*Rho-mask.*master.Rho; 
    t1(:,:,:,TypeIdx)=t1(:,:,:,TypeIdx)+mask.*T1-mask.*master.T1;
    t2(:,:,:,TypeIdx)=t2(:,:,:,TypeIdx)+mask.*T2-mask.*master.T2;
    t2star(:,:,:,TypeIdx)=t2star(:,:,:,TypeIdx)+mask.*T2Star-mask.*master.T2Star;
    massden=massden+mask.*MassDen-mask.*master.MassDen;
    econ(:,:,:,1)=econ(:,:,:,1)+mask.*ECon(1)-mask.*master.ECon(1);
    econ(:,:,:,2)=econ(:,:,:,2)+mask.*ECon(2)-mask.*master.ECon(2);
    econ(:,:,:,3)=econ(:,:,:,3)+mask.*ECon(3)-mask.*master.ECon(3);
    chemshift(TypeIdx)=ChemShift;
    typeflag(TypeIdx)=TypeFlag;
    lineshapeflag(TypeIdx)=LineShapeFlag;
    if ~strcmp(VObj.Model,'Normal')
        for i=1:VObj.TypeNum
            k(:,:,:,(TypeIdx-1)*VObj.TypeNum+i)=k(:,:,:,(TypeIdx-1)*VObj.TypeNum+i)+mask.*K(i);
        end
    end
    p=[];

    VObj.t2star = t2star;
    Phantom_shape_cell{j} = VObj;

end
    %% Crop masks
    for j = 1:length(Phantom_shape_cell)
        VObj = Phantom_shape_cell{j};
        t2star_crop = VObj.t2star(1:VObj.XDim/2, 1:VObj.YDim/2, :);
        %figure();
        %imshow3D(t2star_crop);

        % t2star_crop_interp = interp2(t2star_crop(:,:,1), 5,'spline');
        % figure();
        % imagesc(t2star_crop_interp); colormap gray;
        VObj.t2star = t2star_crop;
        Phantom_shape_cell{j} = VObj;
    end
    save('/Users/jameszhang/Documents/MATLAB/T2star_Resolution_Project/Simulation_Results/3D_MagPurtabation/Cylinder_Phantom.mat', 'Phantom_shape_cell', '-v7.3');
%% Ellipsoid
%=====================================
%Create phantom for MRiLab simulation
%version 1.0
%This script is adapted from MRiLab
%Author: Hongzhang Chen
%=====================================

transmural_array = [0.025, 0.05, 0.10, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8];
a = 0.5; % minor half axis 50 mm
b = 0.5; % major half axis 60 mm
c = 2.4;
theta = pi/2;
Groove = 180;
blood_flow_perc = 0.6;

% r_endo = 0.3; % endo -> 30 mm
% r_epi =  0.5; % epi  -> 50 mm

Phantom_shape_cell = cell(1,length(transmural_array));

for j = 1:length(Phantom_shape_cell) %loop the code 4 times to produce 4 different phantoms
    %parameters of the virtual object are stored in the structure array, VObj
    transmural = transmural_array(j);
    tic
    VObj.Gyro=[267538030.3797];
    AttributeOpt={'Normal','MT','ME','GM'};
    VObj.Model=AttributeOpt{1};
    VObj.Name='VObj_Custom';
    VObj.Notes=' 5 sphere phantom';
    VObj.Type='Muscle';
    VObj.TypeNum=[1];
    VObj.XDim=[256];
    VObj.XDimRes=[0.004]; % cm
    VObj.YDim=[256];
    VObj.YDimRes=[0.004]; % cm
    VObj.ZDim=[24];
    VObj.ZDimRes=[0.1]; 
    [x,y,z]=meshgrid(((-(VObj.XDim-1)/2):((VObj.XDim-1)/2))*VObj.XDimRes,((-(VObj.YDim-1)/2):((VObj.YDim-1)/2))*VObj.YDimRes,((-(VObj.ZDim-1)/2):((VObj.ZDim-1)/2))*VObj.ZDimRes);
    
    %create matrices to store properties (e.g relaxation time t1 and t2)
    rho=zeros([size(x) VObj.TypeNum]);
    t1=zeros([size(x) VObj.TypeNum]);
    t2=zeros([size(x) VObj.TypeNum]);
    t2star=zeros([size(x) VObj.TypeNum]);
    chemshift=zeros(1,VObj.TypeNum);
    typeflag=zeros(1,VObj.TypeNum);
    lineshapeflag=zeros(1,VObj.TypeNum);
    k=zeros([size(x) VObj.TypeNum^2]);
    econ=zeros([size(x) 3]);  
    massden=zeros(size(x));
     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Ellipsoid %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %This part has the parameters of the largest sphere 
    %important parameters are set by using 'sturct master'
    Alpha=[1];
    Color='r';
    Notes='Prolate Ellipsoid (Muscle)';
    p.CenterX=[0];
    p.CenterY=[0];
    p.CenterZ=[-1.2];
    p.FaceNum=[40];
    p.Length = [1]; % XZ
    p.a = a;
    p.b = b;
    p.c = c;
    ChemShift=[0];
    K=[10];
    LineShapeFlag=[0];
    master.ECon=[0.29 0.29 0.29];
    master.MassDen=[1090];
    master.Rho=[0.7];
    master.T1=[1.1];
    master.T2=[0.035];
    master.T2Star=[0.01];
    TypeFlag=[0];
    TypeIdx=[1];

    %create 3D Cylindral virtual object
    [mask, fvc]=VObjEllipsoid(p,x,y,z,flag);
    Obj{1,1}=fvc;
    Obj{2,1}=Color;
    Obj{3,1}=Alpha;
    
    
    %tissue properties assignment---------------------------------
    %element-wise multiplication for matrices
    rho(:,:,:,TypeIdx)=rho(:,:,:,TypeIdx)+mask.*master.Rho; 
    t1(:,:,:,TypeIdx)=t1(:,:,:,TypeIdx)+mask.*master.T1;
    t2(:,:,:,TypeIdx)=t2(:,:,:,TypeIdx)+mask.*master.T2;
    t2star(:,:,:,TypeIdx)=t2star(:,:,:,TypeIdx)+mask.*master.T2Star;
    massden=massden+mask.*master.MassDen;
    econ(:,:,:,1)=econ(:,:,:,1)+mask.*master.ECon(1);
    econ(:,:,:,2)=econ(:,:,:,2)+mask.*master.ECon(2);
    econ(:,:,:,3)=econ(:,:,:,3)+mask.*master.ECon(3);
    %direct assignment for single values
    chemshift(TypeIdx)=ChemShift;
    typeflag(TypeIdx)=TypeFlag;
    lineshapeflag(TypeIdx)=LineShapeFlag;
    %-------------------------------------------------------------
    
    if ~strcmp(VObj.Model,'Normal')
        for i=1:VObj.TypeNum
            k(:,:,:,(TypeIdx-1)*VObj.TypeNum+i)=k(:,:,:,(TypeIdx-1)*VObj.TypeNum+i)+mask.*K(i);
        end
    end
    p=[];


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Ellipsoid 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Sphere 2 is radius-variable 
    Alpha=[1];
    Color='b';
    Notes='Sphere 2';
    p.CenterX=[0];
    p.CenterY=[0];
    p.CenterZ=[-1.2];
    p.FaceNum=[40];
    p.Length = [1]; % XZ
    p.a = a*blood_flow_perc;
    p.b = b*blood_flow_perc;
    p.c = c*blood_flow_perc;
    ChemShift=[0];
    ECon=[0.5 0.5 0.5];
    K=[10];
    LineShapeFlag=[0];
    MassDen=[1000];
    Rho=[1];
    T1=[1];
    T2=[0.1];
    T2Star=[0.0175];
    TypeFlag=[0];
    TypeIdx=[1];
    
    %create 3D sphere virtual object
    [mask, fvc]=VObjEllipsoid(p,x,y,z,flag);
    Obj{1,2}=fvc;
    Obj{2,2}=Color;
    Obj{3,2}=Alpha;

    %tissue properties assignment---------------------------------
    %element-wise multiplication for matrices
    rho(:,:,:,TypeIdx)=rho(:,:,:,TypeIdx)+mask.*Rho-mask.*master.Rho; 
    t1(:,:,:,TypeIdx)=t1(:,:,:,TypeIdx)+mask.*T1-mask.*master.T1;
    t2(:,:,:,TypeIdx)=t2(:,:,:,TypeIdx)+mask.*T2-mask.*master.T2;
    t2star(:,:,:,TypeIdx)=t2star(:,:,:,TypeIdx)+mask.*T2Star-mask.*master.T2Star;
    massden=massden+mask.*MassDen-mask.*master.MassDen;
    econ(:,:,:,1)=econ(:,:,:,1)+mask.*ECon(1)-mask.*master.ECon(1);
    econ(:,:,:,2)=econ(:,:,:,2)+mask.*ECon(2)-mask.*master.ECon(2);
    econ(:,:,:,3)=econ(:,:,:,3)+mask.*ECon(3)-mask.*master.ECon(3);
    %direct assignment for single values
    chemshift(TypeIdx)=ChemShift;
    typeflag(TypeIdx)=TypeFlag;
    lineshapeflag(TypeIdx)=LineShapeFlag;
    %--------------------------------------------------------------
    if ~strcmp(VObj.Model,'Normal')
        for i=1:VObj.TypeNum
            k(:,:,:,(TypeIdx-1)*VObj.TypeNum+i)=k(:,:,:,(TypeIdx-1)*VObj.TypeNum+i)+mask.*K(i);
        end
    end
    p=[];


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Ellipsoid Hemo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Alpha=[1];
    Color='b';
    Notes='Ellipsoid 3';
    p.CenterX=[0.00];
    p.CenterY=[0.00];
    p.CenterZ=[-1.2];
    p.FaceNum=[40];
    p.Length = [1]; % XZ
    p.a1 = a*blood_flow_perc + a*(1-blood_flow_perc)*0.1;
    p.b1 = b*blood_flow_perc + b*(1-blood_flow_perc)*0.1;
    p.c1 = c*blood_flow_perc + c*(1-blood_flow_perc)*0.1;
    p.a2 = a*blood_flow_perc + a*(1-blood_flow_perc)*(transmural+0.1);
    p.b2 = b*blood_flow_perc + b*(1-blood_flow_perc)*(transmural+0.1);
    p.c2 = c*blood_flow_perc + c*(1-blood_flow_perc)*(transmural+0.1);
    ChemShift=[0];
    ECon=[0.5 0.5 0.5];
    K=[10];
    LineShapeFlag=[0];
    MassDen=[1000];
    Rho=[1];
    T1=[1];
    T2=[0.1];
    T2Star=[0.005];
    TypeFlag=[0];
    TypeIdx=[1]; 

    [mask, fvc1, fvc2]=VObjAnnulusEllipsoid(p,x,y,z,theta,Groove,flag);
    Obj{1,2}=fvc1;
    Obj{2,2}=Color;
    Obj{3,2}=Alpha;

    rho(:,:,:,TypeIdx)=rho(:,:,:,TypeIdx)+mask.*Rho-mask.*master.Rho; 
    t1(:,:,:,TypeIdx)=t1(:,:,:,TypeIdx)+mask.*T1-mask.*master.T1;
    t2(:,:,:,TypeIdx)=t2(:,:,:,TypeIdx)+mask.*T2-mask.*master.T2;
    t2star(:,:,:,TypeIdx)=t2star(:,:,:,TypeIdx)+mask.*T2Star-mask.*master.T2Star;
    massden=massden+mask.*MassDen-mask.*master.MassDen;
    econ(:,:,:,1)=econ(:,:,:,1)+mask.*ECon(1)-mask.*master.ECon(1);
    econ(:,:,:,2)=econ(:,:,:,2)+mask.*ECon(2)-mask.*master.ECon(2);
    econ(:,:,:,3)=econ(:,:,:,3)+mask.*ECon(3)-mask.*master.ECon(3);
    chemshift(TypeIdx)=ChemShift;
    typeflag(TypeIdx)=TypeFlag;
    lineshapeflag(TypeIdx)=LineShapeFlag;
    if ~strcmp(VObj.Model,'Normal')
        for i=1:VObj.TypeNum
            k(:,:,:,(TypeIdx-1)*VObj.TypeNum+i)=k(:,:,:,(TypeIdx-1)*VObj.TypeNum+i)+mask.*K(i);
        end
    end
    p=[];

    toc

    VObj.t2star = t2star;
    Phantom_shape_cell{j} = VObj;

end


 %% Crop masks
    for j = 1:length(Phantom_shape_cell)
        VObj = Phantom_shape_cell{j};
        t2star_crop = VObj.t2star(1:VObj.XDim/2, 1:VObj.YDim/2, :);
        %figure();
        %imshow3D(t2star_crop);

        % t2star_crop_interp = interp2(t2star_crop(:,:,1), 5,'spline');
        % figure();
        % imagesc(t2star_crop_interp); colormap gray;
        VObj.t2star = t2star_crop;
        Phantom_shape_cell{j} = VObj;
    end

    save('/Users/jameszhang/Documents/MATLAB/T2star_Resolution_Project/Simulation_Results/3D_MagPurtabation/Ellipsoid_Phantom.mat', 'Phantom_shape_cell', '-v7.3');


    % %%
    % t2(t2<=0) = nan;
    % figure();
    % 
    % 
    % h = slice(t2, [], [], 1:size(t2,3)); 
    % set(h, 'EdgeColor','none', 'FaceColor','interp');
    % alpha 0.1;
