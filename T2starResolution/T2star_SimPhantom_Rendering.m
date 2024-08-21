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


%for j = 1:1:4 %loop the code 4 times to produce 4 different phantoms
    %parameters of the virtual object are stored in the structure array, VObj
    VObj.Gyro=[267538030.3797];
    AttributeOpt={'Normal','MT','ME','GM'};
    VObj.Model=AttributeOpt{1};
    VObj.Name='VObj_Custom';
    VObj.Notes=' 5 sphere phantom';
    VObj.Type='Muscle';
    VObj.TypeNum=[1];
    VObj.XDim=[128];
    VObj.XDimRes=[0.008];
    VObj.YDim=[128];
    VObj.YDimRes=[0.008];
    VObj.ZDim=[128];
    VObj.ZDimRes=[0.008];
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
    Notes='Prolate Cylinder (Muscle)';
    p.CenterX=[0];
    p.CenterY=[0];
    p.CenterZ=[0];
    p.FaceNum=[60];
    p.Length = [1]; % XZ
    p.Radius=[0.4];
    ChemShift=[0];
    K=[10];
    LineShapeFlag=[0];
    master.ECon=[0.29 0.29 0.29];
    master.MassDen=[1090];
    master.Rho=[0.7];
    master.T1=[1.1];
    master.T2=[1];
    master.T2Star=[0.0175];
    TypeFlag=[0];
    TypeIdx=[1];

    %create 3D Cylindral virtual object
    [mask, fvc]=VObjCylinder(p,x,y,z,flag);
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
    Notes='Cylinder 2';
    p.CenterX=[0];
    p.CenterY=[0];
    p.CenterZ=[0];
    p.FaceNum=[60];
    p.Length = [1]; % XZ
    p.Radius=[0.3];% increment of radius between 0.01 to 0.018 at a step of 0.002
    ChemShift=[0];
    ECon=[0.5 0.5 0.5];
    K=[10];
    LineShapeFlag=[0];
    MassDen=[1000];
    Rho=[1];
    T1=[1];
    T2=[2];
    T2Star=[0.01];
    TypeFlag=[0];
    TypeIdx=[1];
    
    %create 3D sphere virtual object
    [mask, fvc]=VObjCylinder(p,x,y,z,flag);
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

    rho_hemo=zeros([size(x) VObj.TypeNum]);
    t1_hemo=zeros([size(x) VObj.TypeNum]);
    t2_hemo=zeros([size(x) VObj.TypeNum]);
    t2star_hemo=zeros([size(x) VObj.TypeNum]);
    chemshift_hemo=zeros(1,VObj.TypeNum);
    typeflag_hemo=zeros(1,VObj.TypeNum);
    lineshapeflag_hemo=zeros(1,VObj.TypeNum);
    k_hemo=zeros([size(x) VObj.TypeNum^2]);
    econ_hemo=zeros([size(x) 3]);
    massden_hemo=zeros(size(x));

    Alpha=[1];
    Color='b';
    Notes='Sphere 3';
    p.CenterX=[0.00];
    p.CenterY=[0.00];
    p.CenterZ=[0.00];
    p.FaceNum=[60];
    p.Length = [1]; % XZ
    p.Radius1=[0.32];
    p.Radius2=[0.36];
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

    theta = pi/2;
    Groove = 180;
    [mask, fvc1, fvc2]=VObjAnnulusCylinder(p,x,y,z,theta,Groove,flag);
    Obj{1,2}=fvc1;
    Obj{2,2}=Color;
    Obj{3,2}=Alpha;

    rho_hemo(:,:,:,TypeIdx)=rho_hemo(:,:,:,TypeIdx)+mask.*Rho-mask.*master.Rho; 
    t1_hemo(:,:,:,TypeIdx)=t1_hemo(:,:,:,TypeIdx)+mask.*T1-mask.*master.T1;
    t2_hemo(:,:,:,TypeIdx)=t2_hemo(:,:,:,TypeIdx)+mask.*T2-mask.*master.T2;
    t2star_hemo(:,:,:,TypeIdx)=t2star_hemo(:,:,:,TypeIdx)+mask.*T2Star-mask.*master.T2Star;
    massden_hemo=massden_hemo+mask.*MassDen-mask.*master.MassDen;
    econ_hemo(:,:,:,1)=econ_hemo(:,:,:,1)+mask.*ECon(1)-mask.*master.ECon(1);
    econ_hemo(:,:,:,2)=econ_hemo(:,:,:,2)+mask.*ECon(2)-mask.*master.ECon(2);
    econ_hemo(:,:,:,3)=econ_hemo(:,:,:,3)+mask.*ECon(3)-mask.*master.ECon(3);
    chemshift_hemo(TypeIdx)=ChemShift;
    typeflag_hemo(TypeIdx)=TypeFlag;
    lineshapeflag_hemo(TypeIdx)=LineShapeFlag;
    if ~strcmp(VObj.Model,'Normal')
        for i=1:VObj.TypeNum
            k_hemo(:,:,:,(TypeIdx-1)*VObj.TypeNum+i)=k(:,:,:,(TypeIdx-1)*VObj.TypeNum+i)+mask.*K(i);
        end
    end
    p=[];

    %%
    % clear data
    % figure();
    % tiledlayout(1,2)
    % nexttile
    % [xq,yq,zq]=meshgrid(((-(VObj.XDim-1)/2):0.5:((VObj.XDim-1)/2))*VObj.XDimRes,((-(VObj.YDim-1)/2):0.5:((VObj.YDim-1)/2))*VObj.YDimRes,((-(VObj.ZDim-1)/2):0.5:((VObj.ZDim-1)/2))*VObj.ZDimRes);
    % mask_interp = interp3(x,y,z,mask,xq,yq,zq,'spline');
    % fv = isosurface(xq,yq,zq,mask_interp,0.5);
    % 
    % p1 = patch(fv,'FaceColor','red','EdgeColor','none');
    % view(3)
    % daspect([1,1,1])
    % axis tight
    % camlight
    % camlight(-80,-10)
    % lighting gouraud
    % title('Triangle Normals')
    % 
    % nexttile
    % fv = isosurface(xq,yq,zq,mask_interp,0.5);
    % p2 = patch(fv,'FaceColor','red','EdgeColor','none');
    % isonormals(xq,yq,zq,mask_interp,p2)
    % view(3)
    % daspect([1 1 1])
    % axis tight
    % camlight
    % camlight(-80,-10)
    % lighting gouraud
    % title('Data Normals')
%% 
t2(t2 == 0) = nan;
t2_hemo(t2_hemo == 0) = nan;

figure();
ax1 = axes;
x_beforemeshgrid = ((-(VObj.XDim-1)/2):((VObj.XDim-1)/2))*VObj.XDimRes;
y_beforemeshgrid = ((-(VObj.YDim-1)/2):((VObj.YDim-1)/2))*VObj.YDimRes;
z_beforemeshgrid = ((-(VObj.ZDim-1)/2):((VObj.ZDim-1)/2))*VObj.ZDimRes;

h = slice(x_beforemeshgrid, y_beforemeshgrid, z_beforemeshgrid, t2, x_beforemeshgrid, y_beforemeshgrid, z_beforemeshgrid); caxis([0 2.5]);
set(h, 'EdgeColor','none', 'FaceColor','interp');
alpha 0.05;
%view(3)
view(47.5441, 39.0131);
daspect([1,1,1]);
%axis off;


ax2 = axes;
% h2 = slice(t2_hemo, [], [], 1:size(t2_hemo,3)); caxis([0 2]);
% set(h2, 'EdgeColor','none', 'FaceColor','interp');

% TODO They are not in the same coordinate

fv = isosurface(xq,yq,zq,mask_interp,0.5);
p1 = patch(fv,'FaceColor','red','EdgeColor','none');
%view(3)
view(47.5441, 39.0131);
daspect([1,1,1])
axis tight
camlight
camlight(-80,-10)
lighting gouraud
alpha 0.3;

linkaxes([ax1, ax2])
ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];
colormap(ax1, gray);
%colormap(ax2, hsv);
%% 3D surface rendering
% figure();
% tiledlayout(1,2)
% nexttile
% data = interp3(mask,3,'spline');
% 
% fv = isosurface(data,.5);
% p1 = patch(fv,'FaceColor','red','EdgeColor','none');
% view(3)
% daspect([1,1,1])
% axis tight
% camlight
% camlight(-80,-10)
% lighting gouraud
% title('Triangle Normals')
% %%
% nexttile
% fv = isosurface(data,.5);
% p2 = patch(fv,'FaceColor','red','EdgeColor','none');
% isonormals(data,p2)
% view(3) 
% daspect([1 1 1])
% axis tight
% camlight 
% camlight(-80,-10) 
% lighting gouraud
% title('Data Normals')


