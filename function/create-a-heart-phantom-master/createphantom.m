%=====================================
%Create phantom for MRiLab simulation
%version 1.0
%This script is adapted from MRiLab
%Author: Hongzhang Chen
%=====================================
for j = 1:1:1 %loop the code 4 times to produce 4 different phantoms
    %parameters of the virtual object are stored in the structure array, VObj
    VObj.Gyro=[267538030.3797];
    AttributeOpt={'Normal','MT','ME','GM'};
    VObj.Model=AttributeOpt{1};
    VObj.Name='VObj_Custom';
    VObj.Notes=' 5 sphere phantom';
    VObj.Type='Muscle';
    VObj.TypeNum=[1];
    VObj.XDim=[64];
    VObj.XDimRes=[0.002];
    VObj.YDim=[64];
    VObj.YDimRes=[0.002];
    VObj.ZDim=[64];
    VObj.ZDimRes=[0.002];
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
     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Sphere 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %This part has the parameters of the largest sphere 
    %important parameters are set by using 'sturct master'
    Alpha=[1];
    Color='r';
    Notes='Sphere 1(Muscle)';
    p.CenterX=[0];
    p.CenterY=[0];
    p.CenterZ=[0];
    p.FaceNum=[30];
    p.Radius=[0.06];
    ChemShift=[0];
    K=[10];
    LineShapeFlag=[0];
    master.ECon=[0.29 0.29 0.29];
    master.MassDen=[1090];
    master.Rho=[0.7];
    master.T1=[1.1];
    master.T2=[0.035];
    master.T2Star=[0.0175];
    TypeFlag=[0];
    TypeIdx=[1];

    %create 3D sphere virtual object
    [mask, fvc]=VObjSphere(p,x,y,z,flag);
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
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Sphere 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Sphere 2 is radius-variable 
    Alpha=[1];
    Color='b';
    Notes='Sphere 2';
    p.CenterX=[0.00];
    p.CenterY=[0.02];
    p.CenterZ=[0.02];
    p.FaceNum=[15];
    p.Radius=[0.01+ j/500 ];% increment of radius between 0.01 to 0.018 at a step of 0.002
    ChemShift=[0];
    ECon=[0.5 0.5 0.5];
    K=[10];
    LineShapeFlag=[0];
    MassDen=[1000];
    Rho=[1];
    T1=[1];
    T2=[0.1];
    T2Star=[0.01];
    TypeFlag=[0];
    TypeIdx=[1];
    
    %create 3D sphere virtual object
    [mask, fvc]=VObjSphere(p,x,y,z,flag);
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
    
    %--------------------------------------------------------------------
    %Sphere 3,4,5 have the same dimensions and remain unchanged
    %--------------------------------------------------------------------
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Sphere 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Alpha=[1];
    Color='b';
    Notes='Sphere 3';
    p.CenterX=[0.00];
    p.CenterY=[-0.02];
    p.CenterZ=[0.02];
    p.FaceNum=[15];
    p.Radius=[0.018];
    ChemShift=[0];
    ECon=[0.5 0.5 0.5];
    K=[10];
    LineShapeFlag=[0];
    MassDen=[1000];
    Rho=[1];
    T1=[1];
    T2=[0.1];
    T2Star=[0.01];
    TypeFlag=[0];
    TypeIdx=[1];
    [mask, fvc]=VObjSphere(p,x,y,z,flag);
    Obj{1,2}=fvc;
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
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Sphere 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Alpha=[1];
    Color='b';
    Notes='Sphere 4';
    p.CenterX=[0.00];
    p.CenterY=[-0.02];
    p.CenterZ=[-0.02];
    p.FaceNum=[15];
    p.Radius=[0.018];
    ChemShift=[0];
    ECon=[0.5 0.5 0.5];
    K=[10];
    LineShapeFlag=[0];
    MassDen=[1000];
    Rho=[1];
    T1=[1];
    T2=[0.1];
    T2Star=[0.01];
    TypeFlag=[0];
    TypeIdx=[1];
    [mask, fvc]=VObjSphere(p,x,y,z,flag);
    Obj{1,2}=fvc;
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

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Sphere 5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Alpha=[1];
    Color='b';
    Notes='Sphere 5';
    p.CenterX=[0.00];
    p.CenterY=[0.02];
    p.CenterZ=[-0.02];
    p.FaceNum=[15];
    p.Radius=[0.018];
    ChemShift=[0];
    ECon=[0.5 0.5 0.5];
    K=[10];
    LineShapeFlag=[0];
    MassDen=[1000];
    Rho=[1];
    T1=[1]; 
    T2=[0.1];
    T2Star=[0.01];
    TypeFlag=[0];
    TypeIdx=[1];
    [mask, fvc]=VObjSphere(p,x,y,z,flag);
    Obj{1,2}=fvc;
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
    %----------------------------------------------------------------------------
    
    
    %Assigning the properties to the struct, VObj
    VObj.Rho=rho;
    VObj.T1=t1;
    VObj.T2=t2;
    VObj.T2Star=t2star;
    VObj.ECon=econ;
    VObj.MassDen=massden;
    VObj.ChemShift=chemshift;
    if strcmp(VObj.Model,'GM')
        k=reshape(k,[size(x), VObj.TypeNum, VObj.TypeNum]);
        VObj.TypeFlag=typeflag;
        VObj.LineShapeFlag=lineshapeflag;
    end
    if ~strcmp(VObj.Model,'Normal')
        VObj.K=k;
    end

    %write the virtual object to .mat file
    phantom = [j];
    outfile = ['phantom' num2str(phantom) '.mat'];
    save(outfile, 'VObj');


end
