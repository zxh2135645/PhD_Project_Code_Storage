%% InVivo Run
%-----------------------------
% Cedars Sinai Randy Yang 2018
%HJrandy.Yang@gmail.com
%-----------------------------
clear all
%% Load Bz;
Bz_name=['Bz3m1_30sliceMEID.mat'];
load([pwd,filesep,'Bz',filesep,Bz_name]);
if strcmp(Bz_name,['Bz5m3_30sliceMEID.mat']);
Bz_30slice=Bz_30slice/2;
end
addpath(genpath(pwd))

%Unic_BzMap=UNIC_B0Map;

%% Bz into 42 channels
%% Set Recon Parameters
Mask_exist=0;
Scan_Type='Cardiac';
IS_Slice=0;
ShimSlice=[14 15 16];

%% Load B0map data & Calculate B0map


B0mapFromRaw_MEDI
%%

%ResolutionMatch
Bz=Bz_30slice;
%% Save the map with the Mask
UNIC_FreqMap_org=UNIC_FreqMap;
UNIC_B0Map_org=UNIC_B0Map;
eval( strcat('B0Freq_org =UNIC_FreqMap_org.Map;'));
if IS_Slice
    SliceMask=zeros(size(Mask));
    SliceMask(:,:,ShimSlice)=1;
    Mask=Mask.*SliceMask;
end
B0Freq_org=B0Freq_org.*Mask;


%% Pre shim Shimm stats

Stats_org=B0_Stats(B0Freq_org,UNIC_FreqMap_org.Parameters.Mask)

%% Prep B0, Bz and Mask
%Mask need to close holes
eval( strcat('B0 =UNIC_B0Map_org.Map;'));
Prep_B_forSolve

%% Solve DC
DC_limit=5;
DC=SolveDC(B0f,Bzf,DC_limit)
DCMat=reshape(round(100*-DC)/100,[7,6]);






%% Theoretical Shim field
clear Shimfield
for c=1:size(Bz,4)
Shimfield(:,:,:,c)=Mask.*Bz(:,:,:,c)*-DC(c);
end
Shimf=sum(Shimfield,4)*UNIC_FreqMap_org.Parameters.CentralFeq;
%Shimf(UNIC_FreqMap_org.Parameters.Mask==0)=nan;




%% Stats Comparison

slicer((Shimf+100)/200,(B0Freq_org+100)/200)
slicer((Shimf+B0Freq_org+100)/200,(B0Freq_org+100)/200)
