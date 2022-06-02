function [q,i]=IBSRread(am,ai,n)
% get address of the files downloaded from IBSR website
% www.bic.mni.mcgill.ca/brainweb/
% and return matrix representing image and ground truth segmentation for
% nth image from 20 normal images
% am is address of manual segmented image. For example 'C:\20Normals_T1_seg\'
% ai is address of input image. For example 'C:\20Normals_T1_brain\'
% nth image from 20 normal images 
% q is manual segmented image
% i is input image
b=[1	24	2	65	57	1	55
2	4	-6	65	56	7	62
5	8	-1	60	57	2	58
4	8	-2	61	55	3	57
6	10	0	63	61	1	61
7	8	0	60	56	1	56
8	4	-2	63	55	3	57
11	3	-1	63	58	2	59
12	3	0	63	59	1	59
13	3	0	63	57	1	57
15	3	0	60	55	1	55
16	3	0	60	51	1	51
17	3	-3	63	57	4	60
100	23	0	63	53	1	53
110	3	-2	64	55	3	57
111	2	-3	64	58	4	61
112	2	-1	63	55	2	56
191	3	-3	63	56	4	59
202	3	-3	63	58	4	61
205	3	-2	63	58	3	60
];
f=fopen([ai int2str(b(n,1)) '_' int2str(b(n,2)) '.buchar']);
a=fread(f);
i=reshape(a,256,256,b(n,4));   
f=fopen([am int2str(b(n,1)) '_' int2str(b(n,2)) '.buchar']);
a=fread(f);
q=reshape(a,256,256,b(n,5));
k=b(n,6)+1:b(n,7)-1;
q=q(:,:,k+b(n,3));
i=i(:,:,k); 
