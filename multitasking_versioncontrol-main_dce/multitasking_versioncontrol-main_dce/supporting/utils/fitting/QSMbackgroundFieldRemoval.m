function QSMparams = QSMbackgroundFieldRemoval(QSMparams)

methodBackgroundFieldRemoval = 'VSHARP';
radius = 5:-1:1;

extractVarFromStruct(QSMparams);

switch methodBackgroundFieldRemoval
    case 'PDF'  
        QSMparams.RDF = PDF(totalField,N_std,Mask,matrix_size,voxel_size,B0_dir);
        mask = imerode(QSMparams.Mask,strel('disk',1));
        QSMparams.RDF = QSMparams.RDF.*mask;
    case 'VSHARP'
        [QSMparams.RDF,mask] = QSMvSHARP(totalField,Mask,radius);
end