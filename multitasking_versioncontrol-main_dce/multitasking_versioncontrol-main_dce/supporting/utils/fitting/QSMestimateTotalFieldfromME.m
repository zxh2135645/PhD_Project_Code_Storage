function QSMparams = QSMestimateTotalFieldfromME(reconME,QSMparams)
%
%   Input
%               reconME     : multi-echo image, Ny*Nx*Nz*Necho(*Ncoils)
%               QSMparams   : image parameters. 
%                             (TEarray,matrix_size,voxel_size)
%
%   Output
%               totalField  : phase at (TEarray(2)-TEarray(1)), in radian
%

if ~isfield(QSMparams,'methodUnwrap') < 3
    QSMparams.methodUnwrap = 'Laplacian';
end

QSMparams.delta_TE = QSMparams.TEarray(2)-QSMparams.TEarray(1);

switch QSMparams.methodUnwrap
    case 'Laplacian'    
        [QSMparams.iFreq_raw, QSMparams.N_std, QSMparams.relres, QSMparams.p0, QSMparams.iter] = Fit_ppm_complex_TE(reconME,QSMparams.TEarray);
%         [QSMparams.iFreq_raw, QSMparams.N_std, QSMparams.relres, QSMparams.p0, QSMparams.iter] = QSMfitTotalFieldRad(reconME,QSMparams.TEarray);
        QSMparams.totalField = unwrapLaplacian(QSMparams.iFreq_raw,QSMparams.matrix_size,QSMparams.voxel_size);
        QSMparams.pFreqw = QSMparams.totalField/QSMparams.delta_TE;
    otherwise
        [QSMparams.iFreq_raw, QSMparams.N_std, QSMparams.relres, QSMparams.p0, QSMparams.iter] = Fit_ppm_complex_TE(reconME,QSMparams.TEarray);
        QSMparams.totalField = unwrapLaplacian(QSMparams.iFreq_raw,QSMparams.matrix_size,QSMparams.voxel_size);
        QSMparams.pFreqw = QSMparams.totalField/QSMparams.delta_TE;
end
