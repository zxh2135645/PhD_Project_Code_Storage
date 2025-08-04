function QSM = QSMestimateChi(QSMparams)

methodQSM = 'MEDI';
extractVarFromStruct(QSMparams);

switch methodQSM
    case 'MEDI'  
        iFreq = totalField;
        save RDF.mat RDF iFreq iFreq_raw iMag N_std Mask matrix_size...
             voxel_size delta_TE CF B0_dir Mask_CSF;
        
        % Morphology enabled dipole inversion with zero reference using CSF (MEDI+0)
        QSM = MEDI_L1('lambda',1000,'lambda_CSF',100,'merit','smv',5);
    case 'iSWIM'
        QSM2 = QSMiSWIM(QSMparams);
end