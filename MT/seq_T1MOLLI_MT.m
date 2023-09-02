function [t_total, Mzmt_total_total, t_readout_array, Mxy_readout_array] = seq_T1MOLLI_MT(TI_array, TD, npulse,...
    alpha, TR, MT_para, MT_prep, num_rampup, M0, restore_pulse, trigger, trigger2)
sham_prep = struct;
sham_prep.flip = d2r(0);
sham_prep.t_delay = 0;
sham_prep.B1SqrdTau = 0;

t_total = [];
Mzmt_total_total = [];
Mxy_readout_array = zeros(1,8);
t_readout_array = zeros(1,8);

for i = 1:length(TI_array)
    if i == 1
        [t, Mz_total, mxys_array, mzs_array, Mzmt_bound_total, Mxy_readout, t_readout] = ...
            seq_block_MT(TD, npulse, alpha, TR, MT_para, MT_prep, M0, num_rampup, restore_pulse);
        
        t_readout_array(i) = t_readout;
        t_total = t;
        Mzmt_total_total = Mz_total;
    elseif i == 4
        TD = trigger2-t_total(end);
        MT_prep.delay = TI_array(i);
        [t, Mz_total, mxys_array, mzs_array, Mzmt_bound_total, Mxy_readout, t_readout] = ...
            seq_block_MT(TD, npulse, alpha, TR, MT_para, MT_prep, M0, num_rampup, restore_pulse);
        idx = find(Mz_total == 1);
        Mz_total(idx) = Mz_total(idx+1);
        
        t_readout_array(i) = t_readout+t_total(end);
        t_total = [t_total t+t_total(end)];
        Mzmt_total_total = [Mzmt_total_total Mz_total];
    elseif i == 2 || i == 3
        TD = trigger + TI_array(i) - t_total(end);
        % fprintf('TD is: %4.2f\n', TD);
        [t, Mz_total, mxys_array, mzs_array, Mzmt_bound_total, Mxy_readout, t_readout] = ...
            seq_block_MT(TD, npulse, alpha, TR, MT_para, sham_prep, M0, num_rampup, restore_pulse);
        idx = find(Mz_total == 1);
        Mz_total(idx) = Mz_total(idx+1);
        
        t_readout_array(i) = t_readout+t_total(end);
        t_total = [t_total t+t_total(end)];
        Mzmt_total_total = [Mzmt_total_total Mz_total];
    else
        TD = trigger2+TI_array(i)-t_total(end);
        % fprintf('TD is: %4.2f\n', TD);
        [t, Mz_total, mxys_array, mzs_array, Mzmt_bound_total, Mxy_readout, t_readout] = ...
            seq_block_MT(TD, npulse, alpha, TR, MT_para, sham_prep, M0, num_rampup, restore_pulse);
        idx = find(Mz_total == 1);
        Mz_total(idx) = Mz_total(idx+1);
        
        t_readout_array(i) = t_readout+t_total(end);
        t_total = [t_total t+t_total(end)];
        Mzmt_total_total = [Mzmt_total_total Mz_total];
    end
    % fprintf('%4.2f\n', t_total(end));
    M0 = [0 0 Mz_total(end) Mzmt_bound_total(end)]';
    Mxy_readout_array(i) = Mxy_readout;
end
end