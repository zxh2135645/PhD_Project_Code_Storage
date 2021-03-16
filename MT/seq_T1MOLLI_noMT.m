function [t_total, Mz_total_total, t_readout_array, Mxy_readout_array] = seq_T1MOLLI_noMT(TI_array, TD, npulse, T1, T2, alpha, TR, prep, num_rampup, Mz0, restore_pulse, trigger, trigger2)

sham_prep = struct;
sham_prep.flip = d2r(0);
sham_prep.t_delay = 0;

t_total = [];
Mz_total_total = [];
Mxy_readout_array = zeros(1,8);
t_readout_array = zeros(1,8);

for i = 1:length(TI_array)
%for i = 1:1
    if i == 1
        [t, Mz_total, mxys_array, mzs_array, Mxy_readout, t_readout] = ...
            seq_block_noMT(TD, npulse, T1, T2, alpha, TR, prep, num_rampup, Mz0, restore_pulse);
        
        t_readout_array(i) = t_readout;
        t_total = t;
        Mz_total_total = Mz_total;
    elseif i == 4
        TD = trigger2-t_total(end);
        prep.delay = TI_array(i);
        [t, Mz_total, mxys_array, mzs_array, Mxy_readout, t_readout] = ...
            seq_block_noMT(TD, npulse, T1, T2, alpha, TR, prep, num_rampup, Mz0, restore_pulse);
        idx = find(Mz_total == 1);
        Mz_total(idx) = Mz_total(idx+1);
        
        t_readout_array(i) = t_readout;
        t_total = [t_total t+t_total(end)];
        Mz_total_total = [Mz_total_total Mz_total];
    elseif i == 2 || i == 3
        TD = trigger + TI_array(i) - t_total(end);
        % fprintf('TD is: %4.2f\n', TD);
        [t, Mz_total, mxys_array, mzs_array, Mxy_readout, t_readout] = ...
            seq_block_noMT(TD, npulse, T1, T2, alpha, TR, sham_prep, num_rampup, Mz0, restore_pulse);
        idx = find(Mz_total == 1);
        Mz_total(idx) = Mz_total(idx+1);
        
        t_readout_array(i) = t_readout;
        t_total = [t_total t+t_total(end)];
        Mz_total_total = [Mz_total_total Mz_total];
    else
        TD = trigger2+TI_array(i)-t_total(end);
        [t, Mz_total, mxys_array, mzs_array, Mxy_readout, t_readout] = ...
            seq_block_noMT(TD, npulse, T1, T2, alpha, TR, sham_prep, num_rampup, Mz0, restore_pulse);
        idx = find(Mz_total == 1);
        Mz_total(idx) = Mz_total(idx+1);
        
        t_readout_array(i) = t_readout;
        t_total = [t_total t+t_total(end)];
        Mz_total_total = [Mz_total_total Mz_total];
    end
    % fprintf('%4.2f\n', t_total(end));
    Mz0 = Mz_total(end);
    Mxy_readout_array(i) = Mxy_readout;
end