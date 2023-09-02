function [t_total, Mz_total_total, t_readout_array, Mxy_readout_array] = seq_T1MOLLI_noMT_bloch(TI_array, TD, npulse, T1, T2, alpha, TR, prep, Mz0, trigger, trigger2, df, varargin)

% extract if RAMP_DOWN is input
if isempty(varargin)
    RAMP_DOWN = 0;
else
    for ii=1:length(varargin)
        if varargin{ii} == 1
            RAMP_DOWN = 1;
        elseif varargin{ii} == 0
            RAMP_DOWN = 0;
        end
    end
end

sham_prep = struct;
sham_prep.flip = d2r(0);
sham_prep.t_delay = 0;

t_total = [];
Mz_total_total = [];
Mxy_readout_array = zeros(1,8);
t_readout_array = zeros(1,8);

%for i = 1:length(TI_array)
 for i = 1:numel(TI_array)
    if i == 1
        [t, Mz_total, Mxy_readout, t_readout] = seq_block_noMT2(TD, npulse, T1, T2, alpha, TR, prep, Mz0, df, RAMP_DOWN);
        
        t_readout_array(i) = t_readout;
        t_total = t;
        Mz_total_total = Mz_total;
    elseif i == 4
        TD = trigger2-t_total(end);
        prep.t_delay = TI_array(i);
        [t, Mz_total, Mxy_readout, t_readout] = seq_block_noMT2(TD, npulse, T1, T2, alpha, TR, prep, Mz0, df, RAMP_DOWN);
        idx = find(Mz_total == 1);
        Mz_total(idx) = Mz_total(idx+1);
        
        t_readout_array(i) = t_readout;
        t_total = [t_total t+t_total(end)];
        Mz_total_total = [Mz_total_total Mz_total];
    elseif i == 2 || i == 3
        TD = trigger + TI_array(i) - t_total(end);
        % fprintf('TD is: %4.2f\n', TD);
        [t, Mz_total, Mxy_readout, t_readout] = seq_block_noMT2(TD, npulse, T1, T2, alpha, TR, sham_prep, Mz0, df, RAMP_DOWN);
        
        idx = find(Mz_total == 1);
        Mz_total(idx) = Mz_total(idx+1);
        
        t_readout_array(i) = t_readout;
        t_total = [t_total t+t_total(end)];
        Mz_total_total = [Mz_total_total Mz_total];
    else
        TD = trigger2+TI_array(i)-t_total(end);
        [t, Mz_total, Mxy_readout, t_readout] = seq_block_noMT2(TD, npulse, T1, T2, alpha, TR, sham_prep, Mz0, df, RAMP_DOWN);
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