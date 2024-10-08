function [t_total, M_total_total, t_readout_array, Mxy_readout_array] = seq_T1MOLLI_noMT_bloch3(TI_array, TD, npulse, T1, T2, alpha, TR, prep, M0, trigger, trigger2, df, adiabatic, RAMP_DOWN, ddt)

% seq_block_noMT3 which adds adiabatic pulse 
% And and Mxy encoding over time
% seq_block_noMT3 which adds ddt as temporal resolution

sham_prep = struct;
sham_prep.flip = d2r(0);
sham_prep.t_delay = 0;

sham_adiabatic = adiabatic;
sham_adiabatic.pulseWidth = 0;

t_total = [];
M_total_total = [];
Mxy_readout_array = zeros(1,8);
t_readout_array = zeros(1,8);

%for i = 1:length(TI_array)
 for i = 1:numel(TI_array)
    if i == 1
        [t, M_total, Mxy_readout, t_readout] = seq_block_noMT4(TD, npulse, T1, T2, alpha, TR, prep, M0, df, RAMP_DOWN, adiabatic, ddt);
        
        t_readout_array(i) = t_readout;
        t_total = t;
        M_total_total = M_total;
        
        % Looking at the spin evolution in Mz direction
        % figure(); plot(t, M_total(3,:));
        % title(cat(2, num2str(i)));
    elseif i == 4
        TD = trigger2-t_total(end);
        prep.t_delay = TI_array(i);
        [t, M_total, Mxy_readout, t_readout] = seq_block_noMT4(TD, npulse, T1, T2, alpha, TR, prep, M0, df, RAMP_DOWN, adiabatic, ddt);
        
        % Not sure if this is working (This helps removing spikes which Mz = 1)
        % Probably due to round in TD/dt
        idx = find(M_total(3,:) == 1);
        M_total(:,idx) = M_total(:,idx+1);
        
        t_readout_array(i) = t_readout;
        t_total = [t_total t+t_total(end)];
        M_total_total = [M_total_total M_total];
        
        % Looking at the spin evolution in Mz direction
        % figure(); plot(t, M_total(3,:));
        % title(cat(2, num2str(i)));
    elseif i == 2 || i == 3
        TD = trigger + TI_array(i) - t_total(end);
        % fprintf('TD is: %4.2f\n', TD);
        [t, M_total, Mxy_readout, t_readout] = seq_block_noMT4(TD, npulse, T1, T2, alpha, TR, sham_prep, M0, df, RAMP_DOWN, sham_adiabatic, ddt);
        
        % Not sure if this is working (This helps removing spikes which Mz = 1)
        % Probably due to round in TD/dt
        idx = find(M_total(3,:) == 1);
        M_total(:,idx) = M_total(:,idx+1);
        
        t_readout_array(i) = t_readout;
        t_total = [t_total t+t_total(end)];
        M_total_total = [M_total_total M_total];
        
        % Looking at the spin evolution in Mz direction
        % figure(); plot(t, M_total(3,:));
        % title(cat(2, num2str(i)));
    else
        TD = trigger2+TI_array(i)-t_total(end);
        [t, M_total, Mxy_readout, t_readout] = seq_block_noMT4(TD, npulse, T1, T2, alpha, TR, sham_prep, M0, df, RAMP_DOWN, sham_adiabatic, ddt);
        
        % Not sure if this is working (This helps removing spikes which Mz = 1)
        % Probably due to round in TD/dt
        idx = find(M_total(3,:) == 1);
        M_total(:,idx) = M_total(:,idx+1);
        
        t_readout_array(i) = t_readout;
        t_total = [t_total t+t_total(end)];
        M_total_total = [M_total_total M_total];
        
        % Looking at the spin evolution in Mz direction
        % figure(); plot(t, M_total(3,:));
        % title(cat(2, num2str(i)));
    end
    % fprintf('%4.2f\n', t_total(end));
    M0 = M_total(:,end);
    Mxy_readout_array(i) = Mxy_readout;
end