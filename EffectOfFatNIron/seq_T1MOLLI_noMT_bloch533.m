function [t_total, M_total_total, t_readout_array, Mxy_readout_array] = seq_T1MOLLI_noMT_bloch533(TI_array, TD, npulse, T1, T2, alpha, TR, prep, M0, trigger, trigger2, df, adiabatic, varargin)

% seq_block_noMT3 which adds adiabatic pulse 
% And and Mxy encoding over time

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

sham_adiabatic = adiabatic;
sham_adiabatic.pulseWidth = 0;

t_total = [];
M_total_total = [];
Mxy_readout_array = zeros(1,8);
t_readout_array = zeros(1,8);

%for i = 1:length(TI_array)
 for i = 1:numel(TI_array)
    if i == 1
        [t, M_total, Mxy_readout, t_readout] = seq_block_noMT3(TD, npulse, T1, T2, alpha, TR, prep, M0, df, RAMP_DOWN, adiabatic);
        
        t_readout_array(i) = t_readout;
        t_total = t;
        M_total_total = M_total;
        
        % Looking at the spin evolution in Mz direction
        % figure(); plot(t, M_total(3,:));
        % title(cat(2, num2str(i)));
    elseif i == 6
        TD = trigger2-t_total(end);
        prep.t_delay = TI_array(i);
        [t, M_total, Mxy_readout, t_readout] = seq_block_noMT3(TD, npulse, T1, T2, alpha, TR, prep, M0, df, RAMP_DOWN, adiabatic);
        
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
    elseif i == 2 || i == 3 || i == 4 || i == 5
        TD = trigger + TI_array(i) - t_total(end);
        % fprintf('TD is: %4.2f\n', TD);
        [t, M_total, Mxy_readout, t_readout] = seq_block_noMT3(TD, npulse, T1, T2, alpha, TR, sham_prep, M0, df, RAMP_DOWN, sham_adiabatic);
        
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
        [t, M_total, Mxy_readout, t_readout] = seq_block_noMT3(TD, npulse, T1, T2, alpha, TR, sham_prep, M0, df, RAMP_DOWN, sham_adiabatic);
        
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