function [trigger, trigger2, half_readout] = MOLLI533_SchemeFunc(TI_array, b1, TR, PhaseEnc, num_rampup)

figure();
HR = 60*1000/(((TI_array(8) - TI_array(7)) + (TI_array(7) - TI_array(6)) + (TI_array(6) - TI_array(4)) + ...
    (TI_array(4) - TI_array(2)) + (TI_array(5) - TI_array(3)) + (TI_array(3) - TI_array(1))) / 6);
window = round(60*1000 / HR);
acq_win = TR*(PhaseEnc+num_rampup);
trigger = window-TI_array(1)-acq_win;
line([0,12],[0,0],'Color', 'black')
hold on;
plot([trigger trigger]/1000, [-b1 b1], 'LineWidth', 1.5)

% 3 bSSFP readout
rectangle('Position', [(window-acq_win)/1000 -b1 acq_win/1000 b1*2], 'FaceColor',[0.2 0.2 0.2])
rectangle('Position', [(trigger+TI_array(3))/1000 -b1 acq_win/1000 b1*2], 'FaceColor',[0.2 0.2 0.2])
rectangle('Position', [(trigger+TI_array(5))/1000 -b1 acq_win/1000 b1*2], 'FaceColor',[0.2 0.2 0.2])
rectangle('Position', [(trigger+TI_array(7))/1000 -b1 acq_win/1000 b1*2], 'FaceColor',[0.2 0.2 0.2])
rectangle('Position', [(trigger+TI_array(8))/1000 -b1 acq_win/1000 b1*2], 'FaceColor',[0.2 0.2 0.2])

% 3 recovery time
% 5 bSSFP readout
trigger2 = 9 * window - TI_array(2) - acq_win;
plot([trigger2 trigger2]/1000, [-b1 b1], 'LineWidth', 1.5)
rectangle('Position', [(9*window-acq_win)/1000 -b1 acq_win/1000 b1*2], 'FaceColor',[0.2 0.2 0.2])
rectangle('Position', [(trigger2+TI_array(4))/1000 -b1 acq_win/1000 b1*2], 'FaceColor',[0.2 0.2 0.2])
rectangle('Position', [(trigger2+TI_array(6))/1000 -b1 acq_win/1000 b1*2], 'FaceColor',[0.2 0.2 0.2])
xlabel('Time (s)')
ylabel('B_1 (Hz)')

half_readout = zeros(8, 1);
center_k = TR * (num_rampup + PhaseEnc/2) / 1000;
half_readout(1) = (window-acq_win)/1000 + center_k;
half_readout(2) = (trigger+TI_array(3))/1000 + center_k;
half_readout(3) = (trigger+TI_array(5))/1000 + center_k;
half_readout(4) = (trigger+TI_array(7))/1000 + center_k;
half_readout(5) = (trigger+TI_array(8))/1000 + center_k;
half_readout(6) = (9*window-acq_win)/1000 + center_k;
half_readout(7) = (trigger2+TI_array(4))/1000 + center_k;
half_readout(8) = (trigger2+TI_array(6))/1000 + center_k;
end