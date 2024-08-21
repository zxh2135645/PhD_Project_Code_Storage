function [Mxy_readout_array_PhaseEnc] = Mxy_PhaseEnc(Mxy_readout_array)
    idx_array = imag(Mxy_readout_array) > 0;
    Mxy_readout_array_PhaseEnc = zeros(size(Mxy_readout_array));
    for i = 1:numel(idx_array)
        if idx_array(i) == 1
            Mxy_readout_array_PhaseEnc(i) =  -abs(Mxy_readout_array(i));
        else
            Mxy_readout_array_PhaseEnc(i) =  abs(Mxy_readout_array(i));
        end
    end
end