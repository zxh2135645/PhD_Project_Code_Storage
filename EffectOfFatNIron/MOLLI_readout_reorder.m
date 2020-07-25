function [Mxy_readout_array] = MOLLI_readout_reorder(Mxy_readout)
Mxy_readout_array = Mxy_readout;
Mxy_readout_array(1) = Mxy_readout(1);
Mxy_readout_array(2) = Mxy_readout(4);
Mxy_readout_array(3) = Mxy_readout(2);
Mxy_readout_array(4) = Mxy_readout(5);
Mxy_readout_array(5) = Mxy_readout(3);
end

