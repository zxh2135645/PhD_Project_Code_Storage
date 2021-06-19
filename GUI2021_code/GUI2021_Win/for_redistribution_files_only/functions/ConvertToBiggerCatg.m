function [labelo] = ConvertToBiggerCatg(label)
    switch label
        case {'MAG', 'PSIR'}
            labelo = 'LGE';
        case {'T1Map'}
            labelo = 'T1';
        case {'T1Map_E'}
            labelo = 'T1_E';
    end
            
end