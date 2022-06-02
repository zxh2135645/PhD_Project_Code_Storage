function [corres_glob] = CVI_NamingConvert(cvi_glob, label, sequence_label)
% label is MAG, PSIR, T1Map originally
% sequence_label = {'MAG', 'PSIR', 'T1Map'}
[check_flg, pos] = find(strcmp(label, sequence_label));

if check_flg
    contain_status = contains(cvi_glob, label);
    if pos <= 2
        sig = sequence_label{3-pos};
        if ~any(contain_status)
            sig_contain_status = contains(cvi_glob, sig);
            corres_glob = cvi_glob(sig_contain_status);
        else
            corres_glob = cvi_glob(contain_status);
        end
    
    elseif pos == 3
       t1_label_new = [label(1:3) upper(label(4:5))];
       contain_status = contains(cvi_glob, t1_label_new);
       corres_glob = cvi_glob(contain_status);
    end

end

