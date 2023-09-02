function Func_plot_chord_analysis_EpiEndo2(MI_Chord_Analysis2, tp_dir2)
% Plot for Chord_Analysis2

% Bad structure
Mipix_mean_epi_hemo_p = MI_Chord_Analysis2(end).Mipix_mean_epi_hemo_p;
Mipix_mean2_epi_hemo_p = MI_Chord_Analysis2(end).Mipix_mean2_epi_hemo_p; 
Mipix_mean3_epi_hemo_p = MI_Chord_Analysis2(end).Mipix_mean3_epi_hemo_p;
Mipix_mean_epi_hemo_n = MI_Chord_Analysis2(end).Mipix_mean_epi_hemo_n;
Mipix_mean2_epi_hemo_n = MI_Chord_Analysis2(end).Mipix_mean2_epi_hemo_n; 
Mipix_mean3_epi_hemo_n = MI_Chord_Analysis2(end).Mipix_mean3_epi_hemo_n;
Mipix_mean_epi_hemo_n_ff_p = MI_Chord_Analysis2(end).Mipix_mean_epi_hemo_n_ff_p;
Mipix_mean2_epi_hemo_n_ff_p = MI_Chord_Analysis2(end).Mipix_mean2_epi_hemo_n_ff_p; 
Mipix_mean3_epi_hemo_n_ff_p = MI_Chord_Analysis2(end).Mipix_mean3_epi_hemo_n_ff_p;
Mipix_mean_epi_hemo_n_ff_n = MI_Chord_Analysis2(end).Mipix_mean_epi_hemo_n_ff_n;
Mipix_mean2_epi_hemo_n_ff_n = MI_Chord_Analysis2(end).Mipix_mean2_epi_hemo_n_ff_n; 
Mipix_mean3_epi_hemo_n_ff_n = MI_Chord_Analysis2(end).Mipix_mean3_epi_hemo_n_ff_n;

Mipix_mean_endo_hemo_p = MI_Chord_Analysis2(end).Mipix_mean_endo_hemo_p;
Mipix_mean2_endo_hemo_p = MI_Chord_Analysis2(end).Mipix_mean2_endo_hemo_p; 
Mipix_mean3_endo_hemo_p = MI_Chord_Analysis2(end).Mipix_mean3_endo_hemo_p;
Mipix_mean_endo_hemo_n = MI_Chord_Analysis2(end).Mipix_mean_endo_hemo_n;
Mipix_mean2_endo_hemo_n = MI_Chord_Analysis2(end).Mipix_mean2_endo_hemo_n; 
Mipix_mean3_endo_hemo_n = MI_Chord_Analysis2(end).Mipix_mean3_endo_hemo_n;
Mipix_mean_endo_hemo_n_ff_p = MI_Chord_Analysis2(end).Mipix_mean_endo_hemo_n_ff_p;
Mipix_mean2_endo_hemo_n_ff_p = MI_Chord_Analysis2(end).Mipix_mean2_endo_hemo_n_ff_p; 
Mipix_mean3_endo_hemo_n_ff_p = MI_Chord_Analysis2(end).Mipix_mean3_endo_hemo_n_ff_p;
Mipix_mean_endo_hemo_n_ff_n = MI_Chord_Analysis2(end).Mipix_mean_endo_hemo_n_ff_n;
Mipix_mean2_endo_hemo_n_ff_n = MI_Chord_Analysis2(end).Mipix_mean2_endo_hemo_n_ff_n; 
Mipix_mean3_endo_hemo_n_ff_n = MI_Chord_Analysis2(end).Mipix_mean3_endo_hemo_n_ff_n;

    Func_plot_chord_analysis_EpiEndo2_Module(Mipix_mean_epi_hemo_p, Mipix_mean2_epi_hemo_p, Mipix_mean3_epi_hemo_p,...
    Mipix_mean_endo_hemo_p, Mipix_mean2_endo_hemo_p, Mipix_mean3_endo_hemo_p, tp_dir2, 'hemo_p');

    Func_plot_chord_analysis_EpiEndo2_Module(Mipix_mean_epi_hemo_n, Mipix_mean2_epi_hemo_n, Mipix_mean3_epi_hemo_n,...
    Mipix_mean_endo_hemo_n, Mipix_mean2_endo_hemo_n, Mipix_mean3_endo_hemo_n, tp_dir2, 'hemo_n');

    Func_plot_chord_analysis_EpiEndo2_Module(Mipix_mean_epi_hemo_n_ff_p, Mipix_mean2_epi_hemo_n_ff_p, Mipix_mean3_epi_hemo_n_ff_p,...
    Mipix_mean_endo_hemo_n_ff_p, Mipix_mean2_endo_hemo_n_ff_p, Mipix_mean3_endo_hemo_n_ff_p, tp_dir2, 'hemo_n_ff_p');

    Func_plot_chord_analysis_EpiEndo2_Module(Mipix_mean_epi_hemo_n_ff_n, Mipix_mean2_epi_hemo_n_ff_n, Mipix_mean3_epi_hemo_n_ff_n,...
    Mipix_mean_endo_hemo_n_ff_n, Mipix_mean2_endo_hemo_n_ff_n, Mipix_mean3_endo_hemo_n_ff_n, tp_dir2, 'hemo_n_ff_n');
end