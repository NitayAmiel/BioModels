% Main script to compute SS NANOG vs B_minus using parameter struct
clear; clc;

for version = 1:3
    params = init_stem_cell_params(version); 
    
    % the original code was executed with different values
    B_radius = 3;
    num_repeats = 2;
    maximal_time_for_stochastic_simulation = 100;

    % Run simulation
    [NANOG_SS, threshold_B, SS_with_B] = run_stem_cell_sweep_det(params, version, 1, B_radius);
    % Call stochastic simulation and plot
    stochastic_plot_from_deterministic(params, threshold_B, B_radius, num_repeats, version, SS_with_B, maximal_time_for_stochastic_simulation);

end