function [NANOG_SS, threshold_B, SS_with_B] = run_stem_cell_sweep_det(params, version, plotting_results, B_radius)
    % Sweep B_minus values from 0 to 100 and compute SS NANOG levels
    epsilon = 3e-1;

    % B_minus sweep setup
    B_range = 0:1:100;
    n_steps = length(B_range);
    % Initial condition: [OCT4, SOX2, OS, NANOG]
    y0 = getONInitialStateValues();
   

    % Storage for NANOG steady-state values
    NANOG_SS = zeros(n_steps, 1);
    SS_values = zeros(n_steps, 4);  % Store [O, S, OS, N]
    % Integration time span
    tspan = [0 500];

    
    ode_func = @(t, y, p) stem_cell_box_ode_struct(t, y, p);

    % Run sweep
    for i = 1:n_steps
        params.B_minus = B_range(i);
        [~, y] = ode45(@(t, y) ode_func(t, y, params), tspan, y0);

        y_ss = y(end, :);
        NANOG_SS(i) = y_ss(4);
        SS_values(i, :) = y_ss;
        % In a case of starting with the initial states of previous value
        % y0 = y_ss';
    end

    if plotting_results == 1
        % Plot results
        figure;
        plot(B_range, NANOG_SS, 'r-', 'LineWidth', 2);
        xlabel('B⁻ [nM]');
        ylabel('Steady-state [NANOG] [nM]');

        title(sprintf('Steady-state NANOG vs B⁻ (version %d – %s)', version, getModelLabel(version)));
        grid on;

        %filename = sprintf('detB_version%d.png', version);
        %saveas(gcf, filename);
    end

    % Determine threshold B⁻ where NANOG becomes 0 (or effectively zero)
    zero_idx = find(NANOG_SS < epsilon, 1, 'first');
    if isempty(zero_idx)
        threshold_B = NaN;
        SS_with_B = [];
    else
        threshold_B = B_range(zero_idx);
        fprintf('Threshold value in version %d: %d\n', version, threshold_B);
        SS_with_B = [];
        for b = (threshold_B - B_radius):(threshold_B + B_radius)
            if b - 1 >= 0
                idx = b;            % since B_range(i) = i-1
                ss_prev = SS_values(idx, :);  % SS at B⁻−1
                SS_with_B(end+1, :) = [ss_prev, b];  % [O S OS N B⁻]
            end
        end
    end
end