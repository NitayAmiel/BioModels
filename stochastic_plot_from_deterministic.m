function stochastic_plot_from_deterministic(params, threshold_B, B_radius, num_repeats, version, init_conditions, maximal_time)
    

    B_vals = init_conditions(:, 5);
    n_B = length(B_vals);
    
    % Run Gillespie simulation
    results = gillespie_stem_cell_box(B_radius, init_conditions, params, num_repeats, version, maximal_time);

    % Interpolate and average traces for NANOG
    time_grid = 0:1:maximal_time;

    NANOG_avg = zeros(n_B, length(time_grid));

    for i = 1:n_B
        NANOG_matrix = zeros(num_repeats, length(time_grid));
        for r = 1:num_repeats
            sim = results{i, r};
            time = sim.time;
            trace = sim.trace(:, 4);  % NANOG = 4th column
            NANOG_interp = interp1(time, trace, time_grid, 'previous', 'extrap');
            NANOG_matrix(r, :) = NANOG_interp;
        end
        NANOG_avg(i, :) = mean(NANOG_matrix, 1);
    end

    % Plot results
    figure;
    hold on;
    colors = lines(n_B);

    for i = 1:n_B
        plot(time_grid, NANOG_avg(i, :), 'LineWidth', 2, 'Color', colors(i, :));
    end

    xlabel('Time [sec]');
    ylabel('[NANOG] [nM] (average of repeats)');
    title(sprintf('Steady-state NANOG vs B⁻ (version %d – %s)', version, getModelLabel(version)));
    legend(arrayfun(@(b) sprintf('B⁻ = %d', b), B_vals, 'UniformOutput', false), ...
       'Location', 'eastoutside');

    grid on;
    %filename = sprintf('nanog_vs_B_version%d.png', version);
    %saveas(gcf, filename);
end


function results = gillespie_stem_cell_box(B_radius, init_conditions, params, num_repeats, version,max_time)
    % Performs Gillespie simulation for each B_minus value around the threshold
    % Inputs:
    %   B_radius - radius around threshold to simulate (e.g., 5 means [threshold-5, ..., threshold+5])
    %   init_conditions - Nx4 array of initial [O, S, OS, N] values for each B_minus value
    %   params - parameter struct
    %   num_repeats - number of simulations per B_minus


    % Extract B values and initial states
    %B_vals = init_conditions(:, 5);
    B_vals = init_conditions(:, 5);
    init_states = init_conditions(:, 1:4);
    n_B = length(B_vals);
    init_states_modified = getONInitialStateValues();


    % Results structure
    results = cell(n_B, num_repeats);

    % Loop over B values
    for i = 1:n_B
        B = B_vals(i);
        params.B_minus = B;

        for r = 1:num_repeats
            % Set initial state for this run
            %In the case of running with previous SS calues 
            %state = init_states(i, :);
            state = init_states_modified;
            t = 0;
            time_points = t;
            state_trace = state;

            while t < max_time
                [tau, new_state] = gillespie_step(state, params, version);
                if isempty(tau), break; end
                t = t + tau;
                state = new_state;
                time_points(end+1) = t; %#ok<AGROW>
                state_trace(end+1, :) = state; %#ok<AGROW>
            end

            % Store result
            results{i, r} = struct('B_minus', B, 'time', time_points, 'trace', state_trace);
        end
    end
end

function [tau, new_state] = gillespie_step(state, p, version)
    % Reactions:
    % [O, S, OS, N]

    O = state(1);
    S = state(2);
    OS = state(3);
    N = state(4);

    A = p.A_plus;
    B = p.B_minus;

    % Calculate propensities (Shea-Ackers-style)
    % Transcription propensities
    a_O = (p.eta_1 + p.a1*A + p.a2*OS + p.a3*OS*N) / (1 + p.eta_2 + p.b1*A + p.b2*OS + p.b3*OS*N);
    a_S = (p.eta_3 + p.c1*A + p.c2*OS + p.c3*OS*N) / (1 + p.eta_4 + p.d1*A + p.d2*OS + p.d3*OS*N);
    a_N = (p.eta_5 + p.e1*OS + p.e2*OS*N) / (1 + p.eta_6 + p.f1*OS + p.f2*OS*N + p.f3*B);

    % All reaction propensities
    a = [ a_O, a_S, a_N, ...                        % transcription
          p.gamma_O * O, p.gamma_S * S, p.gamma_N * N, ...  % degradation
          p.k1c * O * S, ...                        % OS formation
          p.k2c * OS , p.k3c*OS];                   % OS dissociation, OS degragations

    a0 = sum(a);
    if a0 <= 0
        tau = [];
        new_state = state;
        return;
    end

    % Time until next reaction
    tau = -log(rand) / a0;

    % Determine which reaction occurs
    r = find(cumsum(a) >= rand * a0, 1);

    % Update state based on reaction
    new_state = state;
    switch r
        case 1  % O transcription
            new_state(1) = O + 1;
        case 2  % S transcription
            new_state(2) = S + 1;
        case 3  % N transcription
            new_state(4) = N + 1;
        case 4  % O degradation
            if O > 0, new_state(1) = O - 1; end
        case 5  % S degradation
            if S > 0, new_state(2) = S - 1; end
        case 6  % N degradation
            if N > 0, new_state(4) = N - 1; end
        case 7  % OS formation
            if O > 0 && S > 0
                new_state(1) = O - 1;
                new_state(2) = S - 1;
                new_state(3) = OS + 1;
            end
        case 8  % OS dissociation
            if OS > 0
                new_state(1) = O + 1;
                new_state(2) = S + 1;
                new_state(3) = OS - 1;
            end
        case 9  % OS dissociation
            if OS > 0
                new_state(3) = OS - 1;
            end
    end
end

