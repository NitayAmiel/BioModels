function dydt = stem_cell_box_ode_struct(t, y, params)
    % Unpack state variables
    O  = y(1);  % OCT4
    S  = y(2);  % SOX2
    OS = y(3);  % OCT4-SOX2 complex
    N  = y(4);  % NANOG

    % Unpack parameters from struct
    A = params.A_plus;
    B = params.B_minus;

    % ODEs from Equation 1 (Shea-Ackers style logic)
    dO = (params.eta_1 + params.a1*A + params.a2*OS + params.a3*OS*N) / ...
         (1 + params.eta_2 + params.b1*A + params.b2*OS + params.b3*OS*N) ...
         - params.gamma_O*O - params.k1c*O*S + params.k2c*OS;

    dS = (params.eta_3 + params.c1*A + params.c2*OS + params.c3*OS*N) / ...
         (1 + params.eta_4 + params.d1*A + params.d2*OS + params.d3*OS*N) ...
         - params.gamma_S * S - params.k1c * O * S + params.k2c * OS;

    dOS = params.k1c * O * S - params.k2c * OS - params.k3c * OS;

    dN = (params.eta_5 + params.e1 * OS + params.e2 * OS * N) / ...
         (1 + params.eta_6 + params.f1 * OS + params.f2 * OS * N + params.f3 * B) ...
         - params.gamma_N * N;

    % Return derivative vector
    dydt = [dO; dS; dOS; dN];
end
