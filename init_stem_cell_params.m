function params = init_stem_cell_params(version)
    % Complex binding kinetics
    params.k1c = 0.05;    % nM^-1·sec^-1 (association rate of OCT4 and SOX2)
    params.k2c = 0.001;   % sec^-1 (dissociation rate of OS complex)
    if version == 3
        params.k1c = params.k1c*0.75;
        params.k2c = params.k2c*1.5;
    end
    params.k3c = 5;       % sec^-1 (degradation of OS complex)

    % Degradation rates
    params.gamma_O = 1;   % sec^-1 (OCT4 degradation)
    if version == 2
        params.gamma_O = params.gamma_O + 0.25*2.5; %[WWP2]*k_ub
    end
    params.gamma_S = 1;   % sec^-1 (SOX2 degradation)
    params.gamma_N = 1;   % sec^-1 (NANOG degradation)

    % Basal expression parameters
    params.eta_1 = 1e-4;  % unitless (basal OCT4 transcription)
    params.eta_2 = 1e-7;  % unitless (OCT4 transcriptional denominator)
    params.eta_3 = 1e-4;  % unitless (basal SOX2 transcription)
    params.eta_4 = 1e-7;  % unitless (SOX2 transcriptional denominator)
    params.eta_5 = 1e-4;  % unitless (basal NANOG transcription)
    params.eta_6 = 1e-7;  % unitless (NANOG transcriptional denominator)

    % Shea-Ackers numerator coefficients (OCT4)
    params.a1 = 1;        % nM (A⁺ activation of OCT4)
    params.a2 = 0.01;     % nM (OS activation of OCT4)
    params.a3 = 0.2;      % nM^2 (OS*N activation of OCT4)

    % Shea-Ackers denominator coefficients (OCT4)
    params.b1 = 0.0011;   % nM (A⁺ binding for OCT4 repression)
    params.b2 = 0.001;    % nM (OS binding)
    params.b3 = 0.0007;   % nM^2 (OS*N binding)

    % NANOG production regulation (SOX2)
    params.c1 = 1;        % nM (A⁺ activation of SOX2)
    params.c2 = 0.01;     % nM (OS activation)
    params.c3 = 0.2;      % nM^2 (OS*N activation)

    % NANOG repression terms (SOX2)
    params.d1 = 0.0011;   % nM (A⁺ binding for repression)
    params.d2 = 0.001;    % nM (OS binding)
    params.d3 = 0.0007;   % nM^2 (OS*N binding)

    % NANOG activation by OS and OS*N
    params.e1 = 0.005;    % nM (OS activation)
    params.e2 = 0.1;      % nM^2 (OS*N activation)

    % NANOG repression by B⁻
    params.f1 = 0.001;    % nM (OS repression)
    params.f2 = 9.95e-4;  % nM^2 (OS*N repression)
    params.f3 = 0.01;     % nM (B⁻ repression)

    % Set fixed A⁺ value
    params.A_plus = 100;  % nM (external activation signal)
end



%{
function params = init_stem_cell_params(version)
    % Complex binding kinetics
    params.k1c = 0.05;   
    params.k2c = 0.001;
    if version == 3
        params.k1c = params.k1c*0.75;
        params.k2c = params.k2c*1.5;
    end
    params.k3c = 5;      

    % Degradation rates
    params.gamma_O = 1;
    if version == 2
        params.gamma_O = params.gamma_O + 0.25*2.5;
    end
    params.gamma_S = 1;
    params.gamma_N = 1;

    % Basal expression parameters
    params.eta_1 = 1e-4;
    params.eta_2 = 1e-7;
    params.eta_3 = 1e-4;
    params.eta_4 = 1e-7;
    params.eta_5 = 1e-4;
    params.eta_6 = 1e-7;





    % Shea-Ackers numerator coefficients
    params.a1 = 1;
    params.a2 = 0.01;
    params.a3 = 0.2;

    % Shea-Ackers denominator coefficients
    params.b1 = 0.0011;
    params.b2 = 0.001;
    params.b3 = 0.0007;


    % NANOG production regulation
    params.c1 = 1;
    params.c2 = 0.01;
    params.c3 = 0.2;

    % NANOG repression terms
    params.d1 = 0.0011;
    params.d2 = 0.001;
    params.d3 = 0.0007;

    % NANOG activation by OS and OS*N
    params.e1 = 0.005;
    params.e2 = 0.1;


    % NANOG repression by B−
    params.f1 = 0.001;
    params.f2 = 9.95e-4;
    params.f3 = 0.01;

    % Set fixed A⁺ value
    params.A_plus = 100;
end
%}