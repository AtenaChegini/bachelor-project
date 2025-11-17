%% Water-Gas Shift Reactor Simulation for three F_H2O0 values
% CO + H2O ---> CO2 + H2

clc;
clear;

%% --- Parameters ---
R = 8.314;                    % J/mol.K
dp = 2.8e-3;                  % m
epsilon = 0.4;
rho_c = 3730;                 % kg/m^3
rho_B = (1 - epsilon) * rho_c;
D = 0.09;                    % m
L = 2.2;                     % m
S_single = pi * (D/2)^2;     % m^2
num_tubes = 6000;
S_total = num_tubes * S_single;  % m^2
M = S_total * rho_B * 2 ;

%% Inlet flowrates (mol/s) except H2O which will vary
F_CO0   = 23.28;
F_CO2_0 = 94.19;
F_H2_0  = 364.149;
F_N2_0  = 134.354;

F_H2O0_values = [228.93, 116.4, 69.84]; % mol/s
T0 = 590;                      % K fixed temperature
P0 = 1.12 * 101325;            % Pa
XCO0 = 0;

%% Reactor domain
dz = 0.001;
z = 0:dz:2.2;
N = length(z);

figure('Name','CO Conversion vs Reactor Length for Different F_{H2O0}','NumberTitle','off');
hold on;
colors = {'k', 'r', 'b'};          % black, red, blue
lineTypes = {'-', '--', ':'};      % solid, dashed, dotted

for idx = 1:length(F_H2O0_values)
    F_H2O0 = F_H2O0_values(idx);
    % Initialization
    XCO = zeros(1,N);
    XCO(1) = XCO0;

    for i = 1:N-1
        X = XCO(i);
        T_i = T0; % Fixed temperature, no T integration here
        P_i = P0; % Fixed pressure, no P integration here
        P_atm = P_i / 101325;

        % Cp functions (unused since T is fixed, but keep for rate calc)
        Cp_CO  = @(T_i) 27.62 + 5.02e-3 * T_i;
        Cp_H2O = @(T_i) 30.13 + 10.46e-3 * T_i;
        Cp_CO2 = @(T_i) 32.22 + 22.18e-3 * T_i - 3.35e-6 * T_i.^2;
        Cp_H2  = @(T_i) 29.3 - 0.84e-3 * T_i + 2.09e-6 * T_i.^2;
        Cp_N2  = @(T_i) 27.62 + 4.19e-3 * T_i;

        %% Flow rates (mol/s)
        F_CO   = F_CO0 * (1 - X);
        F_H2O  = F_H2O0 - (F_CO0 * X);
        F_H2   = F_H2_0 + F_CO0 * X;
        F_CO2  = F_CO2_0 + F_CO0 * X;
        F_N2   = F_N2_0;
        F_total = F_CO + F_H2O + F_H2 + F_CO2 + F_N2;

        %% Mole fractions
        y_CO   = F_CO / F_total;
        y_H2O  = F_H2O / F_total;
        y_H2   = F_H2 / F_total;
        y_CO2  = F_CO2 / F_total;
        y_N2   = F_N2 / F_total;

        %% Concentrations mol/dm^3 = mol/L
        C_CO   = y_CO  * P_i / (R * T_i) / 1000;
        C_H2O  = y_H2O * P_i / (R * T_i) / 1000;
        C_H2   = y_H2  * P_i / (R * T_i) / 1000;
        C_CO2  = y_CO2 * P_i / (R * T_i) / 1000;

        %% Kinetics
        A = 2623447;
        Ea = 79759;
        Keq = exp(4577.8 / T_i - 4.33);
        beta = (C_CO2 * C_H2) / (C_CO * C_H2O * Keq);

        rate = -A * exp(-Ea/(R*T_i)) * (C_CO^0.74) * (C_H2O^0.47) * (C_CO2^-0.18) * (1 - beta);

        %% dX/dz
        dXdz = (-rate * rho_B * S_total) / F_CO0;

        % Explicit Euler step (simpler here)
        XCO(i+1) = X + dz * dXdz;

        % Prevent unphysical values
        if XCO(i+1) > 1
            XCO(i+1) = 1;
        elseif XCO(i+1) < 0
            XCO(i+1) = 0;
        end
    end

        plot(z, XCO, ...
        'Color', colors{idx}, ...
        'LineStyle', lineTypes{idx}, ...
        'LineWidth', 2, ...
        'MarkerSize', 6);
end

grid on;
xlabel('Reactor length z (m)');
ylabel('CO Conversion X_{CO}');
title('CO Conversion along Reactor Length at T=590 K for different F_{H2O0}');
legend({'H_{2}O/CO = 10', 'H_{2}O/CO = 5', 'H_{2}O/CO = 3'}, 'Location', 'best');
hold off;