%% Water-Gas Shift Reactor Simulation
% CO + H2O ---> CO2 + H2

clc;
clear;

%% --- Parameters ---
R = 8.314;                    % J/mol.K
dp = 2.8e-3;                  % m
epsilon = 0.4;
rho_c = 3730;                % kg/m^3
rho_B = (1 - epsilon) * rho_c;
D = 0.09;                    % m
L = 2.2;                     % m
S_single = pi * (D/2)^2;     % m^2
num_tubes = 6000;
S_total = num_tubes * S_single;  % m^2
M = S_total * rho_B * 2 ;

%% Inlet flowrates (mol/s)
F_CO0   = 23.28;
F_H2O0  = 228.93;
F_CO2_0 = 94.19;
F_H2_0  = 364.149;
F_N2_0  = 134.354;
F_total0 = F_CO0 + F_H2O0 + F_CO2_0 + F_H2_0 + F_N2_0;

%% Initial conditions
T0 = 590;                      % K
P0 = 1.12 * 101325;            % Pa
XCO0 = 0;

%% Reactor domain
dz = 0.001;
z = 0:dz:2;
N = length(z);

%% Initialization
XCO = zeros(1,N); T = zeros(1,N); P = zeros(1,N);
XCO(1) = XCO0; T(1) = T0; P(1) = P0;

for i = 1:N-1
    % Current values
    X = XCO(i); T_i = T(i); P_i = P(i); P_atm = P_i / 101325;
    %% Cp functions (J/mol.K)
    Cp_CO  = @(T_i) 27.62 + 5.02e-3 * T_i;
    Cp_H2O = @(T_i) 30.13 + 10.46e-3 * T_i;
    Cp_CO2 = @(T_i) 32.22 + 22.18e-3 * T_i - 3.35e-6 * T_i.^2;
    Cp_H2  = @(T_i) 29.3 - 0.84e-3 * T_i + 2.09e-6 * T_i.^2;
    Cp_N2  = @(T_i) 27.62 + 4.19e-3 * T_i;

    %% Flow rates  (mol/s)
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

    %% Concentrations in mol/dm^3 = mol/L
    C_CO   = y_CO  * P_i / (R * T_i) / 1000;
    C_H2O  = y_H2O * P_i / (R * T_i) / 1000;
    C_H2   = y_H2  * P_i / (R * T_i) / 1000;
    C_CO2  = y_CO2 * P_i / (R * T_i) / 1000;

    %% Kinetics
    A = 2623447;
    Ea = 79759;
    Keq = @(T) exp(4577.8 / T_i - 4.33);
    P_atm = P_i / 101325;
    F_Pres = @(P_atm) P_atm.^(0.5 - (P_atm/250));

    rCO = @(Cco, Ch2o, Cco2, beta, T) ...
    -A * exp(-Ea./(R*T_i)) .* (Cco.^0.74) .* (Ch2o.^0.47) .* (Cco2.^-0.18) .* (1 - beta);
   
    %% Beta & rate (mol/kg_cat/h)
    Ke = Keq(T_i);
    beta = (C_CO2 * C_H2) / (C_CO * C_H2O * Ke);
    rate = rCO(C_CO, C_H2O, C_CO2, beta, T_i);
    
    %% dX/dz
    dXdz = @(X,T,P,z) (-rate * F_Pres(P/101325) * rho_B * S_total) / F_CO0;

    %% dT/dz
    dH_R = @(T) -4.12e4 + integral(@(T_i) Cp_CO2(T_i) + Cp_H2(T_i) - Cp_CO(T_i) - Cp_H2O(T_i),590,T-i); %(J/mol)
    sigmaFiCp = F_CO * Cp_CO(T_i) + F_H2O * Cp_H2O(T_i) + ...
                F_CO2 * Cp_CO2(T_i) + F_H2 * Cp_H2(T_i) + F_N2 * Cp_N2(T_i);
    dTdz = @(X,T,P,z) dH_R(T_i) * (rate) * F_Pres(P/101325) * rho_B * S_total / sigmaFiCp;
    

    %% dP/dz
    M_vec = [28.01, 18.02, 44.01, 2.016, 28.013];  % g/mol
    mu_CO  = (1.113e-7 * T_i^0.534) / (1 + (94.7 / T_i));
    mu_H2O = (6.1839e-7 * T_i^0.678) / (1 + (847.23/T_i) + (-73930/T_i^2));
    mu_CO2 = (2.148e-7 * T_i^0.46) / (1 + (290 / T_i));
    mu_H2  = (1.797e-7 * T_i^0.685) / (1 + ((-0.59)/T_i) + (140 / T_i^2));
    mu_N2  = (6.56e-7 * T_i^0.608) / (1 + (54.71 / T_i));
    mu_vec = [mu_CO, mu_H2O, mu_CO2, mu_H2, mu_N2];
    y_vec = [y_CO, y_H2O, y_CO2, y_H2, y_N2];

    mu = 0;
    for i_gas = 1:5
        denom = 0;
        for j_gas = 1:5
            if j_gas ~= i_gas
                phi_ij = (1 + sqrt(mu_vec(i_gas)/mu_vec(j_gas)) * (M_vec(j_gas)/M_vec(i_gas))^0.25)^2 / ...
                         ((2*(2^0.5)) * (1 + M_vec(i_gas)/M_vec(j_gas)^0.5));
                denom = denom + y_vec(j_gas) * phi_ij;
            end
        end
        mu = mu + mu_vec(i_gas) / (1 + denom / y_vec(i_gas));
    end

    M_mix = (y_CO*28.01 + y_H2O*18.02 + y_H2*2.016 + y_CO2*44.01 + y_N2*28.013) / 1000;
    rho_g = P_i * M_mix / (R * T_i);
    G = ((F_CO*28.01 + F_H2O*18.02 + F_H2*2.016 + F_CO2*44.01 + F_N2*28.013)/1000) / S_total;
    u_s = G / rho_g;
    Re = G * dp / mu;
    f = ((1 - epsilon)/(epsilon)^3)*(1.75 + 150 * (1 - epsilon)/Re);
    dPdz = @(X,T,P) -f * rho_g * u_s^2 / dp;

    %% -------------------- Runge-Kutta 4th Order --------------------

    % === RK4 for XCO ===
    k1_X = dXdz(X,T_i,P_i);
    k2_X = dXdz(X + 0.5*dz*k1_X, T_i, P_i);
    k3_X = dXdz(X + 0.5*dz*k2_X, T_i, P_i);
    k4_X = dXdz(X + dz*k3_X, T_i, P_i);
    XCO(i+1) = X + dz/6 * (k1_X + 2*k2_X + 2*k3_X + k4_X);

    % === RK4 for T ===
    k1_T = dTdz(X, T_i, P_i);
    k2_T = dTdz(X, T_i + 0.5*dz*k1_T, P_i);
    k3_T = dTdz(X, T_i + 0.5*dz*k2_T, P_i);
    k4_T = dTdz(X, T_i + dz*k3_T, P_i);
    T(i+1) = T_i + dz/6 * (k1_T + 2*k2_T + 2*k3_T + k4_T);

    % === RK4 for P ===
    k1_P = dPdz(X, T_i, P_i);
    k2_P = dPdz(X, T_i, P_i + 0.5*dz*k1_P);
    k3_P = dPdz(X, T_i, P_i + 0.5*dz*k2_P);
    k4_P = dPdz(X, T_i, P_i + dz*k3_P);
    P(i+1) = P_i + dz/6 * (k1_P + 2*k2_P + 2*k3_P + k4_P);

end

%% Plots
figure('Name','Water-Gas Shift Reactor Results','NumberTitle','off','Position',[100 100 900 600]);


subplot(2,2,1);
plot(z, XCO, 'r', 'LineWidth', 2); 
grid on;
xlabel('z (m)'); 
ylabel('X_{CO}');
title('CO Conversion Profile');


subplot(2,2,2);
plot(z, T, 'b', 'LineWidth', 2); 
grid on;
xlabel('z (m)'); 
ylabel('T (K)'); 
title('Temperature Profile');


subplot(2,2,3);
plot(z, P/101325, 'g', 'LineWidth', 2); 
grid on;
xlabel('z (m)');
ylabel('P (atm)');
title('Pressure Profile');

subplot(2,2,4);
axis off;
final_XCO = XCO(end);
final_T = T(end);
final_P = P(end) / 101325;
text(0.2, 0.8, sprintf('Catalyst Mass: %.2f t', M/1000), 'FontSize', 12, 'FontWeight', 'bold');
text(0.2, 0.6, sprintf('Final CO Conversion: %.2f%%', final_XCO*100), 'FontSize', 12);
text(0.2, 0.4, sprintf('Final Temperature: %.2f K', final_T), 'FontSize', 12);
text(0.2, 0.2, sprintf('Final Pressure: %.2f atm', final_P), 'FontSize', 12);

