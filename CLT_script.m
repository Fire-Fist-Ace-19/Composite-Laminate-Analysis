% =========================================================================
% LAMINATE ANALYSIS using CLASSICAL LAMINATION THEORY(CLT)
% Objectives:- 
% 1.ABD matrix generation
% 2.To calculate Effective Young's Moduli in EXtuension and Bending
% 3.To implement successive ply failure analysis with stiffness degradation in composite laminates using different failure criteria (Maximum Stress, Maximum Strain, and Tsai-Hill).
% 4.To determine the First Ply Failure (FPF) under both complete and partial degradation
% 5.To generate carpet plot for Uniaxial tensile Strength of a [0/+ or - 45/90]s Laminate with varying proportion of 0, ±45 and 90 laminae.
% =========================================================================
clear; clc; close all;

%% 1. Input Section
% Laminate stacking sequence 
theta = [0, 45, -45, 90, 90, -45, 45, 0];   % order is top -> bottom
NN    = length(theta);
t_ply = 0.000125;                           % assuming all plies have same thickness in m
h     = NN * t_ply;                         % Total laminate thickness (m)

% Lamina engineering constants (same for all plies here)
E1   = 38.6e9;    % Fibre-direction modulus (in Pa)
E2   = 8.27e9;    % Transverse modulus (in Pa)
v12  = 0.28;      % Major Poisson ratio
G12  = 4.14e9;    % In-plane shear modulus (in Pa)
v21  = v12 * E2 / E1;   % Minor Poisson ratio (by using reciprocity relation)

% Strength parameters (in Pa, all stored as positive magnitudes) 
Xtu = 1062e6;  Xcu = 610e6;   % Fibre tensile / compressive
Ytu =   31e6;  Ycu = 118e6;   % Transverse tensile / compressive
S  =   72e6;                % In-plane shear strength

% Ultimate strain values 
eps_Xtu  = Xtu / E1;   eps_Xcu  = Xcu / E1;
eps_Ytu  = Ytu / E2;   eps_Ycu  = Ycu / E2;
gamma_S = S  / G12;

% Hygrothermal properties 
% if laminate is free → these won't create stresses, only deformation
alpha1  =  8.6e-6;   % CTE along fibre (1/°C)
alpha2  = 22.1e-6;   % CTE transverse  (1/°C)
beta1   =  0.0;      % CME along fibre
beta2   =  0.0;      % CME transverse
Delta_T =  50;       % Temperature change (°C)
Delta_C =  0.0;      % Moisture change (kg/kg or %)

% Loading conditions 
Nx  = 100e3;  Ny  = 0;  Nxy = 0;   % N/m
Mx  = 0;      My  = 0;  Mxy = 0;   % N·m/m


%% 2. Ply z-coordinate Assembly
% building thickness coordinate from bottom (-h/2) to top (+h/2)
% needed for strain variation through thickness (z-dependence)
z = zeros(1, NN+1);
z(1) = -h/2;
for k = 1:NN
    z(k+1) = z(k) + t_ply;   
end

%% 3. Generation of ABBD matrix using the self made ABD function
[A, B, D, Q_cell, Qbar_cell] = fn_ABD(E1, E2, v12, v21, G12, theta, z, NN);
ABD = [A, B; B, D];
% A → in-plane
% B → bending-extension coupling
% D → bending stiffness
fprintf('\n[A], [B], [D] Matrices\n');
fprintf('A Matrix (N/m):\n');      disp(A);
fprintf('B Matrix (N):\n');        disp(B);
fprintf('D Matrix (N·m):\n');      disp(D);
fprintf('ABD Matrix (N·m):\n');      disp(ABD);

%% 4. Hygrothermal Resultant Forces & Moments

NT = zeros(3,1);  MT = zeros(3,1);
NC = zeros(3,1);  MC = zeros(3,1);

for k = 1:NN
    c = cosd(theta(k));  s = sind(theta(k));
    T_inv = [c^2,   s^2,    -2*c*s; 
             s^2,   c^2,   2*c*s; 
            c*s, -c*s, c^2-s^2];
    
    alpha_g = T_inv * [alpha1; alpha2; 0]; % transforming expansion coefficients from local → global
    beta_g  = T_inv * [beta1;  beta2;  0];
    Qb = Qbar_cell{k};
    dz  = z(k+1) - z(k);
    dz2 = z(k+1)^2 - z(k)^2;
    
    NT = NT + Delta_T * Qb * alpha_g * dz;
    NC = NC + Delta_C * Qb * beta_g  * dz;
    MT = MT + 0.5 * Delta_T * Qb * alpha_g * dz2;
    MC = MC + 0.5 * Delta_C * Qb * beta_g  * dz2;
end
N_hygro = NT + NC;
M_hygro = MT + MC;

%% 5. solving [ABD]{ε0, κ} = {N, M} gives mid-plane strains and curvature

N_total = [Nx; Ny; Nxy] + N_hygro;
M_total = [Mx; My; Mxy] + M_hygro;
midplane = ABD \ [N_total; M_total];
eps0  = midplane(1:3);   % Mid-plane strains
kappa = midplane(4:6);   % Curvatures

fprintf('Mid-plane Strains & Curvatures\n');
fprintf('ε_x=%+.4e  ε_y=%+.4e  γ_xy=%+.4e\n', eps0(1), eps0(2), eps0(3));
fprintf('κ_x=%+.4e  κ_y=%+.4e  κ_xy=%+.4e\n\n', kappa(1), kappa(2), kappa(3));

%% 6. Effective Engineering Constants

A_inv = inv(A);
D_inv = inv(D);

Ex_eXtu  = 1 / (h * A_inv(1,1));
Ey_eXtu  = 1 / (h * A_inv(2,2));
Gxy_eXtu = 1 / (h * A_inv(3,3));
vxy_eXtu = -A_inv(1,2) / A_inv(1,1);
vyx_eXtu = -A_inv(1,2) / A_inv(2,2);

Ex_bend  = 12 / (h^3 * D_inv(1,1));
Ey_bend  = 12 / (h^3 * D_inv(2,2));
Gxy_bend = 12 / (h^3 * D_inv(3,3));

fprintf('Effective Engineering Constants\n');
fprintf('--- Extension ---\n');
fprintf('  Ex  = %8.3f GPa\n', Ex_eXtu /1e9);
fprintf('  Ey  = %8.3f GPa\n', Ey_eXtu /1e9);
fprintf('  Gxy = %8.3f GPa\n', Gxy_eXtu/1e9);
fprintf('  νxy = %8.5f\n',     vxy_eXtu);
fprintf('  νyx = %8.5f\n',     vyx_eXtu);
fprintf('--- Bending ---\n');
fprintf('  Ex_b  = %8.3f GPa\n', Ex_bend /1e9);
fprintf('  Ey_b  = %8.3f GPa\n', Ey_bend /1e9);
fprintf('  Gxy_b = %8.3f GPa\n\n', Gxy_bend/1e9);

%% 7. Ply-by-Ply Stress / Strain State 

fprintf('Ply-by-Ply Stresses (mid-plane of each ply)\n');
fprintf('%4s %6s %12s %12s %12s\n','Ply','Angle','σ1 (MPa)','σ2 (MPa)','τ12 (MPa)');
for k = 1:NN
    z_mid = (z(k) + z(k+1)) / 2;
    eps_g = eps0 + z_mid * kappa;             
    eps_l = fn_strain_transform(theta(k), eps_g);   
    eps_l_mech = eps_l - Delta_T*[alpha1;alpha2;0] - Delta_C*[beta1;beta2;0];
    
    sig_l = Q_cell{k} * eps_l_mech;           
    fprintf('%4d %6.1f %12.4f %12.4f %12.4f\n', k, theta(k), ...
            sig_l(1)/1e6, sig_l(2)/1e6, sig_l(3)/1e6);
end
fprintf('\n');

%% User options 
% user selects the type of failure criteria
% Failure criterion: 1=Max Stress, 2=Max Strain, 3=Tsai-Hill
fail_criterion = input('Enter failure criterion (1=Max Stress, 2=Max Strain, 3=Tsai-Hill): ');

%% 8. Progressive Failure Analysis (Optimized Ratio Method)

fprintf('Progressive Failure Analysis\n');

for deg_mode = 1:2   % 1=Complete, 2=Partial
    
    fprintf('\nCriterion: %s | Degradation: %s\n', ...
        criterion_name(fail_criterion), deg_name(deg_mode));
    fprintf('%s\n', repmat('-',1,55));

    ply_failed = zeros(NN,1);
    FPF_load   = 0;
    LPF_load   = 0;

    % Reset material properties
    E1_d  = repmat(E1,  NN, 1);
    E2_d  = repmat(E2,  NN, 1);
    G12_d = repmat(G12, NN, 1);
    v21_d = repmat(v21, NN, 1);

    % Loop until all plies have failed
    while sum(ply_failed) < NN
        
        % 1. Rebuild degraded ABD
        [A_d, B_d, D_d, Q_d, ~] = fn_ABD_degraded(E1_d, E2_d, v12, v21_d, G12_d, theta, z, NN);
        ABD_d = [A_d, B_d; B_d, D_d];

        % Check for instability
        if rcond(ABD_d) < 1e-12
            fprintf('>> LPF (instability detected).\n');
            break;
        end

        % 2. Apply loads
        N_app = [Nx; Ny; Nxy] + N_hygro;
        M_app = [Mx; My; Mxy] + M_hygro;

        mp_d  = ABD_d \ [N_app; M_app];
        eps0_d  = mp_d(1:3);
        kappa_d = mp_d(4:6);

        FI_plies = zeros(NN,1);

        % 3. Failure index calculation
        for k = 1:NN
            if ply_failed(k) == 1
                continue;
            end

            z_mid = (z(k) + z(k+1))/2;
            eps_g = eps0_d + z_mid * kappa_d;
            eps_l = fn_strain_transform(theta(k), eps_g);
            eps_l_mech = eps_l - Delta_T*[alpha1;alpha2;0] - Delta_C*[beta1;beta2;0];

            sig_l  = Q_d{k} * eps_l_mech;
            sig1   = sig_l(1);  
            sig2   = sig_l(2);  
            tau12  = sig_l(3);

            eps1 = eps_l_mech(1); 
            eps2 = eps_l_mech(2); 
            gam12 = eps_l_mech(3);

            X_allow = select_strength(sig1, Xtu, Xcu);
            Y_allow = select_strength(sig2, Ytu, Ycu);

            switch fail_criterion
                case 1 % Max Stress
                    FI_plies(k) = max([abs(sig1/X_allow), abs(sig2/Y_allow), abs(tau12/S)]);
                    
                case 2 % Max Strain
                    X_eps = eps_Xtu; if eps1 < 0, X_eps = eps_Xcu; end
                    Y_eps = eps_Ytu; if eps2 < 0, Y_eps = eps_Ycu; end
                    FI_plies(k) = max([abs(eps1/X_eps), abs(eps2/Y_eps), abs(gam12/gamma_S)]);
                    
                case 3 % Tsai-Hill
                    TH_val = (sig1/X_allow)^2 + (sig2/Y_allow)^2 ...
                           - (sig1*sig2)/(X_allow^2) + (tau12/S)^2;
                    FI_plies(k) = sqrt(abs(TH_val));
            end
        end

        % 4. Find critical ply
        [max_FI, critical_ply] = max(FI_plies);

        if max_FI == 0 || isnan(max_FI)
            break;
        end

        current_load_factor = 1 / max_FI;
        ply_failed(critical_ply) = 1;

        % 5. Logging
        if FPF_load == 0
            FPF_load = current_load_factor * Nx;
            fprintf('>> FPF: Ply %d (%+g°) -> Nx = %.2f kN/m\n', ...
                critical_ply, theta(critical_ply), FPF_load/1e3);
        else
            fprintf('   Ply %d (%+g°) failed -> Nx = %.2f kN/m\n', ...
                critical_ply, theta(critical_ply), current_load_factor*Nx/1e3);
        end

        % 6. Degradation
        if deg_mode == 1
            % Complete
            E1_d(critical_ply)  = 1e-6;
            E2_d(critical_ply)  = 1e-6;
            G12_d(critical_ply) = 1e-6;
        else
            % Partial
            E2_d(critical_ply)  = 1e-6;
            G12_d(critical_ply) = 1e-6;
        end

        v21_d(critical_ply) = v12 * E2_d(critical_ply) / E1_d(critical_ply);
        LPF_load = current_load_factor * Nx;
    end

    fprintf('\nFPF = %.2f kN/m | LPF = %.2f kN/m | LPF/FPF = %.3f\n\n', ...
        FPF_load/1e3, LPF_load/1e3, LPF_load/max(FPF_load,1));
end

%% 9. Carpet Plot — Uniaxial Tensile Ultimate Strength

fprintf('Generating Carpet Plot\n');
pct_0  = 0:2:100;
pct_45 = 0:2:100;
np0 = length(pct_0);  np45 = length(pct_45);
strength_MPa = NaN(np0, np45);

[~, Qbar0]   = fn_reduced_stiffness(E1, E2, v12, v21, G12, 0);
[~, Qbar45]  = fn_reduced_stiffness(E1, E2, v12, v21, G12, 45);
[~, Qbarm45] = fn_reduced_stiffness(E1, E2, v12, v21, G12, -45);
[~, Qbar90]  = fn_reduced_stiffness(E1, E2, v12, v21, G12, 90);

eps_ult = Xtu / E1; 

for i = 1:np0
    for j = 1:np45
        p0  = pct_0(i);
        p45 = pct_45(j);
        p90 = 100 - p0 - p45;
        
        if p90 < -1e-9, continue; end 
        p90 = max(p90, 0);
        
        v0  = p0 / 100;
        v45 = p45 / 100;
        v90 = p90 / 100;
        
        Q_eff = v0*Qbar0 + (v45/2)*Qbar45 + (v45/2)*Qbarm45 + v90*Qbar90;
        S_eff = inv(Q_eff);     
        Ex = 1 / S_eff(1,1);    
        strength_MPa(i,j) = (Ex * eps_ult) / 1e6; 
    end
end

[X_mesh, Y_mesh] = meshgrid(pct_0, pct_45);
figure('Name','Carpet Plot','Color','w', 'Position', [100, 100, 700, 500]);
surf(X_mesh, Y_mesh, strength_MPa', 'EdgeColor', 'none'); 
shading interp; 
xlabel('% 0° Plies', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('% \pm45° Plies', 'FontSize', 12, 'FontWeight', 'bold');
zlabel('Ultimate Strength (MPa)', 'FontSize', 12, 'FontWeight', 'bold');
title('Carpet Plot — Uniaxial Tensile Ultimate Strength', 'FontSize', 13);
colorbar; view(-45, 30); grid on; colormap jet;

figure('Name','Carpet Plot (Contour)','Color','w', 'Position', [850, 100, 600, 500]);
contourf(X_mesh, Y_mesh, strength_MPa', 30, 'LineColor', 'none');
hold on;
contour(X_mesh, Y_mesh, strength_MPa', 15, 'k', 'LineStyle', ':', 'LineWidth', 0.5);
xlabel('% 0° Plies', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('% \pm45° Plies', 'FontSize', 12, 'FontWeight', 'bold');
title('Carpet Plot — Uniaxial Tensile Ultimate Strength', 'FontSize', 13);
colorbar; colormap jet; hold off;
fprintf('Carpet plot complete.\n\n');


%% LOCAL FUNCTIONS 

function [A, B, D, Q_cell, Qbar_cell] = fn_ABD(E1, E2, v12, v21, G12, theta, z, NN)
    A = zeros(3); B = zeros(3); D = zeros(3);
    Q_cell = cell(NN,1);  Qbar_cell = cell(NN,1);
    for k = 1:NN
        [Q, Qbar] = fn_reduced_stiffness(E1, E2, v12, v21, G12, theta(k));
        Q_cell{k} = Q;  Qbar_cell{k} = Qbar;
        zk = z(k+1);  zk1 = z(k);
        A = A + Qbar * (zk - zk1);
        B = B + 0.5 * Qbar * (zk^2 - zk1^2);
        D = D + (1/3) * Qbar * (zk^3 - zk1^3);
    end
end

function [A, B, D, Q_cell, Qbar_cell] = fn_ABD_degraded(E1_arr, E2_arr, v12, v21_arr, G12_arr, theta, z, NN)
    A = zeros(3); B = zeros(3); D = zeros(3);
    Q_cell = cell(NN,1);  Qbar_cell = cell(NN,1);
    for k = 1:NN
        [Q, Qbar] = fn_reduced_stiffness(E1_arr(k), E2_arr(k), v12, v21_arr(k), G12_arr(k), theta(k));
        Q_cell{k} = Q;  Qbar_cell{k} = Qbar;
        zk = z(k+1);  zk1 = z(k);
        A = A + Qbar * (zk - zk1);
        B = B + 0.5 * Qbar * (zk^2 - zk1^2);
        D = D + (1/3) * Qbar * (zk^3 - zk1^3);
    end
end

function [Q, Qbar] = fn_reduced_stiffness(E1, E2, v12, v21, G12, angle_deg)
    denom = 1 - v12*v21;
    Q11 = E1 / denom;   Q22 = E2 / denom;
    Q12 = v12*E2 / denom;   Q66 = G12;
    Q = [Q11, Q12, 0; Q12, Q22, 0; 0, 0, Q66];
    
    c = cosd(angle_deg);  s = sind(angle_deg);
    c2 = c^2; s2 = s^2; cs = c*s;
    
    Qbar11 = Q11*c2^2 + 2*(Q12 + 2*Q66)*s2*c2 + Q22*s2^2;
    Qbar12 = (Q11 + Q22 - 4*Q66)*c2*s2 + Q12*(s2^2 + c2^2);
    Qbar22 = Q11*s2^2 + 2*(Q12 + 2*Q66)*c2*s2 + Q22*c2^2;
    Qbar16 = (Q11 - Q12 - 2*Q66)*cs*c2 - (Q22 - Q12 - 2*Q66)*cs*s2;
    Qbar26 = (Q11 - Q12 - 2*Q66)*cs*s2 - (Q22 - Q12 - 2*Q66)*cs*c2;
    Qbar66 = (Q11 + Q22 - 2*Q12 - 2*Q66)*c2*s2 + Q66*(s2^2 + c2^2);
    
    Qbar = [Qbar11, Qbar12, Qbar16;
            Qbar12, Qbar22, Qbar26;
            Qbar16, Qbar26, Qbar66];
end

function eps_local = fn_strain_transform(angle_deg, eps_global)
    c = cosd(angle_deg);  s = sind(angle_deg);
    T_eps = [c^2,    s^2,   2*c*s; ...
             s^2,    c^2,   -2*c*s; ...
             -c*s,   c*s,   c^2-s^2];
    eps_local = T_eps * eps_global;
end
function X = select_strength(sigma, X_tension, X_compression)
    if sigma >= 0
        X = X_tension;
    else
        X = X_compression;
    end
end
function s = criterion_name(c)
    names = {'Maximum Stress','Maximum Strain','Tsai-Hill'};
    s = names{c};
end

function s = deg_name(d)
    names = {'Complete','Partial'};
    s = names{d};
end