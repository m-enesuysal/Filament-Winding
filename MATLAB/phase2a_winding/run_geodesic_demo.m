%% RUN_GEODESIC_DEMO  Phase-2a S1 demo — tek devre geodesic yol
%
% Parametreler:
%   Elipsoidal dome, R_eq=152.4, r0=45, L_cyl=300, k=0.7, BW_eff=10mm
%
% Kullanim: MATLAB'da bu scripti calistirin.
%           phase1a_geometry/ ve phase2a_winding/utils/ path'te olmalidir.

clear; clc; close all;

%% --- Path ayarlari ---
this_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(this_dir, 'utils'));
addpath(fullfile(this_dir, '..', 'phase1a_geometry'));

%% --- Mandrel parametreleri ---
R_eq   = 152.4;   % Ekvator yaricapi [mm]
r0     = 45;      % Polar aciklik yaricapi [mm]
L_cyl  = 300;     % Silindir uzunlugu [mm]
k      = 0.7;     % Elipsoidal aspect ratio [-]
BW_eff = 10;      % Efektif bant genisligi [mm]

%% --- Phase-1a: Dome profili uret ---
fprintf('Phase-1a: Elipsoidal dome profili uretiliyor...\n');
dome = ellipsoidal_dome_profile(R_eq, r0, k, 1000);
fprintf('  s_total = %.4f mm, h_dome = %.4f mm\n', dome.s_total, dome.h_dome);
fprintf('  kappa_eq = %.6f 1/mm, kappa_pol = %.6f 1/mm\n', ...
        dome.kappa_eq, dome.kappa_pol);

%% --- Phase-2a: Geodesic tek devre ---
fprintf('\nPhase-2a: Geodesic tek devre hesaplaniyor...\n');
circuit = geodesic_single_circuit(dome, R_eq, r0, L_cyl, BW_eff, 2000);

%% --- Sonuc ozeti ---
fprintf('\n========== SONUC OZETI ==========\n');
fprintf('R_E           = %.4f mm\n', circuit.R_E);
fprintf('alpha_eq      = %.4f rad (%.2f deg)\n', circuit.alpha_eq, rad2deg(circuit.alpha_eq));
fprintf('phi_dome      = %.6f rad (%.2f deg)\n', circuit.phi_dome, rad2deg(circuit.phi_dome));
fprintf('phi_cyl       = %.6f rad (%.2f deg)\n', circuit.phi_cyl, rad2deg(circuit.phi_cyl));
fprintf('delta_phi     = %.6f rad (%.2f deg, %.4f tur)\n', ...
        circuit.delta_phi_circuit, ...
        rad2deg(circuit.delta_phi_circuit), ...
        circuit.delta_phi_circuit / (2*pi));
fprintf('Clairaut err  = %.2e mm\n', circuit.clairaut_max_err);
fprintf('lambda_max    = %.2e (S-WIND-05)\n', circuit.lambda_max);
fprintf('Bridging risk = %d nokta\n', numel(circuit.bridging_risk_indices));
fprintf('Toplam nokta  = %d\n', circuit.N_total);
fprintf('=================================\n');
