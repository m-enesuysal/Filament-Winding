%% VERIFY_GEODESIC_PATH  Phase-2a GATE-2a dogrulama — T1-T11 testleri.
%
% 3 dome tipi (hemispherical, ellipsoidal, isotensoid) x 4 senaryo = 12 konfig.
% Her konfigurasyonda T1-T11 testleri calistirilir.
%
% Referans: docs/phase2a_geodesic_math.md Bolum 15.3
% Tarih: 2026-03-10
% Faz: Phase-2a

clear; clc; close all;

%% --- Path ayarlari ---
this_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(this_dir, 'utils'));
addpath(fullfile(this_dir, '..', 'phase1a_geometry'));

%% --- Toleranslar (Karar-11) ---
TOL_CLAIRAUT   = 1e-4;    % T1, T11: Clairaut invariant [mm]
TOL_TURN_ALPHA = 1e-3;    % T2: |alpha_turn - pi/2| [rad]
TOL_TRANSITION = 1e-4;    % T3: silindir-dome gecis [rad] (spline sinir etkisi)
TOL_HEMI_PHI   = 0.01;    % T4: hemispherical phi_dome bagil hata
TOL_EPS_CONV   = 1e-3;    % T5: epsilon yakinsama bagil hata
TOL_CROSS_VAL  = 0.02;    % T6: capraz dogrulama bagil hata (%2)
MU_WET_CARBON  = 0.11;    % T10: islak sarim C/E surtunme katsayisi

total_pass = 0;
total_fail = 0;
test_log = {};

fprintf('=============================================================\n');
fprintf(' PHASE-2a: GEODESIC YOL DOGRULAMA RAPORU (GATE-2a)\n');
fprintf(' Tarih: %s\n', datestr(now, 'yyyy-mm-dd HH:MM'));
fprintf('=============================================================\n\n');

%% =====================================================================
%  TEST SENARYOLARI TANIMI
%  =====================================================================
% 4 senaryo — paylasilmis geometri parametreleri
scenarios = struct( ...
    'name',   {'S1-COPV-Standart', 'S2-Buyuk',      'S3-Kucuk',      'S4-GenisAciklik'}, ...
    'R_eq',   {152.4,               200,              100,              120}, ...
    'r0',     {45,                   60,               20,               50}, ...
    'L_cyl',  {300,                  400,              150,              250}, ...
    'BW_eff', {10,                   12,               6,                8}, ...
    'k_ell',  {0.7,                  0.5,              1.2,              0.8} ...
);

dome_types = {'hemispherical', 'ellipsoidal', 'isotensoid'};

%% =====================================================================
%  ANA TEST DONGUSU: 3 dome x 4 senaryo
%  =====================================================================
for dt = 1:numel(dome_types)
    dtype = dome_types{dt};
    fprintf('\n###########################################################\n');
    fprintf(' DOME TIPI: %s\n', upper(dtype));
    fprintf('###########################################################\n');

    for sc = 1:numel(scenarios)
        S = scenarios(sc);
        fprintf('\n--- %s | %s ---\n', dtype, S.name);
        fprintf('    R_eq=%.1f, r0=%.1f, L_cyl=%.0f, BW_eff=%.1f', ...
                S.R_eq, S.r0, S.L_cyl, S.BW_eff);

        R_E = S.r0 + S.BW_eff / 2;

        % R_E < R_eq kontrolu
        if R_E >= S.R_eq
            fprintf(' [SKIP: R_E >= R_eq]\n');
            continue;
        end
        fprintf(', R_E=%.1f\n', R_E);

        %% --- Dome profil uret ---
        try
            switch dtype
                case 'hemispherical'
                    dome = hemispherical_dome_profile(S.R_eq, S.r0, 1000);
                case 'ellipsoidal'
                    dome = ellipsoidal_dome_profile(S.R_eq, S.r0, S.k_ell, 1000);
                case 'isotensoid'
                    dome = isotensoid_dome_profile(S.R_eq, S.r0, 1000);
            end
        catch ME
            fprintf('  [ERROR] Profil uretimi basarisiz: %s\n', ME.message);
            continue;
        end

        %% --- Geodesic yol uret ---
        try
            circuit = geodesic_single_circuit(dome, S.R_eq, S.r0, ...
                                               S.L_cyl, S.BW_eff, 2000);
            close all;  % Figur pencerelerini kapat
        catch ME
            fprintf('  [ERROR] Yol uretimi basarisiz: %s\n', ME.message);
            continue;
        end

        %% --- T1: Clairaut invariant ---
        err_t1 = circuit.clairaut_max_err;
        [total_pass, total_fail, test_log] = log_test(total_pass, total_fail, test_log, ...
            'T1', dtype, S.name, err_t1 < TOL_CLAIRAUT, ...
            sprintf('max|rho*sin(a)-R_E| = %.2e mm (limit: %.0e)', err_t1, TOL_CLAIRAUT));

        %% --- T2: Dome donus noktasi alpha -> pi/2 ---
        % Donus noktasi: minimum rho (R_E'ye en yakin)
        [rho_min, idx_min] = min(circuit.rho);
        alpha_turn = circuit.alpha(idx_min);
        err_t2 = abs(alpha_turn - pi/2);
        [total_pass, total_fail, test_log] = log_test(total_pass, total_fail, test_log, ...
            'T2', dtype, S.name, err_t2 < TOL_TURN_ALPHA, ...
            sprintf('|alpha_turn - pi/2| = %.2e rad (rho_min=%.2f)', err_t2, rho_min));

        %% --- T3: Silindir-dome gecis surekliligi ---
        % Ekvator noktalarinda (x~0 ve x~L_cyl) alpha ~ alpha_eq olmali
        alpha_eq = circuit.alpha_eq;
        eq_mask = abs(circuit.x) < 1 | abs(circuit.x - S.L_cyl) < 1;
        if any(eq_mask)
            alpha_at_eq = circuit.alpha(eq_mask);
            err_t3 = max(abs(alpha_at_eq - alpha_eq));
        else
            err_t3 = 0;
        end
        [total_pass, total_fail, test_log] = log_test(total_pass, total_fail, test_log, ...
            'T3', dtype, S.name, err_t3 < TOL_TRANSITION, ...
            sprintf('max|alpha_eq_gecis - alpha_eq| = %.2e rad', err_t3));

        %% --- T4: Hemispherical phi_dome ~ pi/2 ---
        if strcmp(dtype, 'hemispherical')
            err_t4 = abs(circuit.phi_dome - pi/2) / (pi/2);
            [total_pass, total_fail, test_log] = log_test(total_pass, total_fail, test_log, ...
                'T4', dtype, S.name, err_t4 < TOL_HEMI_PHI, ...
                sprintf('phi_dome=%.6f, |bagil hata|=%.2e (limit: %.0e)', ...
                        circuit.phi_dome, err_t4, TOL_HEMI_PHI));
        end

        %% --- T5: Epsilon yakinsama ---
        epsilons = [1e-4, 1e-3, 1e-2, 1e-1];
        phi_eps = zeros(size(epsilons));
        for ei = 1:numel(epsilons)
            try
                res_eps = dome_phi_integration(dome, R_E, epsilons(ei));
                phi_eps(ei) = res_eps.phi_dome;
            catch
                phi_eps(ei) = NaN;
            end
        end
        % En kucuk epsilon referans, diger 3 onunla karsilastirilir
        if ~isnan(phi_eps(1))
            phi_ref = phi_eps(1);
            max_conv_err = max(abs(phi_eps(2:end) - phi_ref) / abs(phi_ref));
            pass_t5 = max_conv_err < TOL_EPS_CONV;
        else
            max_conv_err = NaN;
            pass_t5 = false;
        end
        [total_pass, total_fail, test_log] = log_test(total_pass, total_fail, test_log, ...
            'T5', dtype, S.name, pass_t5, ...
            sprintf('eps yakinsama: phi=[%.6f, %.6f, %.6f, %.6f], max_err=%.2e', ...
                    phi_eps(1), phi_eps(2), phi_eps(3), phi_eps(4), max_conv_err));

        %% --- T6: Capraz dogrulama ---
        % Hemispherical: phi_dome -> pi/2 (analitik)
        % Ellipsoidal/Isotensoid: delta_phi_circuit = 4*phi_dome + 2*phi_cyl
        %   formul tutarliligi — bagimsiz yeniden hesap
        phi_dome_indep = dome_phi_integration(dome, R_E);
        phi_cyl_indep  = S.L_cyl * tan(asin(R_E / S.R_eq)) / S.R_eq;
        dphi_indep     = 4 * phi_dome_indep.phi_dome + 2 * phi_cyl_indep;
        err_t6 = abs(circuit.delta_phi_circuit - dphi_indep) / abs(dphi_indep);
        [total_pass, total_fail, test_log] = log_test(total_pass, total_fail, test_log, ...
            'T6', dtype, S.name, err_t6 < TOL_CROSS_VAL, ...
            sprintf('delta_phi: circuit=%.6f, bagimsiz=%.6f, bagil_err=%.2e', ...
                    circuit.delta_phi_circuit, dphi_indep, err_t6));

        %% --- T7: Yol 3D surekliligi ---
        % Ardisik noktalar arasi mesafe — sicrama kontrolu
        dx = diff(circuit.x);
        drho = diff(circuit.rho);
        dphi = diff(circuit.phi);
        rho_mid = (circuit.rho(1:end-1) + circuit.rho(2:end)) / 2;
        ds_3d = sqrt(drho.^2 + (rho_mid .* dphi).^2 + dx.^2);
        ds_median = median(ds_3d);
        ds_max = max(ds_3d);
        jump_ratio = ds_max / (ds_median + 1e-30);
        pass_t7 = jump_ratio < 50;  % max 50x median (adaptif ODE adim dagilimi)
        [total_pass, total_fail, test_log] = log_test(total_pass, total_fail, test_log, ...
            'T7', dtype, S.name, pass_t7, ...
            sprintf('ds_max/ds_median = %.2f (limit: 50)', jump_ratio));

        %% --- T9: Konveksite kontrolu (S-WIND-04) ---
        n_bridging = numel(circuit.bridging_risk_indices);
        % Isotensoid'de km isaret degistirir — beklenen durum
        if strcmp(dtype, 'isotensoid') && n_bridging > 0
            pass_t9 = true;  % isotensoid'de beklenen
            note_t9 = sprintf('%d nokta (isotensoid km<0 beklenir)', n_bridging);
        else
            pass_t9 = (n_bridging == 0);
            note_t9 = sprintf('%d nokta', n_bridging);
        end
        [total_pass, total_fail, test_log] = log_test(total_pass, total_fail, test_log, ...
            'T9', dtype, S.name, pass_t9, ...
            sprintf('bridging risk: %s', note_t9));

        %% --- T10: Kayma emniyet marji (S-WIND-05) ---
        pass_t10 = circuit.lambda_max < MU_WET_CARBON;
        [total_pass, total_fail, test_log] = log_test(total_pass, total_fail, test_log, ...
            'T10', dtype, S.name, pass_t10, ...
            sprintf('lambda_max = %.2e (limit: %.2f)', circuit.lambda_max, MU_WET_CARBON));

        %% --- T11: Yeniden ornekleme sonrasi Clairaut ---
        % T1 ile ayni — yol zaten yeniden orneklenmis cikti
        err_t11 = max(abs(circuit.rho .* sin(circuit.alpha) - R_E));
        [total_pass, total_fail, test_log] = log_test(total_pass, total_fail, test_log, ...
            'T11', dtype, S.name, err_t11 < TOL_CLAIRAUT, ...
            sprintf('resampled max|rho*sin(a)-R_E| = %.2e mm', err_t11));

    end  % scenarios
end  % dome_types

%% =====================================================================
%  T8: Hoop Winding Kisit Kontrolu (ozel test)
%  =====================================================================
fprintf('\n###########################################################\n');
fprintf(' T8: HOOP WINDING KISIT KONTROLU\n');
fprintf('###########################################################\n');

% Hoop durumu: R_E ~ R_eq -> alpha_eq > 85 deg
R_eq_hoop = 100;
r0_hoop   = 95;      % cok buyuk polar aciklik
BW_hoop   = 8;       % R_E = 95 + 4 = 99 -> alpha_eq = arcsin(99/100) = 81.9 deg
% Daha agresif: r0=98, BW=3 -> R_E = 99.5 -> alpha_eq ~84.3 deg
% alpha > 85 icin: R_E/R_eq > sin(85) = 0.9962 -> R_E > 99.62
r0_hoop   = 97;
BW_hoop   = 6;       % R_E = 97 + 3 = 100 >= R_eq -> gecersiz (R_E >= R_eq)
% R_E = R_eq icin dome_phi_integration hata vermeli
dome_hoop = hemispherical_dome_profile(R_eq_hoop, r0_hoop, 200);

hoop_threw_error = false;
try
    geodesic_single_circuit(dome_hoop, R_eq_hoop, r0_hoop, 100, BW_hoop, 500);
    close all;
catch ME
    hoop_threw_error = true;
    fprintf('  Beklenen hata yakalandi: %s\n', ME.message);
end
[total_pass, total_fail, test_log] = log_test(total_pass, total_fail, test_log, ...
    'T8', 'hoop', 'kisit', hoop_threw_error, ...
    sprintf('Hoop kisiti: R_E=%.1f >= R_eq=%.1f -> hata bekleniyor', ...
            r0_hoop + BW_hoop/2, R_eq_hoop));

% alpha_eq > 85 deg testi
r0_hoop2  = 96;
BW_hoop2  = 7;       % R_E = 96 + 3.5 = 99.5 -> alpha_eq = arcsin(99.5/100) = 84.3 deg
% Bu sinirda, 85 altinda — bakmamiz lazim gercek aci
alpha_hoop2 = asind(99.5/100);
fprintf('  alpha_eq = %.1f deg (esik: 85 deg)\n', alpha_hoop2);
if alpha_hoop2 < 85
    % 85 deg ustu icin daha agresif
    r0_hoop2  = 96.5;
    BW_hoop2  = 7;    % R_E = 100 -> yine >= R_eq
end

% Kesin > 85 deg: R_E = 99.7, alpha = arcsin(0.997) = 85.6 deg
r0_hoop3  = 96;
BW_hoop3  = 7.4;     % R_E = 96 + 3.7 = 99.7 -> alpha = 85.6 deg
dome_hoop3 = hemispherical_dome_profile(R_eq_hoop, r0_hoop3, 200);
hoop_threw_error2 = false;
try
    geodesic_single_circuit(dome_hoop3, R_eq_hoop, r0_hoop3, 100, BW_hoop3, 500);
    close all;
catch ME
    hoop_threw_error2 = true;
    fprintf('  Beklenen hata yakalandi: %s\n', ME.message);
end
[total_pass, total_fail, test_log] = log_test(total_pass, total_fail, test_log, ...
    'T8', 'hoop', 'alpha>85', hoop_threw_error2, ...
    sprintf('alpha_eq=%.1f deg > 85 -> hata bekleniyor', ...
            asind((r0_hoop3 + BW_hoop3/2) / R_eq_hoop)));

%% =====================================================================
%  SONUC OZETI
%  =====================================================================
fprintf('\n=============================================================\n');
fprintf(' GATE-2a SONUC OZETI\n');
fprintf('=============================================================\n');
fprintf(' Toplam PASS : %d\n', total_pass);
fprintf(' Toplam FAIL : %d\n', total_fail);
fprintf(' Toplam test : %d\n', total_pass + total_fail);
fprintf('=============================================================\n\n');

% Detayli tablo
fprintf('%-6s %-15s %-20s %-6s %s\n', ...
        'Test', 'Dome', 'Senaryo', 'Sonuc', 'Detay');
fprintf('%s\n', repmat('-', 1, 90));
for i = 1:numel(test_log)
    t = test_log{i};
    if t.pass
        status = 'PASS';
    else
        status = 'FAIL';
    end
    fprintf('%-6s %-15s %-20s %-6s %s\n', ...
            t.test_id, t.dome, t.scenario, status, t.detail);
end

fprintf('\n');
if total_fail == 0
    fprintf('*** GATE-2a: TUMU PASS ***\n');
else
    fprintf('*** GATE-2a: %d FAIL — inceleme gerekli ***\n', total_fail);
end

%% =====================================================================
%  YARDIMCI FONKSIYON
%  =====================================================================
function [tp, tf, lg] = log_test(tp, tf, lg, test_id, dome, scenario, pass, detail)
    entry.test_id  = test_id;
    entry.dome     = dome;
    entry.scenario = scenario;
    entry.pass     = pass;
    entry.detail   = detail;
    lg{end+1} = entry;

    if pass
        tp = tp + 1;
        fprintf('  [PASS] %s: %s\n', test_id, detail);
    else
        tf = tf + 1;
        fprintf('  [FAIL] %s: %s\n', test_id, detail);
    end
end
