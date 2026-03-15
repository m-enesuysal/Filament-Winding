%% PATTERN_SOLVER  Phase-2b S1: Angle-driven pattern çözücü.
%
% Ellipsoidal dome (TEST-02 geometrisi) ile uyumlu p/q çiftleri bulur.
%
% Parametreler:
%   R_eq = 152.4 mm, r0 = 45 mm, L_cyl = 300 mm
%   k = 0.7 (ellipsoidal aspect ratio)
%   BW_eff = 10 mm, d = 1
%   α_eq ≈ 19.15° (R_E = 50 mm'den türetilir)
%   coverage_range = [100, 150]%
%
% Bağımlılıklar:
%   ../phase1a_geometry/ellipsoidal_dome_profile.m
%   ../phase2a_winding/geodesic_single_circuit.m
%   find_compatible_patterns.m
%
% Referans: Karar-13 (angle-driven mod), S-PAT-01, S-PAT-02
% Tarih: 2026-03-15
% Faz: Phase-2b S1

clear; clc;
fprintf('=== Phase-2b S1: Angle-Driven Pattern Solver ===\n\n');

%% --- Parametreler ---
R_eq   = 152.4;    % Ekvator yarıçapı [mm]
r0     = 45.0;     % Polar açıklık yarıçapı [mm]
L_cyl  = 300.0;    % Silindir uzunluğu [mm]
k      = 0.7;      % Ellipsoidal dome aspect ratio [-]
BW_eff = 10.0;     % Efektif bant genişliği [mm]
d      = 1;        % Winding pattern skip index [-]
coverage_range = [100, 150];   % Kaplama aralığı [%]

% Yol bağımlılıkları ekle
script_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, '..', 'phase1a_geometry'));
addpath(fullfile(script_dir, '..', 'phase2a_winding'));

%% --- Girdi doğrulama ---
fprintf('--- Girdi Parametreleri ---\n');
fprintf('  R_eq   = %.1f mm\n', R_eq);
fprintf('  r0     = %.1f mm\n', r0);
fprintf('  L_cyl  = %.1f mm\n', L_cyl);
fprintf('  k      = %.2f\n', k);
fprintf('  BW_eff = %.1f mm\n', BW_eff);
fprintf('  d      = %d\n', d);
fprintf('  Coverage range = [%.0f, %.0f]%%\n\n', coverage_range(1), coverage_range(2));

% L_cyl < 0 kontrolü
if L_cyl < 0
    error('pattern_solver:negativeLcyl', ...
          'L_cyl >= 0 olmalıdır. Girilen: %g mm', L_cyl);
end

% Türetilen parametreler
R_E = r0 + BW_eff / 2;
alpha_eq = asin(R_E / R_eq);

fprintf('--- Türetilen Parametreler ---\n');
fprintf('  R_E      = r0 + BW_eff/2 = %.2f mm\n', R_E);
fprintf('  α_eq     = arcsin(R_E/R_eq) = %.4f rad = %.2f°\n', alpha_eq, alpha_eq*180/pi);
fprintf('  R_E/R_eq = %.4f\n', R_E / R_eq);

% Kaplama referans değerleri
n_max = floor(2 * pi * R_eq * d * cos(alpha_eq) / BW_eff);
p_for_100 = pi * R_eq * d * cos(alpha_eq) / BW_eff;
fprintf('  n_max    = floor(2π·R_eq·d·cos(α_eq)/BW_eff) = %d\n', n_max);
fprintf('  p_100%%   = π·R_eq·d·cos(α_eq)/BW_eff = %.2f\n\n', p_for_100);

%% === ADIM 1: Dome profili oluştur ===
fprintf('--- Adım 1: Ellipsoidal Dome Profili ---\n');
N_dome = 1000;   % Dome profil çözünürlüğü
dome = ellipsoidal_dome_profile(R_eq, r0, k, N_dome);

fprintf('  Dome yüksekliği  = %.2f mm\n', dome.h_dome);
fprintf('  Dome yay uzunluğu = %.2f mm\n', dome.s_total);
fprintf('  θ_p = %.4f rad = %.2f°\n', dome.theta_p, dome.theta_p*180/pi);
fprintf('  κ_m(ekvator) = %.6f 1/mm\n', dome.kappa_m(1));
fprintf('  κ_m(polar)   = %.6f 1/mm\n\n', dome.kappa_m(end));

%% === ADIM 2: Geodesic single circuit ===
fprintf('--- Adım 2: Geodesic Single Circuit ---\n');
geo = geodesic_single_circuit(dome, R_eq, r0, L_cyl, BW_eff);

fprintf('  φ_dome (numerical) = %.6f rad = %.4f°\n', geo.phi_dome_numerical, geo.phi_dome_numerical*180/pi);
fprintf('  φ_tail (analytic)  = %.6f rad = %.4f°\n', geo.phi_tail, geo.phi_tail*180/pi);
fprintf('  φ_dome (total)     = %.6f rad = %.4f°\n', geo.phi_dome, geo.phi_dome*180/pi);
fprintf('  φ_cyl              = %.6f rad = %.4f°\n', geo.phi_cyl, geo.phi_cyl*180/pi);
fprintf('  ε (singülarite)    = %.1e mm\n', geo.epsilon);
fprintf('  Δφ_circuit = 4·φ_dome + 2·φ_cyl\n');
fprintf('             = %.6f rad = %.4f°\n', geo.delta_phi_circuit, geo.delta_phi_circuit*180/pi);
fprintf('  Δφ/(2π)    = %.6f (ideal q/p oranı)\n\n', geo.delta_phi_circuit / (2*pi));

%% === ADIM 3: Uyumlu p/q pattern arama ===
fprintf('--- Adım 3: Pattern Arama (Angle-Driven Mod) ---\n');
patterns = find_compatible_patterns(geo.delta_phi_circuit, alpha_eq, R_eq, BW_eff, d, coverage_range);

if isempty(patterns)
    fprintf('  *** UYARI: Uyumlu pattern bulunamadı! ***\n');
    fprintf('  Coverage aralığını genişletin veya L_cyl ince ayar yapın.\n');
    return;
end

%% === ADIM 4: Sonuçları göster ===
fprintf('\n  Bulunan %d uyumlu pattern (angular_error sıralı):\n\n', numel(patterns));
fprintf('  %4s  %4s  %4s  %10s  %10s  %12s  %12s\n', ...
        'p', 'q', 'n', 'Coverage%', 'Overlap%', 'AngErr[rad]', 'AngErr[°]');
fprintf('  %s\n', repmat('-', 1, 68));

for i = 1:numel(patterns)
    pat = patterns(i);
    fprintf('  %4d  %4d  %4d  %10.2f  %10.2f  %12.6f  %12.4f\n', ...
            pat.p, pat.q, pat.n, pat.coverage_pct, pat.overlap_pct, ...
            pat.angular_error, pat.angular_error_deg);
end

%% === ADIM 5: En iyi pattern detayları ===
best = patterns(1);
fprintf('\n--- En İyi Pattern (minimum açısal hata) ---\n');
fprintf('  p = %d, q = %d, n = %d\n', best.p, best.q, best.n);
fprintf('  gcd(p,q) = %d (S-PAT-02: %s)\n', gcd(best.p, best.q), ...
        iftrue(gcd(best.p, best.q) == 1, 'PASS', 'FAIL'));
fprintf('  Δφ_ideal = 2π·%d/%d = %.6f rad = %.4f°\n', ...
        best.q, best.p, best.delta_phi_ideal, best.delta_phi_ideal*180/pi);
fprintf('  Δφ_actual = %.6f rad\n', geo.delta_phi_circuit);
fprintf('  Açısal hata = %.6f rad = %.4f°\n', best.angular_error, best.angular_error_deg);
fprintf('  Kaplama = %.2f%%\n', best.coverage_pct);
fprintf('  Örtüşme = %.2f%%\n', best.overlap_pct);

% Pattern özet bilgi
fprintf('\n--- Pattern Özet ---\n');
fprintf('  Bir tam pattern %d devreden oluşur.\n', best.p);
fprintf('  Her devre mandrel etrafında ~%.2f tur atar.\n', best.q);
fprintf('  Pattern tamamlandığında fiber başlangıç noktasına döner.\n');

% Bant yerleşim görseli (açısal adım)
angular_step = 2 * pi / best.p;   % Komşu bantlar arası açısal mesafe [rad]
fprintf('  Komşu bantlar arası açısal mesafe = 2π/%d = %.4f rad = %.2f°\n', ...
        best.p, angular_step, angular_step*180/pi);
fprintf('  Çevresel bant aralığı (ekvator) = %.2f mm\n', R_eq * angular_step);
fprintf('  Bant genişliği / cos(α_eq) = %.2f mm (dik footprint)\n', ...
        BW_eff / cos(alpha_eq));
fprintf('  n_max (tam kaplama) = %d bant\n', n_max);

fprintf('\n=== Phase-2b S1 TAMAMLANDI ===\n');


%% === Yardımcı fonksiyon ===
function str = iftrue(cond, true_str, false_str)
    if cond
        str = true_str;
    else
        str = false_str;
    end
end
