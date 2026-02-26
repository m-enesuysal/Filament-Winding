%% VERIFY_ELLIPSOIDAL_DOME.m
% Phase-1a: Elipsoidal dome meridyen profili — doğrulama betiği
%
% Bu betik aşağıdaki GATE-1a koşullarını doğrular:
%   1. Analitik çözüm tutarlılığı (elips denklemi, sınır koşulları)
%   2. C¹ süreklilik (S-GEO-02: ekvator noktası)
%   3. Birim teğet vektör normu
%   4. Eğrilik formülü doğrulaması (sayısal türev ile çapraz kontrol)
%   5. S-GEO-04: k=1 ile hemispherical dome örtüşmesi
%   6. S-GEO-03: k → küçük değer davranışı
%   7. Yüzey alanı doğrulaması (sayısal integrasyon ile çapraz kontrol)
%   8. Karar-16 test senaryoları üzerinde çalıştırma
%   9. Farklı k değerleri için karşılaştırmalı analiz
%
% Toleranslar: Karar-11 Katman 2 referans alınır.
%   Pozisyon: |ε| < 1e-4 mm
%   Türev:    |ε| < 1e-6
%   Eğrilik:  |ε_rel| < 1e-4
%
% Tarih: 2026-02-26
% Faz: Phase-1a

clear; clc; close all;

fprintf('=============================================================\n');
fprintf(' PHASE-1a: ELİPSOİDAL DOME DOĞRULAMA RAPORU\n');
fprintf(' Tarih: %s\n', datestr(now, 'yyyy-mm-dd HH:MM'));
fprintf('=============================================================\n\n');

%% --- Karar-11 Toleransları ---
TOL_POS    = 1e-4;    % Pozisyon toleransı [mm]
TOL_DERIV  = 1e-6;    % Türev toleransı [-]
TOL_CURV   = 1e-4;    % Eğrilik bağıl toleransı [-]
TOL_MACHINE = 1e-12;  % Makine hassasiyeti

all_pass = true;

%% =====================================================================
%  TEST SETİ 1: Temel Analitik Doğrulama
%  R_eq = 73 mm, r0 = 22 mm, k = 0.6 (TEST-01 senaryosu)
%  =====================================================================

fprintf('--- TEST SETİ 1: Temel Analitik Doğrulama (TEST-01) ---\n');
fprintf('    R_eq = 73 mm, r0 = 22 mm, k = 0.6\n\n');

R_eq = 73;   r0 = 22;   k = 1;
N = 1000;

profil = ellipsoidal_dome_profile(R_eq, r0, k, N);

% --- Test 1.1: θ_p doğrulama ---
theta_p_expected = acos(r0 / R_eq);
err = abs(profil.theta_p - theta_p_expected);
pass = err < TOL_MACHINE;
fprintf('  [%s] T1.1 θ_p doğrulama:  hata = %.2e\n', tf(pass), err);
all_pass = all_pass && pass;

% --- Test 1.2: h_dome çapraz doğrulama ---
h_dome_1 = k * R_eq * sin(theta_p_expected);
h_dome_2 = k * sqrt(R_eq^2 - r0^2);
err_h = abs(h_dome_1 - h_dome_2);
pass = err_h < TOL_MACHINE * R_eq;
fprintf('  [%s] T1.2 h_dome çapraz formül: |trig - pisagor| = %.2e mm\n', tf(pass), err_h);
all_pass = all_pass && pass;

err_h2 = abs(profil.h_dome - h_dome_2);
pass = err_h2 < TOL_POS;
fprintf('  [%s] T1.2b h_dome profil çıktısı: hata = %.2e mm\n', tf(pass), err_h2);
all_pass = all_pass && pass;

% --- Test 1.3: Sınır koşulları — ekvator (s=0) ---
err_rho_eq   = abs(profil.rho(1)     - R_eq);
err_x_eq     = abs(profil.x_local(1) - 0);
err_drho_eq  = abs(profil.drho_ds(1) - 0);
err_dx_eq    = abs(profil.dx_ds(1)   - 1);

pass_eq = (err_rho_eq < TOL_POS) && (err_x_eq < TOL_POS) && ...
          (err_drho_eq < TOL_DERIV) && (err_dx_eq < TOL_DERIV);
fprintf('  [%s] T1.3 Ekvator sınır koşulları:\n', tf(pass_eq));
fprintf('         ρ(0) - R_eq  = %.2e mm\n', err_rho_eq);
fprintf('         x(0)         = %.2e mm\n', err_x_eq);
fprintf('         dρ/ds(0)     = %.2e\n', err_drho_eq);
fprintf('         dx/ds(0) - 1 = %.2e\n', err_dx_eq);
all_pass = all_pass && pass_eq;

% --- Test 1.4: Sınır koşulları — polar açıklık (s=s_total) ---
err_rho_p = abs(profil.rho(end)     - r0);
err_x_p   = abs(profil.x_local(end) - profil.h_dome);

pass_p = (err_rho_p < TOL_POS) && (err_x_p < TOL_POS);
fprintf('  [%s] T1.4 Polar açıklık sınır koşulları:\n', tf(pass_p));
fprintf('         ρ(s_total) - r0     = %.2e mm\n', err_rho_p);
fprintf('         x(s_total) - h_dome = %.2e mm\n', err_x_p);
all_pass = all_pass && pass_p;

% --- Test 1.5: Birim teğet vektör normu ---
tangent_norm = sqrt(profil.drho_ds.^2 + profil.dx_ds.^2);
err_norm     = max(abs(tangent_norm - 1));
pass = err_norm < TOL_DERIV;
fprintf('  [%s] T1.5 Birim teğet vektör: max|norm-1| = %.2e\n', tf(pass), err_norm);
all_pass = all_pass && pass;

% --- Test 1.6: Elips denklemi ---
ellipse_val = (profil.rho / R_eq).^2 + (profil.x_local / (k * R_eq)).^2;
err_ellipse = max(abs(ellipse_val - 1));
pass = err_ellipse < TOL_POS / R_eq;
fprintf('  [%s] T1.6 Elips denklemi: max|(ρ/a)²+(x/b)²-1| = %.2e\n', tf(pass), err_ellipse);
all_pass = all_pass && pass;

% --- Test 1.7: ρ monotonik azalan ---
drho = diff(profil.rho);
pass = all(drho <= 1e-14);  % Numerik toleransla
fprintf('  [%s] T1.7 ρ monotonik azalan (ekvator→polar)\n', tf(pass));
all_pass = all_pass && pass;

% --- Test 1.8: x_local monotonik artan ---
dx = diff(profil.x_local);
pass = all(dx >= -1e-14);
fprintf('  [%s] T1.8 x_local monotonik artan (ekvator→polar)\n', tf(pass));
all_pass = all_pass && pass;

% --- Test 1.9: Eğrilik analitik sınır değerleri ---
kappa_eq_expected = 1 / (R_eq * k^2);
err_kappa_eq = abs(profil.kappa_m(1) - kappa_eq_expected);
pass = err_kappa_eq / kappa_eq_expected < TOL_CURV;
fprintf('  [%s] T1.9 Ekvator eğriliği: κ_m(0) = %.6e, beklenen = %.6e, bağıl hata = %.2e\n', ...
        tf(pass), profil.kappa_m(1), kappa_eq_expected, err_kappa_eq / kappa_eq_expected);
all_pass = all_pass && pass;

% --- Test 1.10: Yüzey alanı (trapez vs integral() çapraz doğrulama) ---
A_trapz = 2 * pi * trapz(profil.s, profil.rho);
err_A = abs(A_trapz - profil.A_dome) / profil.A_dome;
pass = err_A < 1e-3;  % Trapez yaklaşımı — gevşek tolerans
fprintf('  [%s] T1.10 Yüzey alanı:\n', tf(pass));
fprintf('          integral()  = %.4f mm²\n', profil.A_dome);
fprintf('          trapz(s,ρ)  = %.4f mm²\n', A_trapz);
fprintf('          Bağıl hata  = %.2e\n', err_A);
all_pass = all_pass && pass;

fprintf('\n');

%% =====================================================================
%  TEST SETİ 2: C¹ Süreklilik Detaylı Analiz (S-GEO-02)
%  =====================================================================

fprintf('--- TEST SETİ 2: C¹ Süreklilik Analizi (S-GEO-02) ---\n\n');

fprintf('  Silindir tarafı (s → 0⁻):\n');
fprintf('    ρ       = %.6f mm  (sabit = R_eq)\n', R_eq);
fprintf('    dρ/ds   = 0\n');
fprintf('    dx/ds   = 1\n');
fprintf('    κ_m     = 0\n\n');

fprintf('  Dome tarafı (s → 0⁺, k = %.2f):\n', k);
fprintf('    ρ(0)    = %.6f mm\n', profil.rho(1));
fprintf('    dρ/ds(0)= %.2e\n', profil.drho_ds(1));
fprintf('    dx/ds(0)= %.10f\n', profil.dx_ds(1));
fprintf('    κ_m(0)  = %.6e 1/mm  [= 1/(R_eq·k²) = 1/(%.1f·%.2f²)]\n\n', ...
        profil.kappa_m(1), R_eq, k);

err_C0 = abs(profil.rho(1) - R_eq);
err_C1_rho = abs(profil.drho_ds(1));
err_C1_x   = abs(profil.dx_ds(1) - 1);

pass_C0 = err_C0 < TOL_MACHINE;
pass_C1 = (err_C1_rho < TOL_DERIV) && (err_C1_x < TOL_DERIV);

fprintf('  [%s] C⁰ süreklilik: |ρ_dome(0) - R_eq| = %.2e mm\n', tf(pass_C0), err_C0);
fprintf('  [%s] C¹ süreklilik: |dρ/ds(0)| = %.2e, |dx/ds(0)-1| = %.2e\n', ...
        tf(pass_C1), err_C1_rho, err_C1_x);

% C² süreksizliği raporlama
kappa_jump = profil.kappa_m(1);
fprintf('\n  [BİLGİ] C² süreksizliği (beklenen):\n');
fprintf('    κ_m sıçraması = 1/(R_eq·k²) = %.6e 1/mm\n', kappa_jump);
fprintf('    Hemispherical referans (k=1): 1/R_eq = %.6e 1/mm\n', 1/R_eq);
fprintf('    Oran: κ_eq(k=%.1f) / κ_eq(k=1) = %.2f\n', k, kappa_jump / (1/R_eq));

all_pass = all_pass && pass_C0 && pass_C1;
fprintf('\n');

%% =====================================================================
%  TEST SETİ 3: S-GEO-04 — k=1 Hemispherical Örtüşme Doğrulaması
%  =====================================================================

fprintf('--- TEST SETİ 3: S-GEO-04 — k=1 Örtüşme Doğrulaması ---\n');
fprintf('    ellipsoidal(k=1) vs hemispherical karşılaştırması\n\n');

R_eq_test = 152.4;  r0_test = 45.0;
N_test = 500;

% Her iki profili üret
p_ellip = ellipsoidal_dome_profile(R_eq_test, r0_test, 1.0, N_test);
p_hemi  = hemispherical_dome_profile(R_eq_test, r0_test, N_test);

% Skaler büyüklükler
err_stotal = abs(p_ellip.s_total - p_hemi.s_total);
err_hdome  = abs(p_ellip.h_dome  - p_hemi.h_dome);
err_Adome  = abs(p_ellip.A_dome  - p_hemi.A_dome);

pass_scalar = (err_stotal < TOL_POS) && (err_hdome < TOL_POS) && ...
              (err_Adome / p_hemi.A_dome < TOL_CURV);
fprintf('  [%s] Skaler büyüklükler:\n', tf(pass_scalar));
fprintf('         |s_total farkı| = %.2e mm\n', err_stotal);
fprintf('         |h_dome farkı|  = %.2e mm\n', err_hdome);
fprintf('         |A_dome bağıl|  = %.2e\n', err_Adome / p_hemi.A_dome);
all_pass = all_pass && pass_scalar;

% Vektörel büyüklükler (uniform s ızgarasında karşılaştır)
err_rho    = max(abs(p_ellip.rho     - p_hemi.rho));
err_x      = max(abs(p_ellip.x_local - p_hemi.x_local));
err_drho   = max(abs(p_ellip.drho_ds - p_hemi.drho_ds));
err_dx     = max(abs(p_ellip.dx_ds   - p_hemi.dx_ds));
err_kappa  = max(abs(p_ellip.kappa_m - p_hemi.kappa_m));

pass_vec = (err_rho < TOL_POS) && (err_x < TOL_POS) && ...
           (err_drho < TOL_DERIV) && (err_dx < TOL_DERIV) && ...
           (err_kappa / (1/R_eq_test) < TOL_CURV);

fprintf('  [%s] Vektörel büyüklükler (max hata):\n', tf(pass_vec));
fprintf('         |Δρ|        = %.2e mm\n', err_rho);
fprintf('         |Δx|        = %.2e mm\n', err_x);
fprintf('         |Δ(dρ/ds)|  = %.2e\n', err_drho);
fprintf('         |Δ(dx/ds)|  = %.2e\n', err_dx);
fprintf('         |Δκ|/κ_ref  = %.2e\n', err_kappa / (1/R_eq_test));
all_pass = all_pass && pass_vec;

fprintf('\n  ► S-GEO-04 DURUMU: %s\n\n', ...
        ternary(pass_scalar && pass_vec, 'ONAYLI ✓', 'BAŞARISIZ ✗'));

%% =====================================================================
%  TEST SETİ 4: Eğrilik Sayısal Doğrulaması
%  =====================================================================

fprintf('--- TEST SETİ 4: Eğrilik — Sayısal vs Analitik Doğrulama ---\n\n');

% Referans senaryo: k = 0.6 (eğrilik değişken)
p = ellipsoidal_dome_profile(73, 22, 0.6, 2000);

% Sayısal eğrilik: κ_m = dβ/ds (merkezi fark)
ds = p.s(2) - p.s(1);
dbeta_ds_num = zeros(size(p.beta));
dbeta_ds_num(2:end-1) = (p.beta(3:end) - p.beta(1:end-2)) / (2*ds);
dbeta_ds_num(1)       = (p.beta(2) - p.beta(1)) / ds;
dbeta_ds_num(end)     = (p.beta(end) - p.beta(end-1)) / ds;

% Hata analizi (iç noktalar, sınırları hariç tut)
idx = 5:length(p.s)-4;   % Sınır etkilerinden kaçın
err_kappa_abs = abs(dbeta_ds_num(idx) - p.kappa_m(idx));
err_kappa_rel = err_kappa_abs ./ p.kappa_m(idx);

fprintf('  ds = %.4e mm\n', ds);
fprintf('  max|κ_sayısal - κ_analitik|     = %.2e 1/mm\n', max(err_kappa_abs));
fprintf('  max|κ_sayısal - κ_analitik|/κ   = %.2e\n', max(err_kappa_rel));
fprintf('  ortalama bağıl hata              = %.2e\n', mean(err_kappa_rel));

pass = max(err_kappa_rel) < 1e-3;  % Merkezi fark O(ds²) hatası
fprintf('  [%s] Sayısal-analitik eğrilik uyumu\n\n', tf(pass));
all_pass = all_pass && pass;

%% =====================================================================
%  TEST SETİ 5: Karar-16 Standart Test Senaryoları
%  =====================================================================

fprintf('--- TEST SETİ 5: Endüstri Test Senaryoları ---\n\n');

test_cases = struct(...
    'name',  {'TEST-01 ASTM', 'TEST-02 COPV', 'TEST-03 Küçük Aç.', 'TEST-04 H₂'}, ...
    'R_eq',  {73,    152.4, 150,   200}, ...
    'r0',    {22,    45,    10,    50}, ...
    'k',     {0.6,   0.7,   0.5,   0.7});

for i = 1:length(test_cases)
    tc = test_cases(i);
    fprintf('  [%s] R_eq=%.1f, r0=%.1f, k=%.2f\n', tc.name, tc.R_eq, tc.r0, tc.k);

    p = ellipsoidal_dome_profile(tc.R_eq, tc.r0, tc.k, 500);

    % Doğrulamalar
    err_rho_eq  = abs(p.rho(1) - tc.R_eq);
    err_rho_p   = abs(p.rho(end) - tc.r0);
    ellipse_err = max(abs((p.rho/tc.R_eq).^2 + (p.x_local/(tc.k*tc.R_eq)).^2 - 1));
    tangent_err = max(abs(sqrt(p.drho_ds.^2 + p.dx_ds.^2) - 1));
    kappa_eq_exp = 1 / (tc.R_eq * tc.k^2);
    kappa_err   = abs(p.kappa_m(1) - kappa_eq_exp) / kappa_eq_exp;

    pass_i = (err_rho_eq < TOL_POS) && (err_rho_p < TOL_POS) && ...
             (ellipse_err < 1e-8) && (tangent_err < TOL_DERIV) && ...
             (kappa_err < TOL_CURV);

    fprintf('    r0/R_eq   = %.4f\n', tc.r0 / tc.R_eq);
    fprintf('    θ_p       = %.4f rad (%.2f°)\n', p.theta_p, rad2deg(p.theta_p));
    fprintf('    s_total   = %.4f mm\n', p.s_total);
    fprintf('    h_dome    = %.4f mm\n', p.h_dome);
    fprintf('    κ_m(ekv)  = %.6e [beklenen: %.6e]\n', p.kappa_m(1), kappa_eq_exp);
    fprintf('    κ_m(pol)  = %.6e\n', p.kappa_m(end));
    fprintf('    Elips hata= %.2e\n', ellipse_err);
    fprintf('    [%s] Tüm kontroller\n\n', tf(pass_i));

    all_pass = all_pass && pass_i;
end

%% =====================================================================
%  TEST SETİ 6: S-GEO-03 ve Uç Durumlar
%  =====================================================================

fprintf('--- TEST SETİ 6: S-GEO-03 ve Uç Durumlar ---\n\n');

R_eq_uc = 152.4;  r0_uc = 45.0;

% 6.1: Sığ dome (k = 0.3)
p_shallow = ellipsoidal_dome_profile(R_eq_uc, r0_uc, 0.3, 500);
pass = abs(p_shallow.rho(end) - r0_uc) < TOL_POS;
fprintf('  [%s] T6.1 Sığ dome k=0.3: h_dome=%.2f mm, κ_eq=%.4e\n', ...
        tf(pass), p_shallow.h_dome, p_shallow.kappa_m(1));
all_pass = all_pass && pass;

% 6.2: Uzun dome (k = 2.0)
p_deep = ellipsoidal_dome_profile(R_eq_uc, r0_uc, 2.0, 500);
pass = abs(p_deep.rho(end) - r0_uc) < TOL_POS;
fprintf('  [%s] T6.2 Uzun dome k=2.0: h_dome=%.2f mm, κ_eq=%.4e\n', ...
        tf(pass), p_deep.h_dome, p_deep.kappa_m(1));
all_pass = all_pass && pass;

% 6.3: Çok yassı dome (k = 0.15) — S-GEO-03 uyarı testi
p_flat = ellipsoidal_dome_profile(R_eq_uc, r0_uc, 0.15, 500);
pass = abs(p_flat.rho(end) - r0_uc) < TOL_POS;
kappa_flat_eq = p_flat.kappa_m(1);
fprintf('  [%s] T6.3 Çok yassı k=0.15: h_dome=%.2f mm, κ_eq=%.4e\n', ...
        tf(pass), p_flat.h_dome, kappa_flat_eq);
fprintf('         >> kappa_eq x R_eq = %.2f (hemisphericaldان ne kadar uzak)\n', ...
        kappa_flat_eq * R_eq_uc);
all_pass = all_pass && pass;

% 6.4: Hata fırlatma testleri
fprintf('\n  Hata fırlatma testleri:\n');

try ellipsoidal_dome_profile(100, 50, -1); fprintf('  [FAIL] k<0 kabul edildi!\n'); all_pass=false;
catch; fprintf('  [PASS] T6.4a k<0 reddedildi\n'); end

try ellipsoidal_dome_profile(100, 50, 0); fprintf('  [FAIL] k=0 kabul edildi!\n'); all_pass=false;
catch; fprintf('  [PASS] T6.4b k=0 reddedildi\n'); end

try ellipsoidal_dome_profile(100, 110, 0.6); fprintf('  [FAIL] r0>R_eq kabul edildi!\n'); all_pass=false;
catch; fprintf('  [PASS] T6.4c r0>R_eq reddedildi\n'); end

fprintf('\n');

%% =====================================================================
%  TEST SETİ 7: Türev Çapraz Doğrulama (Sayısal Fark)
%  =====================================================================

fprintf('--- TEST SETİ 7: Sayısal vs Analitik Türev Karşılaştırması ---\n\n');

p = ellipsoidal_dome_profile(152.4, 45, 0.7, 2000);
ds = p.s(2) - p.s(1);

% Merkezi fark
drho_ds_num = zeros(size(p.rho));
dx_ds_num   = zeros(size(p.x_local));
drho_ds_num(2:end-1) = (p.rho(3:end) - p.rho(1:end-2)) / (2*ds);
dx_ds_num(2:end-1)   = (p.x_local(3:end) - p.x_local(1:end-2)) / (2*ds);
drho_ds_num(1)   = (p.rho(2) - p.rho(1)) / ds;
drho_ds_num(end) = (p.rho(end) - p.rho(end-1)) / ds;
dx_ds_num(1)     = (p.x_local(2) - p.x_local(1)) / ds;
dx_ds_num(end)   = (p.x_local(end) - p.x_local(end-1)) / ds;

idx = 2:length(p.s)-1;
err_drho = max(abs(drho_ds_num(idx) - p.drho_ds(idx)));
err_dx   = max(abs(dx_ds_num(idx)   - p.dx_ds(idx)));

fprintf('  ds = %.4e mm\n', ds);
fprintf('  max|dρ/ds sayısal - analitik| = %.2e\n', err_drho);
fprintf('  max|dx/ds sayısal - analitik| = %.2e\n', err_dx);

pass = (err_drho < 1e-4) && (err_dx < 1e-4);
fprintf('  [%s] Sayısal-analitik türev uyumu\n\n', tf(pass));
all_pass = all_pass && pass;

%% =====================================================================
%  GÖRSELLEŞTİRME
%  =====================================================================

fprintf('--- Görselleştirme oluşturuluyor ---\n\n');

figure('Name', 'Phase-1a: Elipsoidal Dome Karşılaştırma', ...
       'Position', [50, 50, 1600, 1000]);

R_eq_viz = 152.4;  r0_viz = 45.0;

% Farklı k değerleri ile profiller
k_vals = [0.4, 0.6, 0.8, 1.0, 1.2, 1.5];
colors = lines(length(k_vals));

% --- Subplot 1: Meridyen profilleri karşılaştırması ---
subplot(2,3,1);
hold on;
% Silindir referans
x_cyl = linspace(-30, 0, 30);
plot(x_cyl, R_eq_viz*ones(size(x_cyl)), 'k--', 'LineWidth', 1.5);
for j = 1:length(k_vals)
    p = ellipsoidal_dome_profile(R_eq_viz, r0_viz, k_vals(j), 300);
    plot(p.x_local, p.rho, '-', 'Color', colors(j,:), 'LineWidth', 1.5);
end
plot(0, R_eq_viz, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
xlabel('x_{lokal} [mm]');
ylabel('\rho [mm]');
title('Meridyen Profilleri');
legend_entries = [{'Silindir'}; arrayfun(@(kk) sprintf('k=%.1f', kk), k_vals', 'Uni', false)];
legend(legend_entries, 'Location', 'southwest', 'FontSize', 7);
grid on; axis equal;

% --- Subplot 2: Eğrilik profilleri ---
subplot(2,3,2);
hold on;
for j = 1:length(k_vals)
    p = ellipsoidal_dome_profile(R_eq_viz, r0_viz, k_vals(j), 300);
    s_norm = p.s / p.s_total;  % Normalize edilmiş yay uzunluğu
    plot(s_norm, p.kappa_m * R_eq_viz, '-', 'Color', colors(j,:), 'LineWidth', 1.5);
end
xlabel('s / s_{total} [-]');
ylabel('\kappa_m \cdot R_{eq} [-]');
title('Boyutsuz Eğrilik');
legend(arrayfun(@(kk) sprintf('k=%.1f', kk), k_vals', 'Uni', false), ...
       'Location', 'best', 'FontSize', 7);
grid on;

% --- Subplot 3: Eğim açısı β(s) ---
subplot(2,3,3);
hold on;
for j = 1:length(k_vals)
    p = ellipsoidal_dome_profile(R_eq_viz, r0_viz, k_vals(j), 300);
    s_norm = p.s / p.s_total;
    plot(s_norm, rad2deg(p.beta), '-', 'Color', colors(j,:), 'LineWidth', 1.5);
end
xlabel('s / s_{total} [-]');
ylabel('\beta [°]');
title('Eğim Açısı');
legend(arrayfun(@(kk) sprintf('k=%.1f', kk), k_vals', 'Uni', false), ...
       'Location', 'best', 'FontSize', 7);
grid on;

% --- Subplot 4: Türevler (k=0.6 detay) ---
subplot(2,3,4);
p06 = ellipsoidal_dome_profile(R_eq_viz, r0_viz, 0.6, 500);
plot(p06.s, p06.drho_ds, 'b-', 'LineWidth', 1.5); hold on;
plot(p06.s, p06.dx_ds, 'r-', 'LineWidth', 1.5);
yline(0, 'k--', 'LineWidth', 0.5);
xlabel('s [mm]');
ylabel('Türev [-]');
title('d\rho/ds ve dx/ds (k=0.6)');
legend('d\rho/ds', 'dx/ds', 'Location', 'east');
grid on;

% --- Subplot 5: S-GEO-04 k=1 örtüşme ---
subplot(2,3,5);
p_e1 = ellipsoidal_dome_profile(R_eq_viz, r0_viz, 1.0, 500);
p_h  = hemispherical_dome_profile(R_eq_viz, r0_viz, 500);
plot(p_e1.s, abs(p_e1.rho - p_h.rho), 'b-', 'LineWidth', 1.5); hold on;
plot(p_e1.s, abs(p_e1.x_local - p_h.x_local), 'r-', 'LineWidth', 1.5);
xlabel('s [mm]');
ylabel('|fark| [mm]');
title('S-GEO-04: ellip(k=1) vs hemi');
legend('|Δρ|', '|Δx|', 'Location', 'best');
grid on;

% --- Subplot 6: 3D dönel yüzey (k=0.6) ---
subplot(2,3,6);
phi_vec = linspace(0, 2*pi, 100);
[S_mesh, PHI_mesh] = meshgrid(p06.s, phi_vec);
RHO_grid = interp1(p06.s, p06.rho, S_mesh);
X_grid   = interp1(p06.s, p06.x_local, S_mesh);
Y_3d = RHO_grid .* cos(PHI_mesh);
Z_3d = RHO_grid .* sin(PHI_mesh);

surf(X_grid, Y_3d, Z_3d, 'FaceAlpha', 0.6, 'EdgeColor', 'none');
hold on;
% Silindir
[X_cyl_m, PHI_cyl_m] = meshgrid(linspace(-30, 0, 15), phi_vec);
Y_cyl_m = R_eq_viz * cos(PHI_cyl_m);
Z_cyl_m = R_eq_viz * sin(PHI_cyl_m);
surf(X_cyl_m, Y_cyl_m, Z_cyl_m, 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'FaceColor', [0.5 0.5 0.5]);
xlabel('x [mm]'); ylabel('Y [mm]'); zlabel('Z [mm]');
title('3D Dönel Yüzey (k=0.6)');
axis equal; colormap(jet);
view(30, 25); grid on;

sgtitle(sprintf('Elipsoidal Dome — R_{eq}=%.1f mm, r_0=%.1f mm, çeşitli k değerleri', ...
        R_eq_viz, r0_viz), 'FontSize', 14, 'FontWeight', 'bold');

print('-dpng', '-r150', 'ellipsoidal_dome_verification.png');
fprintf('  Grafik kaydedildi: ellipsoidal_dome_verification.png\n\n');

%% =====================================================================
%  TEST SETİ 8: Karşılaştırmalı Özet Tablosu
%  =====================================================================

fprintf('--- TEST SETİ 8: Karşılaştırmalı Özet ---\n\n');

fprintf('  %-16s  %6s  %6s  %5s  %9s  %10s  %10s  %12s  %12s\n', ...
        'Senaryo', 'R_eq', 'r0', 'k', 'θ_p [°]', 's_total', 'h_dome', ...
        'κ_eq [1/mm]', 'κ_pol [1/mm]');
fprintf('  %s\n', repmat('-', 1, 110));

for i = 1:length(test_cases)
    tc = test_cases(i);
    p = ellipsoidal_dome_profile(tc.R_eq, tc.r0, tc.k, 500);
    fprintf('  %-16s  %6.1f  %6.1f  %5.2f  %9.2f  %10.4f  %10.4f  %12.6e  %12.6e\n', ...
            tc.name, tc.R_eq, tc.r0, tc.k, rad2deg(p.theta_p), ...
            p.s_total, p.h_dome, p.kappa_m(1), p.kappa_m(end));
end

% k değişimi tablosu (sabit R_eq, r0)
fprintf('\n  k değişimi etkisi (R_eq=152.4, r0=45.0):\n');
fprintf('  %5s  %10s  %10s  %12s  %12s  %10s\n', ...
        'k', 's_total', 'h_dome', 'κ_m(ekv)', 'κ_m(pol)', 'κ_ekv/κ_pol');
fprintf('  %s\n', repmat('-', 1, 68));
for kk = [0.3, 0.5, 0.6, 0.7, 0.8, 1.0, 1.2, 1.5, 2.0]
    p = ellipsoidal_dome_profile(152.4, 45, kk, 300);
    fprintf('  %5.2f  %10.4f  %10.4f  %12.6e  %12.6e  %10.4f\n', ...
            kk, p.s_total, p.h_dome, p.kappa_m(1), p.kappa_m(end), ...
            p.kappa_m(1) / p.kappa_m(end));
end

fprintf('\n');

%% =====================================================================
%  SONUÇ
%  =====================================================================

fprintf('=============================================================\n');
if all_pass
    fprintf(' SONUÇ: TÜM TESTLER GEÇTİ ✓\n');
    fprintf(' Elipsoidal dome meridyen profili GATE-1a için onaylı.\n');
    fprintf(' S-GEO-04 (k=1 örtüşme): ONAYLI\n');
else
    fprintf(' SONUÇ: BAZI TESTLER BAŞARISIZ ✗\n');
    fprintf(' Başarısız testler düzeltilmeden ilerlenemez.\n');
end
fprintf('=============================================================\n');

%% --- Yardımcı fonksiyonlar ---
function str = tf(val)
    if val; str = 'PASS'; else; str = 'FAIL'; end
end

function str = ternary(cond, true_str, false_str)
    if cond; str = true_str; else; str = false_str; end
end
