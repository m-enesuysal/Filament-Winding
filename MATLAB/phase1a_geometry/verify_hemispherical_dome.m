%% VERIFY_HEMISPHERICAL_DOME.m
% Phase-1a: Hemispherical dome meridyen profili — doğrulama betiği
%
% Bu betik aşağıdaki GATE-1a koşullarını doğrular:
%   1. Analitik çözüm tutarlılığı (kapalı form doğrulaması)
%   2. C¹ süreklilik (S-GEO-02: ekvator noktası)
%   3. Birim teğet vektör normu
%   4. Eğrilik sabitliği (küre özelliği)
%   5. Sınır koşulları (ekvator ve polar açıklık)
%   6. Yüzey alanı (analitik vs sayısal integrasyon)
%   7. Dome yüksekliği (iki bağımsız formül)
%   8. S-GEO-04: Ellipsoidal k=1 örtüşme doğrulaması (placeholder)
%   9. Karar-16 test senaryoları üzerinde çalıştırma
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
fprintf(' PHASE-1a: HEMISPHERICAL DOME DOĞRULAMA RAPORU\n');
fprintf(' Tarih: %s\n', datestr(now, 'yyyy-mm-dd HH:MM'));
fprintf('=============================================================\n\n');

%% --- Karar-11 Toleransları ---
TOL_POS    = 1e-4;    % Pozisyon toleransı [mm]
TOL_DERIV  = 1e-6;    % Türev toleransı [-]
TOL_CURV   = 1e-4;    % Eğrilik bağıl toleransı [-]
TOL_MACHINE = 1e-12;  % Makine hassasiyeti (analitik kapalı form)

all_pass = true;

%% =====================================================================
%  TEST SETİ 1: Temel Analitik Doğrulama (Referans senaryo)
%  R_eq = 152.4 mm, r0 = 45 mm (TEST-02 benzeri)
%  =====================================================================

fprintf('--- TEST SETİ 1: Temel Analitik Doğrulama ---\n');
fprintf('    R_eq = 152.4 mm, r0 = 45.0 mm\n\n');

R_eq = 152.4;
r0   = 45.0;
N    = 1000;

profil = hemispherical_dome_profile(R_eq, r0, N);

% --- Test 1.1: θ_p analitik doğrulama ---
theta_p_expected = acos(r0 / R_eq);
err = abs(profil.theta_p - theta_p_expected);
pass = err < TOL_MACHINE;
fprintf('  [%s] T1.1 θ_p doğrulama:  hata = %.2e\n', tf(pass), err);
all_pass = all_pass && pass;

% --- Test 1.2: s_total analitik doğrulama ---
s_total_expected = R_eq * theta_p_expected;
err = abs(profil.s_total - s_total_expected);
pass = err < TOL_MACHINE;
fprintf('  [%s] T1.2 s_total doğrulama: hata = %.2e mm\n', tf(pass), err);
all_pass = all_pass && pass;

% --- Test 1.3: h_dome iki bağımsız formül ---
h_dome_1 = R_eq * sin(theta_p_expected);         % Trigonometrik
h_dome_2 = sqrt(R_eq^2 - r0^2);                   % Pisagor
err = abs(h_dome_1 - h_dome_2);
pass = err < TOL_MACHINE * R_eq;
fprintf('  [%s] T1.3 h_dome çapraz doğrulama: |trig - pisagor| = %.2e mm\n', ...
        tf(pass), err);
all_pass = all_pass && pass;

err = abs(profil.h_dome - h_dome_2);
pass = err < TOL_POS;
fprintf('  [%s] T1.3b h_dome profil çıktısı: hata = %.2e mm\n', tf(pass), err);
all_pass = all_pass && pass;

% --- Test 1.4: Sınır koşulları — ekvator (s=0) ---
err_rho_eq   = abs(profil.rho(1)     - R_eq);
err_x_eq     = abs(profil.x_local(1) - 0);
err_drho_eq  = abs(profil.drho_ds(1) - 0);
err_dx_eq    = abs(profil.dx_ds(1)   - 1);

pass_eq = (err_rho_eq < TOL_POS) && (err_x_eq < TOL_POS) && ...
          (err_drho_eq < TOL_DERIV) && (err_dx_eq < TOL_DERIV);
fprintf('  [%s] T1.4 Ekvator sınır koşulları:\n', tf(pass_eq));
fprintf('         ρ(0) - R_eq = %.2e mm\n', err_rho_eq);
fprintf('         x(0)       = %.2e mm\n', err_x_eq);
fprintf('         dρ/ds(0)   = %.2e\n', err_drho_eq);
fprintf('         dx/ds(0)   = %.2e\n', err_dx_eq);
all_pass = all_pass && pass_eq;

% --- Test 1.5: Sınır koşulları — polar açıklık (s=s_total) ---
err_rho_p = abs(profil.rho(end)     - r0);
err_x_p   = abs(profil.x_local(end) - profil.h_dome);

pass_p = (err_rho_p < TOL_POS) && (err_x_p < TOL_POS);
fprintf('  [%s] T1.5 Polar açıklık sınır koşulları:\n', tf(pass_p));
fprintf('         ρ(s_total) - r0     = %.2e mm\n', err_rho_p);
fprintf('         x(s_total) - h_dome = %.2e mm\n', err_x_p);
all_pass = all_pass && pass_p;

% --- Test 1.6: Birim teğet vektör normu ---
tangent_norm = sqrt(profil.drho_ds.^2 + profil.dx_ds.^2);
err_norm     = max(abs(tangent_norm - 1));
pass = err_norm < TOL_DERIV;
fprintf('  [%s] T1.6 Birim teğet vektör: max|norm-1| = %.2e\n', tf(pass), err_norm);
all_pass = all_pass && pass;

% --- Test 1.7: Eğrilik sabitliği ---
kappa_expected = 1 / R_eq;
err_kappa_max  = max(abs(profil.kappa_m - kappa_expected));
err_kappa_rel  = err_kappa_max / kappa_expected;
pass = err_kappa_rel < TOL_CURV;
fprintf('  [%s] T1.7 Eğrilik sabitliği: max|κ - 1/R_eq| = %.2e, bağıl = %.2e\n', ...
        tf(pass), err_kappa_max, err_kappa_rel);
all_pass = all_pass && pass;

% --- Test 1.8: ρ monotonik azalan ---
drho = diff(profil.rho);
pass = all(drho <= 0);
fprintf('  [%s] T1.8 ρ monotonik azalan (ekvator→polar)\n', tf(pass));
all_pass = all_pass && pass;

% --- Test 1.9: x_local monotonik artan ---
dx = diff(profil.x_local);
pass = all(dx >= 0);
fprintf('  [%s] T1.9 x_local monotonik artan (ekvator→polar)\n', tf(pass));
all_pass = all_pass && pass;

% --- Test 1.10: Yüzey alanı (sayısal integrasyon vs analitik) ---
A_numeric = 2 * pi * trapz(profil.s, profil.rho);   % Trapez kuralı
A_analytic = profil.A_dome;
err_A_rel = abs(A_numeric - A_analytic) / A_analytic;
pass = err_A_rel < 1e-4;   % Trapez yaklaşımı — daha gevşek tolerans
fprintf('  [%s] T1.10 Yüzey alanı doğrulama:\n', tf(pass));
fprintf('          Analitik  = %.4f mm²\n', A_analytic);
fprintf('          Sayısal   = %.4f mm²\n', A_numeric);
fprintf('          Bağıl hata = %.2e\n', err_A_rel);
all_pass = all_pass && pass;

% --- Test 1.11: Küre denklemi doğrulama ---
% ρ² + x_local² = R_eq² olmalı (tüm noktalarda)
sphere_check = profil.rho.^2 + profil.x_local.^2;
err_sphere = max(abs(sphere_check - R_eq^2));
pass = err_sphere < TOL_POS * R_eq;
fprintf('  [%s] T1.11 Küre denklemi: max|ρ²+x²-R²| = %.2e mm²\n', ...
        tf(pass), err_sphere);
all_pass = all_pass && pass;

fprintf('\n');

%% =====================================================================
%  TEST SETİ 2: C¹ Süreklilik Detaylı Analiz (S-GEO-02)
%  =====================================================================

fprintf('--- TEST SETİ 2: C¹ Süreklilik Analizi (S-GEO-02) ---\n\n');

% Silindir meridyen profili (s < 0 bölgesi, ρ = R_eq sabit)
% Dome meridyen profili (s ≥ 0 bölgesi)
% Geçiş noktası: s = 0

% Silindir tarafından yaklaşım (sol limit)
fprintf('  Silindir tarafı (s → 0⁻):\n');
fprintf('    ρ       = %.6f mm  (sabit = R_eq)\n', R_eq);
fprintf('    dρ/ds   = 0          (sabit yarıçap)\n');
fprintf('    dx/ds   = 1          (birim hız)\n');
fprintf('    κ_m     = 0          (düz çizgi)\n\n');

% Dome tarafından yaklaşım (sağ limit)
fprintf('  Dome tarafı (s → 0⁺):\n');
fprintf('    ρ(0)    = %.6f mm\n', profil.rho(1));
fprintf('    dρ/ds(0)= %.2e\n', profil.drho_ds(1));
fprintf('    dx/ds(0)= %.10f\n', profil.dx_ds(1));
fprintf('    κ_m(0)  = %.6e 1/mm\n\n', profil.kappa_m(1));

% C⁰ süreksizliği (yok — ρ eşleşiyor)
err_C0 = abs(profil.rho(1) - R_eq);
pass_C0 = err_C0 < TOL_MACHINE;
fprintf('  [%s] C⁰ süreklilik: |ρ_dome(0) - ρ_cyl| = %.2e mm\n', tf(pass_C0), err_C0);

% C¹ süreksizliği (yok — türevler eşleşiyor)
err_C1_rho = abs(profil.drho_ds(1) - 0);
err_C1_x   = abs(profil.dx_ds(1)   - 1);
pass_C1 = (err_C1_rho < TOL_DERIV) && (err_C1_x < TOL_DERIV);
fprintf('  [%s] C¹ süreklilik: |dρ/ds_dome(0) - 0| = %.2e\n', tf(pass_C1), err_C1_rho);
fprintf('  [%s] C¹ süreklilik: |dx/ds_dome(0) - 1| = %.2e\n', tf(pass_C1), err_C1_x);

% C² süreksizliği (beklenen — eğrilik sıçraması)
kappa_jump = abs(profil.kappa_m(1) - 0);
fprintf('\n  [BİLGİ] C² süreksizliği (beklenen):\n');
fprintf('    κ_m sıçraması = |1/R_eq - 0| = %.6e 1/mm\n', kappa_jump);
fprintf('    Bu süreksizlik kabul edilir (Karar-10: C¹ yeterli).\n');

all_pass = all_pass && pass_C0 && pass_C1;
fprintf('\n');

%% =====================================================================
%  TEST SETİ 3: Karar-16 Standart Test Senaryoları
%  =====================================================================

fprintf('--- TEST SETİ 3: Endüstri Test Senaryoları ---\n');
fprintf('    (Hemispherical dome tüm senaryolarda çalıştırılır)\n\n');

% Test senaryoları (Karar-16'dan)
test_cases = struct(...
    'name',  {'TEST-01 ASTM Subscale', 'TEST-02 Endüstriyel COPV', ...
              'TEST-03 Küçük Açıklık', 'TEST-04 H₂ Aerospace'}, ...
    'R_eq',  {73,    152.4, 150,   200}, ...
    'r0',    {22,    45,    10,    50});

for i = 1:length(test_cases)
    tc = test_cases(i);
    fprintf('  [%s] R_eq=%.1f, r0=%.1f\n', tc.name, tc.R_eq, tc.r0);

    p = hemispherical_dome_profile(tc.R_eq, tc.r0, 500);

    % Temel doğrulamalar
    ratio = tc.r0 / tc.R_eq;
    err_rho_eq = abs(p.rho(1) - tc.R_eq);
    err_rho_p  = abs(p.rho(end) - tc.r0);
    err_sphere = max(abs(p.rho.^2 + p.x_local.^2 - tc.R_eq^2));
    err_kappa  = max(abs(p.kappa_m - 1/tc.R_eq)) / (1/tc.R_eq);
    tangent_ok = max(abs(sqrt(p.drho_ds.^2 + p.dx_ds.^2) - 1));

    pass_i = (err_rho_eq < TOL_POS) && (err_rho_p < TOL_POS) && ...
             (err_sphere < TOL_POS * tc.R_eq) && (err_kappa < TOL_CURV) && ...
             (tangent_ok < TOL_DERIV);

    fprintf('    r0/R_eq = %.4f\n', ratio);
    fprintf('    θ_p     = %.4f rad (%.2f°)\n', p.theta_p, rad2deg(p.theta_p));
    fprintf('    s_total = %.4f mm\n', p.s_total);
    fprintf('    h_dome  = %.4f mm\n', p.h_dome);
    fprintf('    A_dome  = %.2f mm²\n', p.A_dome);
    fprintf('    Küre denklem hatası = %.2e mm²\n', err_sphere);
    fprintf('    Eğrilik bağıl hata  = %.2e\n', err_kappa);
    fprintf('    [%s] Tüm kontroller\n\n', tf(pass_i));

    all_pass = all_pass && pass_i;
end

%% =====================================================================
%  TEST SETİ 4: Uç Durum Senaryoları
%  =====================================================================

fprintf('--- TEST SETİ 4: Uç Durumlar ---\n\n');

% 4.1: r0 → R_eq yakın (sığ dome)
r0_shallow = 0.95 * R_eq;
p_shallow = hemispherical_dome_profile(R_eq, r0_shallow, 500);
pass = abs(p_shallow.rho(end) - r0_shallow) < TOL_POS;
fprintf('  [%s] T4.1 Sığ dome (r0/R_eq=0.95): θ_p = %.4f rad (%.2f°), h = %.2f mm\n', ...
        tf(pass), p_shallow.theta_p, rad2deg(p_shallow.theta_p), p_shallow.h_dome);
all_pass = all_pass && pass;

% 4.2: r0 çok küçük (derin dome)
r0_deep = 5;  % 5 mm
p_deep = hemispherical_dome_profile(R_eq, r0_deep, 500);
pass = abs(p_deep.rho(end) - r0_deep) < TOL_POS;
fprintf('  [%s] T4.2 Derin dome (r0=5mm): θ_p = %.4f rad (%.2f°), h = %.2f mm\n', ...
        tf(pass), p_deep.theta_p, rad2deg(p_deep.theta_p), p_deep.h_dome);
all_pass = all_pass && pass;

% 4.3: Hata fırlatma testleri
fprintf('\n  Hata fırlatma testleri:\n');

try hemispherical_dome_profile(-1, 10); fprintf('  [FAIL] Negatif R_eq kabul edildi!\n'); all_pass = false;
catch; fprintf('  [PASS] T4.3a Negatif R_eq reddedildi\n'); end

try hemispherical_dome_profile(100, 0); fprintf('  [FAIL] r0=0 kabul edildi!\n'); all_pass = false;
catch; fprintf('  [PASS] T4.3b r0=0 reddedildi\n'); end

try hemispherical_dome_profile(100, 100); fprintf('  [FAIL] r0=R_eq kabul edildi!\n'); all_pass = false;
catch; fprintf('  [PASS] T4.3c r0=R_eq reddedildi\n'); end

try hemispherical_dome_profile(100, 150); fprintf('  [FAIL] r0>R_eq kabul edildi!\n'); all_pass = false;
catch; fprintf('  [PASS] T4.3d r0>R_eq reddedildi\n'); end

fprintf('\n');

%% =====================================================================
%  TEST SETİ 5: Sayısal Türev Karşılaştırması
%  =====================================================================

fprintf('--- TEST SETİ 5: Sayısal vs Analitik Türev ---\n');
fprintf('    (Merkezi fark yöntemi ile çapraz doğrulama)\n\n');

% Yüksek çözünürlüklü profil
p = hemispherical_dome_profile(R_eq, r0, 2000);
ds = p.s(2) - p.s(1);

% Merkezi fark ile sayısal türev (iç noktalar)
drho_ds_num = zeros(size(p.rho));
dx_ds_num   = zeros(size(p.x_local));

% Merkezi fark (2. mertebe doğru)
drho_ds_num(2:end-1) = (p.rho(3:end) - p.rho(1:end-2)) / (2*ds);
dx_ds_num(2:end-1)   = (p.x_local(3:end) - p.x_local(1:end-2)) / (2*ds);

% Tek taraflı fark (sınırlar)
drho_ds_num(1)   = (p.rho(2) - p.rho(1)) / ds;
drho_ds_num(end) = (p.rho(end) - p.rho(end-1)) / ds;
dx_ds_num(1)     = (p.x_local(2) - p.x_local(1)) / ds;
dx_ds_num(end)   = (p.x_local(end) - p.x_local(end-1)) / ds;

% Hata analizi (sadece iç noktalar — sınırlar tek taraflı)
idx = 2:length(p.s)-1;
err_drho = max(abs(drho_ds_num(idx) - p.drho_ds(idx)));
err_dx   = max(abs(dx_ds_num(idx)   - p.dx_ds(idx)));

% Merkezi fark hatası O(ds²) olmalı
fprintf('  ds = %.4e mm\n', ds);
fprintf('  Beklenen merkezi fark hatası: O(ds²) ≈ %.2e\n', ds^2);
fprintf('  max|dρ/ds_sayısal - dρ/ds_analitik| = %.2e\n', err_drho);
fprintf('  max|dx/ds_sayısal - dx/ds_analitik| = %.2e\n', err_dx);

pass = (err_drho < 1e-4) && (err_dx < 1e-4);  % Gevşek tolerans (sayısal fark)
fprintf('  [%s] Sayısal-analitik türev uyumu\n\n', tf(pass));
all_pass = all_pass && pass;

%% =====================================================================
%  GÖRSELLEŞTİRME
%  =====================================================================

fprintf('--- Görselleştirme oluşturuluyor ---\n\n');

% Referans senaryo kullan
p = hemispherical_dome_profile(152.4, 45, 500);

figure('Name', 'Phase-1a: Hemispherical Dome Meridyen Profili', ...
       'Position', [100, 100, 1400, 900]);

% --- Subplot 1: Meridyen profili (ρ, x) düzlemi ---
subplot(2,3,1);
plot(p.x_local, p.rho, 'b-', 'LineWidth', 2);
hold on;
% Silindir bölgesi (referans)
x_cyl = linspace(-50, 0, 50);
rho_cyl = R_eq * ones(size(x_cyl));
plot(x_cyl, rho_cyl, 'k--', 'LineWidth', 1.5);
% İşaret noktaları
plot(0, R_eq, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');  % Ekvator
plot(p.x_local(end), p.rho(end), 'gs', 'MarkerSize', 8, 'MarkerFaceColor', 'g');  % Polar
xlabel('x_{lokal} [mm]');
ylabel('\rho [mm]');
title('Meridyen Profili');
legend('Dome', 'Silindir', 'Ekvator', 'Polar açıklık', 'Location', 'southwest');
grid on; axis equal;

% --- Subplot 2: ρ(s) ve x(s) ---
subplot(2,3,2);
yyaxis left;
plot(p.s, p.rho, 'b-', 'LineWidth', 1.5);
ylabel('\rho [mm]');
yyaxis right;
plot(p.s, p.x_local, 'r-', 'LineWidth', 1.5);
ylabel('x_{lokal} [mm]');
xlabel('s [mm]');
title('\rho(s) ve x(s)');
legend('\rho', 'x', 'Location', 'east');
grid on;

% --- Subplot 3: Türevler ---
subplot(2,3,3);
plot(p.s, p.drho_ds, 'b-', 'LineWidth', 1.5); hold on;
plot(p.s, p.dx_ds, 'r-', 'LineWidth', 1.5);
yline(0, 'k--', 'LineWidth', 0.5);
xlabel('s [mm]');
ylabel('Türev [-]');
title('d\rho/ds ve dx/ds');
legend('d\rho/ds', 'dx/ds', 'Location', 'east');
grid on;

% --- Subplot 4: Eğim açısı β ---
subplot(2,3,4);
plot(p.s, rad2deg(p.beta), 'b-', 'LineWidth', 1.5);
xlabel('s [mm]');
ylabel('\beta [°]');
title('Meridyen Eğim Açısı');
grid on;

% --- Subplot 5: Eğrilik ---
subplot(2,3,5);
plot(p.s, p.kappa_m, 'b-', 'LineWidth', 2);
yline(1/R_eq, 'r--', '1/R_{eq}', 'LineWidth', 1, 'LabelHorizontalAlignment', 'left');
xlabel('s [mm]');
ylabel('\kappa_m [1/mm]');
title('Meridyen Eğriliği (sabit)');
grid on;
ylim([0, 2/R_eq]);

% --- Subplot 6: Dönel yüzey (3D) ---
subplot(2,3,6);
phi_vec = linspace(0, 2*pi, 100);
[S, PHI] = meshgrid(p.s, phi_vec);
RHO_grid = interp1(p.s, p.rho, S);
X_grid   = interp1(p.s, p.x_local, S);

Y_3d = RHO_grid .* cos(PHI);
Z_3d = RHO_grid .* sin(PHI);

surf(X_grid, Y_3d, Z_3d, 'FaceAlpha', 0.6, 'EdgeColor', 'none');
hold on;
% Silindir bölgesi
[X_cyl, PHI_cyl] = meshgrid(linspace(-50, 0, 20), phi_vec);
Y_cyl = R_eq * cos(PHI_cyl);
Z_cyl = R_eq * sin(PHI_cyl);
surf(X_cyl, Y_cyl, Z_cyl, 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'FaceColor', [0.5 0.5 0.5]);
xlabel('x [mm]'); ylabel('Y [mm]'); zlabel('Z [mm]');
title('3D Dönel Yüzey');
axis equal; colormap(jet);
view(30, 25);
grid on;

sgtitle(sprintf('Hemispherical Dome — R_{eq}=%.1f mm, r_0=%.1f mm', R_eq, r0), ...
        'FontSize', 14, 'FontWeight', 'bold');

% Grafik kaydet
print('-dpng', '-r150', 'hemispherical_dome_verification.png');
fprintf('  Grafik kaydedildi: hemispherical_dome_verification.png\n\n');

%% =====================================================================
%  TEST SETİ 6: Tüm Test Senaryoları — Karşılaştırmalı Tablo
%  =====================================================================

fprintf('--- TEST SETİ 6: Karşılaştırmalı Özet Tablosu ---\n\n');

fprintf('  %-20s  %8s  %8s  %8s  %10s  %10s  %12s\n', ...
        'Senaryo', 'R_eq', 'r0', 'r0/R_eq', 'θ_p [°]', 's_total', 'h_dome');
fprintf('  %s\n', repmat('-', 1, 82));

for i = 1:length(test_cases)
    tc = test_cases(i);
    p = hemispherical_dome_profile(tc.R_eq, tc.r0, 500);
    fprintf('  %-20s  %8.1f  %8.1f  %8.4f  %10.2f  %10.4f  %12.4f\n', ...
            tc.name, tc.R_eq, tc.r0, tc.r0/tc.R_eq, ...
            rad2deg(p.theta_p), p.s_total, p.h_dome);
end

fprintf('\n');

%% =====================================================================
%  SONUÇ
%  =====================================================================

fprintf('=============================================================\n');
if all_pass
    fprintf(' SONUÇ: TÜM TESTLER GEÇTİ ✓\n');
    fprintf(' Hemispherical dome meridyen profili GATE-1a için onaylı.\n');
else
    fprintf(' SONUÇ: BAZI TESTLER BAŞARISIZ ✗\n');
    fprintf(' Başarısız testler düzeltilmeden ilerlenemez.\n');
end
fprintf('=============================================================\n');

%% --- Yardımcı fonksiyon ---
function str = tf(val)
    if val
        str = 'PASS';
    else
        str = 'FAIL';
    end
end
