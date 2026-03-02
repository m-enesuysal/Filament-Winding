%% VERIFY_ISOTENSOID_DOME  İzotensoid dome doğrulama test paketi (v2).
%
% Phase-1a: İzotensoid dome profil fonksiyonu doğrulaması.
% 8 test seti, GATE-1a-01 çapraz doğrulama dahil.
%
% v2 DÜZELTME NOTU (2026-02-27):
%   Test Seti 3: ODE düzeltildi — doğru diskriminant yapısı
%     dZ/dY = Y·(1+2q-Y²) / √[(2+2q-Y²)·(Y²-1)·(q+1-Y²)]
%     1/√(Y²-1) → integre edilebilir karekök singülarite
%   Test Seti 4: Eğrilik karşılaştırma — parametrik formülden sayısal fark
%   Test Seti 6: Aspect ratio referans düzeltildi
%     h/R_eq → √2·E(1/2) ≈ 1.9101 (Koussios Fig.4.3 → N_hoop/N_polar)
%
% BAĞIMLILIKLAR:
%   - isotensoid_dome_profile.m (v2, aynı klasörde)
%   - hemispherical_dome_profile.m (Test Seti 8 karşılaştırma)
%   - ellipsoidal_dome_profile.m (Test Seti 8 karşılaştırma)
%
% Tarih: 2026-02-27 (v2)
% Faz: Phase-1a

clear; clc; close all;

fprintf('=============================================================\n');
fprintf(' PHASE-1a: İZOTENSOİD DOME DOĞRULAMA RAPORU (v2)\n');
fprintf(' Tarih: %s\n', datestr(now, 'yyyy-mm-dd HH:MM'));
fprintf('=============================================================\n\n');

% Toleranslar (Karar-11)
TOL_MACHINE  = 1e-12;    % Makine hassasiyeti (analitik hesaplamalar)
TOL_POS      = 1e-4;     % Pozisyon hatası [mm]
TOL_DERIV    = 1e-6;     % Türev hatası [-]
TOL_CURV_REL = 5e-3;     % Eğrilik bağıl hata (sayısal fark sınırı)
TOL_ODE_REL  = 1e-4;     % ODE çapraz doğrulama bağıl toleransı
TOL_AREA_REL = 1e-6;     % Yüzey alanı bağıl hatası

total_pass = 0;
total_fail = 0;

%% =====================================================================
%  TEST SETİ 1: Temel Analitik Doğrulama (TEST-01: R_eq=73, r0=22)
%  =====================================================================
fprintf('--- TEST SETİ 1: Temel Analitik Doğrulama (TEST-01) ---\n');
fprintf('    R_eq = 73 mm, r0 = 22 mm\n\n');

R_eq_t1 = 73;  r0_t1 = 22;
profil1 = isotensoid_dome_profile(R_eq_t1, r0_t1, 1000);

% T1.1: q, m parametreleri
Y_eq_t1   = R_eq_t1 / r0_t1;
q_ref     = Y_eq_t1^2 - 1;
m_ref     = q_ref / (1 + 2*q_ref);
err_q     = abs(profil1.q - q_ref);
err_m     = abs(profil1.m_ell - m_ref);
if err_q < TOL_MACHINE && err_m < TOL_MACHINE
    fprintf('  [PASS] T1.1 q, m doğrulama:  |Δq| = %.2e, |Δm| = %.2e\n', err_q, err_m);
    fprintf('         Y_eq = %f, q = %f, m = %f\n', Y_eq_t1, q_ref, m_ref);
    total_pass = total_pass + 1;
else
    fprintf('  [FAIL] T1.1 q, m doğrulama:  |Δq| = %.2e, |Δm| = %.2e\n', err_q, err_m);
    total_fail = total_fail + 1;
end

% T1.2: h_dome (eliptik integral kapalı form)
sqrt_1p2q = sqrt(1 + 2*q_ref);
K_ref = ellipticF(pi/2, m_ref);
E_ref = ellipticE(pi/2, m_ref);
h_dome_ref = r0_t1 * sqrt_1p2q * (E_ref - K_ref / (1 + 2*q_ref));
err_h = abs(profil1.h_dome - h_dome_ref);
if err_h < TOL_MACHINE
    fprintf('  [PASS] T1.2 h_dome:  profil = %f mm, analitik = %f mm, |Δ| = %.2e mm\n', ...
            profil1.h_dome, h_dome_ref, err_h);
    total_pass = total_pass + 1;
else
    fprintf('  [FAIL] T1.2 h_dome:  |Δ| = %.2e mm\n', err_h);
    total_fail = total_fail + 1;
end

% T1.3: Ekvator sınır koşulları
err_rho_eq  = abs(profil1.rho(1) - R_eq_t1);
err_x_eq    = abs(profil1.x_local(1));
err_drho_eq = abs(profil1.drho_ds(1));
err_dx_eq   = abs(profil1.dx_ds(1) - 1);
if err_rho_eq < TOL_POS && err_x_eq < TOL_POS && ...
   err_drho_eq < TOL_DERIV && err_dx_eq < TOL_DERIV
    fprintf('  [PASS] T1.3 Ekvator sınır koşulları:\n');
    fprintf('         ρ(0) - R_eq  = %.2e mm\n', err_rho_eq);
    fprintf('         x(0)         = %.2e mm\n', err_x_eq);
    fprintf('         dρ/ds(0)     = %.2e\n', err_drho_eq);
    fprintf('         dx/ds(0) - 1 = %.2e\n', err_dx_eq);
    total_pass = total_pass + 1;
else
    fprintf('  [FAIL] T1.3 Ekvator sınır koşulları\n');
    total_fail = total_fail + 1;
end

% T1.4: Polar açıklık sınır koşulları
err_rho_pol = abs(profil1.rho(end) - r0_t1);
err_x_pol   = abs(profil1.x_local(end) - profil1.h_dome);
if err_rho_pol < TOL_POS && err_x_pol < TOL_POS
    fprintf('  [PASS] T1.4 Polar açıklık sınır koşulları:\n');
    fprintf('         ρ(s_total) - r0     = %.2e mm\n', err_rho_pol);
    fprintf('         x(s_total) - h_dome = %.2e mm\n', err_x_pol);
    total_pass = total_pass + 1;
else
    fprintf('  [FAIL] T1.4 Polar açıklık sınır koşulları\n');
    total_fail = total_fail + 1;
end

% T1.5: Birim teğet vektör normu (iç noktalar)
tangent_norm = sqrt(profil1.drho_ds.^2 + profil1.dx_ds.^2);
idx_inner    = 2:(length(profil1.s)-1);
max_norm_err = max(abs(tangent_norm(idx_inner) - 1));
if max_norm_err < TOL_DERIV
    fprintf('  [PASS] T1.5 Birim teğet vektör: max|norm-1| = %.2e (iç noktalar)\n', max_norm_err);
    total_pass = total_pass + 1;
else
    fprintf('  [FAIL] T1.5 Birim teğet vektör: max|norm-1| = %.2e\n', max_norm_err);
    total_fail = total_fail + 1;
end

% T1.6: ρ monotonik azalan (ekvator → polar)
drho = diff(profil1.rho);
if all(drho <= 1e-10)
    fprintf('  [PASS] T1.6 ρ monotonik azalan (ekvator→polar)\n');
    total_pass = total_pass + 1;
else
    fprintf('  [FAIL] T1.6 ρ monotonik BOZULDU\n');
    total_fail = total_fail + 1;
end

% T1.7: x_local monotonik artan (ekvator → polar)
dx = diff(profil1.x_local);
if all(dx >= -1e-10)
    fprintf('  [PASS] T1.7 x_local monotonik artan (ekvator→polar)\n');
    total_pass = total_pass + 1;
else
    fprintf('  [FAIL] T1.7 x_local monotonik BOZULDU\n');
    total_fail = total_fail + 1;
end

% T1.8: Ekvator eğriliği (analitik referans)
kappa_eq_t1  = (1 + q_ref) / (r0_t1 * q_ref * Y_eq_t1);
err_kappa_eq = abs(profil1.kappa_eq - kappa_eq_t1) / abs(kappa_eq_t1);
if err_kappa_eq < TOL_MACHINE
    fprintf('  [PASS] T1.8 Ekvator eğriliği: κ_m(0) = %.6e, beklenen = %.6e, bağıl hata = %.2e\n', ...
            profil1.kappa_eq, kappa_eq_t1, err_kappa_eq);
    total_pass = total_pass + 1;
else
    fprintf('  [FAIL] T1.8 Ekvator eğriliği bağıl hata = %.2e\n', err_kappa_eq);
    total_fail = total_fail + 1;
end

% T1.9: Yüzey alanı (integral vs trapz çapraz doğrulama)
A_trapz = 2*pi * trapz(profil1.s, profil1.rho);
err_A   = abs(A_trapz - profil1.A_dome) / abs(profil1.A_dome);
if err_A < TOL_AREA_REL
    fprintf('  [PASS] T1.9 Yüzey alanı:\n');
    fprintf('          integral()  = %.4f mm²\n', profil1.A_dome);
    fprintf('          trapz(s,ρ)  = %.4f mm²\n', A_trapz);
    fprintf('          Bağıl hata  = %.2e\n', err_A);
    total_pass = total_pass + 1;
else
    fprintf('  [FAIL] T1.9 Yüzey alanı bağıl hatası = %.2e\n', err_A);
    total_fail = total_fail + 1;
end

% T1.10: Clairaut doğrulama
clairaut_err = max(abs(sin(profil1.alpha_w) - r0_t1 ./ profil1.rho));
if clairaut_err < TOL_DERIV
    fprintf('  [PASS] T1.10 Clairaut doğrulama: max|sin(α) - r0/ρ| = %.2e\n', clairaut_err);
    total_pass = total_pass + 1;
else
    fprintf('  [FAIL] T1.10 Clairaut hatası = %.2e\n', clairaut_err);
    total_fail = total_fail + 1;
end

%% =====================================================================
%  TEST SETİ 2: C¹ Süreklilik Analizi (S-GEO-02)
%  =====================================================================
fprintf('\n--- TEST SETİ 2: C¹ Süreklilik Analizi (S-GEO-02) ---\n\n');

fprintf('  Silindir tarafı (s → 0⁻):\n');
fprintf('    ρ       = %f mm  (sabit = R_eq)\n', R_eq_t1);
fprintf('    dρ/ds   = 0\n');
fprintf('    dx/ds   = 1\n');
fprintf('    κ_m     = 0\n\n');

fprintf('  İzotensoid dome tarafı (s → 0⁺):\n');
fprintf('    ρ(0)    = %f mm\n', profil1.rho(1));
fprintf('    dρ/ds(0)= %.2e\n', profil1.drho_ds(1));
fprintf('    dx/ds(0)= %.10f\n', profil1.dx_ds(1));
fprintf('    κ_m(0)  = %.6e 1/mm\n\n', profil1.kappa_m(1));

% C⁰
err_C0 = abs(profil1.rho(1) - R_eq_t1);
if err_C0 < TOL_POS
    fprintf('  [PASS] C⁰ süreklilik: |ρ_dome(0) - R_eq| = %.2e mm\n', err_C0);
    total_pass = total_pass + 1;
else
    fprintf('  [FAIL] C⁰ süreklilik\n');
    total_fail = total_fail + 1;
end

% C¹
err_C1_drho = abs(profil1.drho_ds(1));
err_C1_dx   = abs(profil1.dx_ds(1) - 1);
if err_C1_drho < TOL_DERIV && err_C1_dx < TOL_DERIV
    fprintf('  [PASS] C¹ süreklilik: |dρ/ds(0)| = %.2e, |dx/ds(0)-1| = %.2e\n', ...
            err_C1_drho, err_C1_dx);
    total_pass = total_pass + 1;
else
    fprintf('  [FAIL] C¹ süreklilik\n');
    total_fail = total_fail + 1;
end

% C² bilgi
kappa_eq_hemi = 1/R_eq_t1;
fprintf('\n  [BİLGİ] C² süreksizliği (beklenen):\n');
fprintf('    κ_m sıçraması = %.6e 1/mm (silindir κ=0 → dome κ=%.4e)\n', ...
        profil1.kappa_m(1), profil1.kappa_m(1));
fprintf('    Hemispherical referans (aynı R_eq): 1/R_eq = %.6e 1/mm\n', kappa_eq_hemi);
fprintf('    Oran: κ_eq(isotensoid) / κ_eq(hemispherical) = %.4f\n', ...
        profil1.kappa_m(1) / kappa_eq_hemi);

%% =====================================================================
%  TEST SETİ 3: Eliptik İntegral vs ODE Çapraz Doğrulama (GATE-1a-01)
%  =====================================================================
fprintf('\n--- TEST SETİ 3: Eliptik İntegral vs ODE Çapraz Doğrulama ---\n');
fprintf('    KRİTİK TEST — GATE-1a-01 gereksinimi\n');
fprintf('    v2: DOĞRU ODE — dZ/dY = Y(1+2q-Y²)/√[(2+2q-Y²)(Y²-1)(q+1-Y²)]\n\n');

R_eq_t3 = 152.4;  r0_t3 = 45;
profil3 = isotensoid_dome_profile(R_eq_t3, r0_t3, 2000);

fprintf('  Eliptik integral çözümü (birincil):\n');
fprintf('    s_total = %.8f mm\n', profil3.s_total);
fprintf('    h_dome  = %.8f mm\n', profil3.h_dome);
fprintf('    A_dome  = %.4f mm²\n\n', profil3.A_dome);

% ODE çözümü (bağımsız doğrulama)
q3    = profil3.q;
Y_eq3 = R_eq_t3 / r0_t3;

% Doğru ODE: dZ/dY = Y·(1+2q-Y²) / √[(2+2q-Y²)·(Y²-1)·(q+1-Y²)]
% Alan: Y ∈ [1, Y_eq], paydaki üç terim hep > 0
% Singülariteler: Y=1'de √(Y²-1)→0 (integre edilebilir),
%                 Y=Y_eq'de √(q+1-Y²)→0 (integre edilebilir)
ode_func = @(Y, Z) Y .* (1 + 2*q3 - Y.^2) ./ ...
    sqrt( max((2 + 2*q3 - Y.^2) .* (Y.^2 - 1) .* (q3 + 1 - Y.^2), 1e-60) );

% Başlangıç noktası: Y=1+ε, Z ≈ 2√(2qε/(1+2q)) (lokal seri açılımı)
eps_start = 1e-8;
Y_start   = 1 + eps_start;
Z_start   = 2 * sqrt(2 * q3 * eps_start / (1 + 2*q3));

% Bitiş noktası: Y=Y_eq-ε (ekvator singülaritesinden kaçın)
eps_end = 1e-8;
Y_end   = Y_eq3 - eps_end;

% ODE çözümü
opts_ode = odeset('RelTol', 1e-10, 'AbsTol', 1e-12, 'MaxStep', 0.01);
[Y_ode, Z_ode] = ode45(ode_func, [Y_start, Y_end], Z_start, opts_ode);

fprintf('  ODE çözümü (ikincil, ode45):\n');
fprintf('    Y aralığı: [%.10f, %.10f]\n', Y_ode(1), Y_ode(end));
fprintf('    Z(Y_end)  = %.8f (boyutsuz)\n', Z_ode(end));
fprintf('    N_adım    = %d\n', length(Y_ode));

% h_dome karşılaştırması: Z(Y_eq)·r0
% ODE, Y_eq-ε'da bitiyor → ekvator yakını lokal açılım düzeltmesi
%   Z(Y_eq) - Z(Y_eq-ε) ≈ √(Y_eq·q3/(2(1+q3))) · 2√ε
Z_eq_correction = 2 * sqrt(Y_eq3 * q3 * eps_end / (2*(1+q3)));
h_dome_ode = r0_t3 * (Z_ode(end) + Z_eq_correction);

fprintf('    h_dome(ODE) = %.8f mm (singülarite düzeltmeli)\n\n', h_dome_ode);

% Skaler karşılaştırma
err_h_ode     = abs(h_dome_ode - profil3.h_dome);
err_h_ode_rel = err_h_ode / profil3.h_dome;

if err_h_ode_rel < TOL_ODE_REL
    fprintf('  Skaler karşılaştırma:\n');
    fprintf('  [PASS] h_dome: |Δ| = %.2e mm, bağıl = %.2e\n', err_h_ode, err_h_ode_rel);
    total_pass = total_pass + 1;
else
    fprintf('  Skaler karşılaştırma:\n');
    fprintf('  [FAIL] h_dome: |Δ| = %.2e mm, bağıl = %.2e\n', err_h_ode, err_h_ode_rel);
    total_fail = total_fail + 1;
end

% Vektörel karşılaştırma: ODE çözümünü eliptik integral ile karşılaştır
% Her ODE Y noktasında, eliptik integral ile bağımsız Z hesapla
m3 = profil3.m_ell;
sqrt_1p2q3 = sqrt(1 + 2*q3);

% θ'dan Y'ye: Y = √(1+q·sin²θ) → sinθ = √((Y²-1)/q), θ = arcsin(√((Y²-1)/q))
% Güvenli aralık seç: Y ∈ [1.01·r0/r0, 0.99·Y_eq]
Y_safe_lo = 1 + 0.01*(Y_eq3 - 1);
Y_safe_hi = 1 + 0.99*(Y_eq3 - 1);
idx_safe  = (Y_ode > Y_safe_lo) & (Y_ode < Y_safe_hi);
Y_safe    = Y_ode(idx_safe);
Z_ode_safe = Z_ode(idx_safe);

Z_ell_safe = zeros(size(Y_safe));
for i = 1:length(Y_safe)
    sin2th = (Y_safe(i)^2 - 1) / q3;
    sin2th = min(max(sin2th, 0), 1);
    th_i   = asin(sqrt(sin2th));
    F_i = ellipticF(th_i, m3);
    E_i = ellipticE(th_i, m3);
    Z_ell_safe(i) = sqrt_1p2q3 * (E_i - F_i / (1 + 2*q3));
end

x_ode_safe = r0_t3 * (profil3.h_dome/r0_t3 - Z_ode_safe);
x_ell_safe = r0_t3 * (profil3.h_dome/r0_t3 - Z_ell_safe);

max_dx = max(abs(x_ode_safe - x_ell_safe));
mean_dx = mean(abs(x_ode_safe - x_ell_safe));

fprintf('  Vektörel karşılaştırma (güvenli iç noktalar, N=%d):\n', sum(idx_safe));
if max_dx < TOL_POS
    fprintf('  [PASS] max|Δx_local| = %.2e mm\n', max_dx);
    fprintf('         ort|Δx_local| = %.2e mm\n', mean_dx);
    total_pass = total_pass + 1;
else
    fprintf('  [FAIL] max|Δx_local| = %.2e mm\n', max_dx);
    fprintf('         ort|Δx_local| = %.2e mm\n', mean_dx);
    total_fail = total_fail + 1;
end

% GATE-1a-01 sonucu
gate_pass = (err_h_ode_rel < TOL_ODE_REL) && (max_dx < TOL_POS);
if gate_pass
    fprintf('\n  ► GATE-1a-01 DURUMU: ONAYLI ✓\n');
else
    fprintf('\n  ► GATE-1a-01 DURUMU: BAŞARISIZ ✗\n');
end

%% =====================================================================
%  TEST SETİ 4: Eğrilik — Sayısal vs Analitik Doğrulama (v2)
%  =====================================================================
fprintf('\n--- TEST SETİ 4: Eğrilik — Sayısal (dβ/ds) vs Profil Doğrulama ---\n\n');

% Yüksek çözünürlüklü profil
profil4 = isotensoid_dome_profile(R_eq_t3, r0_t3, 4000);
ds4     = profil4.s(2) - profil4.s(1);
N4      = length(profil4.s);

% Sayısal eğrilik: κ = dβ/ds (merkezi fark)
kappa_num = zeros(N4, 1);
for i = 2:N4-1
    kappa_num(i) = (profil4.beta(i+1) - profil4.beta(i-1)) / (2*ds4);
end

% Analiz aralığı: s/s_total ∈ [0.02, 0.90] (her iki sınır hariç)
s_norm      = profil4.s / profil4.s_total;
idx_analyze = (s_norm > 0.02) & (s_norm < 0.90);

kappa_prof = profil4.kappa_m(idx_analyze);
kappa_n    = kappa_num(idx_analyze);

% Bağıl hata (|κ| > küçük eşik olan noktalarda)
kappa_abs_min = 1e-6;
idx_valid     = abs(kappa_prof) > kappa_abs_min;
if sum(idx_valid) > 0
    rel_err = abs(kappa_n(idx_valid) - kappa_prof(idx_valid)) ./ abs(kappa_prof(idx_valid));
    max_rel = max(rel_err);
    mean_rel = mean(rel_err);
else
    max_rel = 0;
    mean_rel = 0;
end

max_abs_err = max(abs(kappa_n - kappa_prof));

fprintf('  ds = %.4e mm, analiz aralığı: s/s_total ∈ [0.02, 0.90]\n', ds4);
fprintf('  max|κ_sayısal - κ_profil|   = %.2e 1/mm\n', max_abs_err);
fprintf('  max bağıl hata              = %.2e\n', max_rel);
fprintf('  ortalama bağıl hata         = %.2e\n', mean_rel);

if max_rel < TOL_CURV_REL
    fprintf('  [PASS] Sayısal-analitik eğrilik uyumu\n');
    total_pass = total_pass + 1;
else
    fprintf('  [FAIL] Sayısal-analitik eğrilik uyumu\n');
    total_fail = total_fail + 1;
end

% İşaret değişimi (bükülme noktası) kontrolü
kappa_sign_change = any(kappa_prof > 0) && any(kappa_prof < 0);
kappa_sign_change_num = any(kappa_n > 0) && any(kappa_n < 0);
fprintf('\n  [BİLGİ] κ_m işaret değişimi (bükülme noktası):\n');
fprintf('    Profil fonksiyonu: %s\n', mat2str(kappa_sign_change));
fprintf('    Sayısal (dβ/ds):  %s\n', mat2str(kappa_sign_change_num));
fprintf('    κ_eq = %.6e (pozitif — konveks)\n', profil4.kappa_eq);
fprintf('    κ_pol = %.6e (negatif — konkav)\n', profil4.kappa_pol);

%% =====================================================================
%  TEST SETİ 5: Karar-16 Standart Test Senaryoları
%  =====================================================================
fprintf('\n--- TEST SETİ 5: Endüstri Test Senaryoları ---\n\n');

test_cases = {
    'TEST-01 ASTM',    73.0,   22.0;
    'TEST-02 COPV',   152.4,   45.0;
    'TEST-03 Küçük Aç.', 150.0, 10.0;
    'TEST-04 H₂',     200.0,   50.0;
};

for tc = 1:size(test_cases, 1)
    name   = test_cases{tc, 1};
    R_eq_i = test_cases{tc, 2};
    r0_i   = test_cases{tc, 3};

    p = isotensoid_dome_profile(R_eq_i, r0_i);
    Y_eq_i = R_eq_i / r0_i;
    q_i    = Y_eq_i^2 - 1;

    % Ekvator eğriliği beklenen
    kappa_eq_exp = (1 + q_i) / (r0_i * q_i * Y_eq_i);

    % Clairaut doğrulama
    cl_err = max(abs(sin(p.alpha_w) - r0_i ./ p.rho));

    % Sınır koşulları
    bc_ok = abs(p.rho(1) - R_eq_i) < TOL_POS && ...
            abs(p.rho(end) - r0_i) < TOL_POS && ...
            abs(p.x_local(1)) < TOL_POS && ...
            abs(p.x_local(end) - p.h_dome) < TOL_POS;

    all_ok = bc_ok && (cl_err < TOL_DERIV) && ...
             (abs(p.kappa_eq - kappa_eq_exp)/abs(kappa_eq_exp) < TOL_MACHINE);

    fprintf('  [%s] R_eq=%.1f, r0=%.1f\n', name, R_eq_i, r0_i);
    fprintf('    r0/R_eq     = %.4f\n', r0_i/R_eq_i);
    fprintf('    q           = %.4f\n', p.q);
    fprintf('    m           = %.6f\n', p.m_ell);
    fprintf('    s_total     = %.4f mm\n', p.s_total);
    fprintf('    h_dome      = %.4f mm\n', p.h_dome);
    fprintf('    aspect_r    = %.4f (h/R_eq)\n', p.aspect_r);
    fprintf('    κ_m(ekv)    = %.6e [beklenen: %.6e]\n', p.kappa_eq, kappa_eq_exp);
    fprintf('    Clairaut    = %.2e (max hata)\n', cl_err);

    if all_ok
        fprintf('    [PASS] Tüm kontroller\n\n');
        total_pass = total_pass + 1;
    else
        fprintf('    [FAIL] Bir veya daha fazla kontrol başarısız\n\n');
        total_fail = total_fail + 1;
    end
end

%% =====================================================================
%  TEST SETİ 6: S-GEO-01, S-GEO-03 ve Uç Durumlar
%  =====================================================================
fprintf('--- TEST SETİ 6: S-GEO-01, S-GEO-03 ve Uç Durumlar ---\n\n');

% T6.1: Büyük Y_eq (stiff ODE testi)
R_eq_61 = 150;  r0_61 = 10;  % R_eq/r0 = 15
p61 = isotensoid_dome_profile(R_eq_61, r0_61);
fprintf('  T6.1: Küçük polar açıklık (R_eq/r0 = 15, S-GEO-01 stiff ODE)\n');
fprintf('    q = %.2f, aspect_r = %.4f\n', p61.q, p61.aspect_r);
fprintf('    Polar ρ hatası = %.2e mm\n', abs(p61.rho(end) - r0_61));
if abs(p61.rho(end) - r0_61) < TOL_POS
    fprintf('    [PASS] S-GEO-01 (stiff ODE) profil oluşturuldu\n\n');
    total_pass = total_pass + 1;
else
    fprintf('    [FAIL] S-GEO-01\n\n');
    total_fail = total_fail + 1;
end

% T6.2: Küçük Y_eq (sığ dome)
R_eq_62 = 150;  r0_62 = 100;  % R_eq/r0 = 1.5
p62 = isotensoid_dome_profile(R_eq_62, r0_62);
fprintf('  T6.2: Büyük polar açıklık (R_eq/r0 = 1.5)\n');
fprintf('    q = %.4f, aspect_r = %.4f\n', p62.q, p62.aspect_r);
fprintf('    Polar ρ hatası = %.2e mm\n', abs(p62.rho(end) - r0_62));
if abs(p62.rho(end) - r0_62) < TOL_POS
    fprintf('    [PASS] Sığ dome profil oluşturuldu\n\n');
    total_pass = total_pass + 1;
else
    fprintf('    [FAIL] Sığ dome\n\n');
    total_fail = total_fail + 1;
end

% T6.3: Aspect ratio tablosu + asimptotik davranış (DÜZELTİLMİŞ)
fprintf('  T6.3: Aspect ratio (h_dome/R_eq) vs q tablosu\n');
fprintf('         DOĞRU asimptotik limit: √2·E(1/2) ≈ 1.9101\n');
fprintf('         (Koussios Fig.4.3 → N_hoop/N_polar oranı, aspect ratio DEĞİL)\n\n');

R_eq_test = 100;
Y_eq_list = [1.2, 1.5, 2.0, 3.0, 4.0, 5.0, 7.0, 10.0, 15.0, 20.0];
fprintf('     R_eq/r0         q    h/R_eq      h_dome     s_total\n');
fprintf('    ------------------------------------------------------\n');

aspect_list = zeros(size(Y_eq_list));
for i = 1:length(Y_eq_list)
    r0_i = R_eq_test / Y_eq_list(i);
    p_i  = isotensoid_dome_profile(R_eq_test, r0_i);
    aspect_list(i) = p_i.aspect_r;
    fprintf('    %5.1f  %8.2f    %.4f   %10.4f   %10.4f\n', ...
            Y_eq_list(i), p_i.q, p_i.aspect_r, p_i.h_dome, p_i.s_total);
end

% Asimptotik limit: q → ∞ iken h/R_eq → √2·E(π/2, 1/2)
% Türetme: m = q/(1+2q) → 1/2,  √(1+2q)/(Y_eq) → √2·r0/R_eq·Y_eq → √2
% h/(r0·Y_eq) = √(1+2q)/Y_eq · [E(m) - K(m)/(1+2q)]
% q→∞: √(1+2q)/Y_eq → √2,  K/(1+2q) → 0
% → h/R_eq → √2 · E(1/2)
asymp_ref = sqrt(2) * ellipticE(pi/2, 0.5);
fprintf('\n    Asimptotik limit (teorik): √2·E(1/2) = %.4f\n', asymp_ref);
fprintf('    h/R_eq(q=399):                         %.4f\n', aspect_list(end));
err_asymp = abs(aspect_list(end) - asymp_ref) / asymp_ref;
fprintf('    Bağıl fark:                             %.4e\n', err_asymp);

if err_asymp < 0.01   % Sonlu q'da %1'den az fark beklenir
    fprintf('    [PASS] Asimptotik davranış doğru (h/R_eq → %.4f)\n\n', asymp_ref);
    total_pass = total_pass + 1;
else
    fprintf('    [FAIL] Asimptotik davranış hatası = %.2e\n\n', err_asymp);
    total_fail = total_fail + 1;
end

% T6.4: Hata fırlatma testleri
fprintf('  T6.4: Hata fırlatma testleri:\n');
err_tests = {
    'T6.4a r0=0',       @() isotensoid_dome_profile(100, 0);
    'T6.4b r0<0',       @() isotensoid_dome_profile(100, -10);
    'T6.4c r0>R_eq',    @() isotensoid_dome_profile(100, 150);
    'T6.4d r0=R_eq',    @() isotensoid_dome_profile(100, 100);
    'T6.4e R_eq<0',     @() isotensoid_dome_profile(-100, 50);
};
for i = 1:size(err_tests, 1)
    try
        err_tests{i, 2}();
        fprintf('    [FAIL] %s kabul EDİLDİ (hata beklendi)\n', err_tests{i, 1});
        total_fail = total_fail + 1;
    catch
        fprintf('    [PASS] %s reddedildi\n', err_tests{i, 1});
        total_pass = total_pass + 1;
    end
end

%% =====================================================================
%  TEST SETİ 7: Sayısal vs Analitik Türev Karşılaştırması
%  =====================================================================
fprintf('\n--- TEST SETİ 7: Sayısal vs Analitik Türev Karşılaştırması ---\n\n');

p7  = isotensoid_dome_profile(R_eq_t3, r0_t3, 2000);
ds7 = p7.s(2) - p7.s(1);
N7  = length(p7.s);

drho_num = zeros(N7, 1);
dx_num   = zeros(N7, 1);
for i = 2:N7-1
    drho_num(i) = (p7.rho(i+1) - p7.rho(i-1)) / (2*ds7);
    dx_num(i)   = (p7.x_local(i+1) - p7.x_local(i-1)) / (2*ds7);
end

s7_norm    = p7.s / p7.s_total;
idx_inner7 = (s7_norm > 0.01) & (s7_norm < 0.90);
max_drho_err = max(abs(drho_num(idx_inner7) - p7.drho_ds(idx_inner7)));
max_dx_err   = max(abs(dx_num(idx_inner7) - p7.dx_ds(idx_inner7)));

fprintf('  ds = %.4e mm, analiz aralığı: s/s_total ∈ [0.01, 0.90]\n', ds7);
fprintf('  max|dρ/ds sayısal - analitik| = %.2e\n', max_drho_err);
fprintf('  max|dx/ds sayısal - analitik| = %.2e\n', max_dx_err);
if max_drho_err < TOL_DERIV && max_dx_err < TOL_DERIV
    fprintf('  [PASS] Sayısal-analitik türev uyumu\n');
    total_pass = total_pass + 1;
else
    fprintf('  [FAIL] Sayısal-analitik türev uyumu\n');
    total_fail = total_fail + 1;
end

%% =====================================================================
%  Görselleştirme
%  =====================================================================
fprintf('\n--- Görselleştirme oluşturuluyor ---\n\n');

% Ana profil
p_vis = isotensoid_dome_profile(R_eq_t3, r0_t3, 2000);

fig = figure('Position', [50 50 1600 900], 'Color', 'w');
sgtitle(sprintf('İzotensoid Dome (v2) — R_{eq}=%.1f mm, r_0=%.1f mm', R_eq_t3, r0_t3), ...
        'FontSize', 14, 'FontWeight', 'bold');

% --- Subplot 1: Meridyen profilleri karşılaştırması ---
subplot(2,3,1); hold on; grid on;
title('Meridyen Profilleri Karşılaştırması');
xlabel('x_{lokal} [mm]'); ylabel('ρ [mm]');

% Silindir referans
x_cyl = [0, 0];
rho_cyl = [R_eq_t3, R_eq_t3];
plot(x_cyl, rho_cyl, 'k--', 'LineWidth', 1.5);
plot(x_cyl, -rho_cyl, 'k--', 'LineWidth', 1.5);

% Hemispherical
try
    p_hemi = hemispherical_dome_profile(R_eq_t3, r0_t3);
    plot(p_hemi.x_local, p_hemi.rho, 'b-', 'LineWidth', 1.5);
    plot(p_hemi.x_local, -p_hemi.rho, 'b-', 'LineWidth', 1.5);
    has_hemi = true;
catch
    has_hemi = false;
end

% Elipsoidal
try
    p_ell = ellipsoidal_dome_profile(R_eq_t3, r0_t3, 0.7);
    plot(p_ell.x_local, p_ell.rho, 'g-', 'LineWidth', 1.5);
    plot(p_ell.x_local, -p_ell.rho, 'g-', 'LineWidth', 1.5);
    has_ell = true;
catch
    has_ell = false;
end

% İzotensoid
plot(p_vis.x_local, p_vis.rho, 'r-', 'LineWidth', 2);
plot(p_vis.x_local, -p_vis.rho, 'r-', 'LineWidth', 2);

% İşaretler
plot(0, R_eq_t3, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k');
plot(p_vis.h_dome, r0_t3, 'rs', 'MarkerSize', 8, 'MarkerFaceColor', 'r');

legend_entries = {'Silindir'};
if has_hemi, legend_entries{end+1} = 'Hemispherical'; end
if has_ell,  legend_entries{end+1} = 'Elipsoidal (k=0.7)'; end
legend_entries{end+1} = 'İzotensoid';
legend_entries{end+1} = 'Ekvator';
legend_entries{end+1} = 'Polar aç.';
% Legend sadece üst profil çizgileri için
legend(legend_entries, 'Location', 'southwest', 'FontSize', 7);
axis equal;

% --- Subplot 2: Eğrilik (v2 — işaret değişimi gösteriliyor) ---
subplot(2,3,2); hold on; grid on;
title('Meridyen Eğriliği κ_m (v2)');
xlabel('s / s_{total} [-]'); ylabel('κ_m · R_{eq} [-]');

s_n = p_vis.s / p_vis.s_total;
kappa_nondim = p_vis.kappa_m * R_eq_t3;

% s/s_total ∈ [0, 0.95] aralığında göster (polar singülarite hariç)
idx_plot = s_n < 0.95;
plot(s_n(idx_plot), kappa_nondim(idx_plot), 'r-', 'LineWidth', 2);
yline(0, 'k--', 'LineWidth', 0.5);

if has_hemi
    s_n_h = p_hemi.s / p_hemi.s_total;
    plot(s_n_h, p_hemi.kappa_m * R_eq_t3, 'b-', 'LineWidth', 1.5);
end
if has_ell
    s_n_e = p_ell.s / p_ell.s_total;
    plot(s_n_e, p_ell.kappa_m * R_eq_t3, 'g-', 'LineWidth', 1.5);
end

legend_k = {};
legend_k{end+1} = 'İzotensoid';
legend_k{end+1} = 'κ_m = 0';
if has_hemi, legend_k{end+1} = 'Hemispherical'; end
if has_ell,  legend_k{end+1} = 'Elipsoidal (k=0.7)'; end
legend(legend_k, 'Location', 'best', 'FontSize', 7);

% --- Subplot 3: Winding açısı ---
subplot(2,3,3); hold on; grid on;
title('Geodesic Winding Açısı (Clairaut)');
xlabel('s / s_{total} [-]'); ylabel('α [°]');
plot(s_n, rad2deg(p_vis.alpha_w), 'r-', 'LineWidth', 2);
plot(s_n, rad2deg(asin(r0_t3 ./ p_vis.rho)), 'k--', 'LineWidth', 1);
legend('α(s)', 'arcsin(r_0/ρ)', 'Location', 'northwest');

% --- Subplot 4: Türevler ---
subplot(2,3,4); hold on; grid on;
title('dρ/ds ve dx/ds (izotensoid)');
xlabel('s [mm]'); ylabel('Türev [-]');
plot(p_vis.s, p_vis.drho_ds, 'b-', 'LineWidth', 1.5);
plot(p_vis.s, p_vis.dx_ds, 'r-', 'LineWidth', 1.5);
legend('dρ/ds', 'dx/ds', 'Location', 'east');

% --- Subplot 5: GATE-1a-01 çapraz doğrulama ---
subplot(2,3,5); hold on; grid on;
title('GATE-1a-01: Eliptik Int. vs ODE (v2)');
xlabel('Y [-]'); ylabel('|ΔZ| (boyutsuz)');

% Y aralığında x_local farkı
if ~isempty(Y_safe)
    plot(Y_safe, abs(Z_ode_safe - Z_ell_safe), 'b-', 'LineWidth', 1.5);
    ylabel('|ΔZ| (boyutsuz)');
end

% --- Subplot 6: 3D dönel yüzey ---
subplot(2,3,6); hold on; grid on;
title('3D Dönel Yüzey (İzotensoid)');
xlabel('x [mm]'); ylabel('Y [mm]'); zlabel('Z [mm]');

phi_3d = linspace(0, 2*pi, 60);
[PHI, S_grid] = meshgrid(phi_3d, p_vis.s(1:4:end));
RHO_grid = interp1(p_vis.s, p_vis.rho, S_grid(:,1), 'pchip');
X_grid   = interp1(p_vis.s, p_vis.x_local, S_grid(:,1), 'pchip');

Y_3d = RHO_grid .* cos(PHI);
Z_3d = RHO_grid .* sin(PHI);
X_3d = repmat(X_grid, 1, size(PHI, 2));

surf(X_3d, Y_3d, Z_3d, 'FaceAlpha', 0.7, 'EdgeAlpha', 0.1);

% Silindir ekleme
phi_cyl = linspace(0, 2*pi, 60);
x_cyl_3d = linspace(-50, 0, 10);
[PHI_c, X_c] = meshgrid(phi_cyl, x_cyl_3d);
Y_c = R_eq_t3 * cos(PHI_c);
Z_c = R_eq_t3 * sin(PHI_c);
surf(X_c, Y_c, Z_c, 'FaceColor', [0.8 0.8 0.8], 'FaceAlpha', 0.3, 'EdgeAlpha', 0.1);

axis equal; view(30, 20);
colormap(gca, 'jet');
lighting gouraud; camlight;

% Kaydet
saveas(fig, 'isotensoid_dome_verification.png');
fprintf('  Grafik kaydedildi: isotensoid_dome_verification.png\n');

%% =====================================================================
%  TEST SETİ 8: Karşılaştırmalı Özet
%  =====================================================================
fprintf('\n--- TEST SETİ 8: Karşılaştırmalı Özet ---\n\n');

fprintf('  TABLO 1: Dome tipleri karşılaştırması (R_eq=%.1f, r0=%.1f)\n', R_eq_t3, r0_t3);
fprintf('  Dome Tipi            s_total      h_dome      h/R_eq   κ_eq [1/mm]    α_eq [°]\n');
fprintf('  ------------------------------------------------------------------------------\n');

% Hemispherical
%   kappa_eq = kappa_m(1) = 1/R_eq (sabit — küre)
%   alpha_w  = asin(r0/rho) — Clairaut bağıntısı (tüm dome tipleri)
%   NOT: hemispherical_dome_profile struct'ında kappa_eq ve alpha_w
%        alanları tanımlı değil — burada uyumluluk hesabı yapılır.
%        Kalıcı düzeltme: profil fonksiyonlarına bu alanları ekleyin.
if has_hemi
    hemi_kappa_eq = p_hemi.kappa_m(1);              % 1/R_eq
    hemi_alpha_eq = asin(r0_t3 / p_hemi.rho(1));    % Clairaut @ ekvator
    fprintf('  Hemispherical     %10.4f  %10.4f      %.4f  %.6e       %.2f\n', ...
            p_hemi.s_total, p_hemi.h_dome, p_hemi.h_dome/R_eq_t3, ...
            hemi_kappa_eq, rad2deg(hemi_alpha_eq));
end

% Elipsoidal
%   kappa_eq = kappa_m(1) = 1/(R_eq·k²)
%   alpha_w  = asin(r0/rho) — Clairaut bağıntısı
if has_ell
    ell_kappa_eq = p_ell.kappa_m(1);                % 1/(R_eq·k²)
    ell_alpha_eq = asin(r0_t3 / p_ell.rho(1));      % Clairaut @ ekvator
    fprintf('  Elipsoidal k=0.7  %10.4f  %10.4f      %.4f  %.6e       %.2f\n', ...
            p_ell.s_total, p_ell.h_dome, p_ell.h_dome/R_eq_t3, ...
            ell_kappa_eq, rad2deg(ell_alpha_eq));
end

% İzotensoid
fprintf('  İzotensoid        %10.4f  %10.4f      %.4f  %.6e       %.2f\n', ...
        p_vis.s_total, p_vis.h_dome, p_vis.h_dome/R_eq_t3, ...
        p_vis.kappa_eq, rad2deg(p_vis.alpha_w(1)));

fprintf('\n  Fiziksel farklar:\n');
fprintf('    - Hemispherical: κ_eq = 1/R_eq (sabit eğrilik)\n');
fprintf('    - Elipsoidal:    κ_eq = 1/(R_eq·k²) (k ile ayarlanabilir)\n');
fprintf('    - İzotensoid:    κ_eq fizik tarafından belirlenir (netting theory)\n');
fprintf('    - İzotensoid h/R_eq = %.4f (değiştirilemez — optimal dome)\n', p_vis.aspect_r);
fprintf('    - İzotensoid dome''da κ_m İŞARET DEĞİŞTİRİR (bükülme noktası)\n');

% Tablo 2: q değişimi etkisi
fprintf('\n  TABLO 2: q değişimi etkisi (İzotensoid, R_eq=%.1f)\n', R_eq_t3);
fprintf('      r0      Y_eq         q     s_total      h/R_eq   κ_eq [1/mm]\n');
fprintf('  --------------------------------------------------------------\n');
r0_list = [100, 80, 60, 45, 30, 20];
for i = 1:length(r0_list)
    p_t2 = isotensoid_dome_profile(R_eq_t3, r0_list(i));
    fprintf('  %6.1f    %6.2f    %6.2f  %10.4f      %.4f  %.6e\n', ...
            r0_list(i), R_eq_t3/r0_list(i), p_t2.q, ...
            p_t2.s_total, p_t2.aspect_r, p_t2.kappa_eq);
end

%% =====================================================================
%  SONUÇ ÖZETİ
%  =====================================================================
fprintf('\n=============================================================\n');
fprintf(' SONUÇ ÖZETİ\n');
fprintf('=============================================================\n');
fprintf('  Toplam PASS: %d\n', total_pass);
fprintf('  Toplam FAIL: %d\n', total_fail);
if total_fail == 0
    fprintf('\n  ► TÜM TESTLER GEÇTİ ✓\n');
    if gate_pass
        fprintf('  ► GATE-1a-01 ONAYLI ✓\n');
    end
else
    fprintf('\n  ► %d TEST BAŞARISIZ ✗\n', total_fail);
    if ~gate_pass
        fprintf('  ► GATE-1a-01 BAŞARISIZ ✗\n');
    end
end
fprintf('=============================================================\n');
