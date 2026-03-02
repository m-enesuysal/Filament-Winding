function [profil] = isotensoid_dome_profile(R_eq, r0, N_pts)
% ISOTENSOID_DOME_PROFILE  İzotensoid dome meridyen profili üretir (v2).
%
% Phase-1a: Mandrel geometri tanımı — izotensoid dome tipi.
% Karar-3, Karar-4, Karar-5, Karar-8, Karar-9, Karar-11 uyumlu.
%
% v2 DÜZELTME NOTU (2026-02-27):
%   1. d²Z/dθ² formülü düzeltildi:
%      ESKİ (hatalı): -q²·sθcθ·(1+c²θ) / [(1+2q)^(3/2)·Δ³]
%      YENİ (doğru):  -q·sθcθ·[2(1+q)-q·s²θ] / [(1+2q)^(3/2)·Δ³]
%   2. Eğrilik κ_m: Parametrik çözümden analitik türetme
%      κ = (dY'·d²Z - dZ'·d²Y) / {r0·[(dY')²+(dZ')²]^(3/2)}
%      ESKİ (hatalı): 2qY/[r0(1+2q)(Y²-1)] — 1/(Y²-1) singülaritesi
%   3. İzotensoid dome'da κ_m İŞARET DEĞİŞTİRİR (bükülme noktası)
%      Ekvator: κ_m > 0  (konveks),  Polar: κ_m < 0  (konkav)
%   4. Aspect ratio asimptotik limit: √2·E(1/2) ≈ 1.91 (0.6 DEĞİL)
%
% KULLANIM:
%   profil = isotensoid_dome_profile(R_eq, r0)
%   profil = isotensoid_dome_profile(R_eq, r0, N_pts)
%
% GİRDİLER:
%   R_eq  : Ekvator (silindir) yarıçapı [mm]    (> 0)
%   r0    : Polar açıklık yarıçapı [mm]          (0 < r0 < R_eq)
%   N_pts : Profil nokta sayısı (varsayılan: 500)
%
% ÇIKTI struct (Karar-9 uyumlu):
%   .s, .rho, .x_local, .drho_ds, .dx_ds, .beta, .kappa_m, .alpha_w
%   .theta, .R_eq, .r0, .q, .m_ell
%   .s_total, .h_dome, .A_dome, .aspect_r, .kappa_eq, .kappa_pol
%
% MATEMATİKSEL MODEL — KOUSSİOS ELİPTİK İNTEGRAL ÇÖZÜMÜ:
%   Y(θ) = √(1 + q·sin²θ),  θ ∈ [0, π/2], θ=0 polar, θ=π/2 ekvator
%   Z(θ) = √(1+2q)·[E(θ,m) - F(θ,m)/(1+2q)]
%   q = Y_eq² - 1,  m = q/(1+2q)
%
% Referans: Koussios — "Filament Winding: a Unified Approach", Bölüm 3-4
% Tarih: 2026-02-27 (v2)
% Faz: Phase-1a

    %% === Girdi doğrulama ===
    if nargin < 2
        error('isotensoid_dome_profile:yetersizGirdi', ...
              'En az R_eq ve r0 parametreleri gereklidir.');
    end
    if nargin < 3 || isempty(N_pts)
        N_pts = 500;
    end
    if ~isscalar(R_eq) || ~isreal(R_eq) || R_eq <= 0
        error('isotensoid_dome_profile:gecersizReq', ...
              'R_eq pozitif reel skaler olmalıdır. Girilen: %g', R_eq);
    end
    if ~isscalar(r0) || ~isreal(r0) || r0 <= 0
        error('isotensoid_dome_profile:gecersizR0', ...
              'r0 pozitif reel skaler olmalıdır. Girilen: %g', r0);
    end
    if r0 >= R_eq
        error('isotensoid_dome_profile:r0BuyukEsitReq', ...
              'r0 < R_eq olmalıdır. r0 = %g, R_eq = %g', r0, R_eq);
    end
    if ~isscalar(N_pts) || N_pts < 10 || mod(N_pts, 1) ~= 0
        error('isotensoid_dome_profile:gecersizNpts', ...
              'N_pts >= 10 tamsayı olmalıdır. Girilen: %g', N_pts);
    end

    Y_eq = R_eq / r0;
    if Y_eq > 15
        warning('isotensoid_dome_profile:buyukYeq', ...
                'Y_eq = R_eq/r0 = %.1f > 15. Çok küçük polar açıklık — ODE stiff olabilir (S-GEO-01).', Y_eq);
    end
    if Y_eq < 1.2
        warning('isotensoid_dome_profile:kucukYeq', ...
                'Y_eq = R_eq/r0 = %.2f < 1.2. Çok büyük polar açıklık.', Y_eq);
    end

    %% === Koussios parametreleri ===
    q         = Y_eq^2 - 1;
    m_ell     = q / (1 + 2*q);
    sqrt_1p2q = sqrt(1 + 2*q);

    %% === İnce θ ızgarası (eliptik integral çözümü) ===
    N_fine     = 10000;
    theta_fine = linspace(0, pi/2, N_fine)';
    sin_f      = sin(theta_fine);
    cos_f      = cos(theta_fine);
    sin2_f     = sin_f.^2;
    cos2_f     = cos_f.^2;

    % Y(θ) — boyutsuz yarıçap
    Y_fine = sqrt(1 + q * sin2_f);

    % Z(θ) — boyutsuz aksiyel konum (Koussios eliptik integral çözümü)
    F_fine = arrayfun(@(th) ellipticF(th, m_ell), theta_fine);
    E_fine = arrayfun(@(th) ellipticE(th, m_ell), theta_fine);
    Z_fine = sqrt_1p2q * (E_fine - F_fine / (1 + 2*q));

    % Dome yüksekliği (kapalı form — tam eliptik integraller)
    K_complete = ellipticF(pi/2, m_ell);
    E_complete = ellipticE(pi/2, m_ell);
    h_dome     = r0 * sqrt_1p2q * (E_complete - K_complete / (1 + 2*q));

    % Konvansiyon dönüşümü (Karar-4): x_local = 0 ekvator, x_local = h_dome polar
    x_local_fine = h_dome - r0 * Z_fine;

    %% === Birinci türevler dY/dθ ve dZ/dθ (analitik) ===
    % dY/dθ = q·sinθ·cosθ / Y
    dY_f = q * sin_f .* cos_f ./ Y_fine;

    % dZ/dθ = q·(1+cos²θ) / [√(1+2q)·Δ]
    %   Δ = √(1 - m·sin²θ)
    Delta_f = sqrt(1 - m_ell * sin2_f);
    dZ_f    = q * (1 + cos2_f) ./ (sqrt_1p2q * Delta_f);

    %% === İkinci türevler (analitik — v2 DÜZELTİLMİŞ) ===
    % d²Y/dθ² = q·cos(2θ)/Y - q²·sin²θ·cos²θ/Y³
    d2Y_f = q * cos(2*theta_fine) ./ Y_fine ...
          - q^2 * sin2_f .* cos2_f ./ Y_fine.^3;

    % d²Z/dθ² — v2 DOĞRU TÜRETME:
    %   h(θ) = (2-sin²θ) / Δ
    %   dh/dθ = sinθcosθ · [qsin²θ - 2(1+q)] / [(1+2q)·Δ³]
    %   d²Z/dθ² = q/√(1+2q) · dh/dθ
    %           = -q·sinθcosθ·[2(1+q) - q·sin²θ] / [(1+2q)^(3/2)·Δ³]
    %
    %   v1 HATALI formül: -q²·sinθcosθ·(1+cos²θ) / [(1+2q)^(3/2)·Δ³]
    %   Fark: q²(2-s²θ) vs q(2+2q-qs²θ) → 2q fark (ciddi hata)
    d2Z_f = -q * sin_f .* cos_f .* (2*(1+q) - q*sin2_f) ...
          ./ ((1 + 2*q)^(1.5) * Delta_f.^3);

    %% === Eğrilik κ_m — parametrik formül (ince ızgara) ===
    % κ = (dY'·d²Z - dZ'·d²Y) / {r0 · [(dY')² + (dZ')²]^(3/2)}
    numer_f    = dY_f .* d2Z_f - dZ_f .* d2Y_f;
    denom_f    = (dY_f.^2 + dZ_f.^2).^(1.5);
    denom_safe = max(denom_f, 1e-60);
    kappa_m_fine = numer_f ./ (r0 * denom_safe);

    % Sınır değerleri analitik limitlerle düzelt
    % Ekvator (θ=π/2): κ_eq = (1+q) / (r0·q·Y_eq)
    %   Türetme: θ=π/2'de dY=0, d²Y=-q/Y_eq, dZ=q/√(1+q), d²Z=0
    %   pay = 0·0 - (q/√(1+q))·(-q/Y_eq) = q²/(Y_eq·√(1+q))
    %   payda = r0·(q²/(1+q))^(3/2) = r0·q³/(1+q)^(3/2)
    %   κ = (1+q)/(r0·q·Y_eq)
    kappa_eq_ref = (1 + q) / (r0 * q * Y_eq);
    kappa_m_fine(end) = kappa_eq_ref;

    % Polar (θ=0): κ_pol = -(1+2q) / (4·q·r0)
    %   Türetme: θ=0'da dY=0, d²Y=q, dZ=2q/√(1+2q), d²Z=0
    %   pay = 0·0 - (2q/√(1+2q))·q = -2q²/√(1+2q)
    %   payda = r0·(4q²/(1+2q))^(3/2) = r0·8q³/(1+2q)^(3/2)
    %   κ = -(1+2q)/(4qr0)   ← NEGATİF (konkav bölge)
    kappa_pol_ref = -(1 + 2*q) / (4 * q * r0);
    kappa_m_fine(1) = kappa_pol_ref;

    %% === Yay uzunluğu s(θ): cumtrapz + yüksek doğruluk referans ===
    ds_dtheta_f  = r0 * sqrt(dY_f.^2 + dZ_f.^2);
    s_koussios_f = cumtrapz(theta_fine, ds_dtheta_f);

    integrand_s = @(th) r0 * sqrt( ...
        (q * sin(th) .* cos(th) ./ sqrt(1 + q * sin(th).^2)).^2 + ...
        (q * (1 + cos(th).^2) ./ (sqrt_1p2q * sqrt(1 - m_ell * sin(th).^2))).^2);
    s_total_ref = integral(integrand_s, 0, pi/2, ...
                           'RelTol', 1e-12, 'AbsTol', 1e-14);
    s_total = s_total_ref;

    %% === Konvansiyon dönüşümü: s_our = s_total - s_Koussios ===
    s_our_f        = s_total - s_koussios_f;
    s_our_f(end)   = 0;          % θ=π/2 → ekvator, s=0
    s_our_f(1)     = s_total;    % θ=0   → polar, s=s_total

    % Monoton artan s için ters çevir (interpolasyon gereksinimi)
    s_for_interp     = flipud(s_our_f);
    theta_for_interp = flipud(theta_fine);
    kappa_for_interp = flipud(kappa_m_fine);

    %% === Uniform s ızgarası + θ(s) interpolasyonu ===
    s = linspace(0, s_total, N_pts)';
    theta      = interp1(s_for_interp, theta_for_interp, s, 'pchip');
    theta(1)   = pi/2;     % Ekvator
    theta(end)  = 0;        % Polar açıklık

    %% === Tüm büyüklükleri θ'dan analitik hesapla ===
    Y   = sqrt(1 + q * sin(theta).^2);
    rho = r0 * Y;

    F_vals  = arrayfun(@(th) ellipticF(th, m_ell), theta);
    E_vals  = arrayfun(@(th) ellipticE(th, m_ell), theta);
    Z_vals  = sqrt_1p2q * (E_vals - F_vals / (1 + 2*q));
    x_local = h_dome - r0 * Z_vals;

    % Geodesic winding açısı (Clairaut: ρ·sin(α) = r0)
    alpha_w = asin(1 ./ Y);

    %% === Türevler (s cinsinden, zincir kuralı) ===
    sin_th  = sin(theta);  cos_th = cos(theta);
    sin2_th = sin_th.^2;   cos2_th = cos_th.^2;
    dY_dth  = q * sin_th .* cos_th ./ Y;
    Delta   = sqrt(1 - m_ell * sin2_th);
    dZ_dth  = q * (1 + cos2_th) ./ (sqrt_1p2q * Delta);
    ds_dth  = r0 * sqrt(dY_dth.^2 + dZ_dth.^2);
    ds_dth_safe = max(ds_dth, 1e-30);

    % dρ/ds ve dx/ds (s artınca θ azalır → negatif zincir kuralı)
    drho_ds = -(r0 * dY_dth) ./ ds_dth_safe;
    dx_ds   =  (r0 * dZ_dth) ./ ds_dth_safe;

    % Sınır koşulları: her iki uçta dρ/ds=0, dx/ds=1
    drho_ds(1) = 0;    dx_ds(1) = 1;     % Ekvator (θ=π/2)
    drho_ds(end) = 0;  dx_ds(end) = 1;   % Polar   (θ=0)

    %% === β eğim açısı ===
    beta = atan2(-drho_ds, dx_ds);

    %% === κ_m interpolasyonu ===
    kappa_m      = interp1(s_for_interp, kappa_for_interp, s, 'pchip');
    kappa_m(1)   = kappa_eq_ref;
    kappa_m(end)  = kappa_pol_ref;

    %% === Dome yüzey alanı ===
    integrand_A = @(th) 2*pi * r0 * sqrt(1 + q*sin(th).^2) .* ...
        r0 .* sqrt( ...
            (q*sin(th).*cos(th)./sqrt(1+q*sin(th).^2)).^2 + ...
            (q*(1+cos(th).^2)./(sqrt_1p2q*sqrt(1-m_ell*sin(th).^2))).^2);
    A_dome = integral(integrand_A, 0, pi/2, ...
                      'RelTol', 1e-12, 'AbsTol', 1e-14);

    %% === Doğrulama kontrolleri ===
    tangent_norm = sqrt(drho_ds.^2 + dx_ds.^2);
    norm_error   = max(abs(tangent_norm(2:end-1) - 1));
    if norm_error > 1e-6
        warning('isotensoid_dome_profile:tegetNormu', ...
                'Birim teğet vektör normu hatası: %.2e', norm_error);
    end
    assert(abs(rho(1) - R_eq)         < 1e-6,  'Ekvator ρ hatası');
    assert(abs(x_local(1))            < 1e-6,  'Ekvator x hatası');
    assert(abs(rho(end) - r0)         < 1e-4,  'Polar açıklık ρ hatası');
    assert(abs(x_local(end) - h_dome) < 1e-4,  'Polar açıklık x hatası');

    %% === Çıktı struct (Karar-9 uyumlu) ===
    profil.s         = s;
    profil.rho       = rho;
    profil.x_local   = x_local;
    profil.drho_ds   = drho_ds;
    profil.dx_ds     = dx_ds;
    profil.beta      = beta;
    profil.kappa_m   = kappa_m;
    profil.alpha_w   = alpha_w;
    profil.theta     = theta;

    profil.R_eq      = R_eq;
    profil.r0        = r0;
    profil.q         = q;
    profil.m_ell     = m_ell;
    profil.s_total   = s_total;
    profil.h_dome    = h_dome;
    profil.A_dome    = A_dome;
    profil.aspect_r  = h_dome / R_eq;
    profil.kappa_eq  = kappa_eq_ref;
    profil.kappa_pol = kappa_pol_ref;
end
