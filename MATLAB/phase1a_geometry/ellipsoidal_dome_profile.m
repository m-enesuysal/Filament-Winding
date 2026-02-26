function [profil] = ellipsoidal_dome_profile(R_eq, r0, k, N_pts)
% ELLIPSOIDAL_DOME_PROFILE  Elipsoidal dome meridyen profili üretir.
%
% Phase-1a: Mandrel geometri tanımı — elipsoidal dome tipi.
% Karar-3 uyumlu: parametrik meridyen temsili s(t) -> (rho(t), x(t))
% Karar-4 uyumlu: s = 0 ekvator, s artan polar açıklığa doğru
% Karar-5 uyumlu: k = h_dome / R_eq (dome aspect ratio)
% Karar-8 uyumlu: dahili hesaplamalar radyan, uzunluklar mm
% Karar-9 uyumlu: çıktı lookup table formatında {s, rho, x, drho_ds, dx_ds, kappa}
% Karar-11 toleransları uygulanır.
%
% KULLANIM:
%   profil = ellipsoidal_dome_profile(R_eq, r0, k)
%   profil = ellipsoidal_dome_profile(R_eq, r0, k, N_pts)
%
% GİRDİLER:
%   R_eq  : Ekvator (silindir) yarıçapı [mm]             (> 0)
%   r0    : Polar açıklık yarıçapı [mm]                   (0 < r0 < R_eq)
%   k     : Dome aspect ratio [-]                          (k > 0)
%           k < 1: yassı (oblate) dome
%           k = 1: hemispherical dome (özel hal)
%           k > 1: uzun (prolate) dome
%   N_pts : Profil nokta sayısı (isteğe bağlı, varsayılan: 500)
%
% ÇIKTI:
%   profil : Yapı (struct) — aşağıdaki alanları içerir:
%     .s         [N×1] Meridyen yay uzunluğu [mm]         (s=0 ekvator, uniform)
%     .rho       [N×1] Radyal konum [mm]                  (rho = R_eq..r0)
%     .x_local   [N×1] Aksiyel konum (dome başlangıcından) [mm]
%     .drho_ds   [N×1] dρ/ds [-]                          (boyutsuz)
%     .dx_ds     [N×1] dx/ds [-]                          (boyutsuz)
%     .beta      [N×1] Meridyen eğim açısı [rad]
%     .kappa_m   [N×1] Meridyen eğriliği [1/mm]
%     .theta     [N×1] Yardımcı açı parametresi [rad]
%     .R_eq      Ekvator yarıçapı [mm]
%     .r0        Polar açıklık yarıçapı [mm]
%     .k         Dome aspect ratio [-]
%     .theta_p   Polar açıklık açısı [rad]
%     .s_total   Toplam dome yay uzunluğu [mm]
%     .h_dome    Dome yüksekliği [mm]
%     .A_dome    Dome yüzey alanı [mm²]
%
% MATEMATİKSEL MODEL:
%   Elips yarı eksenleri: a = R_eq (radyal), b = k·R_eq (aksiyel)
%   Elips denklemi: (ρ/R_eq)² + (x/(k·R_eq))² = 1
%
%   Parametrik form (θ açı parametresiyle):
%     ρ(θ)       = R_eq · cos(θ)
%     x_local(θ) = k · R_eq · sin(θ)
%
%   Yardımcı fonksiyon:
%     f(θ) = √(sin²θ + k²·cos²θ)
%
%   Yay uzunluğu (eliptik integral — kapalı form yok):
%     ds/dθ  = R_eq · f(θ)
%     s(θ)   = R_eq · ∫₀^θ f(t) dt
%
%   s cinsinden türevler (zincir kuralı):
%     dρ/ds  = −sin(θ) / f(θ)
%     dx/ds  = k·cos(θ) / f(θ)
%
%   Eğim açısı:
%     β(θ) = atan2(sin(θ), k·cos(θ))
%
%   Meridyen eğriliği:
%     κ_m(θ) = k / (R_eq · f³(θ))
%
%   k=1 özel durumu hemispherical dome'a indirgenir (S-GEO-04).
%
% SAYISAL YÖNTEM:
%   s(θ) hesabı için yüksek çözünürlüklü (N_fine=10000) θ ızgarası
%   üzerinde cumtrapz kullanılır. Ardından θ(s) fonksiyonu monoton
%   interpolasyon (pchip) ile uniform s ızgarasına aktarılır.
%   Tüm diğer büyüklükler θ'dan analitik olarak hesaplanır.
%
% C¹ SÜREKLİLİK (S-GEO-02):
%   Ekvator noktasında (θ=0):
%     dρ/ds(0) = −sin(0)/f(0) = 0   → silindir ile eşleşir
%     dx/ds(0) = k·cos(0)/f(0) = 1   → silindir ile eşleşir
%
% EĞRİLİK SINIR DEĞERLERİ:
%   Ekvator:  κ_m(0)   = 1/(R_eq·k²)
%   k < 1: κ_m ekvator > κ_m polar  (yassı dome, eğrilik ekvatorda max)
%   k > 1: κ_m ekvator < κ_m polar  (uzun dome, eğrilik polarda max)
%   k = 1: κ_m sabit = 1/R_eq       (küre)
%
% Referans: Phase-0 Karar-3, Karar-5, Karar-9, Karar-11
% Tarih: 2026-02-26
% Faz: Phase-1a

    %% --- Girdi doğrulama (Karar-5 kısıtları) ---
    if nargin < 3
        error('ellipsoidal_dome_profile:yetersizGirdi', ...
              'En az R_eq, r0 ve k parametreleri gereklidir.');
    end
    if nargin < 4 || isempty(N_pts)
        N_pts = 500;
    end

    % Parametre kontrolleri
    if ~isscalar(R_eq) || ~isreal(R_eq) || R_eq <= 0
        error('ellipsoidal_dome_profile:gecersizReq', ...
              'R_eq pozitif reel skaler olmalıdır. Girilen: %g', R_eq);
    end
    if ~isscalar(r0) || ~isreal(r0) || r0 <= 0
        error('ellipsoidal_dome_profile:gecersizR0', ...
              'r0 pozitif reel skaler olmalıdır. Girilen: %g', r0);
    end
    if r0 >= R_eq
        error('ellipsoidal_dome_profile:r0BuyukEsitReq', ...
              'r0 < R_eq olmalıdır. r0 = %g, R_eq = %g', r0, R_eq);
    end
    if ~isscalar(k) || ~isreal(k) || k <= 0
        error('ellipsoidal_dome_profile:gecersizK', ...
              'k pozitif reel skaler olmalıdır. Girilen: %g', k);
    end
    if k < 0.1
        warning('ellipsoidal_dome_profile:cokKucukK', ...
                'k = %g çok küçük. k < 0.3 endüstriyel uygulamalar için nadir (S-GEO-03).', k);
    end
    if ~isscalar(N_pts) || N_pts < 10 || mod(N_pts, 1) ~= 0
        error('ellipsoidal_dome_profile:gecersizNpts', ...
              'N_pts >= 10 tamsayı olmalıdır. Girilen: %g', N_pts);
    end

    %% --- Temel büyüklükler ---
    theta_p = acos(r0 / R_eq);                % Polar açıklık açısı [rad]
    h_dome  = k * R_eq * sin(theta_p);         % Dome yüksekliği [mm]
    % Alternatif: h_dome = k * sqrt(R_eq^2 - r0^2)

    %% --- Yardımcı fonksiyon f(θ) ---
    f_func = @(th) sqrt(sin(th).^2 + k^2 * cos(th).^2);

    %% --- Yay uzunluğu s(θ): yüksek çözünürlüklü sayısal integrasyon ---
    N_fine  = 10000;   % İnce ızgara çözünürlüğü
    theta_fine = linspace(0, theta_p, N_fine)';
    f_fine     = f_func(theta_fine);
    ds_dtheta_fine = R_eq * f_fine;

    % Birikimli trapezoidal integrasyon → s(θ)
    s_fine = cumtrapz(theta_fine, ds_dtheta_fine);

    % Toplam yay uzunluğu — yüksek doğruluk referans değeri (MATLAB integral)
    s_total_ref = R_eq * integral(f_func, 0, theta_p, ...
                                  'RelTol', 1e-12, 'AbsTol', 1e-14);
    s_total = s_fine(end);

    % cumtrapz doğrulama
    s_total_err = abs(s_total - s_total_ref) / s_total_ref;
    if s_total_err > 1e-6
        warning('ellipsoidal_dome_profile:cumtrapzHata', ...
                'cumtrapz yay uzunluğu bağıl hatası %.2e > 1e-6. N_fine artırılmalı.', ...
                s_total_err);
    end
    % Referans değeri kullan (daha doğru)
    s_total = s_total_ref;

    %% --- Uniform s ızgarası ---
    s = linspace(0, s_total, N_pts)';

    %% --- θ(s) interpolasyonu: monoton (pchip) ---
    % s_fine monoton artan → θ(s) tek değerli ve monoton artan
    % Sınır noktalarını korumak için s_fine(end)'i s_total_ref ile eşitle
    s_fine(end) = s_total;
    theta = interp1(s_fine, theta_fine, s, 'pchip');

    % Sınır noktalarını zorla (makine hassasiyeti)
    theta(1)   = 0;
    theta(end)  = theta_p;

    %% --- Tüm büyüklükleri θ'dan analitik hesapla ---
    f       = f_func(theta);                       % Yardımcı fonksiyon

    % Pozisyon
    rho     = R_eq * cos(theta);                    % Radyal konum [mm]
    x_local = k * R_eq * sin(theta);                % Aksiyel konum [mm]

    % Birinci türevler (s cinsinden — zincir kuralıyla)
    drho_ds = -sin(theta) ./ f;                     % dρ/ds [-]
    dx_ds   =  k * cos(theta) ./ f;                 % dx/ds [-]

    % Eğim açısı
    beta    = atan2(sin(theta), k * cos(theta));     % β [rad]

    % Meridyen eğriliği
    kappa_m = k ./ (R_eq * f.^3);                   % κ_m [1/mm]

    %% --- Dome yüzey alanı (sayısal integrasyon) ---
    % A = 2π·R_eq² · ∫₀^{θ_p} cos(θ)·f(θ) dθ
    integrand_A = @(th) cos(th) .* f_func(th);
    A_integral  = integral(integrand_A, 0, theta_p, ...
                           'RelTol', 1e-12, 'AbsTol', 1e-14);
    A_dome = 2 * pi * R_eq^2 * A_integral;          % [mm²]

    %% --- Birim teğet vektör doğrulama ---
    tangent_norm = sqrt(drho_ds.^2 + dx_ds.^2);
    norm_error   = max(abs(tangent_norm - 1));
    if norm_error > 1e-10
        warning('ellipsoidal_dome_profile:tegetNormu', ...
                'Birim teğet vektör normu hatası: %.2e (beklenen: < 1e-10)', ...
                norm_error);
    end

    %% --- Sınır koşulları doğrulama ---
    % Ekvator (θ = 0)
    assert(abs(rho(1) - R_eq)    < 1e-10, 'Ekvator ρ hatası');
    assert(abs(x_local(1))       < 1e-10, 'Ekvator x hatası');
    assert(abs(drho_ds(1))       < 1e-10, 'Ekvator dρ/ds hatası (C¹ ihlali)');
    assert(abs(dx_ds(1) - 1)     < 1e-10, 'Ekvator dx/ds hatası');

    % Polar açıklık (θ = θ_p)
    assert(abs(rho(end) - r0)    < 1e-6, 'Polar açıklık ρ hatası');
    assert(abs(x_local(end) - h_dome) < 1e-6, 'Polar açıklık x hatası');

    %% --- Elips denklemi doğrulama ---
    ellipse_check = (rho / R_eq).^2 + (x_local / (k * R_eq)).^2;
    ellipse_err   = max(abs(ellipse_check - 1));
    if ellipse_err > 1e-10
        warning('ellipsoidal_dome_profile:elipsDenklemi', ...
                'Elips denklemi hatası: %.2e', ellipse_err);
    end

    %% --- Çıktı struct ---
    profil.s         = s;
    profil.rho       = rho;
    profil.x_local   = x_local;
    profil.drho_ds   = drho_ds;
    profil.dx_ds     = dx_ds;
    profil.beta      = beta;
    profil.kappa_m   = kappa_m;
    profil.theta     = theta;

    % Skaler parametreler
    profil.R_eq      = R_eq;
    profil.r0        = r0;
    profil.k         = k;
    profil.theta_p   = theta_p;
    profil.s_total   = s_total;
    profil.h_dome    = h_dome;
    profil.A_dome    = A_dome;

end
