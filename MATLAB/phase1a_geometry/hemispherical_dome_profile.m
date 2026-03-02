function [profil] = hemispherical_dome_profile(R_eq, r0, N_pts)
% HEMISPHERICAL_DOME_PROFILE  Hemispherical dome meridyen profili üretir.
%
% Phase-1a: Mandrel geometri tanımı — hemispherical dome tipi.
% Karar-3 uyumlu: parametrik meridyen temsili s(t) -> (rho(t), x(t))
% Karar-4 uyumlu: s = 0 ekvator, s artan polar açıklığa doğru
% Karar-8 uyumlu: dahili hesaplamalar radyan, uzunluklar mm
% Karar-9 uyumlu: çıktı lookup table formatında {s, rho, x, drho_ds, dx_ds, kappa}
% Karar-11 toleransları uygulanır.
%
% KULLANIM:
%   profil = hemispherical_dome_profile(R_eq, r0)
%   profil = hemispherical_dome_profile(R_eq, r0, N_pts)
%
% GİRDİLER:
%   R_eq  : Ekvator (silindir) yarıçapı [mm]             (> 0)
%   r0    : Polar açıklık yarıçapı [mm]                   (0 < r0 < R_eq)
%   N_pts : Profil nokta sayısı (isteğe bağlı, varsayılan: 500)
%
% ÇIKTI:
%   profil : Yapı (struct) — aşağıdaki alanları içerir:
%     .s         [N×1] Meridyen yay uzunluğu [mm]         (s=0 ekvator)
%     .rho       [N×1] Radyal konum [mm]                  (rho = R_eq..r0)
%     .x_local   [N×1] Aksiyel konum (dome başlangıcından) [mm]
%     .drho_ds   [N×1] dρ/ds [-]                          (boyutsuz)
%     .dx_ds     [N×1] dx/ds [-]                          (boyutsuz)
%     .beta      [N×1] Meridyen eğim açısı [rad]
%     .kappa_m   [N×1] Meridyen eğriliği [1/mm]
%     .theta     [N×1] Yardımcı açı parametresi [rad]
%     .R_eq      Ekvator yarıçapı [mm]
%     .r0        Polar açıklık yarıçapı [mm]
%     .theta_p   Polar açıklık açısı [rad]
%     .s_total   Toplam dome yay uzunluğu [mm]
%     .h_dome    Dome yüksekliği [mm]
%     .A_dome    Dome yüzey alanı [mm²]
%
% MATEMATİKSEL MODEL:
%   Küre yarıçapı R = R_eq, merkez ekvator noktasında mandrel ekseni üzerinde.
%
%   Parametrik form (yay uzunluğu s cinsinden):
%     θ(s)     = s / R_eq
%     ρ(s)     = R_eq · cos(s / R_eq)
%     x(s)     = R_eq · sin(s / R_eq)
%     dρ/ds    = −sin(s / R_eq)
%     dx/ds    = +cos(s / R_eq)
%     β(s)     = s / R_eq               (eğim açısı = θ)
%     κ_m(s)   = 1 / R_eq               (sabit — küre özelliği)
%
%   Dome sınırları:
%     s = 0       → ρ = R_eq  (ekvator)
%     s = s_total → ρ = r₀    (polar açıklık)
%     s_total = R_eq · arccos(r₀ / R_eq)
%
% SİNGÜLARİTE NOTU (S-GEO-01):
%   Hemispherical dome'da κ_m = 1/R_eq sabittir. Polar açıklık civarında
%   singülarite YOKTUR. Bu, isotensoid dome'dan önemli bir farktır ve
%   numerik avantaj sağlar.
%
% C¹ SÜREKLİLİK (S-GEO-02):
%   Ekvator noktasında (s=0):
%     Silindir: ρ = R_eq, dρ/ds = 0
%     Dome:     ρ = R_eq, dρ/ds = -sin(0) = 0   → C¹ sağlanır
%
% Referans: Phase-0 Karar-3, Karar-9, Karar-11
% Tarih: 2026-02-26
% Faz: Phase-1a

    %% --- Girdi doğrulama (Karar-5 kısıtları) ---
    if nargin < 2
        error('hemispherical_dome_profile:yetersizGirdi', ...
              'En az R_eq ve r0 parametreleri gereklidir.');
    end
    if nargin < 3 || isempty(N_pts)
        N_pts = 500;    % Varsayılan çözünürlük
    end

    % Parametre kontrolleri
    if ~isscalar(R_eq) || ~isreal(R_eq) || R_eq <= 0
        error('hemispherical_dome_profile:gecersizReq', ...
              'R_eq pozitif reel skaler olmalıdır. Girilen: %g', R_eq);
    end
    if ~isscalar(r0) || ~isreal(r0) || r0 <= 0
        error('hemispherical_dome_profile:gecersizR0', ...
              'r0 pozitif reel skaler olmalıdır. Girilen: %g', r0);
    end
    if r0 >= R_eq
        error('hemispherical_dome_profile:r0BuyukEsitReq', ...
              'r0 < R_eq olmalıdır. r0 = %g, R_eq = %g', r0, R_eq);
    end
    if ~isscalar(N_pts) || N_pts < 10 || mod(N_pts, 1) ~= 0
        error('hemispherical_dome_profile:gecersizNpts', ...
              'N_pts >= 10 tamsayı olmalıdır. Girilen: %g', N_pts);
    end

    %% --- Temel büyüklükler ---
    theta_p = acos(r0 / R_eq);          % Polar açıklık açısı [rad]
    s_total = R_eq * theta_p;            % Toplam dome yay uzunluğu [mm]
    h_dome  = R_eq * sin(theta_p);       % Dome yüksekliği [mm]
    %   Alternatif: h_dome = sqrt(R_eq^2 - r0^2)  (sayısal doğrulama için)

    %% --- Yay uzunluğu vektörü ---
    % Uniform örnekleme (hemispherical dome'da adaptif gerekmiyor — κ sabit)
    s = linspace(0, s_total, N_pts)';    % [N×1] sütun vektörü

    %% --- Parametrik hesaplamalar ---
    theta   = s / R_eq;                   % Yardımcı açı [rad]

    % Pozisyon
    rho     = R_eq * cos(theta);           % Radyal konum [mm]
    x_local = R_eq * sin(theta);           % Aksiyel konum (lokal) [mm]

    % Birinci türevler (analitik — yay uzunluğu parametrizasyonu)
    drho_ds = -sin(theta);                 % dρ/ds [-]
    dx_ds   =  cos(theta);                 % dx/ds [-]

    % Eğim açısı
    beta    = theta;                       % β = θ = s/R_eq [rad]

    % Meridyen eğriliği (sabit)
    kappa_m = ones(N_pts, 1) / R_eq;      % κ_m = 1/R_eq [1/mm]

    %% --- Dome yüzey alanı (analitik) ---
    % A = 2π · R_eq² · sin(θ_p) = 2π · R_eq · h_dome
    A_dome = 2 * pi * R_eq * h_dome;      % [mm²]

    %% --- Birim teğet vektör normu doğrulama ---
    tangent_norm = sqrt(drho_ds.^2 + dx_ds.^2);
    norm_error   = max(abs(tangent_norm - 1));
    if norm_error > 1e-14
        warning('hemispherical_dome_profile:tegetNormu', ...
                'Birim teğet vektör normu hatası: %.2e (beklenen: < 1e-14)', ...
                norm_error);
    end

    %% --- Sınır koşulları doğrulama ---
    % Ekvator (s = 0)
    assert(abs(rho(1) - R_eq)    < 1e-12, 'Ekvator ρ hatası');
    assert(abs(x_local(1))       < 1e-12, 'Ekvator x hatası');
    assert(abs(drho_ds(1))       < 1e-12, 'Ekvator dρ/ds hatası (C¹ ihlali)');
    assert(abs(dx_ds(1) - 1)     < 1e-12, 'Ekvator dx/ds hatası');

    % Polar açıklık (s = s_total)
    assert(abs(rho(end) - r0)    < 1e-10, 'Polar açıklık ρ hatası');
    assert(abs(x_local(end) - h_dome) < 1e-10, 'Polar açıklık x hatası');

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
    profil.theta_p   = theta_p;
    profil.s_total   = s_total;
    profil.h_dome    = h_dome;
    profil.A_dome    = A_dome;
    alpha_w = asin(r0 ./ rho);
    alpha_w(end) = pi/2;  % Polar açıklıkta tanım gereği 90°

    profil.alpha_w   = alpha_w;
    profil.kappa_eq  = 1 / R_eq;             % Sabit eğrilik — küre özelliği
    profil.kappa_pol = 1 / R_eq;             % Polar'da da aynı — küre
    profil.aspect_r  = h_dome / R_eq;        % Dome aspect ratio
    % --- YAMA BAŞLANGIÇ (hemispherical) ---
% Geodesic sarım açısı — Clairaut bağıntısı: ρ·sin(α) = r0
%   Hemispherical dome'da ρ = R_eq·cos(θ), θ = s/R_eq
%   → sin(α) = r0 / (R_eq·cos(θ))
%   → α = asin(r0 / ρ)
%   Ekvator: α_eq = asin(r0 / R_eq)
%   Polar:   α_pol = π/2 (tanım gereği — fiber polarda teğet)

% !! DİKKAT: Bu satırlar profil fonksiyonuna eklenecek, ayrı çalıştırılmayacak !!

    % Geodesic sarım açısı (Clairaut)
    


end
