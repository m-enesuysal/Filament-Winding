function result = geodesic_single_circuit(dome_profile, R_eq, r0, L_cyl, BW_eff, N_points)
% GEODESIC_SINGLE_CIRCUIT  Tek devre geodesic sarım yolunu hesaplar.
%
% Phase-2a: Geodesic path — tek devre toplam açısal ilerleme (Δφ_circuit).
% Phase-2a math doc Bölüm 8, Denklem 8.1: Δφ_circuit = 4·φ_dome + 2·φ_cyl
%
% Clairaut ilişkisi: ρ(s)·sin(α(s)) = R_E = r0 + BW_eff/2
% Dome ODE: dφ/ds = R_E / (ρ(s)·√(ρ(s)² - R_E²))   [Eq. 5.5]
% Singülarite yönetimi: analitik kuyruk tamamlama       [Eq. 6.6]
% Silindirik bölge: φ_cyl = L_cyl·tan(α_eq) / R_eq     [Eq. 7.3]
%
% KULLANIM:
%   result = geodesic_single_circuit(dome_profile, R_eq, r0, L_cyl, BW_eff)
%   result = geodesic_single_circuit(dome_profile, R_eq, r0, L_cyl, BW_eff, N_points)
%
% GİRDİLER:
%   dome_profile — Phase-1a profil struct (.s, .rho, .drho_ds, .kappa_m)
%   R_eq         — Ekvator yarıçapı [mm]
%   r0           — Polar açıklık yarıçapı [mm]
%   L_cyl        — Silindir uzunluğu [mm] (>= 0 zorunlu)
%   BW_eff       — Efektif bant genişliği [mm] (= N_tow × BW)
%   N_points     — Yol nokta sayısı (varsayılan: 2000)
%
% ÇIKTI:
%   result struct:
%     .delta_phi_circuit — Tek devre toplam açısal ilerleme [rad]
%     .phi_dome          — Tek dome yarısı φ katkısı [rad]
%     .phi_cyl           — Tek silindir geçişi φ katkısı [rad]
%     .alpha_eq          — Ekvator winding açısı [rad]
%     .R_E               — Efektif polar açıklık yarıçapı [mm]
%     .phi_dome_numerical — Sayısal entegrasyon katkısı [rad]
%     .phi_tail          — Analitik kuyruk tamamlama katkısı [rad]
%     .epsilon           — Singülarite kesme mesafesi [mm]
%
% Referans: Phase-2a math doc (docs/phase2a_geodesic_math.md)
% Karar-11 toleransları uygulanır.
% Tarih: 2026-03-15
% Faz: Phase-2a

    %% --- Girdi doğrulama ---
    if nargin < 5
        error('geodesic_single_circuit:yetersizGirdi', ...
              'En az 5 parametre gereklidir: dome_profile, R_eq, r0, L_cyl, BW_eff');
    end
    if nargin < 6 || isempty(N_points)
        N_points = 2000;
    end

    % L_cyl < 0 kontrolü
    if L_cyl < 0
        error('geodesic_single_circuit:negativeLcyl', ...
              'L_cyl >= 0 olmalıdır. Girilen: %g mm', L_cyl);
    end

    % Parametre kontrolleri
    if R_eq <= 0
        error('geodesic_single_circuit:gecersizReq', ...
              'R_eq pozitif olmalıdır. Girilen: %g', R_eq);
    end
    if r0 <= 0 || r0 >= R_eq
        error('geodesic_single_circuit:gecersizR0', ...
              'r0 ∈ (0, R_eq) olmalıdır. r0 = %g, R_eq = %g', r0, R_eq);
    end
    if BW_eff <= 0
        error('geodesic_single_circuit:gecersizBWeff', ...
              'BW_eff pozitif olmalıdır. Girilen: %g', BW_eff);
    end

    %% --- Efektif polar açıklık ve Clairaut sabiti ---
    R_E = r0 + BW_eff / 2;

    if R_E >= R_eq
        error('geodesic_single_circuit:RE_buyuk_Req', ...
              'R_E (= r0 + BW_eff/2 = %g) >= R_eq (= %g). Fiziksel olarak imkansız.', ...
              R_E, R_eq);
    end

    %% --- Ekvator winding açısı (Karar-7: Konvansiyon A) ---
    alpha_eq = asin(R_E / R_eq);    % [rad]

    % Hoop winding kontrolü (S-WIND-02)
    if alpha_eq > 85 * pi / 180
        error('geodesic_single_circuit:hoopWinding', ...
              'alpha_eq = %.2f° > 85°: hoop winding bölgesi, dome yolu üretilemez.', ...
              alpha_eq * 180 / pi);
    end

    % Polar winding uyarısı (S-WIND-03)
    if alpha_eq < 5 * pi / 180
        warning('geodesic_single_circuit:polarWinding', ...
                'alpha_eq = %.2f° < 5°: polar winding bölgesi, çok az açısal ilerleme.', ...
                alpha_eq * 180 / pi);
    end

    %% --- Dome φ entegrasyonu ---
    % Dome profil verileri
    s_dome   = dome_profile.s;         % [mm], monoton artan, s=0 ekvator
    rho_dome = dome_profile.rho;       % [mm], monoton azalan, R_eq → r0

    % Profil tutarlılık kontrolü
    if abs(rho_dome(1) - R_eq) > 1e-3
        warning('geodesic_single_circuit:profilUyumsuz', ...
                'dome_profile.rho(1) = %g ≠ R_eq = %g. Fark: %g mm', ...
                rho_dome(1), R_eq, abs(rho_dome(1) - R_eq));
    end

    % Dönüş noktası kontrolü
    if min(rho_dome) > R_E
        error('geodesic_single_circuit:turnaroundYok', ...
              'Dome profili dönüş noktasına ulaşmıyor. min(ρ) = %g > R_E = %g', ...
              min(rho_dome), R_E);
    end

    % --- Singülarite kesme (Eq. 6.7): ε = 1e-3 mm ---
    epsilon = 1e-3;     % [mm]

    % s_epsilon bul: ρ(s_eps) = R_E + ε
    % rho_dome monoton azalan → interp1 için flip gerekli
    [rho_sorted, sort_idx] = sort(rho_dome, 'ascend');
    s_sorted = s_dome(sort_idx);
    % Duplicate kontrolü (pchip duplicate x sorununu önlemek için)
    [rho_sorted, unique_idx] = unique(rho_sorted, 'stable');
    s_sorted = s_sorted(unique_idx);

    rho_target = R_E + epsilon;
    s_eps = interp1(rho_sorted, s_sorted, rho_target, 'pchip');

    % --- Sayısal entegrasyon: MATLAB integral ---
    % İnterpolantlar
    rho_interp = griddedInterpolant(s_dome, rho_dome, 'pchip');

    % Integrand: dφ/ds = R_E / (ρ(s)·√(ρ(s)² - R_E²))  [Eq. 5.5]
    integrand = @(s) dome_phi_integrand(s, rho_interp, R_E);

    % Karar-11 Katman 1 toleransları
    phi_dome_numerical = integral(integrand, 0, s_eps, ...
                                  'RelTol', 1e-10, 'AbsTol', 1e-12);

    % --- Analitik kuyruk tamamlama (Eq. 6.6) ---
    % Δφ_tail ≈ (1 / |dρ/ds_turn|) · √(2ε / R_E)
    drho_interp = griddedInterpolant(s_dome, dome_profile.drho_ds, 'pchip');

    % s_turn: ρ = R_E noktası
    rho_turn_target = R_E;
    if rho_turn_target < min(rho_sorted)
        rho_turn_target = min(rho_sorted);
    end
    s_turn = interp1(rho_sorted, s_sorted, rho_turn_target, 'pchip');

    drho_at_turn = drho_interp(s_turn);
    drho_abs = abs(drho_at_turn);

    % Dejenere dönüş noktası kontrolü
    if drho_abs < 1e-8
        error('geodesic_single_circuit:dejenereTurnaround', ...
              '|dρ/ds| = %g < 1e-8 dönüş noktasında. Dejenere geometri.', drho_abs);
    end

    phi_tail = (1 / drho_abs) * sqrt(2 * epsilon / R_E);

    % --- Toplam dome φ ---
    phi_dome = phi_dome_numerical + phi_tail;

    %% --- Silindirik bölge φ ilerleme (Eq. 7.3) ---
    phi_cyl = L_cyl * tan(alpha_eq) / R_eq;

    %% --- Tek devre toplam açısal ilerleme (Eq. 8.1) ---
    % Δφ_circuit = 4·φ_dome + 2·φ_cyl
    delta_phi_circuit = 4 * phi_dome + 2 * phi_cyl;

    %% --- Çıktı struct ---
    result.delta_phi_circuit = delta_phi_circuit;
    result.phi_dome          = phi_dome;
    result.phi_cyl           = phi_cyl;
    result.alpha_eq          = alpha_eq;
    result.R_E               = R_E;
    result.phi_dome_numerical = phi_dome_numerical;
    result.phi_tail          = phi_tail;
    result.epsilon           = epsilon;
    result.s_turn            = s_turn;

end


function val = dome_phi_integrand(s, rho_interp, R_E)
% Dome φ ODE sağ tarafı: R_E / (ρ·√(ρ² - R_E²))
% Vektörleştirilmiş (integral() birden fazla nokta gönderebilir)
    rho = rho_interp(s(:));
    rho2 = rho.^2;
    RE2  = R_E^2;

    % Güvenlik: ρ² - R_E² < 0 olmamalı (Clairaut ihlali)
    diff_sq = rho2 - RE2;
    diff_sq(diff_sq < 0) = eps;     % makine epsilon ile koruma

    val = R_E ./ (rho .* sqrt(diff_sq));
    val = reshape(val, size(s));     % integral() beklediği boyut
end
