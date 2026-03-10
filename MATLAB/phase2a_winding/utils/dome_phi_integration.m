function [result] = dome_phi_integration(dome_profile, R_E, epsilon)
% DOME_PHI_INTEGRATION  Dome bolgesinde paralel aci (phi) entegrasyonu.
%
% Phase-2a: Geodesic sarim yolu — dome phi ODE cozumu.
% Cekirdek ODE (Eq. 5.5): dphi/ds = R_E / (rho(s) * sqrt(rho(s)^2 - R_E^2))
% Singularite yonetimi (Eq. 6.6): analitik kuyruk tamamlama.
%
% KULLANIM:
%   result = dome_phi_integration(dome_profile, R_E)
%   result = dome_phi_integration(dome_profile, R_E, epsilon)
%
% GIRDILER:
%   dome_profile : Phase-1a profil struct (.s, .rho, .drho_ds, .kappa_m, ...)
%   R_E          : Efektif polar aciklik yaricapi [mm] = r0 + BW_eff/2
%   epsilon      : Singularite tamponu [mm] (varsayilan: 1e-3, Eq. 6.7)
%
% CIKTI:
%   result : Struct
%     .phi_dome       — Toplam dome phi katkisi [rad] (numerik + analitik)
%     .phi_numerical  — Numerik integrasyon katkisi [rad]
%     .phi_tail       — Analitik kuyruk katkisi [rad] (Eq. 6.6)
%     .s_epsilon      — Numerik durdurma noktasi [mm]
%     .drho_turn      — |drho/ds| donus noktasinda
%     .s_ode          — ODE ham cikti s noktalari [M x 1]
%     .phi_ode        — ODE ham cikti phi degerleri [M x 1]
%     .alpha_eq       — Ekvator winding acisi [rad]
%
% REFERANS: docs/phase2a_geodesic_math.md Bolum 5, 6

    %% --- Girdi dogrulama ---
    if nargin < 2
        error('dome_phi_integration:yetersizGirdi', ...
              'dome_profile ve R_E parametreleri gereklidir.');
    end
    if nargin < 3 || isempty(epsilon)
        epsilon = 1e-3;  % Eq. 6.7: 10 * tol_interp
    end

    % Struct alan kontrolleri
    required_fields = {'s', 'rho', 'drho_ds', 'R_eq', 'r0', 's_total'};
    for i = 1:numel(required_fields)
        if ~isfield(dome_profile, required_fields{i})
            error('dome_phi_integration:eksikAlan', ...
                  'dome_profile.%s alani bulunamadi.', required_fields{i});
        end
    end

    R_eq = dome_profile.R_eq;

    % Fiziksel kontroller
    if R_E >= R_eq
        error('dome_phi_integration:gecersizRE', ...
              'R_E (%.4f) < R_eq (%.4f) olmalidir.', R_E, R_eq);
    end
    if R_E <= 0
        error('dome_phi_integration:gecersizRE', 'R_E > 0 olmalidir.');
    end
    if epsilon <= 0 || epsilon > 1
        error('dome_phi_integration:gecersizEpsilon', ...
              'epsilon (0, 1] araliginda olmalidir. Girilen: %g', epsilon);
    end

    %% --- Interpolasyon hazirligi ---
    s_prof   = dome_profile.s;
    rho_prof = dome_profile.rho;
    drho_prof = dome_profile.drho_ds;

    % rho(s) ve drho/ds(s) icin pchip interpolantlari olustur
    rho_interp  = griddedInterpolant(s_prof, rho_prof, 'pchip');
    drho_interp = griddedInterpolant(s_prof, drho_prof, 'pchip');

    %% --- s_epsilon tespiti: rho(s_eps) = R_E + epsilon (ters sorgu) ---
    % Dome bolgesinde rho(s) monoton azalir: rho(0)=R_eq -> rho(s_total)=r0
    % Binary search ile s_epsilon bul
    rho_target = R_E + epsilon;

    % Hedef rho degeri profil araliginda mi?
    if rho_target > R_eq
        error('dome_phi_integration:hedefDisinda', ...
              'R_E + epsilon (%.4f) > R_eq (%.4f). Dome icinde bulunamaz.', ...
              rho_target, R_eq);
    end
    if rho_target < dome_profile.r0
        error('dome_phi_integration:hedefDisinda', ...
              'R_E + epsilon (%.4f) < r0 (%.4f). Dome disinda.', ...
              rho_target, dome_profile.r0);
    end

    % Binary search
    s_lo = 0;
    s_hi = dome_profile.s_total;
    tol_bs = 1e-10;  % Karar-11 Katman 3 toleransi
    max_iter = 100;

    for iter = 1:max_iter
        s_mid = (s_lo + s_hi) / 2;
        rho_mid = rho_interp(s_mid);
        if abs(rho_mid - rho_target) < tol_bs
            break;
        end
        if rho_mid > rho_target
            % rho hala buyuk, daha ileriye git
            s_lo = s_mid;
        else
            s_hi = s_mid;
        end
    end
    s_epsilon = s_mid;

    %% --- ODE entegrasyonu: dphi/ds = R_E / (rho * sqrt(rho^2 - R_E^2)) ---
    % Eq. 5.5
    odefun = @(s, phi) R_E / (rho_interp(s) * sqrt(rho_interp(s)^2 - R_E^2));

    % ode45 ile adaptif cozum
    opts = odeset('RelTol', 1e-10, 'AbsTol', 1e-12, ...
                  'MaxStep', s_epsilon / 100);

    [s_ode, phi_ode] = ode45(odefun, [0, s_epsilon], 0, opts);

    phi_numerical = phi_ode(end);

    %% --- Analitik kuyruk tamamlama (Eq. 6.6) ---
    % drho/ds donus noktasinda (s ~ s_total, rho ~ R_E)
    % drho_ds profil sonuna yakin deger
    drho_at_turn = drho_interp(dome_profile.s_total);
    drho_turn_abs = abs(drho_at_turn);

    % Dejenere donus kontrolu (S-GEO-01)
    if drho_turn_abs < 1e-8
        warning('dome_phi_integration:dejenerDonus', ...
                '|drho/ds| donus noktasinda cok kucuk: %.2e. Sonuc guvenilmez olabilir.', ...
                drho_turn_abs);
        % Ikinci derece yaklasim yerine kucuk bir alt sinir kullan
        drho_turn_abs = max(drho_turn_abs, 1e-8);
    end

    % Eq. 6.6: Delta_phi_tail = (1/|drho_turn|) * sqrt(2*epsilon/R_E)
    phi_tail = (1 / drho_turn_abs) * sqrt(2 * epsilon / R_E);

    % Toplam dome phi
    phi_dome = phi_numerical + phi_tail;

    %% --- Ekvator winding acisi ---
    alpha_eq = asin(R_E / R_eq);

    %% --- Cikti struct ---
    result.phi_dome      = phi_dome;
    result.phi_numerical = phi_numerical;
    result.phi_tail      = phi_tail;
    result.s_epsilon     = s_epsilon;
    result.drho_turn     = drho_turn_abs;
    result.s_ode         = s_ode;
    result.phi_ode       = phi_ode;
    result.alpha_eq      = alpha_eq;

end
