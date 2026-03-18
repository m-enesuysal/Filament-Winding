function [result] = dome_phi_integration(dome_profile, R_E, epsilon)
% DOME_PHI_INTEGRATION  Dome bolgesinde paralel aci (phi) entegrasyonu.
%
% Phase-2a: Geodesic sarim yolu — dome phi ODE cozumu.
% Cekirdek ODE (Eq. 5.5): dphi/ds = R_E / (rho(s) * sqrt(rho(s)^2 - R_E^2))
% Singularite yonetimi: birlesik 1. + 2. derece analitik kuyruk tamamlama.
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
%     .phi_tail       — Analitik kuyruk katkisi [rad]
%     .s_epsilon      — Numerik durdurma noktasi [mm] (rho = R_E + epsilon)
%     .s_turn         — Donus noktasi [mm] (rho = R_E)
%     .delta_eps      — Kuyruk mesafesi [mm] (s_turn - s_epsilon)
%     .drho_turn      — |drho/ds| at s_epsilon
%     .d2rho_turn     — |d2rho/ds2| at s_epsilon
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

    %% --- s_turn tespiti: rho(s_turn) = R_E (donus noktasi) ---
    % Binary search ile s_turn bul (rho monoton azalir)
    rho_turn_target = R_E;

    s_lo2 = 0;
    s_hi2 = dome_profile.s_total;
    for iter = 1:max_iter
        s_mid2 = (s_lo2 + s_hi2) / 2;
        rho_mid2 = rho_interp(s_mid2);
        if abs(rho_mid2 - rho_turn_target) < tol_bs
            break;
        end
        if rho_mid2 > rho_turn_target
            s_lo2 = s_mid2;
        else
            s_hi2 = s_mid2;
        end
    end
    s_turn = s_mid2;

    %% --- Analitik kuyruk tamamlama (birlesik 1. + 2. derece) ---
    % Tureve s_epsilon noktasinda (numerik integrasyon siniri) deger
    % NOT: Eski kod s_total'da degerler aliyordu — izotensoid icin drho=0.
    %      Donus noktasi s_turn'de (rho=R_E), s_total'de (rho=r0) degil.

    % 1. turev: a = |drho/ds| at s_epsilon
    a = abs(drho_interp(s_epsilon));

    % 2. turev: b = |d2rho/ds2| at s_epsilon (merkezi sonlu fark)
    ds_fd = min(epsilon / 2, (dome_profile.s_total - s_epsilon) / 4);
    ds_fd = max(ds_fd, 1e-8);  % alt sinir
    drho_plus  = drho_interp(s_epsilon + ds_fd);
    drho_minus = drho_interp(s_epsilon - ds_fd);
    d2rho = (drho_plus - drho_minus) / (2 * ds_fd);
    b = abs(d2rho);

    % --- Birlesik kuyruk formulu ---
    % Taylor yaklasimi: rho(s) ~ R_E + a*(s_turn - s) + (b/2)*(s_turn - s)^2
    % s_epsilon noktasinda: rho = R_E + epsilon, delta_eps = s_turn - s_epsilon
    %
    % delta_eps cozumu (kararlil versiyon):
    %   a*d + (b/2)*d^2 = epsilon  =>  d = 2*eps / (a + sqrt(a^2 + 2*b*eps))
    %
    % Kuyruk integrali:
    %   phi_tail = integral_0^{delta_eps} R_E / (rho * sqrt(rho^2 - R_E^2)) ds
    %
    % Birinci derece baskın (b << a^2/eps): Eq. 6.6'ya indirger
    %   phi_tail = (1/a) * sqrt(2*epsilon/R_E)
    %
    % Ikinci derece baskın (a << sqrt(b*eps)): izotensoid durumu
    %   phi_tail = (2/sqrt(R_E*b)) * arcsinh(sqrt(b*delta_eps/(2*a)))
    %   a->0 limiti: phi_tail = (1/sqrt(R_E*b)) * log(2*b*delta_eps/epsilon)

    if b < 1e-12 && a > 1e-12
        % Saf birinci derece (standart dome'lar)
        delta_eps = epsilon / a;
        phi_tail = (1 / a) * sqrt(2 * epsilon / R_E);
    elseif a < 1e-12 && b > 1e-12
        % Saf ikinci derece (izotensoid donus noktasinda a~0)
        delta_eps = sqrt(2 * epsilon / b);
        phi_tail = (1 / sqrt(R_E * b)) * acosh(1 + epsilon / R_E);
        % Kucuk epsilon/R_E icin: acosh(1+x) ~ sqrt(2x), dolayisiyla
        % phi_tail ~ sqrt(2*epsilon) / (sqrt(R_E) * sqrt(R_E*b)) = sqrt(2*eps/(R_E^2*b))
    else
        % Genel birlesik formul
        discriminant = a^2 + 2 * b * epsilon;
        delta_eps = 2 * epsilon / (a + sqrt(discriminant));

        % phi_tail = (2/sqrt(R_E*b)) * arcsinh(sqrt(b*delta_eps/(2*a)))
        % ancak a->0 icin kararsiz; argumani farkli yazalim:
        %   b*delta_eps = b * 2*eps / (a + sqrt(a^2+2*b*eps))
        % Numerik integral ile dogrudan hesapla (kuyruk bolgesinde rho ~ R_E + eps)
        if b > 1e-12
            % Arcsinh formulu: delta_eps biliniyor, kuyruk integrali analitik
            % rho(u) = R_E + a*u + (b/2)*u^2, u = s_turn - s, u in [0, delta_eps]
            % dphi = R_E / (rho * sqrt(rho^2 - R_E^2)) ds
            % rho^2 - R_E^2 ~ 2*R_E*(rho - R_E) = 2*R_E*(a*u + b*u^2/2)
            % sqrt(rho^2-R_E^2) ~ sqrt(2*R_E) * sqrt(a*u + b*u^2/2)
            % rho ~ R_E (kuyruk bolgesi)
            % dphi/du ~ 1/sqrt(2*R_E) * 1/sqrt(a*u + b*u^2/2)
            %         = 1/sqrt(2*R_E) * 1/sqrt(u*(a + b*u/2))
            % Substitution t = sqrt(u):
            % integral = 2/sqrt(2*R_E) * integral_0^{sqrt(delta_eps)} dt / sqrt(a + b*t^2/2)
            %          = 2/sqrt(R_E*b) * arcsinh(sqrt(b*delta_eps/(2*a)))   [a>0]
            if a > 1e-14
                arg = sqrt(b * delta_eps / (2 * a));
                phi_tail = (2 / sqrt(R_E * b)) * asinh(arg);
            else
                % a ~ 0: integral = 2/sqrt(R_E*b) * integral dt/t (log sonuc)
                % = (1/sqrt(R_E*b)) * log(2*b*delta_eps / (2*a + b*delta_eps) * ??? )
                % Guvenli: rho(u) = R_E + (b/2)*u^2
                % phi = integral_0^de 1/(sqrt(2*R_E)*sqrt(b/2)*u) du
                %     = 1/sqrt(R_E*b) * log(delta_eps / u_min)
                % u_min -> 0 kuyruk ucu, ama fizikte s_turn'de rho=R_E tam
                % acosh formulu daha dogru:
                phi_tail = (1 / sqrt(R_E * b)) * acosh(1 + epsilon / R_E);
            end
        else
            % b ~ 0, a ~ 0: her ikisi de cok kucuk — dejenere durum
            warning('dome_phi_integration:dejenerDonus', ...
                    '|drho/ds|=%.2e ve |d2rho/ds2|=%.2e donus noktasinda cok kucuk.', a, b);
            % Alt sinir ile birinci derece formul kullan
            a_safe = max(a, 1e-8);
            delta_eps = epsilon / a_safe;
            phi_tail = (1 / a_safe) * sqrt(2 * epsilon / R_E);
        end
    end

    % Toplam dome phi
    phi_dome = phi_numerical + phi_tail;

    %% --- Ekvator winding acisi ---
    alpha_eq = asin(R_E / R_eq);

    %% --- Cikti struct ---
    result.phi_dome      = phi_dome;
    result.phi_numerical = phi_numerical;
    result.phi_tail      = phi_tail;
    result.s_epsilon     = s_epsilon;
    result.s_turn        = s_turn;
    result.delta_eps     = delta_eps;
    result.drho_turn     = a;          % |drho/ds| at s_epsilon
    result.d2rho_turn    = b;          % |d2rho/ds2| at s_epsilon
    result.s_ode         = s_ode;
    result.phi_ode       = phi_ode;
    result.alpha_eq      = alpha_eq;

end
