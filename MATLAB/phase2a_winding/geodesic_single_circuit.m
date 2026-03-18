function [circuit] = geodesic_single_circuit(dome_profile, R_eq, r0, L_cyl, BW_eff, N_points)
% GEODESIC_SINGLE_CIRCUIT  Tek devre geodesic fiber yolu uretimi.
%
% Phase-2a: Tam mandrel uzerinde tek devre (single circuit) geodesic yol.
% Dome-1 donus -> silindir -> dome-2 donus -> geri donus -> dome-1 donus.
%
% Matematiksel temel: docs/phase2a_geodesic_math.md
%   - Clairaut iliskisi (Eq. 4.1): rho*sin(alpha) = R_E
%   - Dome ODE (Eq. 5.5): dphi/ds = R_E / (rho*sqrt(rho^2 - R_E^2))
%   - Singularite yonetimi (Eq. 6.6): analitik kuyruk tamamlama
%   - Silindir heliks (Eq. 7.3): phi_cyl = L_cyl*tan(alpha_eq)/R_eq
%   - Toplam (Eq. 8.1): delta_phi_circuit = 4*phi_dome + 2*phi_cyl
%
% KULLANIM:
%   circuit = geodesic_single_circuit(dome_profile, R_eq, r0, L_cyl, BW_eff)
%   circuit = geodesic_single_circuit(dome_profile, R_eq, r0, L_cyl, BW_eff, N_points)
%
% GIRDILER:
%   dome_profile : Phase-1a profil struct (ellipsoidal/hemispherical/isotensoid)
%   R_eq         : Ekvator yaricapi [mm]
%   r0           : Polar aciklik yaricapi [mm]
%   L_cyl        : Silindir uzunlugu [mm]
%   BW_eff       : Efektif bant genisligi [mm]
%   N_points     : Yol nokta sayisi (varsayilan: 2000)
%
% CIKTI:
%   circuit : Struct — Bolum 15.2 arayuzu
%
% REFERANS: docs/phase2a_geodesic_math.md

    %% --- Varsayilan parametreler ---
    if nargin < 6 || isempty(N_points)
        N_points = 2000;
    end

    %% --- On hesap (Bolum 12.2) ---
    R_E      = r0 + BW_eff / 2;              % Karar-5
    alpha_eq = asin(R_E / R_eq);             % Eq. 4.3

    % Hoop kontrolu (S-WIND-02)
    if alpha_eq > deg2rad(85)
        error('geodesic_single_circuit:hoopKisiti', ...
              'alpha_eq = %.1f deg > 85 deg. Hoop — dome yolu uretilemez.', ...
              rad2deg(alpha_eq));
    end
    % Polar uyari (S-WIND-03)
    if alpha_eq < deg2rad(5)
        warning('geodesic_single_circuit:polarUyari', ...
                'alpha_eq = %.1f deg < 5 deg.', rad2deg(alpha_eq));
    end

    fprintf('=== Geodesic Single Circuit ===\n');
    fprintf('R_eq=%.2f, r0=%.2f, R_E=%.2f, L_cyl=%.1f, BW_eff=%.1f\n', ...
            R_eq, r0, R_E, L_cyl, BW_eff);
    fprintf('alpha_eq = %.4f rad (%.2f deg)\n', alpha_eq, rad2deg(alpha_eq));

    %% --- Dome phi entegrasyonu (Bolum 5-6) ---
    dome_result = dome_phi_integration(dome_profile, R_E);
    phi_dome = dome_result.phi_dome;
    fprintf('phi_dome = %.6f rad (%.2f deg)  [num=%.6f, tail=%.6f]\n', ...
            phi_dome, rad2deg(phi_dome), ...
            dome_result.phi_numerical, dome_result.phi_tail);

    %% --- Silindir phi katkisi (Bolum 7) ---
    phi_cyl = L_cyl * tan(alpha_eq) / R_eq;  % Eq. 7.3
    fprintf('phi_cyl = %.6f rad (%.2f deg)\n', phi_cyl, rad2deg(phi_cyl));

    %% --- Toplam devre (Bolum 8) ---
    delta_phi_circuit = 4 * phi_dome + 2 * phi_cyl;  % Eq. 8.1
    fprintf('delta_phi = %.6f rad (%.2f deg, %.4f tur)\n', ...
            delta_phi_circuit, rad2deg(delta_phi_circuit), ...
            delta_phi_circuit / (2*pi));

    %% --- Dome phi(s) fonksiyonu olustur ---
    % ODE cozumu s in [0, s_eps] verir; s_eps..s_total arasi kuyruk
    s_dome_total = dome_profile.s_total;
    s_eps        = dome_result.s_epsilon;

    % Fiber donus noktasi: rho(s_turn) = R_E (Eq. 4.4, alpha = pi/2)
    % rho monoton azalan -> flip ile ters sorgu
    s_turn = interp1(flipud(dome_profile.rho), flipud(dome_profile.s), R_E, 'pchip');

    % ODE ciktisi: s in [0, s_eps] -> phi in [0, phi_numerical]
    phi_ode_interp = griddedInterpolant(dome_result.s_ode, ...
                                         dome_result.phi_ode, 'pchip');

    % Kuyruk bolgesi: s in [s_eps, s_turn]
    % Eq. 6.6 yaklasimi: phi(s) ~ phi_num + phi_tail * sqrt((s-s_eps)/(s_turn-s_eps))
    % Birlesik phi(s) tablosu olustur
    N_main = 500;
    N_tail_pts = 50;
    s_main = linspace(0, s_eps, N_main)';
    phi_main = phi_ode_interp(s_main);

    s_tail = linspace(s_eps, s_turn, N_tail_pts + 1)';
    s_tail = s_tail(2:end);  % ilki zaten s_eps
    t_norm = (s_tail - s_eps) / (s_turn - s_eps + 1e-30);
    phi_tail_vals = dome_result.phi_numerical + dome_result.phi_tail * sqrt(t_norm);

    s_dome_tab   = [s_main; s_tail];
    phi_dome_tab = [phi_main; phi_tail_vals];

    % phi_dome_func: s -> phi(s) dome icinde [0, s_turn] -> [0, phi_dome]
    phi_dome_func = @(s_q) interp1(s_dome_tab, phi_dome_tab, s_q, 'pchip');

    %% --- Profil interpolantlari ---
    rho_interp   = griddedInterpolant(dome_profile.s, dome_profile.rho, 'pchip');
    x_interp     = griddedInterpolant(dome_profile.s, dome_profile.x_local, 'pchip');
    drho_interp  = griddedInterpolant(dome_profile.s, dome_profile.drho_ds, 'pchip');
    kappa_interp = griddedInterpolant(dome_profile.s, dome_profile.kappa_m, 'pchip');

    %% --- Segment nokta sayilari ---
    s_half = s_turn + L_cyl + s_turn;
    N_half = round(N_points / 2);
    N_d = max(round(N_half * s_turn / s_half), 30);
    N_c = max(N_half - 2 * N_d, 30);

    %% --- Dome uniform gridi ---
    s_d = linspace(0, s_turn, N_d)';
    rho_d = rho_interp(s_d);
    x_d   = x_interp(s_d);
    phi_d = phi_dome_func(s_d);  % 0 -> phi_dome (ekvator -> donus)

    %% ========== TAM DEVRE MONTAJI ==========
    % 6 segment, phi surekli birikiyor (hep artis)
    %
    % FWD-1: Dome-1 donus -> ekvator-1  (s_dome: total->0, phi: 0->phi_dome)
    % FWD-2: Silindir eq1 -> eq2        (x: 0->L_cyl,     phi: ->+phi_cyl)
    % FWD-3: Dome-2 ekvator -> donus    (s_dome: 0->total, phi: ->+phi_dome)
    % RET-4: Dome-2 donus -> ekvator    (s_dome: total->0, phi: ->+phi_dome)
    % RET-5: Silindir eq2 -> eq1        (x: L_cyl->0,     phi: ->+phi_cyl)
    % RET-6: Dome-1 ekvator -> donus    (s_dome: 0->total, phi: ->+phi_dome)

    phi_acc = 0;  % kumulatif phi

    % --- FWD-1: Dome-1 donus -> ekvator ---
    % s_dome goes total->0 (fiber donustan ekvator'a gider)
    % phi(s) dome icinde: s=total'de phi_dome, s=0'da 0
    % Fiber donustan basladiginda phi = 0, ekvator'a vardiginda phi += phi_dome
    seg1_s_dome = flipud(s_d);            % total -> 0
    seg1_rho    = rho_interp(seg1_s_dome);
    seg1_x      = -x_interp(seg1_s_dome); % dome-1: negatif x
    seg1_phi    = phi_acc + (phi_dome - flipud(phi_d));
    % flipud(phi_d) = phi(total), phi(total-ds), ..., phi(0) = phi_dome,..,0
    % phi_dome - flipud(phi_d) = 0, ..., phi_dome  (monoton artan)
    phi_acc = seg1_phi(end);              % = phi_dome

    % --- FWD-2: Silindir ---
    s_c = linspace(0, L_cyl, N_c)';
    seg2_rho = R_eq * ones(N_c, 1);
    seg2_x   = s_c;
    seg2_phi = phi_acc + s_c * tan(alpha_eq) / R_eq;
    phi_acc  = seg2_phi(end);             % phi_dome + phi_cyl

    % --- FWD-3: Dome-2 ekvator -> donus ---
    seg3_s_dome = s_d;                    % 0 -> total
    seg3_rho    = rho_interp(seg3_s_dome);
    seg3_x      = L_cyl + x_interp(seg3_s_dome);
    seg3_phi    = phi_acc + phi_d;
    phi_acc     = seg3_phi(end);          % phi_dome + phi_cyl + phi_dome

    % --- RET-4: Dome-2 donus -> ekvator ---
    seg4_s_dome = flipud(s_d);
    seg4_rho    = rho_interp(seg4_s_dome);
    seg4_x      = L_cyl + x_interp(seg4_s_dome);
    seg4_phi    = phi_acc + (phi_dome - flipud(phi_d));
    phi_acc     = seg4_phi(end);          % 2*phi_dome + phi_cyl + phi_dome

    % --- RET-5: Silindir geri ---
    seg5_rho = R_eq * ones(N_c, 1);
    seg5_x   = flipud(s_c);              % L_cyl -> 0
    seg5_phi = phi_acc + (L_cyl - flipud(s_c)) * tan(alpha_eq) / R_eq;
    phi_acc  = seg5_phi(end);

    % --- RET-6: Dome-1 ekvator -> donus ---
    seg6_s_dome = s_d;
    seg6_rho    = rho_interp(seg6_s_dome);
    seg6_x      = -x_interp(seg6_s_dome);
    seg6_phi    = phi_acc + phi_d;
    phi_acc     = seg6_phi(end);          % should be ~ delta_phi_circuit

    %% --- Birlestir (birlesim noktalarinda tekrar yok) ---
    rho_all = [seg1_rho; seg2_rho(2:end); seg3_rho(2:end); ...
               seg4_rho(2:end); seg5_rho(2:end); seg6_rho(2:end)];
    x_all   = [seg1_x; seg2_x(2:end); seg3_x(2:end); ...
               seg4_x(2:end); seg5_x(2:end); seg6_x(2:end)];
    phi_all = [seg1_phi; seg2_phi(2:end); seg3_phi(2:end); ...
               seg4_phi(2:end); seg5_phi(2:end); seg6_phi(2:end)];

    N_total = numel(rho_all);

    % Global yay uzunlugu (kumulatif)
    drho_d = diff(rho_all);
    dx_d   = diff(x_all);
    dphi_d = diff(phi_all);
    rho_mid = (rho_all(1:end-1) + rho_all(2:end)) / 2;
    ds_inc  = sqrt(drho_d.^2 + (rho_mid .* dphi_d).^2 + dx_d.^2);
    s_all   = [0; cumsum(ds_inc)];

    %% --- Winding acisi (Clairaut, Eq. 4.2) ---
    alpha_all = asin(min(R_E ./ rho_all, 1.0));

    %% --- Normal egrilik (S-WIND-04, Eq. 9.2) ---
    kn_all = compute_normal_curvature(rho_all, x_all, alpha_all, ...
                                       dome_profile, L_cyl, R_eq, ...
                                       rho_interp, drho_interp, kappa_interp);

    bridging_idx = find(kn_all < -1e-10);
    if ~isempty(bridging_idx)
        warning('geodesic_single_circuit:bridging', ...
                'S-WIND-04: %d noktada kn < 0 (bridging riski).', ...
                numel(bridging_idx));
    end

    %% --- Clairaut dogrulama (Karar-11 Katman 2, Eq. 14.1) ---
    delta_c = abs(rho_all .* sin(alpha_all) - R_E);
    clairaut_err = max(delta_c);
    fprintf('Clairaut max sapma: %.2e mm (limit: 1e-4)\n', clairaut_err);

    %% --- Kayma egilimi (S-WIND-05, Eq. 9.4-9.5) ---
    lambda_vals = zeros(N_total, 1);
    valid_kn = kn_all > 1e-15;
    lambda_vals(valid_kn) = delta_c(valid_kn) ./ ...
                            (rho_all(valid_kn) .* kn_all(valid_kn));
    lambda_max = max(lambda_vals);

    %% --- Toplam phi dogrulamasi ---
    phi_check = phi_all(end);
    phi_err   = abs(phi_check - delta_phi_circuit) / delta_phi_circuit;
    fprintf('Phi dogrulama: birikimli=%.6f, teorik=%.6f, bagil hata=%.2e\n', ...
            phi_check, delta_phi_circuit, phi_err);

    %% --- Cikti struct (Bolum 15.2) ---
    circuit.s     = s_all;
    circuit.rho   = rho_all;
    circuit.x     = x_all;
    circuit.phi   = phi_all;
    circuit.alpha = alpha_all;
    circuit.kn    = kn_all;
    circuit.lambda_max         = lambda_max;
    circuit.delta_phi_circuit  = delta_phi_circuit;
    circuit.phi_dome           = phi_dome;
    circuit.phi_cyl            = phi_cyl;
    circuit.bridging_risk_indices = bridging_idx;
    circuit.alpha_eq = alpha_eq;
    circuit.R_E      = R_E;
    circuit.N_total  = N_total;
    circuit.clairaut_max_err = clairaut_err;

    fprintf('Yol uretildi: %d nokta, bridging: %d nokta\n', ...
            N_total, numel(bridging_idx));

    %% ========== 3D GORSELLISTIRME ==========
    plot_geodesic_3d(circuit, dome_profile, R_eq, r0, L_cyl);
end


%% ================================================================
function kn = compute_normal_curvature(rho, x, alpha, dome_profile, ...
                                        L_cyl, R_eq, rho_i, drho_i, kappa_i)
% Normal egrilik hesabi (Eq. 9.2): kn = km*cos^2(a) + kp*sin^2(a)

    N = numel(rho);
    kn = zeros(N, 1);
    s_dome_total = dome_profile.s_total;

    % x_local -> s_dome ters sorgu icin interpolant
    % dome_profile.x_local monoton artan, dome_profile.s monoton artan
    s_from_x = griddedInterpolant(dome_profile.x_local, dome_profile.s, 'pchip');

    for i = 1:N
        aa = alpha(i);

        if x(i) < -1e-6
            % Dome-1 bolgesi (x < 0)
            x_loc = min(abs(x(i)), dome_profile.x_local(end));
            s_loc = s_from_x(x_loc);
            s_loc = max(0, min(s_loc, s_dome_total));
            km = kappa_i(s_loc);
            dr = drho_i(s_loc);
        elseif x(i) > L_cyl + 1e-6
            % Dome-2 bolgesi (x > L_cyl)
            x_loc = min(x(i) - L_cyl, dome_profile.x_local(end));
            s_loc = s_from_x(x_loc);
            s_loc = max(0, min(s_loc, s_dome_total));
            km = kappa_i(s_loc);
            dr = drho_i(s_loc);
        else
            % Silindir bolgesi: km = 0
            km = 0;
            dr = 0;
        end

        % Paralel egrilik (Eq. 9.3)
        rr = rho(i);
        dr2 = dr^2;
        if dr2 < 1 - 1e-12 && rr > 1e-10
            kappa_p = -dr / (rr * sqrt(1 - dr2));
        else
            kappa_p = 0;
        end

        kn(i) = km * cos(aa)^2 + kappa_p * sin(aa)^2;
    end
end


%% ================================================================
function plot_geodesic_3d(circuit, dome_profile, R_eq, r0, L_cyl)
% 3D gorsellistirme: mandrel yuzeyi + tek devre geodesic yolu

    figure('Name', 'Phase-2a: Geodesic Single Circuit', ...
           'Position', [100, 100, 1200, 800], 'Color', 'w');

    %% --- Mandrel yuzeyi ---
    N_phi_surf = 80;
    phi_surf = linspace(0, 2*pi, N_phi_surf);
    N_s = numel(dome_profile.s);
    s_idx = round(linspace(1, N_s, min(N_s, 120)));

    rho_d = dome_profile.rho(s_idx);
    x_d   = dome_profile.x_local(s_idx);

    % Dome-1 yuzeyi (x < 0)
    [PHI_m, RHO_m] = meshgrid(phi_surf, rho_d);
    X_d1 = RHO_m .* cos(PHI_m);
    Y_d1 = RHO_m .* sin(PHI_m);
    Z_d1 = repmat(-x_d, 1, N_phi_surf);

    % Silindir yuzeyi
    N_cs = 20;
    x_cs = linspace(0, L_cyl, N_cs)';
    [PHI_c, ~] = meshgrid(phi_surf, x_cs);
    X_c = R_eq * cos(PHI_c);
    Y_c = R_eq * sin(PHI_c);
    Z_c = repmat(x_cs, 1, N_phi_surf);

    % Dome-2 yuzeyi (x > L_cyl)
    X_d2 = RHO_m .* cos(PHI_m);
    Y_d2 = RHO_m .* sin(PHI_m);
    Z_d2 = repmat(L_cyl + x_d, 1, N_phi_surf);

    hold on;
    s_args = {'FaceAlpha', 0.2, 'EdgeAlpha', 0.08, 'FaceColor', [0.7 0.85 1]};
    surf(X_d1, Y_d1, Z_d1, s_args{:});
    surf(X_c,  Y_c,  Z_c,  s_args{:});
    surf(X_d2, Y_d2, Z_d2, s_args{:});

    %% --- Geodesic yol (Eq. 9.1) ---
    Xp = circuit.rho .* cos(circuit.phi);
    Yp = circuit.rho .* sin(circuit.phi);
    Zp = circuit.x;

    plot3(Xp, Yp, Zp, 'r-', 'LineWidth', 2.0);

    % Baslangic / donus / bitis isaretleri
    plot3(Xp(1), Yp(1), Zp(1), 'go', 'MarkerSize', 12, ...
          'MarkerFaceColor', 'g', 'DisplayName', 'Baslangic (D1 donus)');
    mid = round(numel(Xp)/2);
    plot3(Xp(mid), Yp(mid), Zp(mid), 'bs', 'MarkerSize', 12, ...
          'MarkerFaceColor', 'b', 'DisplayName', 'Dome-2 donus');
    plot3(Xp(end), Yp(end), Zp(end), 'r^', 'MarkerSize', 12, ...
          'MarkerFaceColor', 'r', 'DisplayName', 'Bitis (D1 donus)');

    % Boss daireleri
    th_b = linspace(0, 2*pi, 100);
    xb = r0*cos(th_b); yb = r0*sin(th_b);
    plot3(xb, yb, -dome_profile.x_local(end)*ones(size(th_b)), ...
          'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    plot3(xb, yb, (L_cyl+dome_profile.x_local(end))*ones(size(th_b)), ...
          'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');

    axis equal; grid on;
    xlabel('X [mm]'); ylabel('Y [mm]'); zlabel('Z (mandrel ekseni) [mm]');
    if isfield(dome_profile, 'k')
        dome_label = sprintf('Elipsoidal (k=%.1f)', dome_profile.k);
    else
        dome_label = sprintf('AR=%.2f', dome_profile.aspect_r);
    end
    title(sprintf(['Phase-2a: Geodesic Tek Devre — %s\n' ...
                   'R_{eq}=%.1f, r_0=%.1f, L_{cyl}=%.0f, ' ...
                   '\\alpha_{eq}=%.1f\\circ, ' ...
                   '\\Delta\\phi=%.1f\\circ (%.3f tur)'], ...
                  dome_label, R_eq, r0, L_cyl, ...
                  rad2deg(circuit.alpha_eq), ...
                  rad2deg(circuit.delta_phi_circuit), ...
                  circuit.delta_phi_circuit / (2*pi)));
    legend('Location', 'best');
    view(35, 25);
    hold off;

    %% --- Alt grafik: alpha ve phi profilleri ---
    figure('Name', 'Phase-2a: Winding Profiles', ...
           'Position', [150, 150, 1000, 400], 'Color', 'w');

    subplot(1,2,1);
    plot(circuit.s, rad2deg(circuit.alpha), 'b-', 'LineWidth', 1.2);
    xlabel('s [mm]'); ylabel('\alpha [deg]');
    title('Winding Acisi');
    grid on;
    yline(rad2deg(circuit.alpha_eq), 'r--', ...
          sprintf('\\alpha_{eq}=%.1f\\circ', rad2deg(circuit.alpha_eq)));

    subplot(1,2,2);
    plot(circuit.s, rad2deg(circuit.phi), 'b-', 'LineWidth', 1.2);
    xlabel('s [mm]'); ylabel('\phi [deg]');
    title('Paralel Aci Ilerlemesi');
    grid on;
end
