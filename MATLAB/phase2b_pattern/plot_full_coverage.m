%% PLOT_FULL_COVERAGE  Phase-2b tam kaplama sarim gorsellestirme — TEST-02 x 3 dome tipi.
%
% Her dome tipi icin (ellipsoidal k=0.7, hemispherical, isotensoid):
%   Sol panel: 3D mandrel yuzeyi + tum devreler (renk = devre no)
%   Sag panel: Sol ve sag dome star pattern (sarim sirasi baglantilari)
%
% En kucuk p pattern kullanilir.
% Baslik: dome tipi, devre sayisi, toplam fiber uzunlugu
%
% Cikti: docs/figures/full_coverage_{dome_type}.png (3 dosya)
%
% Parametreler: TEST-02 (Endustriyel COPV)
%   R_eq=152.4, r0=45, L_cyl=300, BW_eff=10, k=0.7
%
% Tarih: 2026-03-16
% Faz: Phase-2b S2

clear; clc; %close all;
fprintf('=== Phase-2b Full Coverage Gorselleri ===\n\n');

%% --- Yol bagimliliklari ---
script_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, '..', 'phase1a_geometry'));
addpath(fullfile(script_dir, '..', 'phase2a_winding'));

fig_dir = fullfile(script_dir, '..', '..', 'docs', 'figures');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

%% --- TEST-02 parametreleri ---
R_eq   = 152.4;
r0     = 45.0;
L_cyl  = 1000.0;
BW_eff = 10.0;
k_ell  = 0.7;
d      = 1;
coverage_range = [100, 150];

R_E      = r0 + BW_eff / 2;
alpha_eq = asin(R_E / R_eq);

%% --- Dome tipleri ---
dome_types  = {'ellipsoidal', 'hemispherical', 'isotensoid'};
dome_labels = {'Ellipsoidal (k=0.7)', 'Hemispherical', 'Isotensoid'};

%% --- Ana dongu ---
for dt = 1:3
    fprintf('--- %s ---\n', dome_labels{dt});

    % Dome profili
    switch dome_types{dt}
        case 'ellipsoidal'
            dome = ellipsoidal_dome_profile(R_eq, r0, k_ell, 1000);
        case 'hemispherical'
            dome = hemispherical_dome_profile(R_eq, r0, 1000);
        case 'isotensoid'
            dome = isotensoid_dome_profile(R_eq, r0, 1000);
    end

    % Geodesic hesabi
    geo = geodesic_single_circuit(dome, R_eq, r0, L_cyl, BW_eff);

    % Pattern arama
    patterns = find_compatible_patterns(geo.delta_phi_circuit, alpha_eq, ...
                                         R_eq, BW_eff, L_cyl, d, coverage_range);

    if isempty(patterns)
        fprintf('  Pattern bulunamadi, atlaniyor.\n\n');
        continue;
    end

    % En kucuk p pattern
    pat = patterns(1);
    p_pat = pat.p;
    q_pat = pat.q;
    L_cyl_adj = pat.L_cyl_adj;

    fprintf('  p=%d, q=%d, coverage=%.1f%%, L_cyl_adj=%.2f mm\n', ...
            p_pat, q_pat, pat.coverage_pct, L_cyl_adj);

    %% === Dome geodesic yol hesabi (Phase-2a composite approach) ===
    % Singularite-free: numerical (0→s_cut) + analytical tail (s_cut→s_turn)
    rho_d  = dome.rho;
    s_d    = dome.s;
    x_d    = dome.x_local;
    h_dome = x_d(end);

    % Donus noktasi s_turn (geo ciktisi — rho = R_E noktasi)
    s_turn = geo.s_turn;
    phi_dome_half = geo.phi_dome;

    % --- (1) Numerical region: s = 0 to s_cut (rho safely above R_E) ---
    epsilon_phi = 1.0;   % mm — singularity safety margin
    rho_cut = R_E + epsilon_phi;

    % s_cut: s noktasi where rho(s_cut) = R_E + epsilon_phi
    [rho_sorted, sort_idx] = sort(rho_d, 'ascend');
    s_sorted = s_d(sort_idx);
    [rho_uniq, uniq_idx] = unique(rho_sorted, 'stable');
    s_uniq = s_sorted(uniq_idx);
    s_cut = interp1(rho_uniq, s_uniq, rho_cut, 'pchip');

    N_num = 500;
    s_num = linspace(0, s_cut, N_num)';
    rho_num = interp1(s_d, rho_d, s_num, 'pchip');
    x_num   = interp1(s_d, x_d,   s_num, 'pchip');
    rho_num(1) = R_eq;  % ekvator sinir

    rho2_diff = rho_num.^2 - R_E^2;    % tumu pozitif (rho > R_E + 1mm)
    integrand_phi = R_E ./ (rho_num .* sqrt(rho2_diff));
    phi_num = cumtrapz(s_num, integrand_phi);
    phi_at_cut = phi_num(end);

    % --- (2) Analytical tail: s_cut to s_turn ---
    % Phase-2a Eq. 6.6: phi(s) ~ phi_cut + phi_remaining * sqrt(t_norm)
    phi_remaining = phi_dome_half - phi_at_cut;

    N_tail = 100;
    s_tail = linspace(s_cut, s_turn, N_tail + 1)';
    s_tail = s_tail(2:end);   % ilki zaten s_num'da
    t_norm = (s_tail - s_cut) / (s_turn - s_cut + 1e-30);
    phi_tail_vals = phi_at_cut + phi_remaining * sqrt(t_norm);

    rho_tail = interp1(s_d, rho_d, s_tail, 'pchip');
    rho_tail = max(rho_tail, R_E);   % clamp: turnaround altina dusme
    x_tail   = interp1(s_d, x_d,   s_tail, 'pchip');

    % --- (3) Combine: s_fwd [0, s_turn], phi_fwd [0, phi_dome_half] ---
    s_fwd   = [s_num; s_tail];
    rho_fwd = [rho_num; rho_tail];
    x_fwd   = [x_num; x_tail];
    phi_fwd = [phi_num; phi_tail_vals];

    % Kesin sinir degerleri
    rho_fwd(1) = R_eq;   x_fwd(1) = 0;
    rho_fwd(end) = R_E;  phi_fwd(end) = phi_dome_half;

    fprintf('  Dome yol nokta: %d (num=%d + tail=%d), s_turn=%.2f mm\n', ...
            length(s_fwd), N_num, N_tail, s_turn);

    % Silindir phi (L_cyl_adj ile)
    phi_cyl = L_cyl_adj * tan(alpha_eq) / R_eq;

    %% === Tek devre referans yolu (6 segment) ===
    % Devre: sag dome -> silindir sol -> sol dome -> silindir sag
    N_cyl_pts = 80;

    % Segment 1: Sag dome forward (ekvator -> donus)
    seg1_x   = L_cyl_adj/2 + x_fwd;
    seg1_rho = rho_fwd;
    seg1_phi = phi_fwd;

    % Segment 2: Sag dome return (donus -> ekvator)
    seg2_x   = flipud(seg1_x);
    seg2_rho = flipud(seg1_rho);
    seg2_phi = 2*phi_dome_half - flipud(phi_fwd);

    % Segment 3: Silindir sol (sag ekvator -> sol ekvator)
    phi_off3 = 2 * phi_dome_half;
    seg3_x   = linspace(L_cyl_adj/2, -L_cyl_adj/2, N_cyl_pts)';
    seg3_rho = R_eq * ones(N_cyl_pts, 1);
    seg3_phi = linspace(phi_off3, phi_off3 + phi_cyl, N_cyl_pts)';

    % Segment 4: Sol dome forward (ekvator -> donus)
    phi_off4 = 2*phi_dome_half + phi_cyl;
    seg4_x   = -(L_cyl_adj/2 + x_fwd);
    seg4_rho = rho_fwd;
    seg4_phi = phi_off4 + phi_fwd;

    % Segment 5: Sol dome return (donus -> ekvator)
    seg5_x   = flipud(seg4_x);
    seg5_rho = flipud(seg4_rho);
    seg5_phi = phi_off4 + 2*phi_dome_half - flipud(phi_fwd);

    % Segment 6: Silindir sag (sol ekvator -> sag ekvator)
    phi_off6 = 4*phi_dome_half + phi_cyl;
    seg6_x   = linspace(-L_cyl_adj/2, L_cyl_adj/2, N_cyl_pts)';
    seg6_rho = R_eq * ones(N_cyl_pts, 1);
    seg6_phi = linspace(phi_off6, phi_off6 + phi_cyl, N_cyl_pts)';

    % Birlestir (tekrar eden uc noktalari atla)
    circuit_x   = [seg1_x;   seg2_x(2:end);   seg3_x(2:end); ...
                   seg4_x(2:end); seg5_x(2:end); seg6_x(2:end)];
    circuit_rho = [seg1_rho; seg2_rho(2:end); seg3_rho(2:end); ...
                   seg4_rho(2:end); seg5_rho(2:end); seg6_rho(2:end)];
    circuit_phi = [seg1_phi; seg2_phi(2:end); seg3_phi(2:end); ...
                   seg4_phi(2:end); seg5_phi(2:end); seg6_phi(2:end)];

    delta_phi_circuit = circuit_phi(end);  % = 4*phi_dome + 2*phi_cyl

    %% === Smooth resampling: 500 esit aralikli noktaya pchip ===
    % Kumulatif yay uzunlugu parametresi
    dx = diff(circuit_x);
    drho = diff(circuit_rho);
    dphi_seg = diff(circuit_phi);
    rho_avg = (circuit_rho(1:end-1) + circuit_rho(2:end)) / 2;
    ds_arc = sqrt(dx.^2 + drho.^2 + (rho_avg .* dphi_seg).^2);
    s_cum = [0; cumsum(ds_arc)];
    s_total_circuit = s_cum(end);

    N_resample = 1500;
    s_uniform = linspace(0, s_total_circuit, N_resample)';

    circuit_x   = interp1(s_cum, circuit_x,   s_uniform, 'pchip');
    circuit_rho = interp1(s_cum, circuit_rho, s_uniform, 'pchip');
    circuit_phi = interp1(s_cum, circuit_phi, s_uniform, 'pchip');

    % Sinir degerlerini koru
    circuit_phi(end) = delta_phi_circuit;

    fprintf('  Circuit resampled: %d -> %d nokta (pchip, s_total=%.1f mm)\n', ...
            length(s_cum), N_resample, s_total_circuit);

    %% === Fiber uzunlugu ===
    cy_ref = circuit_rho .* cos(circuit_phi);
    cz_ref = circuit_rho .* sin(circuit_phi);
    ds_fib = sqrt(diff(circuit_x).^2 + diff(cy_ref).^2 + diff(cz_ref).^2);
    fiber_one = sum(ds_fib);
    fiber_total_m = p_pat * fiber_one / 1000;
    fprintf('  Fiber: %.1f m (tek devre %.1f mm)\n', fiber_total_m, fiber_one);

    %% === DOGRULAMA BLOGU ===
    fprintf('\n  --- Dogrulama ---\n');

    % (A) Segment sinir phi degerleri ve sureklilik
    bnd = struct();
    bnd.seg1_end   = seg1_phi(end);
    bnd.seg2_start = seg2_phi(1);
    bnd.seg2_end   = seg2_phi(end);
    bnd.seg3_start = seg3_phi(1);
    bnd.seg3_end   = seg3_phi(end);
    bnd.seg4_start = seg4_phi(1);
    bnd.seg4_end   = seg4_phi(end);
    bnd.seg5_start = seg5_phi(1);
    bnd.seg5_end   = seg5_phi(end);
    bnd.seg6_start = seg6_phi(1);
    bnd.seg6_end   = seg6_phi(end);

    fprintf('  phi segment sinirlari [rad] / [deg]:\n');
    fprintf('    seg1 end  : %10.6f (%7.2f deg)\n', bnd.seg1_end, rad2deg(bnd.seg1_end));
    fprintf('    seg2 start: %10.6f (%7.2f deg)  gap: %.2e\n', bnd.seg2_start, rad2deg(bnd.seg2_start), abs(bnd.seg2_start - bnd.seg1_end));
    fprintf('    seg2 end  : %10.6f (%7.2f deg)\n', bnd.seg2_end, rad2deg(bnd.seg2_end));
    fprintf('    seg3 start: %10.6f (%7.2f deg)  gap: %.2e\n', bnd.seg3_start, rad2deg(bnd.seg3_start), abs(bnd.seg3_start - bnd.seg2_end));
    fprintf('    seg3 end  : %10.6f (%7.2f deg)\n', bnd.seg3_end, rad2deg(bnd.seg3_end));
    fprintf('    seg4 start: %10.6f (%7.2f deg)  gap: %.2e\n', bnd.seg4_start, rad2deg(bnd.seg4_start), abs(bnd.seg4_start - bnd.seg3_end));
    fprintf('    seg4 end  : %10.6f (%7.2f deg)\n', bnd.seg4_end, rad2deg(bnd.seg4_end));
    fprintf('    seg5 start: %10.6f (%7.2f deg)  gap: %.2e\n', bnd.seg5_start, rad2deg(bnd.seg5_start), abs(bnd.seg5_start - bnd.seg4_end));
    fprintf('    seg5 end  : %10.6f (%7.2f deg)\n', bnd.seg5_end, rad2deg(bnd.seg5_end));
    fprintf('    seg6 start: %10.6f (%7.2f deg)  gap: %.2e\n', bnd.seg6_start, rad2deg(bnd.seg6_start), abs(bnd.seg6_start - bnd.seg5_end));
    fprintf('    seg6 end  : %10.6f (%7.2f deg)\n', bnd.seg6_end, rad2deg(bnd.seg6_end));

    % (B) Phi monoton artan mi?
    dphi = diff(circuit_phi);
    n_neg = sum(dphi < 0);
    n_zero = sum(dphi == 0);
    fprintf('  phi monotonluk: dphi<0 = %d, dphi=0 = %d, toplam = %d\n', ...
            n_neg, n_zero, length(dphi));
    if n_neg > 0
        fprintf('  *** BUG: phi monoton DEGIL! %d azalan segment ***\n', n_neg);
        [min_dphi, min_idx] = min(dphi);
        fprintf('  *** min(dphi) = %.6e rad, index = %d ***\n', min_dphi, min_idx);
    else
        fprintf('  OK: phi kesinlikle monoton artan.\n');
    end

    % (C) delta_phi_circuit vs ideal 2*pi*q/p
    ideal_dphi = 2*pi*q_pat/p_pat;
    dphi_err = abs(delta_phi_circuit - ideal_dphi);
    fprintf('  delta_phi_circuit = %.10f rad (%.4f deg)\n', delta_phi_circuit, rad2deg(delta_phi_circuit));
    fprintf('  ideal 2*pi*q/p   = %.10f rad (%.4f deg)\n', ideal_dphi, rad2deg(ideal_dphi));
    fprintf('  hata             = %.4e rad (%.4e deg)\n', dphi_err, rad2deg(dphi_err));

    % (D) phi_dome ve phi_cyl
    fprintf('  phi_dome_half = %.6f rad (%.2f deg)\n', phi_dome_half, rad2deg(phi_dome_half));
    fprintf('  phi_cyl       = %.6f rad (%.2f deg)\n', phi_cyl, rad2deg(phi_cyl));
    fprintf('  4*phi_dome + 2*phi_cyl = %.10f rad\n', 4*phi_dome_half + 2*phi_cyl);

    % (E) sqrt(Y^2+Z^2) = rho kontrolu (tum devreler icin)
    max_radial_err = 0;
    for i = 1:p_pat
        phi_i = circuit_phi + (i-1) * delta_phi_circuit;
        y_i = circuit_rho .* cos(phi_i);
        z_i = circuit_rho .* sin(phi_i);
        radial_dist = sqrt(y_i.^2 + z_i.^2);
        radial_err = max(abs(radial_dist - circuit_rho));
        if radial_err > max_radial_err
            max_radial_err = radial_err;
        end
    end
    fprintf('  max|sqrt(Y^2+Z^2) - rho| = %.4e mm (tum %d devre)\n', ...
            max_radial_err, p_pat);

    % (F) Polar aciklik kontrolu
    rho_min_circuit = min(circuit_rho);
    fprintf('  rho_min = %.4f mm, R_E = %.1f mm, r0 = %.1f mm\n', ...
            rho_min_circuit, R_E, r0);
    if rho_min_circuit < r0
        fprintf('  *** BUG: Fiber polar acikliga giriyor! rho_min < r0 ***\n');
    elseif rho_min_circuit < R_E - 0.1
        fprintf('  *** UYARI: rho_min < R_E (beklenen donus noktasi) ***\n');
    else
        fprintf('  OK: rho_min >= R_E, polar aciklik icine girmiyor.\n');
    end

    % (G) Dome-silindir gecis noktalari pozisyon uyumu
    % seg2 end -> seg3 start (sag ekvator)
    seg2_end_xyz = [seg2_x(end), seg2_rho(end)*cos(seg2_phi(end)), seg2_rho(end)*sin(seg2_phi(end))];
    seg3_start_xyz = [seg3_x(1), seg3_rho(1)*cos(seg3_phi(1)), seg3_rho(1)*sin(seg3_phi(1))];
    gap23 = norm(seg2_end_xyz - seg3_start_xyz);

    % seg3 end -> seg4 start (sol ekvator)
    seg3_end_xyz = [seg3_x(end), seg3_rho(end)*cos(seg3_phi(end)), seg3_rho(end)*sin(seg3_phi(end))];
    seg4_start_xyz = [seg4_x(1), seg4_rho(1)*cos(seg4_phi(1)), seg4_rho(1)*sin(seg4_phi(1))];
    gap34 = norm(seg3_end_xyz - seg4_start_xyz);

    % seg5 end -> seg6 start (sol ekvator)
    seg5_end_xyz = [seg5_x(end), seg5_rho(end)*cos(seg5_phi(end)), seg5_rho(end)*sin(seg5_phi(end))];
    seg6_start_xyz = [seg6_x(1), seg6_rho(1)*cos(seg6_phi(1)), seg6_rho(1)*sin(seg6_phi(1))];
    gap56 = norm(seg5_end_xyz - seg6_start_xyz);

    fprintf('  Gecis noktasi 3D gap: seg2-3=%.4e, seg3-4=%.4e, seg5-6=%.4e mm\n', ...
            gap23, gap34, gap56);
    if max([gap23, gap34, gap56]) > 0.01
        fprintf('  *** BUG: Segment gecis noktalarinda %f mm boslik! ***\n', max([gap23, gap34, gap56]));
    else
        fprintf('  OK: Tum gecis noktalari surekli (<%0.2e mm).\n', max([gap23, gap34, gap56]));
    end

    fprintf('  --- Dogrulama sonu ---\n\n');

    %% === 3D Mandrel mesh ===
    N_th = 80;
    th_mesh = linspace(0, 2*pi, N_th)';

    % Silindir
    N_cm = 40;
    x_cm = linspace(-L_cyl_adj/2, L_cyl_adj/2, N_cm);
    X_cm = repmat(x_cm, N_th, 1);
    Y_cm = R_eq * cos(th_mesh) * ones(1, N_cm);
    Z_cm = R_eq * sin(th_mesh) * ones(1, N_cm);

    % Sag dome mesh
    N_dm = 50;
    idx_dm = round(linspace(1, length(rho_d), N_dm));
    x_dm   = x_d(idx_dm);
    rho_dm = rho_d(idx_dm);
    X_dr = repmat((L_cyl_adj/2 + x_dm)', N_th, 1);
    Y_dr = cos(th_mesh) * rho_dm';
    Z_dr = sin(th_mesh) * rho_dm';

    % Sol dome mesh (ayna)
    X_dl = fliplr(-X_dr);
    Y_dl = fliplr(Y_dr);
    Z_dl = fliplr(Z_dr);

    %% === FIGURE ===
    fig = figure('Position', [50, 50, 1600, 700], 'Color', 'w');

    %% --- Sol panel: 3D mandrel + sarim yollari ---
    ax1 = subplot(1, 2, 1);
    hold on;

    % Mandrel yuzeyi (yari saydam)
    surf(X_cm, Y_cm, Z_cm, 'FaceColor', [0.85 0.85 0.85], ...
         'FaceAlpha', 0.20, 'EdgeColor', 'none');
    surf(X_dr, Y_dr, Z_dr, 'FaceColor', [0.78 0.83 0.90], ...
         'FaceAlpha', 0.20, 'EdgeColor', 'none');
    surf(X_dl, Y_dl, Z_dl, 'FaceColor', [0.78 0.83 0.90], ...
         'FaceAlpha', 0.20, 'EdgeColor', 'none');

    % Ekvator cemberleri (referans)
    plot3( L_cyl_adj/2 * ones(N_th,1), R_eq*cos(th_mesh), R_eq*sin(th_mesh), ...
          '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5);
    plot3(-L_cyl_adj/2 * ones(N_th,1), R_eq*cos(th_mesh), R_eq*sin(th_mesh), ...
          '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5);

    % Donus cemberleri (rho = R_E)
    x_turn_val = x_fwd(end);   % dome x at turnaround
    plot3( (L_cyl_adj/2 + x_turn_val) * ones(N_th,1), ...
          R_E*cos(th_mesh), R_E*sin(th_mesh), ...
          ':', 'Color', [0.8 0.3 0.3], 'LineWidth', 0.8);
    plot3(-(L_cyl_adj/2 + x_turn_val) * ones(N_th,1), ...
          R_E*cos(th_mesh), R_E*sin(th_mesh), ...
          ':', 'Color', [0.8 0.3 0.3], 'LineWidth', 0.8);

    % Sarim yollari (her devre farkli renk)
    cmap_c = turbo(max(p_pat, 2));
    for i = 1:p_pat
        phi_i = circuit_phi + (i-1) * delta_phi_circuit;
        y_i = circuit_rho .* cos(phi_i);
        z_i = circuit_rho .* sin(phi_i);
        plot3(circuit_x, y_i, z_i, '-', 'Color', [cmap_c(i,:) 0.85], 'LineWidth', 1.2);
    end

    axis equal; grid on; box on;
    xlabel('X [mm]', 'FontSize', 11);
    ylabel('Y [mm]', 'FontSize', 11);
    zlabel('Z [mm]', 'FontSize', 11);
    title(sprintf('%s — %d Circuits — %.1f m fiber\np=%d, q=%d, coverage=%.0f%%', ...
          dome_labels{dt}, p_pat, fiber_total_m, p_pat, q_pat, pat.coverage_pct), ...
          'FontSize', 12);
    view(30, 20);
    camlight('headlight');
    lighting gouraud;
    material dull;

    colormap(ax1, turbo(max(p_pat, 2)));
    cb1 = colorbar(ax1);
    cb1.Label.String = 'Devre No';
    cb1.Label.FontSize = 10;
    clim(ax1, [1 p_pat]);
    set(ax1, 'FontSize', 10);

    %% --- Sag panel: Star pattern (sol + sag dome) ---
    ax2 = subplot(1, 2, 2);
    hold on; axis equal; box on;

    % Ideal delta_phi (tam kapanma icin)
    delta_phi_ideal = 2*pi*q_pat/p_pat;

    % Donus noktasi acilari (her devre icin)
    phi_right = zeros(1, p_pat);
    phi_left  = zeros(1, p_pat);
    for i = 1:p_pat
        phi_start_i = (i-1) * delta_phi_ideal;
        phi_right(i) = mod(phi_start_i + phi_dome_half, 2*pi);
        phi_left(i)  = mod(phi_start_i + 3*phi_dome_half + phi_cyl, 2*pi);
    end

    % Cember merkez konumlari
    x_off_R =  1.4;   % Sag dome
    x_off_L = -1.4;   % Sol dome

    theta_c = linspace(0, 2*pi, 200)';

    % Referans cemberler
    plot(x_off_R + cos(theta_c), sin(theta_c), '-', ...
         'Color', [0.75 0.75 0.75], 'LineWidth', 1.2);
    plot(x_off_L + cos(theta_c), sin(theta_c), '-', ...
         'Color', [0.75 0.75 0.75], 'LineWidth', 1.2);

    % Ic cember (donus capi orani: R_E/R_eq)
    r_inner = R_E / R_eq;
    plot(x_off_R + r_inner*cos(theta_c), r_inner*sin(theta_c), ':', ...
         'Color', [0.85 0.85 0.85], 'LineWidth', 0.8);
    plot(x_off_L + r_inner*cos(theta_c), r_inner*sin(theta_c), ':', ...
         'Color', [0.85 0.85 0.85], 'LineWidth', 0.8);

    % Star renk haritasi
    cmap_s = turbo(max(p_pat, 2));

    % --- Sag dome star ---
    % Baglanti cizgileri (sarim sirasinda)
    for i = 1:p_pat
        j = mod(i, p_pat) + 1;
        x1 = x_off_R + cos(phi_right(i)); y1 = sin(phi_right(i));
        x2 = x_off_R + cos(phi_right(j)); y2 = sin(phi_right(j));
        plot([x1 x2], [y1 y2], '-', 'Color', [cmap_s(i,:) 0.55], 'LineWidth', 1.5);
    end

    % Donus noktasi isaretcileri
    for i = 1:p_pat
        plot(x_off_R + cos(phi_right(i)), sin(phi_right(i)), 'o', ...
             'MarkerSize', 8, 'MarkerFaceColor', cmap_s(i,:), ...
             'MarkerEdgeColor', 'k', 'LineWidth', 0.6);
    end

    % --- Sol dome star ---
    for i = 1:p_pat
        j = mod(i, p_pat) + 1;
        x1 = x_off_L + cos(phi_left(i)); y1 = sin(phi_left(i));
        x2 = x_off_L + cos(phi_left(j)); y2 = sin(phi_left(j));
        plot([x1 x2], [y1 y2], '-', 'Color', [cmap_s(i,:) 0.55], 'LineWidth', 1.5);
    end

    for i = 1:p_pat
        plot(x_off_L + cos(phi_left(i)), sin(phi_left(i)), 'o', ...
             'MarkerSize', 8, 'MarkerFaceColor', cmap_s(i,:), ...
             'MarkerEdgeColor', 'k', 'LineWidth', 0.6);
    end

    % Baslangic noktasi (sari yildiz)
    plot(x_off_R + cos(phi_right(1)), sin(phi_right(1)), 'h', ...
         'MarkerSize', 15, 'MarkerFaceColor', [1 0.8 0], ...
         'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
    plot(x_off_L + cos(phi_left(1)),  sin(phi_left(1)),  'h', ...
         'MarkerSize', 15, 'MarkerFaceColor', [1 0.8 0], ...
         'MarkerEdgeColor', 'k', 'LineWidth', 1.5);

    % Etiketler
    text(x_off_R, -1.35, 'Sag Dome', 'HorizontalAlignment', 'center', ...
         'FontSize', 11, 'FontWeight', 'bold', 'Color', [0.2 0.2 0.6]);
    text(x_off_L, -1.35, 'Sol Dome', 'HorizontalAlignment', 'center', ...
         'FontSize', 11, 'FontWeight', 'bold', 'Color', [0.2 0.2 0.6]);

    % Aci etiketleri (0, 90, 180, 270)
    for ang_deg = [0, 90, 180, 270]
        ang_rad = ang_deg * pi / 180;
        text(x_off_R + 1.18*cos(ang_rad), 1.18*sin(ang_rad), ...
             sprintf('%d', ang_deg), 'HorizontalAlignment', 'center', ...
             'FontSize', 8, 'Color', [0.4 0.4 0.4]);
        text(x_off_L + 1.18*cos(ang_rad), 1.18*sin(ang_rad), ...
             sprintf('%d', ang_deg), 'HorizontalAlignment', 'center', ...
             'FontSize', 8, 'Color', [0.4 0.4 0.4]);
    end

    title(sprintf('Star Pattern (p=%d, q=%d, n=%d)', p_pat, q_pat, 2*p_pat), ...
          'FontSize', 12);
    xlim([-3.1 3.1]); ylim([-1.65 1.65]);
    set(ax2, 'XTick', [], 'YTick', [], 'FontSize', 10);

    colormap(ax2, turbo(max(p_pat, 2)));
    cb2 = colorbar(ax2);
    cb2.Label.String = 'Sarim sirasi';
    cb2.Label.FontSize = 10;
    clim(ax2, [1 p_pat]);

    %% === Kaydet (full coverage) ===
    out_file = fullfile(fig_dir, sprintf('full_coverage_%s.png', dome_types{dt}));
    exportgraphics(fig, out_file, 'Resolution', 200);
    fprintf('  Kaydedildi: %s\n', out_file);
    %close(fig);

    %% ================================================================
    %% === EK FIGURE: Ilk 5 devre (ayri renkler) ===
    %% ================================================================
    n_show = min(5, p_pat);
    colors_5 = [0.0 0.2 0.8;    % 1: mavi
                0.9 0.1 0.1;    % 2: kirmizi
                0.9 0.8 0.0;    % 3: sari
                0.0 0.7 0.2;    % 4: yesil
                0.6 0.1 0.7];   % 5: mor
    color_names = {'Mavi','Kirmizi','Sari','Yesil','Mor'};

    fig5 = figure('Position', [80, 80, 1200, 700], 'Color', 'w');
    hold on;

    % Mandrel yuzeyi
    surf(X_cm, Y_cm, Z_cm, 'FaceColor', [0.88 0.88 0.88], ...
         'FaceAlpha', 0.15, 'EdgeColor', 'none');
    surf(X_dr, Y_dr, Z_dr, 'FaceColor', [0.82 0.86 0.92], ...
         'FaceAlpha', 0.15, 'EdgeColor', 'none');
    surf(X_dl, Y_dl, Z_dl, 'FaceColor', [0.82 0.86 0.92], ...
         'FaceAlpha', 0.15, 'EdgeColor', 'none');

    % Ekvator ve donus cemberleri
    plot3( L_cyl_adj/2 * ones(N_th,1), R_eq*cos(th_mesh), R_eq*sin(th_mesh), ...
          '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.5);
    plot3(-L_cyl_adj/2 * ones(N_th,1), R_eq*cos(th_mesh), R_eq*sin(th_mesh), ...
          '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.5);
    plot3( (L_cyl_adj/2 + x_turn_val) * ones(N_th,1), ...
          R_E*cos(th_mesh), R_E*sin(th_mesh), ...
          ':', 'Color', [0.8 0.3 0.3], 'LineWidth', 0.8);
    plot3(-(L_cyl_adj/2 + x_turn_val) * ones(N_th,1), ...
          R_E*cos(th_mesh), R_E*sin(th_mesh), ...
          ':', 'Color', [0.8 0.3 0.3], 'LineWidth', 0.8);

    % Ilk 5 devre
    leg_entries = gobjects(n_show, 1);
    for i = 1:n_show
        phi_i = circuit_phi + (i-1) * delta_phi_circuit;
        y_i = circuit_rho .* cos(phi_i);
        z_i = circuit_rho .* sin(phi_i);
        leg_entries(i) = plot3(circuit_x, y_i, z_i, '-', ...
            'Color', colors_5(i,:), 'LineWidth', 2.0);
    end

    axis equal; grid on; box on;
    xlabel('X [mm]', 'FontSize', 12);
    ylabel('Y [mm]', 'FontSize', 12);
    zlabel('Z [mm]', 'FontSize', 12);
    title(sprintf('%s — Ilk %d Devre (p=%d, q=%d)', ...
          dome_labels{dt}, n_show, p_pat, q_pat), 'FontSize', 14);
    view(30, 20);
    camlight('headlight');
    lighting gouraud;
    material dull;

    legend(leg_entries, ...
        arrayfun(@(k) sprintf('Devre %d (%s)', k, color_names{k}), ...
                 1:n_show, 'UniformOutput', false), ...
        'Location', 'eastoutside', 'FontSize', 10);

    % Kaydet (first5)
    out5 = fullfile(fig_dir, sprintf('first5_%s.png', dome_types{dt}));
    exportgraphics(fig5, out5, 'Resolution', 200);
    fprintf('  Kaydedildi: %s\n\n', out5);
    %close(fig5);
end

fprintf('=== Tamamlandi: 3 full coverage figure kaydedildi ===\n');
