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

clear; clc; close all;
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
L_cyl  = 300.0;
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

    %% === Dome geodesic yol hesabi ===
    rho_d  = dome.rho;
    s_d    = dome.s;
    x_d    = dome.x_local;
    h_dome = x_d(end);

    % Donus noktasi (rho = R_E)
    idx_turn = find(rho_d <= R_E, 1, 'first');
    if isempty(idx_turn), idx_turn = length(rho_d); end

    % Forward yol (ekvator -> donus)
    rho_fwd = rho_d(1:idx_turn);
    s_fwd   = s_d(1:idx_turn);
    x_fwd   = x_d(1:idx_turn);

    % phi birikimi (sayisal)
    rho2_diff = max(rho_fwd.^2 - R_E^2, 1e-10);
    integrand_phi = R_E ./ (rho_fwd .* sqrt(rho2_diff));
    phi_fwd_raw = cumtrapz(s_fwd, integrand_phi);

    % Normalize: geo.phi_dome daha dogru (tail correction dahil)
    phi_dome_half = geo.phi_dome;
    if phi_fwd_raw(end) > 0
        phi_fwd = phi_fwd_raw * (phi_dome_half / phi_fwd_raw(end));
    else
        phi_fwd = phi_fwd_raw;
    end

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

    %% === Fiber uzunlugu ===
    cy_ref = circuit_rho .* cos(circuit_phi);
    cz_ref = circuit_rho .* sin(circuit_phi);
    ds_fib = sqrt(diff(circuit_x).^2 + diff(cy_ref).^2 + diff(cz_ref).^2);
    fiber_one = sum(ds_fib);
    fiber_total_m = p_pat * fiber_one / 1000;
    fprintf('  Fiber: %.1f m (tek devre %.1f mm)\n', fiber_total_m, fiber_one);

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
    x_turn_val = x_d(idx_turn);
    plot3( (L_cyl_adj/2 + x_turn_val) * ones(N_th,1), ...
          R_E*cos(th_mesh), R_E*sin(th_mesh), ...
          ':', 'Color', [0.8 0.3 0.3], 'LineWidth', 0.8);
    plot3(-(L_cyl_adj/2 + x_turn_val) * ones(N_th,1), ...
          R_E*cos(th_mesh), R_E*sin(th_mesh), ...
          ':', 'Color', [0.8 0.3 0.3], 'LineWidth', 0.8);

    % Polar aciklik cemberleri (rho = r0)
    plot3( (L_cyl_adj/2 + h_dome) * ones(N_th,1), ...
          r0*cos(th_mesh), r0*sin(th_mesh), ...
          '-', 'Color', [0.3 0.3 0.8], 'LineWidth', 0.5);
    plot3(-(L_cyl_adj/2 + h_dome) * ones(N_th,1), ...
          r0*cos(th_mesh), r0*sin(th_mesh), ...
          '-', 'Color', [0.3 0.3 0.8], 'LineWidth', 0.5);

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

    %% === Kaydet ===
    out_file = fullfile(fig_dir, sprintf('full_coverage_%s.png', dome_types{dt}));
    exportgraphics(fig, out_file, 'Resolution', 200);
    fprintf('  Kaydedildi: %s\n\n', out_file);
    close(fig);
end

fprintf('=== Tamamlandi: 3 full coverage figure kaydedildi ===\n');
