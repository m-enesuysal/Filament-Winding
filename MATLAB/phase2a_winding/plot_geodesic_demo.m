%% PLOT_GEODESIC_DEMO  Phase-2a geodesic yol goruntusu — 3 dome tipi.
%
% S1 senaryosu (R_eq=152.4, r0=45, L_cyl=300, BW_eff=10) icin referans
% CSV dosyalarini okur ve her dome tipi icin uc figure olusturur:
%   Figure 1 (2 panel): Sol — rho vs x, Sag — alpha vs phi
%   Figure 2 (3D):      Mandrel yuzeyi + geodesik yol
%
% Cikti:
%   docs/figures/geodesic_demo_{dome}.png   (2 panelli meridyen + sarim acisi)
%   docs/figures/geodesic_3d_{dome}.png     (3D mandrel + geodesik yol)
%
% Tarih: 2026-03-14
% Faz: Phase-2a S6

clear; clc; %close all;

%% --- Yapilandirma ---
this_dir  = fileparts(mfilename('fullpath'));
ref_dir   = fullfile(this_dir, 'reference_data');
fig_dir   = fullfile(this_dir, '..', '..', 'docs', 'figures');

if ~exist(fig_dir, 'dir')
    mkdir(fig_dir);
end

% S1 parametreleri
R_eq   = 152.4;
r0     = 45;
L_cyl  = 300;
BW_eff = 10;
R_E    = r0 + BW_eff / 2;  % = 50

dome_types = {'hemispherical', 'ellipsoidal', 'isotensoid'};
dome_labels = {'Hemispherical (k=1.0)', 'Ellipsoidal (k=0.7)', 'Isotensoid'};

% Renk paleti — 6 segment icin
seg_colors = [
    0.20  0.40  0.80    % dome1_out (mavi)
    0.85  0.33  0.10    % cyl_fwd   (turuncu)
    0.47  0.67  0.19    % dome2_in  (yesil)
    0.49  0.18  0.56    % dome2_out (mor)
    0.93  0.69  0.13    % cyl_ret   (sari)
    0.30  0.75  0.93    % dome1_in  (acik mavi)
];
seg_names = {'Dome-1 out', 'Cyl fwd', 'Dome-2 in', ...
             'Dome-2 out', 'Cyl ret', 'Dome-1 in'};

%% --- Ana dongu ---
for dt = 1:numel(dome_types)
    dtype = dome_types{dt};
    csv_file = fullfile(ref_dir, sprintf('geodesic_ref_%s_S1.csv', dtype));

    if ~isfile(csv_file)
        warning('CSV bulunamadi: %s', csv_file);
        continue;
    end

    % --- CSV oku ---
    data = readtable(csv_file);
    s     = data.s;
    rho   = data.rho;
    x     = data.x;
    phi   = data.phi;
    alpha = data.alpha;
    n     = height(data);

    % --- Segment sinirlari tespit et ---
    drho = diff(rho);
    seg_bounds = [1];

    for i = 2:numel(drho)
        if drho(i-1) > 0 && drho(i) < 0
            seg_bounds(end+1) = i; %#ok<AGROW>
        elseif drho(i-1) < 0 && drho(i) > 0
            seg_bounds(end+1) = i; %#ok<AGROW>
        end
    end
    seg_bounds(end+1) = n;

    if numel(seg_bounds) >= 7
        seg_idx = seg_bounds(1:7);
    else
        seg_idx = round(linspace(1, n, 7));
    end

    % --- Dome-2 donus noktasi indeksi ---
    % rho'nun yerel minimumu (turnaround, baslangic/bitis haric)
    dome2_turn_idx = [];
    for i = 2:n-1
        if rho(i) < rho(i-1) && rho(i) <= rho(i+1)
            dome2_turn_idx = i;
            break;
        end
    end
    if isempty(dome2_turn_idx)
        dome2_turn_idx = round(n/2);
    end

    % =====================================================================
    %  FIGURE 1: 2-panelli meridyen + sarim acisi
    % =====================================================================
    fig1 = figure('Position', [100 100 1200 500], 'Color', 'w');

    % --- Sol panel: rho vs x ---
    subplot(1, 2, 1);
    hold on; grid on; box on;

    for seg = 1:min(6, numel(seg_idx)-1)
        idx = seg_idx(seg):seg_idx(seg+1);
        plot(x(idx), rho(idx), '-', 'Color', seg_colors(seg,:), ...
             'LineWidth', 1.8);
        plot(x(idx), -rho(idx), '-', 'Color', seg_colors(seg,:), ...
             'LineWidth', 1.0, 'HandleVisibility', 'off');
    end

    yline( R_E, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8, ...
           'HandleVisibility', 'off');
    yline(-R_E, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8, ...
           'HandleVisibility', 'off');
    yline( R_eq, ':', 'Color', [0.7 0.3 0.3], 'LineWidth', 0.8, ...
           'HandleVisibility', 'off');
    yline(-R_eq, ':', 'Color', [0.7 0.3 0.3], 'LineWidth', 0.8, ...
           'HandleVisibility', 'off');

    rho_max = max(rho) * 1.15;
    ylim([-rho_max, rho_max]);
    xlabel('x [mm]', 'FontSize', 11);
    ylabel('\rho [mm]', 'FontSize', 11);
    title(sprintf('%s — Meridyen Profili', dome_labels{dt}), 'FontSize', 12);

    n_seg = min(6, numel(seg_idx)-1);
    legend(seg_names(1:n_seg), 'Location', 'southoutside', ...
           'Orientation', 'horizontal', 'FontSize', 8);

    % --- Sag panel: alpha vs phi ---
    subplot(1, 2, 2);
    hold on; grid on; box on;

    for seg = 1:min(6, numel(seg_idx)-1)
        idx = seg_idx(seg):seg_idx(seg+1);
        plot(rad2deg(phi(idx)), rad2deg(alpha(idx)), '-', ...
             'Color', seg_colors(seg,:), 'LineWidth', 1.8);
    end

    alpha_eq = asin(R_E / R_eq);
    yline(rad2deg(alpha_eq), '--', sprintf('\\alpha_{eq}=%.1f°', rad2deg(alpha_eq)), ...
           'Color', [0.5 0.5 0.5], 'LineWidth', 0.8, ...
           'LabelHorizontalAlignment', 'left', 'FontSize', 9);
    yline(90, ':', '\alpha=90°', 'Color', [0.7 0.3 0.3], ...
           'LineWidth', 0.8, 'LabelHorizontalAlignment', 'right', 'FontSize', 9);

    xlabel('\phi [deg]', 'FontSize', 11);
    ylabel('\alpha [deg]', 'FontSize', 11);
    title(sprintf('%s — Sarim Acisi', dome_labels{dt}), 'FontSize', 12);

    delta_phi = phi(end);
    text(0.98, 0.05, sprintf('\\Delta\\phi = %.2f°', rad2deg(delta_phi)), ...
         'Units', 'normalized', 'HorizontalAlignment', 'right', ...
         'FontSize', 10, 'FontWeight', 'bold', ...
         'BackgroundColor', [1 1 0.85], 'EdgeColor', [0.5 0.5 0.5]);

    sgtitle(sprintf('Phase-2a Geodesic Demo — %s S1  (R_{eq}=%.1f, r_0=%g, L_{cyl}=%g, BW_{eff}=%g)', ...
            dome_labels{dt}, R_eq, r0, L_cyl, BW_eff), ...
            'FontSize', 13, 'FontWeight', 'bold');

    out_file1 = fullfile(fig_dir, sprintf('geodesic_demo_%s.png', dtype));
    exportgraphics(fig1, out_file1, 'Resolution', 150);
    fprintf('Kaydedildi: %s\n', out_file1);
   % close(fig1);

    % =====================================================================
    %  FIGURE 2: 3D mandrel yuzeyi + geodesik yol
    % =====================================================================
    fig2 = figure('Position', [100 100 900 700], 'Color', 'w');
    hold on;

    % --- Mandrel yuzeyini olustur (donum yuzeyi) ---
    % Meridyen profili: ilk yari (dome1_turn -> equator -> cyl -> equator -> dome2_turn)
    % Indeksler: 1 .. dome2_turn_idx arasi mandrel'in bir tarafini verir
    half_idx = 1:dome2_turn_idx;

    % Benzersiz x degerleri icin: monoton artan x profili cikart
    x_half = x(half_idx);
    rho_half = rho(half_idx);

    % x monoton artmiyor olabilir (dome1'de x azalir). Siralayalim.
    [x_sorted, sort_idx] = sort(x_half);
    rho_sorted = rho_half(sort_idx);

    % Tekrarlayan x degerlerini kaldir (silindir bolgesi)
    [x_unique, unique_idx] = unique(x_sorted, 'stable');
    rho_unique = rho_sorted(unique_idx);

    % Azimuthal aci dizisi (0 -> 2*pi tam tur)
    n_theta = 80;
    theta_surf = linspace(0, 2*pi, n_theta);

    % Yuzey koordinatlari
    [X_surf, Theta_surf] = meshgrid(x_unique, theta_surf);
    Rho_surf = repmat(rho_unique', n_theta, 1);
    Y_surf = Rho_surf .* cos(Theta_surf);
    Z_surf = Rho_surf .* sin(Theta_surf);

    % Mandrel yuzeyi — yari saydam
    surf(X_surf, Y_surf, Z_surf, ...
         'FaceColor', [0.6 0.75 0.9], ...
         'FaceAlpha', 0.15, ...
         'EdgeColor', [0.7 0.8 0.9], ...
         'EdgeAlpha', 0.08, ...
         'DisplayName', 'Mandrel yuzeyi');

    % --- Geodesik yol (3D) ---
    % Kartezyen: X_geo = x, Y_geo = rho*cos(phi), Z_geo = rho*sin(phi)
    X_geo = x;
    Y_geo = rho .* cos(phi);
    Z_geo = rho .* sin(phi);

    plot3(X_geo, Y_geo, Z_geo, '-', ...
          'Color', [0.85 0.15 0.15], 'LineWidth', 1.6, ...
          'DisplayName', 'Geodesik yol');

    % --- Ozel noktalar ---
    % Baslangic: Dome-1 turnaround (i=1) — yesil daire
    plot3(X_geo(1), Y_geo(1), Z_geo(1), 'o', ...
          'MarkerSize', 10, 'MarkerFaceColor', [0.2 0.7 0.2], ...
          'MarkerEdgeColor', 'k', 'LineWidth', 1.2, ...
          'DisplayName', sprintf('Baslangic (Dome-1, \\rho=%.0f)', rho(1)));

    % Dome-2 turnaround — mavi kare
    plot3(X_geo(dome2_turn_idx), Y_geo(dome2_turn_idx), Z_geo(dome2_turn_idx), 's', ...
          'MarkerSize', 10, 'MarkerFaceColor', [0.2 0.4 0.8], ...
          'MarkerEdgeColor', 'k', 'LineWidth', 1.2, ...
          'DisplayName', sprintf('Dome-2 donus (\\rho=%.0f)', rho(dome2_turn_idx)));

    % Bitis: son nokta — kirmizi ucgen
    plot3(X_geo(end), Y_geo(end), Z_geo(end), '^', ...
          'MarkerSize', 10, 'MarkerFaceColor', [0.85 0.15 0.15], ...
          'MarkerEdgeColor', 'k', 'LineWidth', 1.2, ...
          'DisplayName', sprintf('Bitis (\\Delta\\phi=%.1f°)', rad2deg(delta_phi)));

    % --- Eksen ve gorunum ayarlari ---
    axis equal;
    grid on; box on;
    xlabel('X [mm]', 'FontSize', 11);
    ylabel('Y [mm]', 'FontSize', 11);
    zlabel('Z [mm]', 'FontSize', 11);

    title(sprintf('Phase-2a — %s S1:  3D Geodesik Yol', dome_labels{dt}), ...
          'FontSize', 13, 'FontWeight', 'bold');

    legend('Location', 'southoutside', 'Orientation', 'horizontal', ...
           'FontSize', 9, 'NumColumns', 2);

    % Kamera acisi: hafif yukari ve yandan
    view(135, 25);

    % Isik efekti
    camlight('headlight');
    lighting gouraud;
    material dull;

    out_file2 = fullfile(fig_dir, sprintf('geodesic_3d_%s.png', dtype));
    exportgraphics(fig2, out_file2, 'Resolution', 150);
    fprintf('Kaydedildi: %s\n', out_file2);
    %close(fig2);
end

fprintf('\nTamamlandi: 6 figure docs/figures/ dizinine kaydedildi.\n');
