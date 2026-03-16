%% PLOT_PATTERN_DEMO  Phase-2b pattern görselleştirme — TEST-02 × 3 dome tipi.
%
% Her dome tipi için 2 panelli figure:
%   Sol  — p vs coverage% scatter (en düşük angular_error vurgulu)
%   Sağ  — Touchpoint dağılımı (en küçük p, sarım sırası renkli)
%
% Çıktı: docs/figures/pattern_demo_{dome_type}.png (3 dosya)
%
% Parametreler: TEST-02 (Endüstriyel COPV)
%   R_eq=152.4, r0=45, L_cyl=300, BW_eff=10, k=0.7(ellipsoidal)
%
% Tarih: 2026-03-16
% Faz: Phase-2b S2

clear; clc; %close all;
fprintf('=== Phase-2b Pattern Demo Görselleri ===\n\n');

%% --- Yol bağımlılıkları ---
script_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, '..', 'phase1a_geometry'));
addpath(fullfile(script_dir, '..', 'phase2a_winding'));

fig_dir = fullfile(script_dir, '..', '..', 'docs', 'figures');
if ~exist(fig_dir, 'dir')
    mkdir(fig_dir);
end

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

%% --- Ana döngü ---
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

    % Geodesic ve pattern arama
    geo = geodesic_single_circuit(dome, R_eq, r0, L_cyl, BW_eff);
    patterns = find_compatible_patterns(geo.delta_phi_circuit, alpha_eq, ...
                                         R_eq, BW_eff, L_cyl, d, coverage_range);
    n_pat = numel(patterns);
    fprintf('  %d pattern bulundu\n', n_pat);

    if n_pat == 0
        fprintf('  UYARI: Pattern yok, figure atlanıyor.\n\n');
        continue;
    end

    % Vektörler
    p_vec   = [patterns.p];
    cov_vec = [patterns.coverage_pct];
    err_vec = [patterns.angular_error_deg];
    ovl_vec = [patterns.overlap_pct];

    % En düşük angular error pattern
    [~, best_idx] = min(err_vec);

    % En küçük p pattern (touchpoint demo)
    tp_pat = patterns(1);  % p sıralı, ilk eleman en küçük p

    %% === Figure oluştur ===
    fig = figure('Position', [100, 100, 1400, 600], 'Color', 'w');

    %% --- Sol panel: p vs Coverage scatter ---
    ax1 = subplot(1, 2, 1);
    hold on; grid on; box on;

    % Ana scatter — angular error'a göre renk
    scatter(p_vec, cov_vec, 80, err_vec, 'filled', 'MarkerEdgeColor', [0.3 0.3 0.3], 'LineWidth', 0.8);
    cb = colorbar;
    cb.Label.String = 'Angular Error [deg]';
    cb.Label.FontSize = 11;
    colormap(ax1, flipud(hot(256)));

    % En iyi pattern vurgula
    scatter(p_vec(best_idx), cov_vec(best_idx), 200, 'p', ...
            'MarkerEdgeColor', [0 0.5 0], 'MarkerFaceColor', [0.3 0.9 0.3], 'LineWidth', 2);

    % 100% coverage çizgisi
    yline(100, '--', '100% Coverage', 'Color', [0.2 0.2 0.8], 'LineWidth', 1.2, ...
          'LabelHorizontalAlignment', 'left', 'FontSize', 10);

    % Etiketler
    xlabel('p (devre sayisi)', 'FontSize', 12);
    ylabel('Coverage [%]', 'FontSize', 12);
    title(sprintf('%s — p vs Coverage\nTEST-02: R_{eq}=%.1f, r_0=%.0f, L_{cyl}=%.0f', ...
          dome_labels{dt}, R_eq, r0, L_cyl), 'FontSize', 13);

    % En iyi pattern annotasyon
    text(p_vec(best_idx) + 0.8, cov_vec(best_idx) + 1.5, ...
         sprintf('Best: p=%d,q=%d\n%.4f deg', patterns(best_idx).p, ...
                 patterns(best_idx).q, err_vec(best_idx)), ...
         'FontSize', 9, 'Color', [0 0.4 0], 'FontWeight', 'bold', ...
         'BackgroundColor', [0.95 1 0.95], 'EdgeColor', [0 0.5 0]);

    set(ax1, 'FontSize', 11);

    %% --- Sağ panel: Touchpoint dağılımı (Kartezyen polar benzeri) ---
    ax2 = subplot(1, 2, 2);
    hold on; grid on; box on; axis equal;

    p_tp = tp_pat.p;
    q_tp = tp_pat.q;
    n_tp = tp_pat.n;  % n = 2*p (her devre 2 ekvator geçişi)

    % Touchpoint açıları (Eq. 5.6): TP_i = [(i-1)*q] mod p
    tp_indices = zeros(1, p_tp);
    for j = 1:p_tp
        tp_indices(j) = mod((j-1) * q_tp, p_tp);
    end

    % Ekvator üzerinde açısal konumlar
    angular_step = 2 * pi / p_tp;
    tp_angles_fwd = tp_indices * angular_step;          % Forward geçiş
    tp_angles_ret = tp_angles_fwd + pi;                  % Return geçiş (karşı dome)
    tp_angles_ret = mod(tp_angles_ret, 2*pi);

    % Ekvator çemberi çiz
    theta_circle = linspace(0, 2*pi, 200);
    plot(cos(theta_circle), sin(theta_circle), '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.5);
    plot(0.75*cos(theta_circle), 0.75*sin(theta_circle), ':', 'Color', [0.85 0.85 0.85], 'LineWidth', 0.8);

    % Sarım sırası renklendirme
    cmap = turbo(p_tp);

    % Sarım sırası bağlantı çizgileri (forward'lar arası)
    for j = 1:p_tp-1
        x_line = [cos(tp_angles_fwd(j)), cos(tp_angles_fwd(j+1))];
        y_line = [sin(tp_angles_fwd(j)), sin(tp_angles_fwd(j+1))];
        plot(x_line, y_line, '-', 'Color', [cmap(j,:) 0.35], 'LineWidth', 1.2);
    end

    % Forward touchpoint'ler (dış çember, r=1)
    for j = 1:p_tp
        plot(cos(tp_angles_fwd(j)), sin(tp_angles_fwd(j)), 'o', 'MarkerSize', 9, ...
             'MarkerFaceColor', cmap(j,:), 'MarkerEdgeColor', 'k', 'LineWidth', 0.8);
    end

    % Return touchpoint'ler (iç çember, r=0.75)
    for j = 1:p_tp
        plot(0.75*cos(tp_angles_ret(j)), 0.75*sin(tp_angles_ret(j)), 'o', 'MarkerSize', 6, ...
             'MarkerFaceColor', cmap(j,:)*0.6, 'MarkerEdgeColor', [0.4 0.4 0.4], 'LineWidth', 0.5);
    end

    % İlk noktayı özel işaretle (yıldız)
    plot(cos(tp_angles_fwd(1)), sin(tp_angles_fwd(1)), 'h', 'MarkerSize', 16, ...
         'MarkerFaceColor', [1 0.8 0], 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);

    % Açı etiketleri (0°, 90°, 180°, 270°)
    for ang_deg = [0, 90, 180, 270]
        ang_rad = ang_deg * pi / 180;
        text(1.15*cos(ang_rad), 1.15*sin(ang_rad), sprintf('%d°', ang_deg), ...
             'HorizontalAlignment', 'center', 'FontSize', 9, 'Color', [0.4 0.4 0.4]);
    end

    title(sprintf('Touchpoint Dagilimi (p=%d, q=%d, n=%d)\nDis=forward, ic=return, sari=baslangic', ...
          p_tp, q_tp, n_tp), 'FontSize', 13);

    xlim([-1.35 1.35]); ylim([-1.35 1.35]);
    set(ax2, 'XTick', [], 'YTick', [], 'FontSize', 10);

    % Colorbar — sarım sırası
    cb2 = colorbar(ax2);
    cb2.Label.String = 'Sarim sirasi';
    cb2.Label.FontSize = 11;
    colormap(ax2, turbo(p_tp));
    clim(ax2, [1 p_tp]);

    %% --- Kaydet ---
    out_file = fullfile(fig_dir, sprintf('pattern_demo_%s.png', dome_types{dt}));
    exportgraphics(fig, out_file, 'Resolution', 200);
    fprintf('  Kaydedildi: %s\n\n', out_file);
    %close(fig);
end

fprintf('=== Tamamlandi: 3 figure docs/figures/ altina kaydedildi ===\n');
