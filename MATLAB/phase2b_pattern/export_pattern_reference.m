%% EXPORT_PATTERN_REFERENCE  TEST-01..04 pattern tablolarını CSV olarak dışa aktar.
%
% Her test senaryosu için ellipsoidal dome (k=0.7) ile pattern arama yapılır
% ve sonuçlar reference_data/ dizinine CSV olarak kaydedilir.
%
% CSV sütunları:
%   p, q, n, dk, type, coverage_pct, overlap_pct,
%   delta_L_cyl_mm, L_cyl_adj_mm, angular_error_rad
%
% type: LEADING (+1), LAGGING (-1), NEITHER (0) — Diophantine p·dk − n_ref·d
%
% Referans: Phase-2b S2, Karar-16 test senaryoları
% Tarih: 2026-03-16
% Faz: Phase-2b S2

clear; clc;
fprintf('=== Phase-2b S2: Referans CSV Dışa Aktarım ===\n\n');

%% --- Yol bağımlılıkları ---
script_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, '..', 'phase1a_geometry'));
addpath(fullfile(script_dir, '..', 'phase2a_winding'));

output_dir = fullfile(script_dir, 'reference_data');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

%% --- Ortak parametreler ---
BW_eff = 10.0;
d      = 1;
k      = 0.7;
coverage_range = [100, 150];

%% --- TEST-01..04 parametreleri (Karar-16) ---
test_ids   = {'TEST-01', 'TEST-02', 'TEST-03', 'TEST-04'};
test_names = {'ASTM_Subscale', 'Industrial_COPV', 'Small_Opening', 'H2_Aerospace'};
test_Req   = [73,    152.4, 150,   200  ];
test_r0    = [22,    45,    10,    50   ];
test_Lcyl  = [200,   300,   250,   400  ];

%% --- Her senaryo için pattern arama ve CSV dışa aktarım ---
for t = 1:4
    fprintf('--- %s: %s ---\n', test_ids{t}, test_names{t});
    fprintf('  R_eq=%.1f, r0=%.1f, L_cyl=%.1f, k=%.2f, BW_eff=%.1f\n', ...
            test_Req(t), test_r0(t), test_Lcyl(t), k, BW_eff);

    % Dome profili ve geodesic hesabı
    dome = ellipsoidal_dome_profile(test_Req(t), test_r0(t), k, 1000);
    geo  = geodesic_single_circuit(dome, test_Req(t), test_r0(t), test_Lcyl(t), BW_eff);

    R_E_t     = test_r0(t) + BW_eff / 2;
    alpha_t   = asin(R_E_t / test_Req(t));
    cos_alpha = cos(alpha_t);
    n_ref     = floor(2 * pi * test_Req(t) * d * cos_alpha / BW_eff);

    fprintf('  α_eq=%.4f° R_E=%.2f n_ref=%d Δφ=%.6f rad\n', ...
            alpha_t*180/pi, R_E_t, n_ref, geo.delta_phi_circuit);

    % Pattern arama
    patterns = find_compatible_patterns(geo.delta_phi_circuit, alpha_t, ...
                                         test_Req(t), BW_eff, test_Lcyl(t), d, coverage_range);
    n_pat = numel(patterns);
    fprintf('  Bulunan pattern sayısı: %d\n', n_pat);

    if n_pat == 0
        fprintf('  UYARI: Pattern bulunamadı, CSV atlanıyor.\n\n');
        continue;
    end

    % CSV dosyası oluştur
    csv_file = fullfile(output_dir, sprintf('%s_ellipsoidal_patterns.csv', test_ids{t}));
    fid = fopen(csv_file, 'w');
    if fid == -1
        error('CSV dosyası açılamadı: %s', csv_file);
    end

    % Başlık satırı
    fprintf(fid, 'p,q,n,dk,type,coverage_pct,overlap_pct,delta_L_cyl_mm,L_cyl_adj_mm,angular_error_rad\n');

    for i = 1:n_pat
        pat = patterns(i);

        % dk hesabı: dk = q mod p
        dk = mod(pat.q, pat.p);

        % Diophantine tipi: p·dk − n_ref·d
        dioph = pat.p * dk - n_ref * d;
        if dioph == 1
            type_str = 'LEADING';
        elseif dioph == -1
            type_str = 'LAGGING';
        else
            type_str = 'NEITHER';
        end

        % CSV satırı yaz
        fprintf(fid, '%d,%d,%d,%d,%s,%.6f,%.6f,%.6f,%.6f,%.10e\n', ...
                pat.p, pat.q, pat.n, dk, type_str, ...
                pat.coverage_pct, pat.overlap_pct, ...
                pat.delta_L_cyl, pat.L_cyl_adj, pat.angular_error);
    end

    fclose(fid);
    fprintf('  CSV kaydedildi: %s\n\n', csv_file);
end

fprintf('=== Dışa aktarım tamamlandı ===\n');
