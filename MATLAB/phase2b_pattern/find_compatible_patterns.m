function patterns = find_compatible_patterns(delta_phi_circuit, alpha_eq, R_eq, BW_eff, L_cyl, d, coverage_range)
% FIND_COMPATIBLE_PATTERNS  Verilen Δφ_circuit için uyumlu p/q pattern çiftleri bulur.
%
% Phase-2b: Pattern generation — angle-driven mod.
% Karar-13: p/q notasyonu, gcd(p,q)=1 zorunlu (S-PAT-02).
% K-2b-04: Varsayılan sıralama p küçükten büyüğe.
%
% MATEMATİK:
%   Pattern kapatma koşulu: p · Δφ_circuit = 2π · q
%   Toplam bant sayısı: n = 2·p (her devre 2 ekvator geçişi: forward + return)
%   Tam kaplama bant sayısı: n_max = floor(2π·R_eq·d·cos(α_eq) / BW_eff)
%   Ekvator kaplama: coverage = n · BW_eff / (2π · R_eq · d · cos(α_eq)) × 100%
%   Overlap: coverage - 100%
%   L_cyl ince ayar: ΔL = (2πq/p − Δφ_actual) · R_eq / (2·tan(α_eq))
%   L_cyl_adj = L_cyl + ΔL; L_cyl_adj < 0 → pattern elenir.
%
% KULLANIM:
%   patterns = find_compatible_patterns(delta_phi_circuit, alpha_eq, R_eq, BW_eff, L_cyl)
%   patterns = find_compatible_patterns(..., d, coverage_range)
%
% GİRDİLER:
%   delta_phi_circuit — Tek devre toplam açısal ilerleme [rad]
%   alpha_eq          — Ekvator winding açısı [rad]
%   R_eq              — Ekvator yarıçapı [mm]
%   BW_eff            — Efektif bant genişliği [mm]
%   L_cyl             — Mevcut silindir uzunluğu [mm]
%   d                 — Winding pattern skip index [-] (varsayılan: 1)
%   coverage_range    — [min, max] kaplama yüzdesi (varsayılan: [100, 150])
%
% ÇIKTI:
%   patterns — Struct dizisi (p küçükten büyüğe sıralı), alanlar:
%     .p              — Devre sayısı (circuits per pattern)
%     .q              — Tam tur sayısı (full revolutions per circuit)
%     .n              — Toplam bant sayısı (= 2·p)
%     .coverage_pct   — Ekvator kaplama yüzdesi [%]
%     .overlap_pct    — Örtüşme yüzdesi [%] (= coverage - 100)
%     .delta_phi_ideal — İdeal Δφ_circuit = 2πq/p [rad]
%     .angular_error  — |Δφ_circuit - 2πq/p| [rad]
%     .angular_error_deg — |Δφ_circuit - 2πq/p| [°]
%     .delta_L_cyl    — Silindir ince ayar miktarı [mm]
%     .L_cyl_adj      — Ayarlanmış silindir uzunluğu [mm]
%
% Referans: Karar-13, S-PAT-01, S-PAT-02, K-2b-04
% Tarih: 2026-03-16
% Faz: Phase-2b

    %% --- Varsayılan parametreler ---
    if nargin < 6 || isempty(d)
        d = 1;
    end
    if nargin < 7 || isempty(coverage_range)
        coverage_range = [100, 150];
    end

    %% --- Girdi doğrulama ---
    if delta_phi_circuit <= 0
        error('find_compatible_patterns:gecersizDeltaPhi', ...
              'delta_phi_circuit pozitif olmalıdır. Girilen: %g', delta_phi_circuit);
    end
    if alpha_eq <= 0 || alpha_eq >= pi/2
        error('find_compatible_patterns:gecersizAlpha', ...
              'alpha_eq ∈ (0, π/2) olmalıdır. Girilen: %g rad (%.2f°)', ...
              alpha_eq, alpha_eq * 180 / pi);
    end
    if R_eq <= 0
        error('find_compatible_patterns:gecersizReq', ...
              'R_eq pozitif olmalıdır. Girilen: %g', R_eq);
    end
    if BW_eff <= 0
        error('find_compatible_patterns:gecersizBWeff', ...
              'BW_eff pozitif olmalıdır. Girilen: %g', BW_eff);
    end
    if d < 1 || mod(d, 1) ~= 0
        error('find_compatible_patterns:gecersizD', ...
              'd pozitif tamsayı olmalıdır. Girilen: %g', d);
    end
    if L_cyl < 0
        error('find_compatible_patterns:negativeLcyl', ...
              'L_cyl >= 0 olmalıdır. Girilen: %g', L_cyl);
    end
    if numel(coverage_range) ~= 2 || coverage_range(1) >= coverage_range(2)
        error('find_compatible_patterns:gecersizCoverage', ...
              'coverage_range = [min, max] olmalı ve min < max.');
    end

    coverage_min = coverage_range(1);
    coverage_max = coverage_range(2);
    tan_alpha = tan(alpha_eq);        % ΔL_cyl hesabı için

    %% --- p aralığını coverage kısıtından belirle ---
    % n = 2·p (her devre 2 ekvator geçişi: forward + return)
    % coverage = n · BW_eff / (2π · R_eq · d · cos(α_eq)) × 100
    %          = p · BW_eff / (π · R_eq · d · cos(α_eq)) × 100
    % p = coverage/100 × π · R_eq · d · cos(α_eq) / BW_eff
    cos_alpha = cos(alpha_eq);
    n_max = floor(2 * pi * R_eq * d * cos_alpha / BW_eff);  % tam kaplama bant sayısı
    p_for_100 = pi * R_eq * d * cos_alpha / BW_eff;         % p for 100% coverage

    p_min = ceil(coverage_min / 100 * p_for_100);
    p_max = floor(coverage_max / 100 * p_for_100);

    if p_min < 1
        p_min = 1;
    end
    if p_max < p_min
        warning('find_compatible_patterns:aralikBos', ...
                'Coverage aralığı [%.1f, %.1f]%% için geçerli p yok (p aralığı boş). p_100%% = %.2f', ...
                coverage_min, coverage_max, p_for_100);
        patterns = [];
        return;
    end

    %% --- İdeal oran ---
    ratio = delta_phi_circuit / (2 * pi);   % q/p ideal oranı

    %% --- p taraması: her p için en yakın q bul ---
    results = [];
    count = 0;

    for p = p_min:p_max
        % En yakın tam q
        q = round(p * ratio);

        if q < 1
            continue;
        end

        % S-PAT-02: gcd(p,q) = 1 kontrolü
        if gcd(p, q) ~= 1
            continue;
        end

        % d kontrolü: skip index ile uyumluluk
        % d=1 temel durum, d>1 ileri faz
        if d > 1
            % d-diamond pattern: d·q mod p ilişkisi
            % Temel kontrol: gcd(p, d) = 1 olmalı (d-diamond uyumluluk)
            if gcd(p, d) ~= 1
                continue;
            end
        end

        % İdeal Δφ ve hata
        delta_phi_ideal = 2 * pi * q / p;
        angular_error   = abs(delta_phi_circuit - delta_phi_ideal);

        % L_cyl ince ayar (Karar-13: S-PAT-01)
        % ΔL = (2πq/p − Δφ_actual) · R_eq / (2·tan(α_eq))
        delta_L_cyl = (delta_phi_ideal - delta_phi_circuit) * R_eq / (2 * tan_alpha);
        L_cyl_adj   = L_cyl + delta_L_cyl;

        % L_cyl_adj < 0 → fiziksel olarak imkansız, pattern elenir
        if L_cyl_adj < 0
            continue;
        end

        % Kaplama hesabı
        % n = 2·p, coverage = n·BW_eff / (2π·R_eq·d·cos(α_eq)) × 100
        n = 2 * p;
        coverage_pct = n * BW_eff / (2 * pi * R_eq * d * cos_alpha) * 100;
        overlap_pct  = coverage_pct - 100;

        % Sonuç kaydet
        count = count + 1;
        results(count).p              = p;                          %#ok<AGROW>
        results(count).n              = n;                          %#ok<AGROW>
        results(count).q              = q;                          %#ok<AGROW>
        results(count).coverage_pct   = coverage_pct;               %#ok<AGROW>
        results(count).overlap_pct    = overlap_pct;                %#ok<AGROW>
        results(count).delta_phi_ideal = delta_phi_ideal;           %#ok<AGROW>
        results(count).angular_error   = angular_error;             %#ok<AGROW>
        results(count).angular_error_deg = angular_error * 180/pi;  %#ok<AGROW>
        results(count).delta_L_cyl    = delta_L_cyl;                %#ok<AGROW>
        results(count).L_cyl_adj      = L_cyl_adj;                  %#ok<AGROW>
    end

    %% --- Sonuçları p'ye göre sırala (K-2b-04: küçükten büyüğe) ---
    if isempty(results)
        warning('find_compatible_patterns:sonucYok', ...
                'Coverage aralığı [%.1f, %.1f]%% ve gcd(p,q)=1 koşulunda uyumlu pattern bulunamadı.', ...
                coverage_min, coverage_max);
        patterns = [];
    else
        [~, sort_idx] = sort([results.p], 'ascend');
        patterns = results(sort_idx);
    end

end
