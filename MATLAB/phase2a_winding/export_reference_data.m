%% EXPORT_REFERENCE_DATA  Phase-2a referans CSV dosyalarini disa aktar.
%
% 3 dome tipi x 4 senaryo = 12 CSV dosyasi.
% Dizin: MATLAB/phase2a_winding/reference_data/
% Format: s, rho, x, phi, alpha, kn (6 sutun, baslik satirli)
%
% Dosya adlandirma: geodesic_ref_{dome}_{senaryo}.csv
%   ornek: geodesic_ref_hemispherical_S1.csv
%
% Referans: docs/phase2a_geodesic_math.md
% Faz: Phase-2a S2 — dogrulama referans verisi
% Tarih: 2026-03-11

clear; clc; close all;

%% --- Path ayarlari ---
this_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(this_dir, 'utils'));
addpath(fullfile(this_dir, '..', 'phase1a_geometry'));

out_dir = fullfile(this_dir, 'reference_data');
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

%% --- Senaryolar (verify_geodesic_path.m ile ayni) ---
scenarios = struct( ...
    'name',   {'S1-COPV-Standart', 'S2-Buyuk',      'S3-Kucuk',      'S4-GenisAciklik'}, ...
    'tag',    {'S1',               'S2',             'S3',             'S4'}, ...
    'R_eq',   {152.4,               200,              100,              120}, ...
    'r0',     {45,                   60,               20,               50}, ...
    'L_cyl',  {300,                  400,              150,              250}, ...
    'BW_eff', {10,                   12,               6,                8}, ...
    'k_ell',  {0.7,                  0.5,              1.2,              0.8} ...
);

dome_types = {'hemispherical', 'ellipsoidal', 'isotensoid'};

%% --- Ana dongu: 3 dome x 4 senaryo ---
n_exported = 0;

fprintf('=============================================================\n');
fprintf(' PHASE-2a: REFERANS CSV DISA AKTARIMI\n');
fprintf(' Tarih: %s\n', datestr(now, 'yyyy-mm-dd HH:MM'));
fprintf(' Dizin: %s\n', out_dir);
fprintf('=============================================================\n\n');

for dt = 1:numel(dome_types)
    dtype = dome_types{dt};

    for sc = 1:numel(scenarios)
        S = scenarios(sc);
        R_E = S.r0 + S.BW_eff / 2;

        % R_E < R_eq kontrolu
        if R_E >= S.R_eq
            fprintf('[SKIP] %s | %s : R_E >= R_eq\n', dtype, S.name);
            continue;
        end

        %% --- Dome profil uret ---
        switch dtype
            case 'hemispherical'
                dome = hemispherical_dome_profile(S.R_eq, S.r0, 1000);
            case 'ellipsoidal'
                dome = ellipsoidal_dome_profile(S.R_eq, S.r0, S.k_ell, 1000);
            case 'isotensoid'
                dome = isotensoid_dome_profile(S.R_eq, S.r0, 1000);
        end

        %% --- Geodesic yol uret ---
        circuit = geodesic_single_circuit(dome, S.R_eq, S.r0, ...
                                           S.L_cyl, S.BW_eff, 2000);
        close all;

        %% --- CSV yaz ---
        fname = sprintf('geodesic_ref_%s_%s.csv', dtype, S.tag);
        fpath = fullfile(out_dir, fname);

        % Veri matrisi: [s, rho, x, phi, alpha, kn]
        data = [circuit.s(:), circuit.rho(:), circuit.x(:), ...
                circuit.phi(:), circuit.alpha(:), circuit.kn(:)];

        % Dosya yaz
        fid = fopen(fpath, 'w');
        if fid == -1
            error('Dosya acilamadi: %s', fpath);
        end

        % Baslik satiri
        fprintf(fid, 's,rho,x,phi,alpha,kn\n');

        % Veri satirlari (%.10g — yeterli hassasiyet, kompakt)
        for i = 1:size(data, 1)
            fprintf(fid, '%.10g,%.10g,%.10g,%.10g,%.10g,%.10g\n', data(i, :));
        end

        fclose(fid);
        n_exported = n_exported + 1;

        fprintf('[OK] %s | %-18s -> %s  (%d nokta)\n', ...
                dtype, S.name, fname, size(data, 1));
    end
end

fprintf('\n=============================================================\n');
fprintf(' Toplam: %d CSV dosyasi disa aktarildi.\n', n_exported);
fprintf('=============================================================\n');
