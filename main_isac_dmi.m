%% main_isac_dmi.m
% Fixed guard-band multi-modulation LFM ISAC waveform simulation.
clear; clc; close all;

rootDir = fileparts(mfilename('fullpath'));
addpath(rootDir);
rng(42);

cfg = isac_core('default_config', rootDir);
if cfg.output.saveResults && ~exist(cfg.output.resultsDir, 'dir'), mkdir(cfg.output.resultsDir); end
if cfg.output.saveFigures && ~exist(cfg.output.figuresDir, 'dir'), mkdir(cfg.output.figuresDir); end

model = isac_core('build_model', cfg);
baseline = isac_core('evaluate_lfm_baseline', cfg, model);
mods = cfg.modulation.orders;
allResults = cell(numel(mods), 1);

fprintf('=== Fixed guard-band FA-DMI ISAC simulation ===\n');
fprintf('N=%d, C=%d, K=%d, Fs=%.2f MHz, B_R=%.2f MHz\n', ...
    cfg.signal.N, cfg.signal.C, 2*cfg.signal.C, cfg.signal.Fs/1e6, cfg.signal.B_R/1e6);
isac_core('print_lfm_baseline', baseline);

for im = 1:numel(mods)
    M = mods(im);
    fprintf('\n--- Optimizing M=%d ---\n', M);

    opt = isac_core('optimize_dmi', M, cfg, model);
    evalRes = isac_core('evaluate_solution', M, opt.xBest, cfg, model);
    result = isac_core('merge_structs', opt, evalRes);
    result.M = M;
    result.modulationName = isac_core('modulation_name', M);
    allResults{im} = result;

    if cfg.output.saveResults
        save(fullfile(cfg.output.resultsDir, sprintf('results_M%d.mat', M)), 'cfg', 'model', 'result');
    end

    fprintf(['M=%d %-6s | J=%.4f | Ic=%.3f bit/use | IcNorm=%.4f | IrNorm=%.4f | ' ...
        'pL=%.4f pR=%.4f Pc=%.4f RCR=%.2f dB | PSLR=%.2f dB ISLR=%.2f dB PAPR=%.2f dB | ' ...
        'BER=%.3g SER=%.3g Pd=%.3f feasible=%d\n'], ...
        M, result.modulationName, result.objective, result.Ic, result.IcNorm, result.IrNorm, ...
        result.xBest(1), result.xBest(2), result.Pc, result.RCR_dB, ...
        result.constraints.PSLR_dB, result.constraints.ISLR_dB, result.constraints.PAPR_dB, ...
        result.BER, result.SER, result.P_D, result.isFeasible);
end

isac_core('print_diagnostic_table', allResults);

if cfg.output.saveResults
    save(fullfile(cfg.output.resultsDir, 'results_all.mat'), 'cfg', 'model', 'allResults');
end
if cfg.output.showFigures || cfg.output.saveFigures
    isac_core('plot_results', cfg, model, allResults);
end

fprintf('\nDone. saveResults=%d, saveFigures=%d.\n', cfg.output.saveResults, cfg.output.saveFigures);

