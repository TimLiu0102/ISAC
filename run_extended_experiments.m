%% run_extended_experiments.m
% Extended experiments for fixed guard-band FA-DMI ISAC simulation.
clear; clc; close all;

rootDir = fileparts(mfilename('fullpath'));
addpath(rootDir);
rng(123);

cfg = isac_core('default_config', rootDir);
cfg.output.saveResults = false;
cfg.output.saveFigures = false;
cfg.output.showFigures = true;

fprintf('=== Extended ISAC experiments ===\n');
[baseResults, model] = optimize_all_modulations(cfg);
baseline = isac_core('evaluate_lfm_baseline', cfg, model);
isac_core('print_lfm_baseline', baseline);
isac_core('print_diagnostic_table', baseResults);

isac_core('plot_results', cfg, model, baseResults);
plot_pd_with_lfm_baseline(cfg, baseline, baseResults);
experiment_comm_snr_sweep(cfg, baseResults);
experiment_beta_r_sweep(cfg, baseResults);
experiment_fixed_symbol_vs_bit_rate(cfg, baseResults);
experiment_constraint_sensitivity(cfg);

fprintf('\nExtended experiments finished. No MAT/PNG files were saved.\n');

function [allResults, model] = optimize_all_modulations(cfg)
model = isac_core('build_model', cfg);
mods = cfg.modulation.orders;
allResults = cell(numel(mods), 1);
for im = 1:numel(mods)
    M = mods(im);
    fprintf('Optimizing base Case-A M=%d...\n', M);
    opt = isac_core('optimize_dmi', M, cfg, model);
    evalRes = isac_core('evaluate_solution', M, opt.xBest, cfg, model);
    r = isac_core('merge_structs', opt, evalRes);
    r.M = M;
    r.modulationName = isac_core('modulation_name', M);
    allResults{im} = r;
end
end

function plot_pd_with_lfm_baseline(cfg, baseline, allResults)
figure('Name', 'Pd versus radar SNR with LFM-only baseline');
plot(baseline.P_D_snr_db, baseline.P_D_curve, 'k--', 'LineWidth', 1.6, 'DisplayName', 'LFM-only');
hold on;
for i = 1:numel(allResults)
    r = allResults{i};
    plot(r.P_D_snr_db, r.P_D_curve, 'o-', 'LineWidth', 1.2, 'DisplayName', r.modulationName);
end
grid on; xlabel('Radar detection SNR (dB)'); ylabel('P_D');
ylim([0 1]); legend('Location', 'best');
title(sprintf('P_D vs radar SNR, P_{FA}=%.1e', cfg.eval.P_FA));
end

function experiment_comm_snr_sweep(cfg, baseResults)
fprintf('\n=== Communication SNR sweep, fixed optimized W ===\n');
snrList = cfg.experiment.commSnrDbList;
mods = cellfun(@(r) r.M, baseResults);
Ic = zeros(numel(mods), numel(snrList));
IcNorm = Ic; BER = Ic; SER = Ic;

for im = 1:numel(mods)
    r = baseResults{im};
    const = isac_core('make_constellation', r.M);
    for is = 1:numel(snrList)
        cfgS = cfg;
        cfgS.mi.sigma_c2 = max(r.Pc, 1e-6) / (10^(snrList(is)/10));
        modelS = isac_core('build_model', cfgS);
        W = isac_core('build_W', r.xBest, modelS.C);
        Ic(im, is) = isac_core('estimate_comm_mi_saa', const, W, cfgS, modelS, cfgS.mi.Q, cfgS.mi.L);
        IcNorm(im, is) = Ic(im, is) / (modelS.K * log2(r.M));
        [BER(im, is), SER(im, is)] = isac_core('evaluate_ber_ser', const, W, cfgS, modelS, cfgS.eval.numBerTrials);
    end
end

figure('Name', 'Communication SNR sweep');
subplot(2,2,1); plot_metric_lines(snrList, Ic, baseResults); xlabel('SNR_c (dB)'); ylabel('I_c (bit/use)'); title('Absolute communication MI');
subplot(2,2,2); plot_metric_lines(snrList, IcNorm, baseResults); xlabel('SNR_c (dB)'); ylabel('I_c norm'); title('Normalized communication MI');
subplot(2,2,3); semilogy_metric_lines(snrList, BER, baseResults); xlabel('SNR_c (dB)'); ylabel('BER'); title('BER');
subplot(2,2,4); semilogy_metric_lines(snrList, SER, baseResults); xlabel('SNR_c (dB)'); ylabel('SER'); title('SER');
end

function experiment_beta_r_sweep(cfg, baseResults)
fprintf('\n=== Residual radar-interference beta_r sweep, fixed optimized W ===\n');
betaList = cfg.experiment.betaRList;
mods = cellfun(@(r) r.M, baseResults);
IcNorm = zeros(numel(mods), numel(betaList));
BER = IcNorm; SER = IcNorm;

for im = 1:numel(mods)
    r = baseResults{im};
    const = isac_core('make_constellation', r.M);
    for ib = 1:numel(betaList)
        cfgB = cfg;
        cfgB.mi.beta_r = betaList(ib);
        modelB = isac_core('build_model', cfgB);
        W = isac_core('build_W', r.xBest, modelB.C);
        Ic = isac_core('estimate_comm_mi_saa', const, W, cfgB, modelB, cfgB.mi.Q, cfgB.mi.L);
        IcNorm(im, ib) = Ic / (modelB.K * log2(r.M));
        [BER(im, ib), SER(im, ib)] = isac_core('evaluate_ber_ser', const, W, cfgB, modelB, cfgB.eval.numBerTrials);
    end
end

figure('Name', 'Residual radar interference beta_r sweep');
subplot(1,3,1); semilogx_metric_lines(betaList, IcNorm, baseResults); xlabel('\beta_r'); ylabel('I_c norm'); title('Communication MI');
subplot(1,3,2); loglog_metric_lines(betaList + eps, max(BER, 1e-5), baseResults); xlabel('\beta_r'); ylabel('BER'); title('BER');
subplot(1,3,3); loglog_metric_lines(betaList + eps, max(SER, 1e-5), baseResults); xlabel('\beta_r'); ylabel('SER'); title('SER');
end

function experiment_fixed_symbol_vs_bit_rate(cfg, caseAResults)
fprintf('\n=== Case A fixed symbol rate vs Case B fixed bit rate ===\n');
mods = cfg.modulation.orders;
caseB = cell(numel(mods), 1);
rbOverBr = cfg.signal.SRBR * log2(cfg.experiment.fixedBitRateReferenceM);

for im = 1:numel(mods)
    M = mods(im);
    cfgB = cfg;
    cfgB.signal.SRBR = rbOverBr / log2(M);
    modelB = isac_core('build_model', cfgB);
    fprintf('Optimizing Case-B fixed Rb M=%d, SRBR=%.4f...\n', M, cfgB.signal.SRBR);
    opt = isac_core('optimize_dmi', M, cfgB, modelB);
    evalRes = isac_core('evaluate_solution', M, opt.xBest, cfgB, modelB);
    r = isac_core('merge_structs', opt, evalRes);
    r.M = M;
    r.modulationName = isac_core('modulation_name', M);
    r.SRBR = cfgB.signal.SRBR;
    caseB{im} = r;
end

fprintf('%6s %10s %10s %10s %10s %10s %10s\n', 'Mod', 'Ic_A', 'Ic_B', 'IcN_A', 'IcN_B', 'Pd_A', 'Pd_B');
for im = 1:numel(mods)
    a = caseAResults{im}; b = caseB{im};
    fprintf('%6s %10.3f %10.3f %10.4f %10.4f %10.3f %10.3f\n', ...
        a.modulationName, a.Ic, b.Ic, a.IcNorm, b.IcNorm, a.P_D, b.P_D);
end

modsCat = categorical(cellfun(@(r) r.modulationName, caseAResults, 'UniformOutput', false));
icA = cellfun(@(r) r.Ic, caseAResults); icB = cellfun(@(r) r.Ic, caseB);
icnA = cellfun(@(r) r.IcNorm, caseAResults); icnB = cellfun(@(r) r.IcNorm, caseB);
pdA = cellfun(@(r) r.P_D, caseAResults); pdB = cellfun(@(r) r.P_D, caseB);

figure('Name', 'Fixed symbol rate versus fixed bit rate');
subplot(1,3,1); bar(modsCat, [icA(:), icB(:)]); ylabel('I_c (bit/use)'); legend('Fixed R_s', 'Fixed R_b', 'Location', 'best'); title('Absolute MI');
subplot(1,3,2); bar(modsCat, [icnA(:), icnB(:)]); ylabel('I_c norm'); title('Normalized MI');
subplot(1,3,3); bar(modsCat, [pdA(:), pdB(:)]); ylabel('P_D'); ylim([0 1]); title('Detection probability');
end

function experiment_constraint_sensitivity(cfg)
fprintf('\n=== Constraint-threshold sensitivity, M=%d ===\n', cfg.experiment.constraintSweepM);
M = cfg.experiment.constraintSweepM;
pslrList = cfg.experiment.pslrMaxList;
paprList = cfg.experiment.paprMaxList;
pslrObj = zeros(size(pslrList)); pslrFeas = false(size(pslrList)); pslrIc = pslrObj;
paprObj = zeros(size(paprList)); paprFeas = false(size(paprList)); paprIc = paprObj;

for i = 1:numel(pslrList)
    cfgS = cfg; cfgS.constraints.PSLR_max = pslrList(i);
    modelS = isac_core('build_model', cfgS);
    opt = isac_core('optimize_dmi', M, cfgS, modelS);
    evalRes = isac_core('evaluate_solution', M, opt.xBest, cfgS, modelS);
    r = isac_core('merge_structs', opt, evalRes);
    pslrObj(i) = r.objective; pslrIc(i) = r.IcNorm; pslrFeas(i) = r.isFeasible;
    fprintf('PSLRmax=%6.1f dB | feasible=%d | J=%.4f | IcNorm=%.4f | actual PSLR=%.2f dB\n', ...
        pslrList(i), r.isFeasible, r.objective, r.IcNorm, r.constraints.PSLR_dB);
end

for i = 1:numel(paprList)
    cfgS = cfg; cfgS.constraints.PAPR_max = 10^(paprList(i)/10);
    modelS = isac_core('build_model', cfgS);
    opt = isac_core('optimize_dmi', M, cfgS, modelS);
    evalRes = isac_core('evaluate_solution', M, opt.xBest, cfgS, modelS);
    r = isac_core('merge_structs', opt, evalRes);
    paprObj(i) = r.objective; paprIc(i) = r.IcNorm; paprFeas(i) = r.isFeasible;
    fprintf('PAPRmax=%6.1f dB | feasible=%d | J=%.4f | IcNorm=%.4f | actual PAPR=%.2f dB\n', ...
        paprList(i), r.isFeasible, r.objective, r.IcNorm, r.constraints.PAPR_dB);
end

figure('Name', 'Constraint threshold sensitivity');
subplot(1,2,1);
yyaxis left; plot(pslrList, pslrObj, 'o-', pslrList, pslrIc, 's-', 'LineWidth', 1.2); ylabel('Objective / I_c norm');
yyaxis right; stem(pslrList, double(pslrFeas), 'filled'); ylim([-0.1 1.1]); ylabel('Feasible');
grid on; xlabel('PSLR_{max} (dB)'); title(sprintf('PSLR sensitivity, %s', isac_core('modulation_name', M)));
subplot(1,2,2);
yyaxis left; plot(paprList, paprObj, 'o-', paprList, paprIc, 's-', 'LineWidth', 1.2); ylabel('Objective / I_c norm');
yyaxis right; stem(paprList, double(paprFeas), 'filled'); ylim([-0.1 1.1]); ylabel('Feasible');
grid on; xlabel('PAPR_{max} (dB)'); title(sprintf('PAPR sensitivity, %s', isac_core('modulation_name', M)));
end

function plot_metric_lines(x, Y, results)
hold on;
for i = 1:size(Y, 1)
    plot(x, Y(i, :), 'o-', 'LineWidth', 1.2, 'DisplayName', results{i}.modulationName);
end
grid on; legend('Location', 'best');
end

function semilogy_metric_lines(x, Y, results)
hold on;
for i = 1:size(Y, 1)
    semilogy(x, max(Y(i, :), 1e-5), 'o-', 'LineWidth', 1.2, 'DisplayName', results{i}.modulationName);
end
grid on; legend('Location', 'best');
end

function semilogx_metric_lines(x, Y, results)
hold on;
for i = 1:size(Y, 1)
    semilogx(x + eps, Y(i, :), 'o-', 'LineWidth', 1.2, 'DisplayName', results{i}.modulationName);
end
grid on; legend('Location', 'best');
end

function loglog_metric_lines(x, Y, results)
hold on;
for i = 1:size(Y, 1)
    loglog(x, Y(i, :), 'o-', 'LineWidth', 1.2, 'DisplayName', results{i}.modulationName);
end
grid on; legend('Location', 'best');
end

