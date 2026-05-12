function varargout = isac_core(action, varargin)
%ISAC_CORE Shared dispatcher for the compact ISAC MATLAB project.
switch lower(action)
    case 'default_config'
        varargout{1} = default_config(varargin{:});
    case 'build_model'
        varargout{1} = build_model(varargin{:});
    case 'evaluate_lfm_baseline'
        varargout{1} = evaluate_lfm_baseline(varargin{:});
    case 'print_lfm_baseline'
        print_lfm_baseline(varargin{:});
    case 'optimize_dmi'
        varargout{1} = optimize_dmi(varargin{:});
    case 'evaluate_solution'
        varargout{1} = evaluate_solution(varargin{:});
    case 'evaluate_design'
        varargout{1} = evaluate_design(varargin{:});
    case 'plot_results'
        plot_results(varargin{:});
    case 'print_diagnostic_table'
        print_diagnostic_table(varargin{:});
    case 'merge_structs'
        varargout{1} = merge_structs(varargin{:});
    case 'modulation_name'
        varargout{1} = modulation_name(varargin{:});
    case 'make_constellation'
        varargout{1} = make_constellation(varargin{:});
    case 'build_w'
        varargout{1} = build_W(varargin{:});
    case 'estimate_comm_mi_saa'
        varargout{1} = estimate_comm_mi_saa(varargin{:});
    case 'evaluate_ber_ser'
        [varargout{1}, varargout{2}] = evaluate_ber_ser(varargin{:});
    otherwise
        error('Unknown isac_core action: %s', action);
end
end

function cfg = default_config(rootDir)
cfg.signal.N = 96;
cfg.signal.Fs = 20e6;
cfg.signal.T = cfg.signal.N / cfg.signal.Fs;
cfg.signal.B_R = 8e6;
cfg.signal.f_delta = 0.25;
cfg.signal.SRBR = 0.16;
cfg.signal.beta_p = 0.25;
cfg.signal.C = 3;
cfg.signal.rrcSpanSymbols = 6;
cfg.modulation.orders = [2 4 16 64];
cfg.power.P_r = 1.0;
cfg.power.P_c_max = 0.25;
cfg.power.RCR_min = 6.0;
cfg.constraints.PSLR_max = -7.0;
cfg.constraints.ISLR_max = -2.0;
cfg.constraints.PAPR_max = 8.0;
cfg.constraints.B_main_max = 11;
cfg.constraints.mainlobeHalfWidth = 4;
cfg.mi.Q = 12;
cfg.mi.L = 24;
cfg.mi.omega_c = 0.55;
cfg.mi.omega_r = 0.45;
cfg.mi.sigma_c2 = 2e-2;
cfg.mi.sigma_r2 = 2e-2;
cfg.mi.beta_r = 0.03;
cfg.channel.h_rc_gain = 0.15;
cfg.channel.h_cr_gain = 1.0;
cfg.opt.pGrid = linspace(0, 0.18, 4);
cfg.opt.thetaGrid = linspace(0, 2*pi, 4);
cfg.opt.thetaGrid(end) = [];
cfg.opt.localIters = 8;
cfg.opt.initialStep = [0.025; 0.025; pi/6; pi/6];
cfg.opt.penaltyWeight = 8.0;
cfg.opt.feasTol = 1e-8;
cfg.eval.numBerTrials = 100;
cfg.eval.numPdTrials = 300;
cfg.eval.P_FA = 1e-3;
cfg.eval.radarDetectSnrDb = 8;
cfg.eval.radarDetectSnrDbList = -5:5:15;
cfg.eval.snrDbList = 0:5:20;
cfg.eval.paprSamples = 12;
cfg.experiment.commSnrDbList = -5:5:20;
cfg.experiment.betaRList = [0 1e-3 1e-2 1e-1];
cfg.experiment.fixedBitRateReferenceM = 4;
cfg.experiment.pslrMaxList = [-13 -15 -18 -20];
cfg.experiment.paprMaxList = [6 7 8];
cfg.experiment.constraintSweepM = 4;
cfg.output.resultsDir = fullfile(rootDir, 'results');
cfg.output.figuresDir = fullfile(rootDir, 'figures');
cfg.output.saveResults = false;
cfg.output.saveFigures = false;
cfg.output.showFigures = true;
end

function model = build_model(cfg)
N = cfg.signal.N; Fs = cfg.signal.Fs; B_R = cfg.signal.B_R; C = cfg.signal.C; Ts = 1/Fs;
s_r = generate_lfm_unit(N, Fs, B_R);
P = build_pulse_matrix(N, Fs, cfg.signal.SRBR * B_R, cfg.signal.beta_p, C, cfg.signal.rrcSpanSymbols);
fc1 = -B_R/2 - cfg.signal.f_delta * B_R;
fc2 =  B_R/2 + cfg.signal.f_delta * B_R;
n = (0:N-1).';
G = [diag(exp(1j*2*pi*fc1*Ts*n))*P, diag(exp(1j*2*pi*fc2*Ts*n))*P];
Hc = eye(N);
Hcr = cfg.channel.h_cr_gain * eye(N);
numTargets = min(8, N);
Xr = zeros(N, numTargets);
for k = 1:numTargets
    Xr(:, k) = delay_vector(s_r, k-1);
end
Sigma_g = diag(exp(-0.25 * (0:numTargets-1)));
Hrc = build_radar_comm_coupling(N, cfg.channel.h_rc_gain, G, Xr);
Rv = cfg.mi.sigma_c2 * eye(N) + cfg.mi.beta_r * cfg.power.P_r * (Hcr*s_r) * (Hcr*s_r)';
Sg = Xr * Sigma_g * Xr';
model = struct('s_r', s_r, 'P', P, 'G', G, 'fc1', fc1, 'fc2', fc2, ...
    'Hc', Hc, 'Hcr', Hcr, 'Hrc', Hrc, 'Rv', hermitian_regularize(Rv), ...
    'Xr', Xr, 'Sigma_g', Sigma_g, 'Sg', hermitian_regularize(Sg), ...
    'I_r0', real_logdet(eye(N) + (1/cfg.mi.sigma_r2)*Sg), ...
    'K', 2*C, 'C', C, 'N', N, 'Fs', Fs);
end

function baseline = evaluate_lfm_baseline(cfg, model)
rngState = rng; cleanupObj = onCleanup(@() rng(rngState)); %#ok<NASGU>
const = make_constellation(4);
W0 = zeros(model.K, model.K);
waveMet = compute_waveform_metrics(W0, cfg, model, const, 1);
P_D0 = evaluate_pd(W0, cfg, model, cfg.eval.numPdTrials, cfg.eval.P_FA);
pdCurve = zeros(size(cfg.eval.radarDetectSnrDbList));
cfgPd = cfg;
for k = 1:numel(cfg.eval.radarDetectSnrDbList)
    cfgPd.eval.radarDetectSnrDb = cfg.eval.radarDetectSnrDbList(k);
    pdCurve(k) = evaluate_pd(W0, cfgPd, model, cfg.eval.numPdTrials, cfg.eval.P_FA);
end
baseline = struct('Ir0', model.I_r0, 'P_D0', P_D0, 'P_D_snr_db', cfg.eval.radarDetectSnrDbList, ...
    'P_D_curve', pdCurve, 'PSLR_dB', waveMet.PSLR_dB, 'ISLR_dB', waveMet.ISLR_dB, ...
    'PAPR', waveMet.PAPR, 'PAPR_dB', waveMet.PAPR_dB, 'B_main', waveMet.B_main, ...
    'constraints', waveMet);
end

function print_lfm_baseline(b)
fprintf('\n=== LFM-only radar baseline ===\n');
fprintf('Ir0=%.4f bit/use | Pd0=%.3f | PSLR=%.2f dB | ISLR=%.2f dB | PAPR=%.2f dB | Bmain=%d\n', ...
    b.Ir0, b.P_D0, b.PSLR_dB, b.ISLR_dB, b.PAPR_dB, b.B_main);
end

function opt = optimize_dmi(M, cfg, model)
const = make_constellation(M);
best = struct('x', [], 'score', -inf);
pGrid = cfg.opt.pGrid; tGrid = cfg.opt.thetaGrid;
for iL = 1:numel(pGrid)
    for iR = 1:numel(pGrid)
        for iTL = 1:numel(tGrid)
            for iTR = 1:numel(tGrid)
                x = project_power_rcr([pGrid(iL); pGrid(iR); tGrid(iTL); tGrid(iTR)], cfg, model);
                cand = objective_with_penalty(M, const, x, cfg, model);
                if cand.score > best.score, best = cand; end
            end
        end
    end
end
x = best.x; step = cfg.opt.initialStep; history = zeros(cfg.opt.localIters, 3);
for it = 1:cfg.opt.localIters
    improved = false;
    base = objective_with_penalty(M, const, x, cfg, model);
    dirs = [eye(4), -eye(4)];
    for d = 1:size(dirs, 2)
        xTry = x + step .* dirs(:, d);
        xTry(3:4) = mod(xTry(3:4), 2*pi);
        xTry = project_power_rcr(xTry, cfg, model);
        cand = objective_with_penalty(M, const, xTry, cfg, model);
        if cand.score > base.score
            x = xTry; base = cand; improved = true;
        end
    end
    if ~improved, step = 0.65 * step; end
    history(it, :) = [base.objective, base.violation, base.score];
end
final = objective_with_penalty(M, const, x, cfg, model);
opt = final; opt.xBest = x; opt.gridBest = best; opt.history = history;
end

function cand = objective_with_penalty(M, const, x, cfg, model)
met = evaluate_design(M, const, x, cfg, model, cfg.mi.Q, cfg.mi.L, cfg.eval.paprSamples);
violation = sum(max(0, met.g).^2);
cand = met; cand.x = x; cand.violation = violation;
cand.score = met.objective - cfg.opt.penaltyWeight * violation;
end

function res = evaluate_solution(M, x, cfg, model)
const = make_constellation(M);
W = build_W(x, model.C);
res = evaluate_design(M, const, x, cfg, model, max(2*cfg.mi.Q, 40), max(2*cfg.mi.L, 80), cfg.eval.paprSamples);
[BER, SER] = evaluate_ber_ser(const, W, cfg, model, cfg.eval.numBerTrials);
P_D = evaluate_pd(W, cfg, model, cfg.eval.numPdTrials, cfg.eval.P_FA);
pdCurve = zeros(size(cfg.eval.radarDetectSnrDbList));
cfgPd = cfg;
for k = 1:numel(cfg.eval.radarDetectSnrDbList)
    cfgPd.eval.radarDetectSnrDb = cfg.eval.radarDetectSnrDbList(k);
    pdCurve(k) = evaluate_pd(W, cfgPd, model, cfg.eval.numPdTrials, cfg.eval.P_FA);
end
res.W = W; res.BER = BER; res.SER = SER; res.P_D = P_D;
res.P_D_snr_db = cfg.eval.radarDetectSnrDbList; res.P_D_curve = pdCurve;
end

function met = evaluate_design(M, const, x, cfg, model, Q, L, paprSamples)
W = build_W(x, model.C);
Pc = communication_power(x, model);
if Pc <= 0, RCR = inf; else, RCR = 10*log10(cfg.power.P_r/Pc); end
Ic = estimate_comm_mi_saa(const, W, cfg, model, Q, L);
Ir = compute_radar_mi(W, cfg, model);
IcNorm = Ic / (model.K * log2(M));
IrNorm = Ir / max(model.I_r0, eps);
objective = cfg.mi.omega_c * IcNorm + cfg.mi.omega_r * IrNorm;
waveMet = compute_waveform_metrics(W, cfg, model, const, paprSamples);
g = [Pc-cfg.power.P_c_max; cfg.power.RCR_min-RCR; waveMet.PSLR_dB-cfg.constraints.PSLR_max; ...
     waveMet.ISLR_dB-cfg.constraints.ISLR_max; waveMet.PAPR-cfg.constraints.PAPR_max; ...
     waveMet.B_main-cfg.constraints.B_main_max];
met = struct('Pc', Pc, 'RCR_dB', RCR, 'Ic', Ic, 'Ir', Ir, 'IcNorm', IcNorm, ...
    'IrNorm', IrNorm, 'objective', objective, 'constraints', waveMet, 'g', g, ...
    'isFeasible', all(g <= cfg.opt.feasTol));
end

function plot_results(cfg, model, allResults)
mods = zeros(numel(allResults), 1);
IcNorm = mods; IcAbs = mods; Ir = mods; BER = mods; SER = mods; PD = mods;
PSLR = mods; ISLR = mods; PAPR = mods; Obj = mods;
for i = 1:numel(allResults)
    r = allResults{i};
    mods(i) = r.M; IcNorm(i) = r.IcNorm; IcAbs(i) = r.Ic; Ir(i) = r.IrNorm;
    BER(i) = r.BER; SER(i) = r.SER; PD(i) = r.P_D; PSLR(i) = r.constraints.PSLR_dB;
    ISLR(i) = r.constraints.ISLR_dB; PAPR(i) = r.constraints.PAPR_dB; Obj(i) = r.objective;
end
fig = figure('Name', 'MI summary');
subplot(2,1,1); plot(mods, IcNorm, 'o-', mods, Ir, 's-', mods, Obj, '^-', 'LineWidth', 1.4);
grid on; xlabel('M'); ylabel('Normalized value'); legend('I_c norm', 'I_r norm', 'Objective', 'Location', 'best');
subplot(2,1,2); plot(mods, IcAbs, 'o-', 'LineWidth', 1.4);
grid on; xlabel('M'); ylabel('I_c (bit/use)'); title('Absolute communication mutual information');
save_figure_if_needed(fig, cfg, 'mi_summary.png');
fig = figure('Name', 'BER SER PD');
yyaxis left; semilogy(mods, max(BER,1e-5), 'o-', mods, max(SER,1e-5), 's-', 'LineWidth', 1.4);
ylabel('BER / SER'); ylim([1e-5 1]);
yyaxis right; plot(mods, PD, '^-', 'LineWidth', 1.4); ylabel('P_D'); ylim([0 1]);
grid on; xlabel('M'); legend('BER', 'SER', 'P_D', 'Location', 'best');
save_figure_if_needed(fig, cfg, 'ber_ser_pd.png');
fig = figure('Name', 'Radar waveform metrics');
bar(categorical(cellstr(string(mods))), [PSLR, ISLR, PAPR]);
grid on; ylabel('dB'); legend('PSLR', 'ISLR', 'PAPR', 'Location', 'best');
save_figure_if_needed(fig, cfg, 'waveform_metrics.png');
fig = figure('Name', 'Pd versus radar detection SNR');
hold on;
for i = 1:numel(allResults)
    r = allResults{i};
    plot(r.P_D_snr_db, r.P_D_curve, 'o-', 'LineWidth', 1.2, 'DisplayName', r.modulationName);
end
grid on; xlabel('Radar detection SNR (dB)'); ylabel('P_D'); ylim([0 1]);
legend('Location', 'best'); title('Detection probability under different radar SNRs');
save_figure_if_needed(fig, cfg, 'pd_vs_snr.png');
plot_optimized_spectrum_style(cfg, model, allResults);
for i = 1:numel(allResults)
    r = allResults{i};
    fig = figure('Name', sprintf('Expected autocorrelation M=%d', r.M));
    plot(r.constraints.lags, 20*log10(abs(r.constraints.R)+1e-12), 'LineWidth', 1.2);
    grid on; xlabel('Lag'); ylabel('Expected autocorr (dB)'); ylim([-70 5]);
    title(sprintf('Expected autocorrelation, %s', r.modulationName));
    save_figure_if_needed(fig, cfg, sprintf('autocorr_M%d.png', r.M));
end
end

function plot_optimized_spectrum_style(cfg, model, allResults)
fig = figure('Name', 'Optimized Magnitude Spectrum', 'Color', [0.94 0.92 0.88]);
set(fig, 'Position', [100 80 560 620]);
oldRng = rng;
cleanupObj = onCleanup(@() rng(oldRng)); %#ok<NASGU>
rng(2026);

N = model.N;
fNorm = ((0:N-1).' - floor(N/2)) / N;
lfm = sqrt(cfg.power.P_r) * model.s_r;
S_lfm = abs(fftshift(fft(lfm)));
scale = 500 / max(S_lfm);
lineColors = [0.0000 0.4470 0.7410;
              0.8500 0.3250 0.0980;
              0.9290 0.6940 0.1250;
              0.4940 0.1840 0.5560;
              0.4660 0.6740 0.1880];

ax1 = subplot(2, 1, 1);
plot(ax1, fNorm, scale * S_lfm, 'Color', lineColors(1,:), 'LineWidth', 0.9, ...
    'DisplayName', 'LFM-only');
hold(ax1, 'on');

ax2 = subplot(2, 1, 2);
hold(ax2, 'on');

for i = 1:numel(allResults)
    r = allResults{i};
    const = make_constellation(r.M);
    a = sample_symbols(const, model.K, 1);
    sComm = model.G * r.W * a;
    sJoint = lfm + sComm;
    cidx = 1 + mod(i, size(lineColors, 1) - 1);
    labelText = sprintf('%s opt', r.modulationName);

    S_joint = abs(fftshift(fft(sJoint)));
    S_comm = abs(fftshift(fft(sComm)));
    plot(ax1, fNorm, scale * S_joint, 'Color', lineColors(cidx,:), 'LineWidth', 0.75, ...
        'DisplayName', labelText);
    plot(ax2, fNorm, scale * S_comm, 'Color', lineColors(cidx,:), 'LineWidth', 0.75, ...
        'DisplayName', labelText);
end

format_spectrum_axis(ax1);
title(ax1, sprintf('Magnitude Spectrum; fixed f_\\Delta = %.0f%%, SRBR = %.2f', ...
    100*cfg.signal.f_delta, cfg.signal.SRBR), 'FontWeight', 'bold', 'FontSize', 9);
legend(ax1, 'Location', 'northeast', 'FontSize', 6);

format_spectrum_axis(ax2);
title(ax2, sprintf('Communication Embed Spectrum; fixed f_\\Delta = %.0f%%, SRBR = %.2f', ...
    100*cfg.signal.f_delta, cfg.signal.SRBR), 'FontWeight', 'bold', 'FontSize', 9);
legend(ax2, 'Location', 'northeast', 'FontSize', 6);

if cfg.output.saveFigures
    if ~exist(cfg.output.figuresDir, 'dir'), mkdir(cfg.output.figuresDir); end
    saveas(fig, fullfile(cfg.output.figuresDir, 'optimized_magnitude_spectrum.png'));
end
if ~cfg.output.showFigures, close(fig); end
end

function format_spectrum_axis(ax)
grid(ax, 'off');
box(ax, 'on');
xlim(ax, [-0.5 0.5]);
ylim(ax, [0 600]);
set(ax, 'XTick', [-0.5 0 0.5], 'YTick', 0:100:600, 'FontSize', 9, 'Color', 'w');
xlabel(ax, 'Frequency (Hz), normalized');
ylabel(ax, 'Magnitude');
end

function print_diagnostic_table(allResults)
fprintf('\n=== Diagnostic comparison ===\n');
fprintf('%6s %8s %8s %9s %9s %9s %9s %9s %9s %9s %8s %8s %8s\n', ...
    'Mod','pL','pR','Pc','RCRdB','Ic','IcNorm','IrNorm','PSLRdB','PAPRdB','BER','SER','Pd');
for i = 1:numel(allResults)
    r = allResults{i};
    fprintf('%6s %8.4f %8.4f %9.4f %9.2f %9.3f %9.4f %9.4f %9.2f %9.2f %8.3g %8.3g %8.3f\n', ...
        r.modulationName, r.xBest(1), r.xBest(2), r.Pc, r.RCR_dB, r.Ic, r.IcNorm, r.IrNorm, ...
        r.constraints.PSLR_dB, r.constraints.PAPR_dB, r.BER, r.SER, r.P_D);
end
fprintf('\nNote: Ic is absolute throughput in bit/use; IcNorm is normalized by K*log2(M).\n');
fprintf('Single-point Pd is for the configured radarDetectSnrDb only; use the Pd-vs-SNR figure for scenario trends.\n');
end

function out = merge_structs(a, b)
out = a; names = fieldnames(b);
for k = 1:numel(names), out.(names{k}) = b.(names{k}); end
end

function W = build_W(x, C)
pL = max(real(x(1)),0); pR = max(real(x(2)),0);
W = blkdiag(sqrt(pL)*exp(1j*x(3))*eye(C), sqrt(pR)*exp(1j*x(4))*eye(C));
end

function x = project_power_rcr(x, cfg, model)
x = real(x(:)); x(1:2) = max(x(1:2), 0); x(3:4) = mod(x(3:4), 2*pi);
pcLimit = min(cfg.power.P_c_max, cfg.power.P_r / (10^(cfg.power.RCR_min/10)));
pc = communication_power(x, model);
if pc > pcLimit && pc > 0, x(1:2) = x(1:2) * (pcLimit / pc); end
end

function Pc = communication_power(x, model)
W = build_W(x, model.C);
Pc = max(real(trace(model.G * W * W' * model.G')), 0);
end

function Ic = estimate_comm_mi_saa(const, W, cfg, model, Q, L)
A = model.Hc * model.G * W;
cholRv = chol(model.Rv, 'lower');
K = model.K; vals = zeros(Q,1);
for q = 1:Q
    a = sample_symbols(const, K, 1);
    v = complex_gaussian_cov(cholRv, 1);
    y = A*a + v;
    rTrue = cholRv \ (y - A*a);
    dTrue = real(rTrue'*rTrue);
    d = zeros(L,1);
    for ell = 1:L
        b = sample_symbols(const, K, 1);
        rb = cholRv \ (y - A*b);
        d(ell) = real(rb'*rb);
    end
    vals(q) = log2(L) - logsumexp_real(-(d-dTrue)) / log(2);
end
Ic = max(0, min(K*log2(const.M), mean(vals)));
end

function Ir = compute_radar_mi(W, cfg, model)
Rrc = model.Hrc * model.G * W * W' * model.G' * model.Hrc';
Reff = hermitian_regularize(cfg.mi.sigma_r2*eye(model.N) + Rrc);
L = chol(Reff, 'lower');
Swhite = hermitian_regularize(L \ model.Sg / L');
Ir = max(real(real_logdet(eye(model.N) + Swhite)), 0);
end

function met = compute_waveform_metrics(W, cfg, model, const, paprSamples)
R = expected_autocorr(W, cfg, model);
lags = (-(model.N-1):(model.N-1)).';
idx0 = find(lags == 0, 1);
mainIdx = abs(lags) <= cfg.constraints.mainlobeHalfWidth; sideIdx = ~mainIdx;
mainPeak = abs(R(idx0)) + 1e-12;
PSLR = max(abs(R(sideIdx))) / mainPeak;
ISLR = sum(abs(R(sideIdx)).^2) / (sum(abs(R(mainIdx)).^2) + 1e-12);
halfPowerIdx = abs(R) >= mainPeak/sqrt(2);
Bmain = max(lags(halfPowerIdx)) - min(lags(halfPowerIdx)) + 1;
paprMax = 0;
for q = 1:paprSamples
    a = sample_symbols(const, model.K, 1);
    s = sqrt(cfg.power.P_r)*model.s_r + model.G*W*a;
    paprMax = max(paprMax, max(abs(s).^2)/(mean(abs(s).^2)+1e-12));
end
met = struct('R', R, 'lags', lags, 'PSLR_dB', 20*log10(PSLR+1e-12), ...
    'ISLR_dB', 10*log10(ISLR+1e-12), 'PAPR', paprMax, ...
    'PAPR_dB', 10*log10(paprMax+1e-12), 'B_main', Bmain);
end

function R = expected_autocorr(W, cfg, model)
N = model.N; Qm = model.G * W * W' * model.G';
R = zeros(2*N-1,1);
for idx = 1:(2*N-1)
    lag = idx - N; Js = shift_matrix_sparse(N, lag);
    R(idx) = cfg.power.P_r*(model.s_r'*Js*model.s_r) + trace(Js*Qm);
end
R = R / (abs(R(N)) + 1e-12);
end

function [BER, SER] = evaluate_ber_ser(const, W, cfg, model, numTrials)
A = model.Hc * model.G * W; K = model.K; cholRv = chol(model.Rv, 'lower');
Aw = cholRv \ A; errBits = 0; totBits = 0; errSyms = 0; totSyms = 0;
for tr = 1:numTrials
    [a, idx, bits] = sample_symbols(const, K, 1);
    y = A*a + complex_gaussian_cov(cholRv, 1);
    aHatLin = pinv(Aw) * (cholRv \ y);
    idxHat = nearest_constellation(aHatLin, const.symbols);
    bitsHat = const.bits(idxHat, :); bitsTrue = bits(:,:,1);
    errBits = errBits + sum(bitsHat(:) ~= bitsTrue(:));
    totBits = totBits + numel(bitsTrue);
    errSyms = errSyms + sum(idxHat(:) ~= idx(:)); totSyms = totSyms + numel(idx);
end
BER = errBits/max(totBits,1); SER = errSyms/max(totSyms,1);
end

function P_D = evaluate_pd(W, cfg, model, numTrials, P_FA)
Rrc = model.Hrc * model.G * W * W' * model.G' * model.Hrc';
Reff = hermitian_regularize(cfg.mi.sigma_r2*eye(model.N) + Rrc);
cholR = chol(Reff, 'lower');
st = delay_vector(model.s_r, 2); st = st / max(norm(st), 1e-12);
alpha = sqrt(10^(cfg.eval.radarDetectSnrDb/10) * cfg.mi.sigma_r2);
den = max(real(st' * (Reff \ st)), 1e-12);
T1 = zeros(numTrials,1);
eta = -log(P_FA);
for k = 1:numTrials
    y1 = alpha*st + complex_gaussian_cov(cholR, 1);
    T1(k) = abs(st' * (Reff \ y1)).^2 / den;
end
P_D = mean(T1 > eta);
end

function const = make_constellation(M)
switch M
    case 2
        symbols = [-1; 1]; bits = [0; 1];
    case 4
        re = [-1; 1]; [I,Q] = ndgrid(re,re);
        symbols = (I(:)+1j*Q(:))/sqrt(2); bits = [I(:)>0, Q(:)>0];
    otherwise
        mSide = sqrt(M); levels = -(mSide-1):2:(mSide-1);
        [I,Q] = ndgrid(levels, levels);
        symbols = I(:)+1j*Q(:); symbols = symbols/sqrt(mean(abs(symbols).^2));
        bits = int_bits((0:M-1).', log2(M)); [~,order] = sortrows([I(:), Q(:)]);
        symbols = symbols(order); bits = bits(order,:);
end
const = struct('M', M, 'symbols', symbols(:), 'bits', logical(bits), ...
    'bitsPerSymbol', log2(M), 'name', modulation_name(M));
end

function name = modulation_name(M)
switch M
    case 2, name = 'BPSK';
    case 4, name = 'QPSK';
    otherwise, name = sprintf('%dQAM', M);
end
end

function [a, idx, bits] = sample_symbols(const, K, numCols)
if nargin < 3, numCols = 1; end
idx = randi(numel(const.symbols), K, numCols);
a = const.symbols(idx);
bits = false(K, const.bitsPerSymbol, numCols);
for n = 1:numCols, bits(:,:,n) = const.bits(idx(:,n), :); end
end

function s = generate_lfm_unit(N, Fs, B)
t = ((0:N-1).' - (N-1)/2) / Fs;
s = exp(1j*pi*(B/(N/Fs))*t.^2);
s = s / norm(s);
end

function P = build_pulse_matrix(N, Fs, Rs, beta, C, spanSymbols)
Tsym = 1/Rs; t = ((0:N-1).' - (N-1)/2)/Fs;
centers = ((1:C).' - (C+1)/2) * Tsym; P = zeros(N,C);
for c = 1:C
    tau = (t - centers(c)) / Tsym;
    pulse = raised_cosine_pulse(tau, beta);
    pulse(abs(tau) > spanSymbols/2) = 0;
    if norm(pulse) > 0, pulse = pulse/norm(pulse); end
    P(:,c) = pulse;
end
end

function y = raised_cosine_pulse(t, beta)
y = sinc_local(t) .* cos(pi*beta*t) ./ (1 - (2*beta*t).^2);
sing = abs(1 - (2*beta*t).^2) < 1e-8;
if beta > 0, y(sing) = (pi/4) * sinc_local(1/(2*beta)); end
y(~isfinite(y)) = 0;
end

function H = build_radar_comm_coupling(N, gain, G, Xr)
n = (0:N-1).';
phase1 = exp(1j*2*pi*0.073*n);
phase2 = exp(1j*(2*pi*0.137*n + 0.35*sin(2*pi*n/N)));
H = 0.45*eye(N) + 0.35*diag(phase1)*delay_matrix(N,2) + 0.20*diag(phase2)*delay_matrix(N,-3);
Qg = orth(full(G)); Qx = orth(full(Xr)); r = min(size(Qg,2), size(Qx,2));
if r > 0, H = 0.35*H + 0.65*(Qx(:,1:r)*Qg(:,1:r)'); end
H = gain * H / max(norm(H,'fro')/sqrt(N), 1e-12);
end

function idx = nearest_constellation(z, symbols)
idx = zeros(numel(z),1);
for k = 1:numel(z), [~,idx(k)] = min(abs(z(k)-symbols).^2); end
end

function x = complex_gaussian_cov(cholCov, numCols)
z = (randn(size(cholCov,1), numCols) + 1j*randn(size(cholCov,1), numCols))/sqrt(2);
x = cholCov * z;
end

function J = shift_matrix_sparse(N, lag)
if lag >= 0
    rows = (1+lag):N; cols = 1:(N-lag);
else
    rows = 1:(N+lag); cols = (1-lag):N;
end
J = sparse(rows, cols, 1, N, N);
end

function D = delay_matrix(N, lag)
D = shift_matrix_sparse(N, lag);
end

function y = delay_vector(x, lag)
x = x(:); N = numel(x); y = zeros(N,1);
if lag >= 0, y((1+lag):N) = x(1:(N-lag));
else, lag = abs(lag); y(1:(N-lag)) = x((1+lag):N); end
end

function y = real_logdet(A)
A = hermitian_regularize(A);
[R,p] = chol(A);
if p == 0
    y = 2*sum(log2(abs(diag(R))+realmin));
else
    y = sum(log2(max(real(eig(A)), realmin)));
end
y = real(y);
end

function A = hermitian_regularize(A)
A = (A + A')/2 + 1e-10*eye(size(A));
end

function y = logsumexp_real(x)
x = real(x(:)); m = max(x); y = m + log(sum(exp(x-m)));
end

function y = sinc_local(x)
y = ones(size(x)); idx = abs(x) > 1e-12;
y(idx) = sin(pi*x(idx))./(pi*x(idx));
end

function bits = int_bits(vals, width)
bits = false(numel(vals), width);
for b = 1:width, bits(:, width-b+1) = bitget(vals, b) > 0; end
end

function save_figure_if_needed(fig, cfg, filename)
if cfg.output.saveFigures
    if ~exist(cfg.output.figuresDir, 'dir'), mkdir(cfg.output.figuresDir); end
    saveas(fig, fullfile(cfg.output.figuresDir, filename));
end
if ~cfg.output.showFigures, close(fig); end
end
