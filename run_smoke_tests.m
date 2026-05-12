%% run_smoke_tests.m
% Lightweight checks for the fixed guard-band FA-DMI ISAC code.
clear; clc;

rootDir = fileparts(mfilename('fullpath'));
addpath(rootDir);
rng(7);

cfg = isac_core('default_config', rootDir);
cfg.mi.Q = 4;
cfg.mi.L = 8;
cfg.eval.paprSamples = 4;
model = isac_core('build_model', cfg);

assert(model.N == cfg.signal.N);
assert(model.C == cfg.signal.C);
assert(model.K == 2 * cfg.signal.C);
assert(all(size(model.G) == [cfg.signal.N, model.K]));

for M = cfg.modulation.orders
    const = isac_core('make_constellation', M);
    avgPow = mean(abs(const.symbols).^2);
    assert(abs(avgPow - 1) < 1e-12, 'Constellation is not unit power.');

    x = [0.05; 0.04; 0.1; 0.2];
    W = isac_core('build_W', x, model.C);
    assert(all(size(W) == [model.K, model.K]));

    met = isac_core('evaluate_design', M, const, x, cfg, model, cfg.mi.Q, cfg.mi.L, cfg.eval.paprSamples);
    assert(isfinite(met.objective));
    assert(isfinite(met.Pc));
    assert(isfinite(met.constraints.PSLR_dB));
    assert(isfinite(met.constraints.ISLR_dB));
    assert(isfinite(met.constraints.PAPR));
end

fprintf('Smoke tests passed for M = %s.\n', mat2str(cfg.modulation.orders));

