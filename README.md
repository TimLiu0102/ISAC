# Fixed Guard-Band FA-DMI ISAC MATLAB Simulation

This folder contains a self-contained MATLAB R2020a implementation for multi-modulation LFM dual-function radar-communication waveform simulation under fixed guard-band constraints.

## Five-File Layout

- `main_isac_dmi.m`: runs the full BPSK/QPSK/16QAM/64QAM simulation, optimization, evaluation, result saving and plotting.
- `run_extended_experiments.m`: runs the LFM-only baseline, radar-SNR detection curve, communication-SNR sweep, residual-interference sweep, fixed-symbol/fixed-bit-rate comparison, and constraint-sensitivity checks.
- `run_smoke_tests.m`: checks constellation normalization, matrix dimensions, waveform metrics and one objective evaluation.
- `isac_core.m`: shared implementation library and dispatcher.
- `README.md`: this note.

## Outputs

Running `main_isac_dmi.m` creates:

- `results/results_M*.mat`: per-modulation optimal embedding and metrics.
- `results/results_all.mat`: all modulation results.
- `figures/*.png`: MI, BER/SER/P_D, waveform-metric and spectrum/autocorrelation figures.

## Notes

- The default configuration is intentionally lightweight so the full chain runs quickly.
- Increase `cfg.mi.Q`, `cfg.mi.L`, `cfg.opt.pGrid`, `cfg.eval.numBerTrials`, and `cfg.eval.numPdTrials` in `default_isac_config.m` for publication-grade Monte Carlo runs.
- No Communications Toolbox functions are required.
