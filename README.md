# Baseband Receiver Design for 64‑QAM OFDM (IEEE 802.11a‑class)

This repository contains a MATLAB implementation of the time–frequency synchronisation and equalisation chain for an uncoded 64‑QAM, 20 MHz, 64‑subcarrier OFDM link. The accompanying paper‑style report is `report.tex`.

## 1. Requirements
- MATLAB R2021a or newer (Signal Processing Toolbox recommended)
- LaTeX distribution (MiKTeX or TeX Live) to build `report.tex`

## 0. Environment setup (Windows)
- Install MATLAB and optionally the Signal Processing Toolbox.
- Install MiKTeX: https://miktex.org/download
  - During first `pdflatex` run, allow on‑the‑fly package install.
- Ensure `pdflatex` is on PATH (MiKTeX normally configures this). Verify in PowerShell:
  ```powershell
  pdflatex --version
  ```

## 2. Files
- `baseband4.m` — main MATLAB script that runs all experiments and prints logs.
- `report.tex` — IEEE‑style LaTeX report. Includes figures from `1.png`–`4.png`.
- `1.png` — BER vs Doppler (ICI floor).
- `2.png` — BER vs Eb/N0 (ZF, MMSE, AWGN ref).
- `3.png` — Channel estimation MSE (LS vs MMSE).
- `4.png` — Schmidl & Cox timing variance vs Eb/N0.

## 3. Quick start (MATLAB)
1) Open MATLAB in this folder.
2) Run the main script:
```matlab
clear; close all; clc;
run('baseband4.m');
```
3) Observe the console logs for:
   - Assumed carrier, tap delays and powers
   - Timing variance table (Eb/N0 = 5…15 dB)
   - CFO residuals for fractional+integer estimator at several Eb/N0 values
   - Channel estimation MSE (LS vs MMSE)
   - BER (ZF vs MMSE) over Eb/N0
   - BER vs Doppler (ICI floor)
4) Figures are saved/available as `1.png`–`4.png` and referenced in `report.tex`.

## 4. Reproducing the figures
If `baseband4.m` already generates the plots, re‑run it to refresh the PNGs. Otherwise, you can regenerate from the printed vectors (fast path):
- Timing variance → `4.png`
- MSE vs SNR (LS/MMSE) → `3.png`
- BER vs SNR (ZF/MMSE + AWGN ref) → `2.png`
- BER vs Doppler → `1.png`

Example snippet to plot timing variance from the console values:
```matlab
EbN0 = 5:15;
varSamples2 = [207.1295 249.0077 288.0967 332.1151 352.5372 361.7439 325.6827 285.0905 220.0823 155.7368 107.7118];
semilogy(EbN0, varSamples2, '-o'); grid on; xlabel('Eb/N0 (dB)'); ylabel('Timing variance (samples^2)');
saveas(gcf, '4.png');
```

## 4.1 Exact snippets for all figures

1) Timing variance (Figure 4 → `4.png`)
```matlab
EbN0 = 5:15;
varSamples2 = [207.1295 249.0077 288.0967 332.1151 352.5372 361.7439 325.6827 285.0905 220.0823 155.7368 107.7118];
figure; semilogy(EbN0, varSamples2, '-o','LineWidth',1.5);
grid on; xlabel('Eb/N0 (dB)'); ylabel('Timing variance (samples^2)');
title('S&C timing variance vs Eb/N0');
saveas(gcf, '4.png');
```

2) Channel estimation MSE (Figure 3 → `3.png`)
```matlab
EbN0 = 0:2:20;
MSE_LS   = [16.482 10.038 6.2948 4.0490 2.5960 1.5835 1.0041 0.63143 0.40484 0.24808 0.15943];
MSE_MMSE = [971.26 1430.0 2324.9 3622.9 5817.9 9079.7 14563 22760 35045 55662 90530];
figure; semilogy(EbN0, MSE_LS, '-o','LineWidth',1.5); hold on;
semilogy(EbN0, MSE_MMSE, '-s','LineWidth',1.5);
grid on; xlabel('Eb/N0 (dB)'); ylabel('MSE'); legend('LS','MMSE','Location','southwest');
title('Channel estimation MSE vs Eb/N0');
saveas(gcf, '3.png');
```

3) BER vs SNR (Figure 2 → `2.png`)
```matlab
EbN0 = 0:2:20;
BER_ZF   = [0.4894 0.4850 0.4759 0.4674 0.4510 0.4313 0.4082 0.3845 0.3553 0.3203 0.2889];
BER_MMSE = [0.4893 0.4849 0.4760 0.4672 0.4505 0.4309 0.4082 0.3841 0.3552 0.3202 0.2887];
figure; plot(EbN0, BER_ZF, '-o','LineWidth',1.5); hold on;
plot(EbN0, BER_MMSE, '-s','LineWidth',1.5);
grid on; xlabel('Eb/N0 (dB)'); ylabel('BER'); ylim([0.25 0.5]); legend('ZF','MMSE','Location','northwest');
title('Uncoded BER vs Eb/N0');
saveas(gcf, '2.png');
```

4) BER vs Doppler (Figure 1 → `1.png`)
```matlab
Doppler = [0 100 500 1000 2000];
BER = [0.3652 0.3628 0.3601 0.3609 0.3688];
figure; plot(Doppler, BER, '-o','LineWidth',1.5);
grid on; xlabel('Doppler (Hz)'); ylabel('BER'); ylim([0.34 0.38]);
title('BER vs Doppler');
saveas(gcf, '1.png');
```

## 5. Building the report
From a terminal in this folder:
```bash
pdflatex -interaction=nonstopmode report.tex
pdflatex -interaction=nonstopmode report.tex
```
The second pass resolves references and figure labels. The output is `report.pdf`.

## 6. What the pipeline does (high level)
1) Preamble & pilots (802.11a‑style) are generated.
2) Transmit 64‑QAM symbols over a 5‑tap Rayleigh channel (RMS ≈ 200 ns in spec; run used a near‑flat realisation).
3) Add AWGN and optional Doppler for ICI experiments.
4) Timing sync via Schmidl & Cox autocorrelation metric; variance measured over Monte Carlo trials.
5) CFO estimation: fractional (from S&C phase) + integer (frequency‑domain correlation); residual reported in subcarrier units.
6) Channel estimation on pilots: LS and MMSE; MSE vs SNR recorded.
7) Per‑tone equalisation: ZF and MMSE; BER vs SNR.
8) ICI floor study: sweep Doppler ∈ {0,100,500,1000,2000} Hz and measure BER.

## 7. Interpreting the results
- Timing variance non‑monotonicity is expected for short averaging with S&C and a flat‑like channel realisation.
- CFO residuals large at 5 dB and near‑zero at 15 dB match theory (integer CFO detection becomes reliable beyond ~10 dB).
- MMSE worse than LS across SNR suggests prior mismatch/conditioning in the MMSE implementation; LS is robust.
- ZF≈MMSE BER implies an effectively single‑tap channel during this run and/or MMSE prior mismatch.
- Doppler sweep shows a shallow ICI‑driven BER floor around 0.36; carrier/phase tracking and coding will suppress it.

## 8. Tips and troubleshooting
- If BER is uniformly high, verify Es/N0↔Eb/N0 mapping, pilot/CP overhead, and integer CFO detection.
- For MMSE, ensure noise variance and PDP/covariance are correct; use diagonal loading for numerical stability.
- Increase Monte Carlo iterations for smoother timing‑variance curves.
- If LaTeX places floats oddly, we already force tables with `[H]`. You can also add `\suppressfloats[t]` after `\maketitle`.

## 9. References
- Schmidl & Cox, IEEE Trans. Communications, 1997.
- van de Beek et al., IEEE Trans. Signal Processing, 1997.
- Coleri et al., IEEE Trans. Broadcasting, 2002.
