% MATLAB_OFDM_64QAM_Simulation_CP.m
% Comprehensive simulation of time-frequency synchronization and equalization
% for a 64-QAM, 20 MHz-band, 64-subcarrier OFDM link (CP = 16 samples).
% This version uses an 802.11a-style preamble/pilot layout and includes:
%   - Schmidl & Cox timing sync (variance vs Eb/N0)
%   - Two-stage CFO estimator (fractional from long training + integer via correlation)
%   - 5-tap Rayleigh channel ~200 ns RMS delay
%   - Pilot-assisted LS vs MMSE channel estimation (MSE vs Eb/N0)
%   - Per-tone MMSE and ZF equalizers (BER vs Eb/N0)  <-- with CP correctly handled
%   - BER vs Doppler (ICI) at Eb/N0 = 15 dB
%
% Author: (updated) 2025-08-19

clear; clc; close all; rng(0);

%% System parameters
Fs   = 20e6;      % sampling rate
N    = 64;        % FFT size
Ncp  = 16;        % cyclic prefix
df   = Fs/N;      % subcarrier spacing
Tfft = N/Fs;
Tcp  = Ncp/Fs;
Tsym = Tfft + Tcp;
k    = 6;         % bits per 64-QAM
M    = 64;

% S&C parameters (802.11a-like long training made of two identical halves)
L          = 32;    % half length (two identical 32-sample halves => 64 total)
short_len  = 16;    % short symbol length
num_short  = 10;    % number of short symbols (10 x 16 = 160 samples)

% Pilot & data allocation per 802.11a
active_pos = [-26:-1, 1:26];       % 52 active tones (DC excluded)
idx_active = 33 + active_pos;      % MATLAB indices
pilot_rel  = [-21 -7 7 21];        % pilot tones
pilot_idx  = 33 + pilot_rel;       % pilot indices
pilot_pattern = [1 1 1 -1].';      % BPSK pilot pattern (column)
pilot_symbol_freq = zeros(N,1);
pilot_symbol_freq(pilot_idx) = pilot_pattern;
data_idx = setdiff(idx_active, pilot_idx).';

% Simulation control
EbN0_timing = 5:1:15;
EbN0_ber    = 0:2:20;
nFrames_timing = 2000;
nFrames_ber    = 2000;

% Channel PDP target
tau_rms_target = 200e-9;                 % 200 ns
tap_delays_s   = [0 1 2 4 8]*(1/Fs);     % [0 50ns 100ns 200ns 400ns]

% Compute exponential PDP with target RMS delay
alpha = find_alpha_for_rms(tap_delays_s, tau_rms_target);
p = exp(-alpha*tap_delays_s); p = p/sum(p);

fprintf('Assumed carrier frequency (for Doppler calcs) = 2.4 GHz (GUESS).');
fprintf('Tap delays (s): %s', mat2str(tap_delays_s));
fprintf('Tap powers (normalized): %s', mat2str(p.',6));

%% Helpers
ifft64 = @(X) ifft(X, N);
fft64  = @(x) fft(x, N);

add_cp    = @(x) [x(end-Ncp+1:end); x];
remove_cp = @(x) x(Ncp+1:Ncp+N);

map64qam   = @(bits) qammod(bits, M, 'InputType','bit', 'UnitAveragePower',true);
demap64qam = @(sym)  qamdemod(sym, M, 'OutputType','bit','UnitAveragePower',true);

%% Construct 802.11a-style preamble
% Short training: synthesize a 16-sample QPSK tile and repeat 10x
short_sym_time     = exp(1j*pi/2*(randi([0 3], short_len, 1)));
short_preamble_time= repmat(short_sym_time, num_short, 1); % 160 samples

% Long training: two identical halves (32 + 32 = 64)
long_half       = exp(1j*2*pi*rand(L,1));
long_symbol_time= [long_half; long_half];                   % 64 samples

% Full preamble: short (160) + long (64)
preamble_time = [short_preamble_time; long_symbol_time];

% Block pilot OFDM symbol (frequency)
block_pilot_freq = zeros(N,1);
block_pilot_freq(data_idx)  = 1;               % simple +1 on all data tones
block_pilot_freq(pilot_idx) = pilot_pattern;
pilot_time_no_cp = ifft64(block_pilot_freq);
pilot_time       = add_cp(pilot_time_no_cp);   % *** CP added ***

%% 1) Timing sync variance (Schmidl & Cox) using long training (L=32)
fprintf('Running timing variance experiment (Schmidl & Cox with 802.11a-style preamble)...');
timing_var = zeros(size(EbN0_timing));
for idx=1:length(EbN0_timing)
    EbN0dB = EbN0_timing(idx);
    sigma2 = compute_noise_variance(EbN0dB, k, N, Ncp);
    errors = zeros(nFrames_timing,1);

    for f=1:nFrames_timing
        % one random data OFDM symbol (not used by S&C but keeps frame realistic)
        Xd = zeros(N,1);
        bits_tmp = randi([0 1], numel(data_idx)*k, 1);
        Xd(data_idx) = map64qam(bits_tmp);
        data_time_no_cp = ifft64(Xd);
        data_time       = add_cp(data_time_no_cp);

        tx_frame = [preamble_time; pilot_time; data_time];
        rx = tx_frame + sqrt(sigma2/2)*(randn(size(tx_frame))+1j*randn(size(tx_frame)));

        % S&C metric around expected long start
        search_start = max(1, length(short_preamble_time)-20);
        search_end   = search_start + 200;
        Mmet = zeros(search_end,1);
        for d = search_start:search_end-2*L
            r1 = rx(d:d+L-1);
            r2 = rx(d+L:d+2*L-1);
            P  = sum(conj(r1).*r2);
            R  = sum(abs(r2).^2);
            Mmet(d) = abs(P)^2/(R^2 + eps);
        end
        [~, d_est] = max(Mmet);
        true_long_start = length(short_preamble_time)+1;
        errors(f) = d_est - true_long_start;
    end
    timing_var(idx) = var(errors);
    fprintf('Eb/N0= %2.1f dB: timing variance = %.4f (samples^2)', EbN0dB, timing_var(idx));
end

figure; plot(EbN0_timing, timing_var, '-o'); grid on;
xlabel('Eb/N0 (dB)'); ylabel('Timing estimation variance (samples^2)');
title('Schmidl & Cox timing detection variance vs Eb/N0 (802.11a-style preamble)');

%% 2) CFO estimator test: fractional (S&C long) + integer (pilot corr)
fprintf('Testing CFO estimator (fractional + integer) with 802.11a-style preamble...');
EbN0_test = [5 10 15];
cfo_true_list = [-0.4 -0.1 0.05 0.3]; % in subcarrier units
residuals = zeros(length(EbN0_test), length(cfo_true_list));

for ie=1:length(EbN0_test)
    EbN0dB = EbN0_test(ie);
    sigma2 = compute_noise_variance(EbN0dB, k, N, Ncp);
    for ic=1:length(cfo_true_list)
        eps_true_sc = cfo_true_list(ic);    % subcarrier units
        eps_true_cs = eps_true_sc / N;      % cycles per sample
        ntrials = 500; resid = zeros(ntrials,1);
        for t=1:ntrials
            tx_frame = [preamble_time; pilot_time];
            n = (0:length(tx_frame)-1).';
            rx = tx_frame .* exp(1j*2*pi*eps_true_cs*n);
            rx = rx + sqrt(sigma2/2)*(randn(size(rx))+1j*randn(size(rx)));

            % Assume timing found; long start:
            d = length(short_preamble_time)+1;
            % Fractional CFO from long halves
            r1 = rx(d:d+L-1); r2 = rx(d+L:d+2*L-1);
            P  = sum(conj(r1).*r2);
            eps_frac_hat_cs = angle(P)/(2*pi*L); % cycles/sample

            % Correct fractional
            rx_corr = rx .* exp(-1j*2*pi*eps_frac_hat_cs*(0:length(rx)-1).');

            % Integer CFO via correlation on pilot OFDM (strip CP then FFT)
            pilot_pos = d + 2*L; % pilot starts right after 64-sample long
            pilot_rx_cp = rx_corr(pilot_pos : pilot_pos+Ncp+N-1);
            pilot_rx_td = remove_cp(pilot_rx_cp);
            Y = fft64(pilot_rx_td);
            ref = block_pilot_freq;

            corr_vals = zeros(N,1);
            for shift = 0:N-1
                ref_shift = circshift(ref, shift);
                corr_vals(shift+1) = abs(ref_shift'*conj(Y));
            end
            [~, shift_hat] = max(corr_vals);
            shift_hat = shift_hat - 1;

            eps_hat_total_cs = eps_frac_hat_cs + shift_hat/N;   % cycles/sample
            % Residual in subcarrier units:
            resid(t) = abs(eps_true_sc - eps_hat_total_cs*N);
        end
        residuals(ie,ic) = mean(resid);
        fprintf('EbN0=%d dB, true eps=%.3f subcarrier: residual (avg)=%.4f subcarrier units', ...
                 EbN0dB, eps_true_sc, residuals(ie,ic));
    end
end
fprintf('Fractional+Integer residuals (avg):'); disp(residuals);

%% 3) Channel estimation MSE: LS vs MMSE (block pilot), with CP handled
fprintf('Running channel estimation MSE comparison (LS vs MMSE) with 802.11a pilot layout...');
EbN0_list = 0:2:20;
mse_ls   = zeros(size(EbN0_list));
mse_mmse = zeros(size(EbN0_list));

% Channel covariance across tones from PDP
Rhh = channel_freq_covariance_from_pdp(p, tap_delays_s, N, Fs);

for ii=1:length(EbN0_list)
    EbN0dB = EbN0_list(ii);
    sigma2 = compute_noise_variance(EbN0dB, k, N, Ncp);
    ntrials = 1000; se_ls=0; se_mmse=0;

    for t=1:ntrials
        % draw random channel taps
        h_time = (randn(length(p),1)+1j*randn(length(p),1))/sqrt(2) .* sqrt(p(:));
        H_freq = fft([h_time; zeros(N-length(h_time),1)], N);

        % Transmit single pilot OFDM symbol with CP
        tx = pilot_time;
        % Convolution through channel
        rx = conv(tx, [h_time; zeros(N-length(h_time),1)]);
        rx = rx(1:length(tx));
        % AWGN
        rx = rx + sqrt(sigma2/2)*(randn(size(rx))+1j*randn(size(rx)));

        % Receiver: remove CP then FFT
        y_td = remove_cp(rx(1:Ncp+N));
        Y    = fft64(y_td);

        Xp = block_pilot_freq;           % known pilot symbol (freq)
        H_true = H_freq;

        % LS on pilot positions
        Hls_pilots = Y(pilot_idx)./Xp(pilot_idx);
        % Interpolate (linear) across all tones
        Hls_full = interp1(pilot_idx, Hls_pilots, 1:N, 'linear', 'extrap').';

        % MMSE (Wiener) interpolation from pilots to all tones
        Rp   = Rhh(pilot_idx, pilot_idx);
        Rhp  = Rhh(:, pilot_idx);
        P    = numel(pilot_idx);
        Xp_mat = diag(Xp(pilot_idx));
        A    = Xp_mat*Rp*Xp_mat' + sigma2*eye(P);
        W    = Rhp / A * Xp_mat;              % N x P
        Hmmse = (W * Hls_pilots).';

        se_ls   = se_ls   + mean(abs(H_true - Hls_full).^2);
        se_mmse = se_mmse + mean(abs(H_true - Hmmse.').^2);
    end

    mse_ls(ii)   = se_ls/ntrials;
    mse_mmse(ii) = se_mmse/ntrials;
    fprintf('EbN0=%d dB: MSE_LS=%.4e, MSE_MMSE=%.4e', EbN0dB, mse_ls(ii), mse_mmse(ii));
end

figure; semilogy(EbN0_list, mse_ls, '-o', EbN0_list, mse_mmse, '-x'); grid on;
legend('LS','MMSE'); xlabel('Eb/N0 (dB)'); ylabel('MSE');
title('Channel estimation MSE: LS vs MMSE (802.11a pilots, CP handled)');

%% 4) BER simulations (MMSE vs ZF) — with CP on TX, strip CP at RX
fprintf('Running BER simulations (MMSE equalizer recommended). This may take longer...');
BER_mmse = zeros(size(EbN0_ber));
BER_zf   = zeros(size(EbN0_ber));

for ii=1:length(EbN0_ber)
    EbN0dB = EbN0_ber(ii);
    sigma2 = compute_noise_variance(EbN0dB, k, N, Ncp);
    nErr_mmse=0; nErr_zf=0; nBitsTotal=0;

    for f=1:nFrames_ber
        % Bits & mapping on data tones
        nData = numel(data_idx);
        bits = randi([0 1], nData*k, 1);
        symbols = map64qam(bits);

        X = zeros(N,1);
        X(data_idx)  = symbols;
        X(pilot_idx) = pilot_pattern;

        data_time_no_cp = ifft64(X);
        data_time       = add_cp(data_time_no_cp);

        tx_frame = [preamble_time; pilot_time; data_time];

        % Rayleigh channel (draw new per frame)
        h_time = (randn(length(p),1)+1j*randn(length(p),1))/sqrt(2) .* sqrt(p(:));

        % Apply channel (linear convolution)
        rx = conv(tx_frame, [h_time; zeros(N-length(h_time),1)]);
        rx = rx(1:length(tx_frame));

        % Add small random CFO (subcarrier units)
        eps_sc = (randn*0.02);
        n = (0:length(rx)-1).';
        rx = rx .* exp(1j*2*pi*(eps_sc/N)*n);

        % AWGN
        rx = rx + sqrt(sigma2/2)*(randn(size(rx))+1j*randn(size(rx)));

        % Receiver: assume long start known (or a coarse detector)
        long_start = length(short_preamble_time) + 1;
        pilot_pos  = long_start + 2*L;   % pilot AFTER 64-sample long

        % --- Pilot: strip CP and FFT
        pilot_rx_cp = rx(pilot_pos : pilot_pos+Ncp+N-1);
        pilot_rx_td = remove_cp(pilot_rx_cp);
        Y_pilot     = fft64(pilot_rx_td);

        % LS channel estimate on pilots + interpolation
        Hls_pilots = Y_pilot(pilot_idx)./block_pilot_freq(pilot_idx);
        Hls_full   = interp1(pilot_idx, Hls_pilots, 1:N, 'linear', 'extrap').';

        % Per-tone equalizers
        H_est  = Hls_full;
        G_zf   = 1./(H_est + 1e-12);
        G_mmse = conj(H_est)./(abs(H_est).^2 + sigma2);

        % --- Data: strip CP and FFT
        data_pos     = pilot_pos + (Ncp + N);
        data_rx_cp   = rx(data_pos : data_pos+Ncp+N-1);
        data_rx_td   = remove_cp(data_rx_cp);
        Y_data       = fft64(data_rx_td);

        % Equalize
        eq_zf   = G_zf   .* Y_data;
        eq_mmse = G_mmse .* Y_data;

        rx_symbols_zf   = eq_zf(data_idx);
        rx_symbols_mmse = eq_mmse(data_idx);

        % Demap
        bits_hat_zf   = demap64qam(rx_symbols_zf);
        bits_hat_mmse = demap64qam(rx_symbols_mmse);

        % Count errors
        nErr_zf    = nErr_zf   + sum(bits_hat_zf   ~= bits);
        nErr_mmse  = nErr_mmse + sum(bits_hat_mmse ~= bits);
        nBitsTotal = nBitsTotal + length(bits);
    end

    BER_zf(ii)   = nErr_zf   / nBitsTotal;
    BER_mmse(ii) = nErr_mmse / nBitsTotal;
    fprintf('Eb/N0=%d dB: BER_ZF=%.3e, BER_MMSE=%.3e', EbN0dB, BER_zf(ii), BER_mmse(ii));
end

% AWGN "ballpark" bound for 64-QAM (uncoded)
EbN0_lin = 10.^(EbN0_ber/10);
EsN0_lin = EbN0_lin * k;                % ignoring CP overhead in this bound
ser_awgn = 4*(1-1/sqrt(M)) .* qfunclocal(sqrt(3/(M-1)*EsN0_lin));
ber_awgn_approx = ser_awgn / k;

figure; semilogy(EbN0_ber, BER_mmse,'-o', EbN0_ber, BER_zf,'-x', EbN0_ber, ber_awgn_approx,'--'); grid on;
legend('MMSE (sim)','ZF (sim)','AWGN (approx)');
xlabel('Eb/N0 (dB)'); ylabel('BER');
title('BER: uncoded 64-QAM OFDM (802.11a pilots/preamble, CP handled)');

%% 5) BER vs Doppler (ICI) at Eb/N0 = 15 dB
fprintf('Testing ICI impact at Doppler = 100 Hz and higher values...');
fc = 2.4e9;
fD_list = [0 100 500 1000 2000];
BER_dop  = zeros(size(fD_list));
EbN0dB = 15; sigma2 = compute_noise_variance(EbN0dB, k, N, Ncp);

for ii=1:length(fD_list)
    fD = fD_list(ii);
    nErr=0; nBits=0; nRuns=500;
    rho = besselj(0, 2*pi*fD*Tsym); % symbol-to-symbol correlation

    for r=1:nRuns
        nData = numel(data_idx);
        bits  = randi([0 1], nData*k, 1);
        symbols = map64qam(bits);

        X = zeros(N,1); X(data_idx)=symbols; X(pilot_idx)=pilot_pattern;
        data_time_no_cp = ifft64(X);
        data_time       = add_cp(data_time_no_cp);
        tx_frame = [preamble_time; pilot_time; data_time];

        % Time-varying channel (first-order Gauss-Markov with corr rho)
        h_prev = (randn(length(p),1)+1j*randn(length(p),1))/sqrt(2).*sqrt(p(:));
        h_time = rho*h_prev + sqrt(max(0,1-rho^2))*(randn(length(p),1)+1j*randn(length(p),1))/sqrt(2).*sqrt(p(:));

        rx = conv(tx_frame, [h_time; zeros(N-length(h_time),1)]);
        rx = rx(1:length(tx_frame));

        rx = rx + sqrt(sigma2/2)*(randn(size(rx))+1j*randn(size(rx)));

        long_start = length(short_preamble_time) + 1;
        pilot_pos  = long_start + 2*L;

        % Pilot
        pilot_rx_cp = rx(pilot_pos : pilot_pos+Ncp+N-1);
        pilot_rx_td = remove_cp(pilot_rx_cp);
        Y_pilot     = fft64(pilot_rx_td);
        Hls_pilots  = Y_pilot(pilot_idx)./block_pilot_freq(pilot_idx);
        H_est       = interp1(pilot_idx, Hls_pilots, 1:N, 'linear', 'extrap').';

        % Data
        data_pos   = pilot_pos + (Ncp + N);
        data_rx_cp = rx(data_pos : data_pos+Ncp+N-1);
        data_rx_td = remove_cp(data_rx_cp);
        Y_data     = fft64(data_rx_td);

        G_mmse = conj(H_est)./(abs(H_est).^2 + sigma2);
        eq_mmse = G_mmse .* Y_data;
        rx_symbols_mmse = eq_mmse(data_idx);             % *** fixed indexing ***
        bits_hat = demap64qam(rx_symbols_mmse);

        nErr  = nErr + sum(bits_hat ~= bits);
        nBits = nBits + length(bits);
    end
    BER_dop(ii) = nErr/nBits;
    fprintf('Doppler=%d Hz: BER=%.3e', fD, BER_dop(ii));
end

figure; semilogy(fD_list, BER_dop, '-o'); grid on;
xlabel('Doppler (Hz)'); ylabel('BER at Eb/N0=15 dB');
title('BER vs Doppler (ICI floor behavior)');

%% ---------------- Local helper functions ----------------
function alpha = find_alpha_for_rms(tau_vec, tau_rms_target)
    obj = @(a) (compute_rms_for_alpha(a, tau_vec) - tau_rms_target).^2;
    a_low = 1e2; a_high = 1e10;
    try
        alpha = fminbnd(obj, a_low, a_high, optimset('TolX',1e-6,'Display','off'));
        if ~isfinite(alpha) || alpha<=0, error('bad alpha'); end
    catch
        grid = logspace(log10(a_low), log10(a_high), 200);
        vals = arrayfun(@(a) obj(a), grid);
        [~, idxmin] = min(vals); alpha = grid(idxmin);
    end
end

function rms = compute_rms_for_alpha(alpha, tau_vec)
    p = exp(-alpha * tau_vec);
    p = p / sum(p);
    tau_mean = sum(p .* tau_vec);
    rms = sqrt(sum(p .* (tau_vec - tau_mean).^2));
end

function sigma2 = compute_noise_variance(EbN0dB, k, N, Ncp)
    EbN0 = 10^(EbN0dB/10);
    EsN0 = EbN0 * k;
    useful_frac = N/(N+Ncp);      % CP overhead
    EsN0_eff = EsN0 * useful_frac;
    sigma2 = 1/max(EsN0_eff, 1e-12); % noise variance per complex sample
end

function Rhh = channel_freq_covariance_from_pdp(p, delays, N, Fs)
    % R_HH[i,j] = sum_t p_t * exp(-j*2*pi*(f_i - f_j)*tau_t)
    Rhh = zeros(N,N);
    for i=1:N
        for j=1:N
            df_ij = (i-j)*Fs/N; % (i-j)*Δf
            s = 0;
            for t=1:length(p)
                s = s + p(t)*exp(-1j*2*pi*df_ij*delays(t));
            end
            Rhh(i,j) = s;
        end
    end
end

function y = qfunclocal(x)
    y = 0.5*erfc(x./sqrt(2));
end
