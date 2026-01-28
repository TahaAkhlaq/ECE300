% Taha Akhlaq - Communication Theory Problem Set 6
clc; clear; close all;

q = @(x) 0.5 * erfc(x / sqrt(2));

%% Problem 1

% Parameters
SNR_dB = 4;
gamma_b = 10^(SNR_dB/10); 

% Case 1: Coherent Orthogonal 
P_coh1 = q(sqrt(gamma_b));

% Case 2: Coherent Optimal FSK 
rho = -0.217;
P_coh2 = q(sqrt((1 - rho) * gamma_b));

% Case 1: Noncoherent Orthogonal
P_noncoh1 = 0.5 * exp(-gamma_b / 2);

fprintf('Problem 1:\n');
fprintf('  Coherent Orthogonal (P_coh1): %.4f\n', P_coh1);
fprintf('  Coherent Optimal (P_coh2): %.4f\n', P_coh2);
fprintf('  Noncoherent Orthogonal (P_noncoh1): %.4f\n', P_noncoh1);
fprintf('  Noncoherent Optimal (P_noncoh2): [See Written Work]\n');
fprintf('\n');

%% Problem 2d & 2e

% Coefficients
ak = [2.5, 1.5, 1.0, 1.0, 0.5, 0.5];
bk = [1, 2, 4, 5, 9, 10];

% SNR Vector
gamma_dB_vec = 4:0.1:14;
gamma_lin_vec = 10.^(gamma_dB_vec / 10);

% Pre-allocate
Pe_detailed = zeros(size(gamma_lin_vec));
Pe_simple   = zeros(size(gamma_lin_vec));

% Calculations
for k = 1:length(gamma_lin_vec)
    gb = gamma_lin_vec(k);
    Pe_detailed(k) = sum(ak .* q(sqrt(bk .* gb)));
    Pe_simple(k)   = ak(1) * q(sqrt(bk(1) * gb));
end

% Output at 10 dB
idx_10 = find(abs(gamma_dB_vec - 10) < 1e-5);
fprintf('Problem 2:\n');
fprintf('  Detailed Pe: %.4e\n', Pe_detailed(idx_10));
fprintf('  Simple Pe: %.4e\n', Pe_simple(idx_10));
fprintf('\n');

% Plot
figure('Name', 'Problem 2e');
semilogy(gamma_dB_vec, Pe_detailed, 'b-', 'LineWidth', 2); 
hold on;
semilogy(gamma_dB_vec, Pe_simple, 'r--', 'LineWidth', 2);
grid on;
title('Problem 2e: Union Bound vs SNR');
xlabel('SNR per bit \gamma_b (dB)');
ylabel('Probability of Error P_e');
legend('Detailed', 'Simplified');
axis([4 14 1e-7 1]);

%% Problem 6d

% Parameters
Rs = 48e3;
beta = 0.2;
L = 16;
span = 3;

% Transmit Filter
h_tx = rcosdesign(beta, span, L, 'sqrt');

% Matched Filter Output
g_matched = conv(h_tx, h_tx);


peak_val = max(abs(g_matched));
norm_scale = 1 / sqrt(peak_val);
h_tx = h_tx * norm_scale;
g_matched = conv(h_tx, h_tx);

% Plot
figure('Name', 'Problem 6d');
subplot(2,1,1); 
stem(h_tx, 'filled');
title('Transmit Pulse (\surdRC)');
xlabel('Samples'); ylabel('Amplitude'); 
grid on;

subplot(2,1,2); 
stem(g_matched, 'filled');
title('Matched Filter Output (RC)');
xlabel('Samples'); ylabel('Amplitude'); 
grid on;

%% Problem 6e

% Peak center index
mid_idx = ceil(length(g_matched)/2);

idx_all = 1:length(g_matched);
is_symbol = mod(idx_all - mid_idx, L) == 0;

% Interference indices
idx_interf = idx_all(is_symbol & (idx_all ~= mid_idx));

% Worst Case ISI
vals_interf = g_matched(idx_interf);
I_worst = sum(abs(vals_interf));

% SIR Calculation
P_sig = 1^2;
P_int_worst = I_worst^2;

SIR0_lin = P_sig / P_int_worst;
SIR0_dB  = 10 * log10(SIR0_lin);

fprintf('Problem 6:\n');
fprintf('  e: Theoretical SIR_0: %.2f dB\n', SIR0_dB);

%% Problem 6f

% SNIR = SIR0 - 5dB
target_SNIR_dB = SIR0_dB - 5;
target_SNIR_lin = 10^(target_SNIR_dB/10);

% 1/SNIR = 1/SNR + 1/SIR
inv_SNR = (1/target_SNIR_lin) - (1/SIR0_lin);
req_SNR_lin = 1 / inv_SNR;
req_SNR_dB  = 10 * log10(req_SNR_lin);

fprintf('  f: Target SNIR: %.2f dB\n', target_SNIR_dB);
fprintf('  f: Required SNR: %.2f dB\n', req_SNR_dB);
fprintf('\n');

%% Problem 7a

% Generate bits and map to QPSK
total_bits = 1e5;
raw_bits = randi([0 1], 1, total_bits);

N_sym = total_bits / 2;
qpsk_seq = zeros(1, N_sym);

for k = 1:N_sym
    pair = raw_bits(2*k-1 : 2*k);
    if isequal(pair, [0 0]),     s =  1 + 1j;
    elseif isequal(pair, [0 1]), s = -1 + 1j;
    elseif isequal(pair, [1 1]), s = -1 - 1j;
    else,                        s =  1 - 1j;
    end
    qpsk_seq(k) = s;
end

% Upsample and Filter
seq_up = upsample(qpsk_seq, L);
waveform_tx = conv(seq_up, h_tx);

% Plot Envelope
figure('Name', 'Problem 7a');
plot(abs(waveform_tx(1:1000))); 
title('7a Transmitted Envelope');
xlabel('Sample'); ylabel('Magnitude');
grid on;

% Calculate Transmitted SIR
[peak_h, loc_h] = max(abs(h_tx));

% Identify interference samples
idx_all = 1:length(h_tx);
idx_int = idx_all(mod(idx_all - loc_h, L) == 0 & idx_all ~= loc_h);

% Calculate Powers
P_sig_tx = 2 * peak_h^2;
P_isi_tx = 2 * sum(abs(h_tx(idx_int)).^2);

SIR_tx_dB = 10 * log10(P_sig_tx / P_isi_tx);
fprintf('Problem 7:\n');
fprintf('  a: Transmitted Average SIR: %.2f dB\n', SIR_tx_dB);

%% Problem 7b

% Target Noise Power 
Ps_out = 2; 
Pn_out_req = Ps_out / (10^(req_SNR_dB/10));

% Input Noise Variance
E_filt = sum(abs(h_tx).^2);
Pn_in = Pn_out_req / E_filt;
n_scale = sqrt(Pn_in / 2);

% Generate Noise and Add to Signal
noise = n_scale * (randn(size(waveform_tx)) + 1j*randn(size(waveform_tx)));
rx_sig = waveform_tx + noise;

% Matched Filter
mf_out = conv(rx_sig, h_tx);

% Plot Envelope
figure('Name', 'Problem 7b');
plot(abs(mf_out(1:1000))); 
title('7b Matched Filter Envelope (Noisy)'); 
xlabel('Sample'); ylabel('Magnitude');
grid on;

% Output SNIR
c_idx = ceil(length(g_matched)/2);
k_idx = 1:length(g_matched);
isi_mask = mod(k_idx - c_idx, L) == 0 & (k_idx ~= c_idx);
g_isi = g_matched(isi_mask);

% Average ISI Power
Pi_out = 2 * sum(abs(g_isi).^2);
P_total = Pi_out + Pn_out_req;

SNIR_final = 10 * log10(Ps_out / P_total);
fprintf('  b: Output SNIR (Average): %.2f dB\n', SNIR_final);

%% Problem 7c

delay_samples = span * L;
start_idx = delay_samples + 1;

% Downsample Matched Filter Output
sample_indices = start_idx : L : (start_idx + (N_sym-1)*L);
rx_soft = mf_out(sample_indices);

% Decode Bits
rx_bits = zeros(1, total_bits);

for k = 1:length(rx_soft)
    val = rx_soft(k);
    re = real(val); im = imag(val);
    
    if re >= 0
        if im >= 0, pair = [0 0]; % Q1
        else,       pair = [1 0]; % Q4
        end
    else
        if im >= 0, pair = [0 1]; % Q2
        else,       pair = [1 1]; % Q3
        end
    end
    rx_bits(2*k-1 : 2*k) = pair;
end

% Count Errors
bit_errs = sum(raw_bits ~= rx_bits);
BER_val = bit_errs / total_bits;

% Symbol Errors
tx_pairs = reshape(raw_bits, 2, []).';
rx_pairs = reshape(rx_bits, 2, []).';
sym_errs = sum(any(tx_pairs ~= rx_pairs, 2));
SER_val = sym_errs / N_sym;

fprintf('  c: Bit Errors: %d\n', bit_errs);
fprintf('  c: BER: %.2e\n', BER_val);
fprintf('  c: Symbol Errors: %d\n', sym_errs);
fprintf('  c: SER: %.2e\n', SER_val);