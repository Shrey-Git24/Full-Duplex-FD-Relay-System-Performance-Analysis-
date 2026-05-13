clc;
clear;
close all;

%% =========================================================
%%        FULL DUPLEX DF RELAY SYSTEM : BER vs SNR
%%     Simulation + Theoretical BER + Throughput Analysis
%% =========================================================

%% ================= SYSTEM PARAMETERS =================
c = 3e8;                        % Speed of light
f = 2.4e9;                      % Carrier frequency
lambda = c/f;

% Distances
d1 = 10;                        % Source -> Relay
d2 = 6;                         % RSI loop
d3 = 16;                        % Relay -> Destination

path_loss_exp = 2.5;

%% ================= MODULATION PARAMETERS =================
M_values = [4 16 32];           % QPSK, 16-QAM, 32-QAM
nBits = 1e5;

%% ================= NAKAGAMI PARAMETERS =================
m1 = 1;      omega1 = 4;        % S->R
m2 = 0.5;    omega2 = 0.0625;   % RSI
m3 = 2;      omega3 = 4;        % R->D

%% ================= FIXED RSI =================
alpha_dB = -10;                 % Fixed RSI
alpha = 10^(alpha_dB/10);

%% ================= SNR RANGE =================
SNR_dB = 0:2:30;
SNR_linear = 10.^(SNR_dB./10);

%% ================= PATH LOSS =================
PL1 = d1^(-path_loss_exp);
PL2 = d2^(-path_loss_exp);
PL3 = d3^(-path_loss_exp);

%% ================= PREALLOCATIONS =================
ABER_FD = zeros(length(M_values), length(SNR_dB));
ABER_HD = zeros(length(M_values), length(SNR_dB));

BER_Theory = zeros(length(M_values), length(SNR_dB));

Throughput_FD = zeros(length(M_values), length(SNR_dB));
Throughput_HD = zeros(length(M_values), length(SNR_dB));

%% =========================================================
%%                  MAIN SIMULATION
%% =========================================================

for m_idx = 1:length(M_values)

    M = M_values(m_idx);
    k = log2(M);

    % Ensure bits are multiple of k
    nBits_mod = k * floor(nBits/k);

    fprintf('\n====================================\n');
    fprintf('Running Simulation for %d-QAM\n', M);
    fprintf('====================================\n');

    for i = 1:length(SNR_dB)

        %% ============================================
        %%        SNR AND POWER CALCULATIONS
        %% ============================================

        snr_linear = SNR_linear(i);

        No = 1e-3;

        % Signal Power from SNR
        Power_s = snr_linear * No;
        Power_r = Power_s;

        sigma = sqrt(No/2);

        total_error_FD = 0;
        total_error_HD = 0;

        total_bits = 0;

        %% ============================================
        %%              RANDOM DATA
        %% ============================================

        data = randi([0 1], nBits_mod, 1);

        %% ============================================
        %%              QAM MODULATION
        %% ============================================

        tx_s = qammod(data, M, ...
            'InputType', 'bit', ...
            'UnitAveragePower', true);

        %% ============================================
        %%          NAKAGAMI FADING CHANNELS
        %% ============================================

        h12 = sqrt(PL1) .* ...
            sqrt(gamrnd(m1, omega1/m1, size(tx_s)));

        h22 = sqrt(PL2) .* ...
            sqrt(gamrnd(m2, omega2/m2, size(tx_s)));

        h23 = sqrt(PL3) .* ...
            sqrt(gamrnd(m3, omega3/m3, size(tx_s)));

        %% ============================================
        %%              RSI SIGNAL
        %% ============================================

        s2 = 2*randi([0 1], length(tx_s), 1) - 1;

        %% ============================================
        %%                  AWGN
        %% ============================================

        noise_r = sigma * ...
            (randn(size(tx_s)) + 1j*randn(size(tx_s)));

        noise_d = sigma * ...
            (randn(size(tx_s)) + 1j*randn(size(tx_s)));

        %% =================================================
        %%              FULL DUPLEX RELAY
        %% =================================================

        received_r_fd = ...
            sqrt(Power_s).*h12.*tx_s + ...
            sqrt(Power_r * alpha).*h22.*s2 + ...
            noise_r;

        %% Equalization at Relay
        r_fd = received_r_fd ./ (sqrt(Power_s).*h12);

        %% Relay Detection
        rx_r_bits_fd = qamdemod(r_fd, M, ...
            'OutputType', 'bit', ...
            'UnitAveragePower', true);

        %% Re-Modulation
        tx_d_fd = qammod(rx_r_bits_fd, M, ...
            'InputType', 'bit', ...
            'UnitAveragePower', true);

        %% Destination Reception
        received_d_fd = ...
            sqrt(Power_r).*h23.*tx_d_fd + ...
            noise_d;

        %% Equalization
        r2_fd = received_d_fd ./ (sqrt(Power_r).*h23);

        %% Final Detection
        rx_d_bits_fd = qamdemod(r2_fd, M, ...
            'OutputType', 'bit', ...
            'UnitAveragePower', true);

        %% BER Calculation
        error_fd = sum(data ~= rx_d_bits_fd);

        total_error_FD = total_error_FD + error_fd;

        %% =================================================
        %%              HALF DUPLEX RELAY
        %% =================================================

        % No RSI present in HD mode
        received_r_hd = ...
            sqrt(Power_s).*h12.*tx_s + ...
            noise_r;

        %% Equalization
        r_hd = received_r_hd ./ (sqrt(Power_s).*h12);

        %% Relay Detection
        rx_r_bits_hd = qamdemod(r_hd, M, ...
            'OutputType', 'bit', ...
            'UnitAveragePower', true);

        %% Re-Modulation
        tx_d_hd = qammod(rx_r_bits_hd, M, ...
            'InputType', 'bit', ...
            'UnitAveragePower', true);

        %% Destination Reception
        received_d_hd = ...
            sqrt(Power_r).*h23.*tx_d_hd + ...
            noise_d;

        %% Equalization
        r2_hd = received_d_hd ./ (sqrt(Power_r).*h23);

        %% Final Detection
        rx_d_bits_hd = qamdemod(r2_hd, M, ...
            'OutputType', 'bit', ...
            'UnitAveragePower', true);

        %% BER Calculation
        error_hd = sum(data ~= rx_d_bits_hd);

        total_error_HD = total_error_HD + error_hd;

        total_bits = total_bits + nBits_mod;

        %% ============================================
        %%              AVERAGE BER
        %% ============================================

        ABER_FD(m_idx, i) = total_error_FD / total_bits;

        ABER_HD(m_idx, i) = total_error_HD / total_bits;

        %% ============================================
        %%              THROUGHPUT
        %% ============================================

        Throughput_FD(m_idx, i) = ...
            k * (1 - ABER_FD(m_idx, i));

        Throughput_HD(m_idx, i) = ...
            0.5 * k * (1 - ABER_HD(m_idx, i));

        %% ============================================
        %%          THEORETICAL BER (M-QAM)
        %% ============================================

        BER_Theory(m_idx, i) = ...
            (4/k) * (1 - 1/sqrt(M)) * ...
            qfunc(sqrt((3*k*snr_linear)/(M-1)));

    end
end

%% =========================================================
%%                      BER PLOTS
%% =========================================================

figure;

markers = {'o-', 's-', '^-'};

for m_idx = 1:length(M_values)

    %% FD BER
    semilogy(SNR_dB, ABER_FD(m_idx,:), ...
        markers{m_idx}, ...
        'LineWidth', 1.8);

    hold on;

    %% HD BER
    semilogy(SNR_dB, ABER_HD(m_idx,:), ...
        '--', ...
        'LineWidth', 1.8);

    %% THEORETICAL BER
    semilogy(SNR_dB, BER_Theory(m_idx,:), ...
        ':', ...
        'LineWidth', 2.2);

end

grid on;

xlabel('SNR (dB)', ...
    'FontWeight', 'bold');

ylabel('Bit Error Rate (BER)', ...
    'FontWeight', 'bold');

title('BER vs SNR for FD and HD Relay Systems under Nakagami-m Fading', ...
    'FontWeight', 'bold');

legend( ...
    'FD-QPSK', 'HD-QPSK', 'Theory-QPSK', ...
    'FD-16QAM', 'HD-16QAM', 'Theory-16QAM', ...
    'FD-32QAM', 'HD-32QAM', 'Theory-32QAM', ...
    'Location', 'southwest');

set(gca, 'FontSize', 12);

%% =========================================================
%%                  THROUGHPUT PLOTS
%% =========================================================

figure;

for m_idx = 1:length(M_values)

    %% FD Throughput
    plot(SNR_dB, Throughput_FD(m_idx,:), ...
        markers{m_idx}, ...
        'LineWidth', 1.8);

    hold on;

    %% HD Throughput
    plot(SNR_dB, Throughput_HD(m_idx,:), ...
        '--', ...
        'LineWidth', 1.8);

end

grid on;

xlabel('SNR (dB)', ...
    'FontWeight', 'bold');

ylabel('Normalized Throughput', ...
    'FontWeight', 'bold');

title('Throughput vs SNR for FD and HD Relay Systems', ...
    'FontWeight', 'bold');

legend( ...
    'FD-QPSK', 'HD-QPSK', ...
    'FD-16QAM', 'HD-16QAM', ...
    'FD-32QAM', 'HD-32QAM', ...
    'Location', 'southeast');

set(gca, 'FontSize', 12);

%% =========================================================
%%                  DISPLAY SUMMARY
%% =========================================================

fprintf('\n============================================\n');
fprintf('Simulation Completed Successfully\n');
fprintf('============================================\n');

fprintf('FEATURES INCLUDED:\n');
fprintf('1. Full Duplex DF Relay System\n');
fprintf('2. Half Duplex Relay Comparison\n');
fprintf('3. Nakagami-m Fading Channels\n');
fprintf('4. Path Loss Modeling\n');
fprintf('5. Multiple QAM Schemes\n');
fprintf('6. Monte Carlo BER Simulation\n');
fprintf('7. Theoretical BER Analysis\n');
fprintf('8. Throughput Analysis\n');
fprintf('9. Residual Self-Interference Modeling\n');

fprintf('============================================\n');