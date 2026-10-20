%  Average BLER vs SNR for L-antenna Rayleigh (MRC)
clear; clc; close all;
%% ---------------- System parameters ----------------------------
U        = 1;          % number of active users 
L        = 1;          % conventional antenna num 
M        = 400;       % finite blocklength
sigma    = sqrt(2);    % per-antenna channel std  => E[|g_l|^2] = sigma^2
Bit_num  = 100;        % number of info bits
R_c      = Bit_num / M;   % coding rate

log2e    = 1/log(2);

%% ---------------- SNR range ------------------------------------
snr_dB_vec  = -30:2:30;              % SNR in dB
nSNR        = numel(snr_dB_vec);

epsilon_inst_vec = zeros(1, nSNR);   % instantaneous (no fading)
epsilon_avg_vec  = zeros(1, nSNR);   % average over Rayleigh fading (L antennas)

%% ===============================================================
%   LOOP over SNR
% ===============================================================
for iSNR = 1:nSNR
    
    snr_dB  = snr_dB_vec(iSNR);
    snr_lin = 10^(snr_dB/10);              % linear SNR

    % noise variance: average SNR = E[|g|^2] / sigma_eta2 = sigma^2 / sigma_eta2
    sigma_eta2 = sigma^2 / snr_lin;

    %% ---- instaneous BLER----------------------
    % After MRA, average SNR: E[gamma] = L * snr_lin
    gamma_hat = L * snr_lin;

    C_inst     = 0.5 * log2(1 + gamma_hat);
    V_dis_inst = 0.5 * gamma_hat .* (gamma_hat + 2) ./ (1 + gamma_hat).^2 .* (log2e.^2);

    epsilon_inst_vec(iSNR) = qfunc( (C_inst - R_c) ./ sqrt(V_dis_inst / M) );

    %% ---- average BLER：E[epsilon(t_L)] ----------------------------

    integrand = @(t) avg_bler_integrand(t, sigma, sigma_eta2, M, R_c, L);

    % ∫_0^∞ epsilon(t) f_{t_L}(t) dt
    epsilon_avg_vec(iSNR) = integral(integrand, 0, Inf, ...
                                     'ArrayValued', true); 
end

%% ===============================================================
%   Plot: BLER vs SNR
% ===============================================================
figure;
% semilogy(snr_dB_vec, epsilon_inst_vec, 'ko-','LineWidth',1.8); hold on;
semilogy(snr_dB_vec, epsilon_avg_vec,  'rs--','LineWidth',1.8);

grid on;
xlabel('SNR (dB)');
ylabel('BLER');
% title(sprintf('Instantaneous vs Average BLER over Rayleigh (L = %d)', L));
% legend('Instantaneous BLER (no fading)', ...
%        'Average BLER over Rayleigh (Gamma, MRC)', ...
%        'Location','southwest');

%% ===============================================================
%   Local function for integrand: epsilon(t_L) * f_{t_L}(t_L)
% ===============================================================
function y = avg_bler_integrand(t, sigma, sigma_eta2, M, R_c, L)

    log2e = 1/log(2);

    % t 在这里就是 t_L (合并后的总信道功率)
    % Instantaneous SNR gamma(t_L) = t_L / sigma_eta^2
    gamma_t = t ./ sigma_eta2;

    % Finite-blocklength normal approximation
    C_t = 0.5 .* log2(1 + gamma_t);                % capacity
    V_t = 0.5 .* gamma_t .* (gamma_t + 2) ./ ...
          (1 + gamma_t).^2 .* (log2e.^2);          % dispersion

    Xi_t = (C_t - R_c) ./ sqrt(V_t ./ M);          % Q-argument
    eps_t = qfunc(Xi_t);                           % epsilon(t_L)

    % PDF of t_L ~ Gamma(L, sigma^2)
    % f_{t_L}(t) = t^{L-1} e^{-t/sigma^2} / (Gamma(L) * sigma^{2L}), t >= 0
    coeff = 1 ./ (gamma(L) .* sigma.^(2*L));
    f_t   = coeff .* t.^(L-1) .* exp(-t ./ sigma.^2);

    y = eps_t .* f_t;                              % integrand
end
