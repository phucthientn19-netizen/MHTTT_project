clear; clc; close all;

%% ---------------- System parameters ----------------------------
M        = 400;      % finite blocklength
N        = 10;       % number of FAS ports
W        = 0.5;        % antenna array normalized length
sigma    = sqrt(2);  % per-antenna signal std
mu2      = 0.97;     % correlation coefficient squared
rho_thres = 1;       % eigenvalue threshold for dominant set
Bit_num   = 100;
% SNR range (per-antenna SNR in dB): snr = sigma^2 / sigma_eta^2
snr_dB_vec = -30:2:30;
nSNR       = numel(snr_dB_vec);

%% ----------- Toeplitz correlation matrix Sigma ----------------
Sigma = zeros(N,N);
for row = 1:N
    x = ((1:N) - row) * (2*pi*W/(N-1));
    Sigma(row,:) = sinc(x);  % normalized sinc
end

% Eigen-decomposition
e_v   = eig(Sigma);
S_dom = e_v(e_v > rho_thres);   % dominant eigenvalues > threshold
Num_eig = length(S_dom);

% Block sizes via BlockCorrelation() function
L_vec = BlockCorrelation(N, S_dom, Num_eig, mu2);
B     = length(L_vec);          % number of blocks

fprintf('Number of blocks B = %d\n', B);
fprintf('Block sizes L_vec = ['); fprintf('%d ', L_vec); fprintf(']\n');

%% ================= Monte-Carlo simulation =================
MC = 2e5;   % số mẫu mô phỏng

% Cholesky decomposition để tạo kênh tương quan
L_chol = chol(Sigma,'lower');

gFAS_power = zeros(MC,1);

for i = 1:MC
    
    % Rayleigh channel vector
    g = sigma/sqrt(2) * (randn(N,1) + 1j*randn(N,1));
    
    % apply spatial correlation
    g_corr = L_chol * g;
    
    % FAS chọn antenna tốt nhất
    gFAS_power(i) = max(abs(g_corr).^2);

end

%% Histogram Monte-Carlo
[counts,edges] = histcounts(gFAS_power,60,'Normalization','pdf');

t_hist = (edges(1:end-1)+edges(2:end))/2;

%% ---------------- Gauss-Laguerre nodes / weights --------------
N_GL  = 128;                     % GL order
n     = N_GL;
alpha = 2*(0:n-1) + 1;          % alpha_n = 2n+1
beta  = (1:n-1).^2;             % beta_n = n^2 (n>=1)

J = diag(alpha) + diag(sqrt(beta),1) + diag(sqrt(beta),-1);
[W_GL, D_GL] = eig(J);
x_GL = diag(D_GL);
w_GL = W_GL(1,:).^2;

% sort
[x_GL, idx] = sort(x_GL);
w_GL = w_GL(idx)';

%% ---------------- t = |g_FAS|^2 grid and exact PDF ------------
t_vec = logspace(-4, 3, 800);   % channel power grid
Nt    = numel(t_vec);

% ===== exact PDF via Gauss-Laguerre =====
fas_pdf = zeros(size(t_vec));
for k = 1:Nt
    fas_pdf(k) = fas_pdf_true(t_vec(k), mu2, B, L_vec, x_GL, w_GL);
end
fas_pdf = fas_pdf ./ trapz(t_vec, fas_pdf);   % normalize

% non-uniform grid spacing for integration
dt_vec = [diff(t_vec), t_vec(end)-t_vec(end-1)];

%% ---------------- Gamma Approximation of FAS power -------------
% Compute moments from exact PDF
m1 = sum(t_vec .* fas_pdf .* dt_vec);      % E[t]
m2 = sum(t_vec.^2 .* fas_pdf .* dt_vec);   % E[t^2]
var_t = m2 - m1^2;

% Gamma parameters (shape k, scale theta)
k_gamma     = m1^2 / var_t;
theta_gamma = var_t / m1;

% Gamma approximate PDF on the same grid
fas_pdf_gamma = gampdf(t_vec, k_gamma, theta_gamma);
fas_pdf_gamma = fas_pdf_gamma ./ trapz(t_vec, fas_pdf_gamma);   % renormalize

% Exact CDF & Gamma CDF
fas_cdf_exact = cumtrapz(t_vec, fas_pdf);
fas_cdf_exact = fas_cdf_exact ./ fas_cdf_exact(end);

fas_cdf_gamma = gamcdf(t_vec, k_gamma, theta_gamma);

fas_pdf_GL = fas_pdf;     % GL-PDF
fas_pdf_exact = fas_pdf;  % tạm thời giống nhau

%% ================= Figure 1 =================
figure;
set(gcf,'Position',[200 200 800 500]);

yyaxis left
plot(t_vec, fas_cdf_exact,'b-','LineWidth',2); hold on
ylabel('CDF')
ylim([0 1])

yyaxis right
plot(t_vec, fas_pdf,'r--','LineWidth',2)
ylabel('PDF')
ylim([0 0.2])

xlabel('|g_{FAS}|^2')
xlim([0 25])
grid on
text(8,0.85,'Exact CDF (blue solid)','Color','b','FontSize',11)
yyaxis right
text(6,0.05,'Derived PDF (red dashed)','Color','r','FontSize',11)

%% Inset subplot
axes('Position',[0.55 0.25 0.35 0.35])

histogram(gFAS_power,60,'Normalization','count')
hold on

% Monte-Carlo CDF
[g_sorted,idx] = sort(gFAS_power);
cdf_mc = (1:length(g_sorted))/length(g_sorted);

yyaxis right
plot(g_sorted,cdf_mc,'k--','LineWidth',2)

xlim([0 20])
grid on

% ===== chú thích inset =====
text(7,1500,'Histogram Monte-Carlo samples','FontSize',10)
yyaxis right
text(8,0.7,'Monte-Carlo CDF (black dashed)','FontSize',10)
%% ===============================================================
%   Figure 1: PDF & CDF comparison (dual y-axis)
% ===============================================================
figure;
set(gcf,'Position',[200 200 800 500]);

yyaxis left
plot(t_vec, fas_pdf_GL,'k--','LineWidth',2); hold on
plot(t_vec, fas_pdf_gamma,'b-','LineWidth',2)
plot(t_vec, fas_pdf_exact,'b:','LineWidth',2)

ylabel('PDF')
ylim([0 0.2])

yyaxis right
plot(t_vec, fas_cdf_exact,'r--','LineWidth',2)
plot(t_vec, fas_cdf_gamma,'m-','LineWidth',2)

ylabel('CDF')
ylim([0 1])

xlabel('|g_{FAS}|^2')
xlim([0 20])

grid on

legend({'GL-PDF (N_{GL}=128)', ...
        'Gamma-Approx. PDF', ...
        'Exact-PDF', ...
        'Exact-CDF', ...
        'Gamma-Approx. CDF'}, ...
        'Location','northeast')
%% ===============================================================
%   BLER functions: exact & high-SNR asymptotic
% ===============================================================
log2e = 1/log(2);
Rc    = Bit_num/M;           % coding rate

% Exact finite-blocklength BLER under normal approximation
bler_exact = @(gamma) bler_exact_fun(gamma, M, Rc, log2e);

% High-SNR BLER asymptotic:
bler_asym = @(gamma) bler_asym_fun(gamma, M, Rc, log2e);

%% ===============================================================
%   Figure 2: Instantaneous BLER vs SNR (t fixed at E[t])
% ===============================================================
snr_lin_vec = 10.^(snr_dB_vec/10);
gamma_inst_exact = zeros(size(snr_lin_vec));
bler_inst_exact  = zeros(size(snr_lin_vec));
bler_inst_asym   = zeros(size(snr_lin_vec));

t_eff = 1;   % representative channel power for "instantaneous" plot

for iSNR = 1:nSNR
    snr_lin     = snr_lin_vec(iSNR);
    sigma_eta2  = sigma^2 / snr_lin;   % noise variance

    den = sigma_eta2 / t_eff;
    gamma_val = 1 ./ den;

    gamma_inst_exact(iSNR) = gamma_val;

    bler_inst_exact(iSNR) = bler_exact(gamma_val);
    bler_inst_asym(iSNR)  = bler_asym(gamma_val);
end

figure;
semilogy(snr_dB_vec, bler_inst_exact, 'ko-','LineWidth',1.8); hold on;
semilogy(snr_dB_vec, bler_inst_asym,  'r--','LineWidth',1.8);
grid on;
xlabel('SNR (dB)');
ylabel('Instantaneous BLER');
% title(sprintf('Instantaneous BLER vs SNR (t fixed at E[t]=%.3f)', t_eff));
legend('Exact BLER expression','High-SNR BLER asymptotic', ...
    'Location','southwest');

%% ===============================================================
%   Figure 3: Average BLER vs SNR (4 curves)
% ===============================================================
bler_ExactPDF_ExactBLER = zeros(size(snr_lin_vec));
bler_ExactPDF_AsymBLER  = zeros(size(snr_lin_vec));
bler_GammaPDF_ExactBLER = zeros(size(snr_lin_vec));
bler_GammaPDF_AsymBLER  = zeros(size(snr_lin_vec));

for iSNR = 1:nSNR
    snr_lin    = snr_lin_vec(iSNR);
    sigma_eta2 = sigma^2 / snr_lin;

    % SNR for each t on grid
    den_vec   =  sigma_eta2 ./ t_vec; % SINR denominator
    gamma_vec = 1 ./ den_vec;  % SINR

    % Instantaneous BLER (exact & asymptotic)
    eps_exact_vec = bler_exact(gamma_vec);
    eps_asym_vec  = bler_asym(gamma_vec);

    % 1) Exact PDF + Exact BLER
    bler_ExactPDF_ExactBLER(iSNR) = sum(eps_exact_vec .* fas_pdf .* dt_vec);

    % 2) Exact PDF + Asymptotic BLER
    bler_ExactPDF_AsymBLER(iSNR)  = sum(eps_asym_vec  .* fas_pdf .* dt_vec);

    % 3) Gamma PDF + Exact BLER
    bler_GammaPDF_ExactBLER(iSNR) = sum(eps_exact_vec .* fas_pdf_gamma .* dt_vec);

    % 4) Gamma PDF + Asymptotic BLER
    bler_GammaPDF_AsymBLER(iSNR)  = sum(eps_asym_vec  .* fas_pdf_gamma .* dt_vec);
end

figure;
semilogy(snr_dB_vec, bler_ExactPDF_ExactBLER, 'ko-','LineWidth',1.8); hold on;
semilogy(snr_dB_vec, bler_ExactPDF_AsymBLER,  'rs--','LineWidth',1.8);
semilogy(snr_dB_vec, bler_GammaPDF_ExactBLER, 'b^-.','LineWidth',1.8);
semilogy(snr_dB_vec, bler_GammaPDF_AsymBLER,  'm*:','LineWidth',1.8);

grid on;
xlabel('SNR (dB)');
ylabel('Average BLER');
% title('Average BLER over FAS: Exact vs Gamma Approx & BLER Asymptotic');
legend( ...
    'Exact PDF + Exact BLER', ...
    'Exact PDF + Asymptotic BLER', ...
    'Gamma PDF + Exact BLER', ...
    'Gamma PDF + Asymptotic BLER', ...
    'Location','southwest');
%% ===============================================================
%   Local functions
% ===============================================================

function [L] = BlockCorrelation(N, rho, Num_eig, mu2)
    % Dominant set of eigenvalues
    lambda_eff = rho(1:Num_eig);
    
    % Algorithm 1 from your code
    L = zeros(Num_eig,1);
    counter = 0;
    index = ones(length(L),1);
    
    while counter < N && sum(index) > 0
        L(index==1) = L(index==1) + ones(sum(index),1);
        counter = sum(L);
    
        index(index==1) = abs((L(index==1)-1)*mu2 + 1 - lambda_eff(index==1)) >= ...
                abs((L(index==1))*mu2 + 1 - lambda_eff(index==1));
    end
end

function f = fas_pdf_true(t, mu2, B, L_vec, x_GL, w_GL)
    % FAS channel power PDF via Gauss-Laguerre approximation
    xnorm  = t/(1-mu2);                  % x = t / (1 - mu^2)
    lambda = 2*mu2*x_GL/(1-mu2);         % non-centrality parameters

    Phi = ncx2cdf(xnorm, 2, lambda);     % CDF (df=2)
    phi = ncx2pdf(xnorm, 2, lambda);     % PDF (df=2)

    S = zeros(B,1);
    % S_j(t) = sum_i w_i * Phi_i^{L_j}
    for j = 1:B
        Lj = L_vec(j);
        S(j) = sum( w_GL .* (Phi.^Lj) );
    end

    f = 0;
    for b = 1:B
        Lb  = L_vec(b);
        T_b = sum( w_GL .* (Lb/(1-mu2)) .* (Phi.^(Lb-1)) .* phi );
        prod_other = 1;
        for j = 1:B
            if j == b, continue; end
            prod_other = prod_other * S(j);
        end
        f = f + T_b * prod_other;
    end
end

function epsL = bler_exact_fun(gamma, M, Rc, log2e)
    C = 0.5 * log2(1 + gamma);                         % capacity
    V = 0.5 .* gamma .* (gamma + 2) ./ (gamma + 1).^2 ...
        .* (log2e.^2);                                 % dispersion
    Xi = (C - Rc) ./ sqrt(V ./ M);
    epsL = qfunc(Xi);
end

function epsL_as = bler_asym_fun(gamma, M, Rc, log2e)
    V_inf = 0.5 * (log2e.^2);
    x_as  = sqrt(M./V_inf) .* (0.5*log2(1+gamma) - Rc);

    epsL_as = zeros(size(gamma));
    idx_pos = x_as > 0;
    epsL_as(idx_pos) = 1 ./ (sqrt(2*pi) .* x_as(idx_pos)) .* ...
                       exp(-0.5 .* x_as(idx_pos).^2);
    epsL_as(~idx_pos) = 0.5;
end
