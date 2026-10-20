%% BƯỚC 1: SINH DỮ LIỆU HUẤN LUYỆN CHO FAS
clear; clc;

% 1. Thông số hệ thống
N = 30;                 % Tổng số cổng (ports)
W = 1;                  % Chiều dài chuẩn hóa của ăng-ten
num_samples = 200000;   % Số lượng mẫu dữ liệu (100,000 để train cho tốt)
K_idx = [1, 8, 15, 23, 30]; % Vị trí 5 cổng mà ta sẽ "đo đạc" (Partial observation)
K = length(K_idx);

fprintf('Đang tạo ma trận tương quan không gian...\n');
% 2. Tạo ma trận tương quan không gian (Toeplitz Correlation Matrix)
Sigma = zeros(N,N);
for row = 1:N
    x = ((1:N) - row) * (2*pi*W/(N-1));
    Sigma(row,:) = sinc(x/pi); % Hàm sinc(x) trong Matlab là sin(pi*x)/(pi*x)
end

% Lấy căn bậc 2 của ma trận tương quan để tạo nhiễu có tương quan
R_sqrt = sqrtm(Sigma); 

% 3. Sinh kênh truyền Rayleigh fading có tương quan
fprintf('Đang sinh %d mẫu kênh truyền...\n', num_samples);
% H_indep là kênh Rayleigh độc lập (Gaussian phức)
H_indep = (randn(N, num_samples) + 1i * randn(N, num_samples)) / sqrt(2);

% H_corr là kênh đã có tương quan không gian
H_corr = R_sqrt * H_indep;

% Lấy bình phương biên độ (Channel Power Gain |g|^2)
Channel_Power = abs(H_corr).^2;

% 4. Chuẩn bị dữ liệu cho Machine Learning
% Đầu vào X: Biên độ tại 5 cổng đo được
X_data = Channel_Power(K_idx, :); % Kích thước: K x num_samples

% Đầu ra Y: Biên độ tại toàn bộ N cổng (để model học hình dạng kênh)
Y_data = Channel_Power;           % Kích thước: N x num_samples

% Lưu dữ liệu lại để dùng cho bước sau
save('FAS_ML_Data.mat', 'X_data', 'Y_data', 'N', 'K_idx', 'num_samples');
fprintf('Đã sinh xong dữ liệu và lưu vào FAS_ML_Data.mat!\n');