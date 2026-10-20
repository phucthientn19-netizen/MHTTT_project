%% BƯỚC 2 (CẬP NHẬT): HUẤN LUYỆN MẠNG NƠ-RON DỰ ĐOÁN KÊNH
clear; clc;

% 1. Tải dữ liệu
fprintf('Đang tải dữ liệu...\n');
load('FAS_ML_Data.mat');

% Chia dữ liệu: 80% Train, 20% Test
train_ratio = 0.8;
num_train = floor(num_samples * train_ratio);

X_train = X_data(:, 1:num_train);
Y_train = Y_data(:, 1:num_train);

X_test = X_data(:, num_train+1:end);
Y_test = Y_data(:, num_train+1:end);

% 2. Khởi tạo mạng Neural Network (Regression)
% ĐIỀU CHỈNH 1: Giảm nhẹ số nơ-ron xuống [64, 32] để train nhanh hơn
hiddenLayerSize = [256, 128]; 

% ĐIỀU CHỈNH 2: Đổi sang thuật toán 'trainscg' (phù hợp cho Dữ liệu lớn)
net = fitnet(hiddenLayerSize, 'trainscg'); 

% Thiết lập các thông số huấn luyện
net.trainParam.epochs = 500;      % Tăng số epoch vì SCG cần nhiều bước lặp hơn LM
net.trainParam.max_fail = 20;     % Dừng sớm nếu validation không cải thiện (Early Stopping)
net.trainParam.showWindow = true; % Hiện cửa sổ tiến trình train

% 3. Tiến hành huấn luyện
fprintf('Bắt đầu huấn luyện mô hình bằng thuật toán SCG...\n');
tic; % Bắt đầu đếm thời gian
[net, tr] = train(net, X_train, Y_train);
thoigian = toc; % Kết thúc đếm thời gian

fprintf('Thời gian huấn luyện: %.2f giây.\n', thoigian);

% Lưu mô hình đã huấn luyện
save('FAS_Trained_Net.mat', 'net', 'X_test', 'Y_test', 'N');
fprintf('Huấn luyện hoàn tất và đã lưu model vào FAS_Trained_Net.mat!\n');