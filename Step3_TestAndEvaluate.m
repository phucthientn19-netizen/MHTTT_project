%% BƯỚC 3: KIỂM THỬ VÀ ĐÁNH GIÁ (CHỌN CỔNG TỐT NHẤT)
clear; clc;

load('FAS_Trained_Net.mat');
num_test = size(X_test, 2);

fprintf('Đang dự đoán tập test (%d mẫu)...\n', num_test);
% 1. Chạy model ML để dự đoán toàn bộ kênh từ vài điểm
Y_pred = net(X_test);

% 2. Thuật toán chọn cổng
% - Cổng tốt nhất thực tế (Ground Truth - Quét cạn toàn bộ N cổng)
[Max_Gain_True, Best_Port_True] = max(Y_test, [], 1);

% - Cổng tốt nhất do ML dự đoán
[~, Best_Port_ML] = max(Y_pred, [], 1);

% Trích xuất giá trị kênh thực tế tại cổng mà ML đã chọn
% (Lưu ý: Ta phải tính trên kênh THỰC TẾ Y_test, chứ ko phải giá trị ảo Y_pred)
Max_Gain_ML_Choice = zeros(1, num_test);
for i = 1:num_test
    idx_chosen_by_ML = Best_Port_ML(i);
    Max_Gain_ML_Choice(i) = Y_test(idx_chosen_by_ML, i);
end

% 3. Đánh giá độ chính xác (Accuracy của việc chọn đúng vị trí)
% Đo lường xem bao nhiêu % ML chọn trúng đúng y chang cổng tốt nhất
Acc = sum(Best_Port_ML == Best_Port_True) / num_test * 100;
fprintf('Độ chính xác chọn ĐÚNG vị trí cổng tốt nhất tuyệt đối: %.2f%%\n', Acc);

% 4. Đánh giá hao hụt hiệu năng kênh (Performance Loss)
% Dù không chọn đúng vị trí, nhưng nếu cổng ML chọn có giá trị kênh gần bằng 
% cổng tốt nhất thì vẫn là xuất sắc.
Ratio_Gain = mean(Max_Gain_ML_Choice ./ Max_Gain_True) * 100;
fprintf('Tỉ lệ cường độ kênh ML giữ lại được so với Quét cạn: %.2f%%\n', Ratio_Gain);

% 5. So sánh BLER (Giả lập đơn giản)
% Giả sử công thức Outage Probability / BLER xấp xỉ tỉ lệ nghịch với SNR
snr_threshold = 1.0; % Ngưỡng mô phỏng
OP_True = sum(Max_Gain_True < snr_threshold) / num_test;
OP_ML = sum(Max_Gain_ML_Choice < snr_threshold) / num_test;

fprintf('\n=== SO SÁNH OUTAGE PROBABILITY (OP) ===\n');
fprintf('OP khi dùng Quét cạn (Optimal)  : %.5f\n', OP_True);
fprintf('OP khi dùng ML dự đoán (5 cổng) : %.5f\n', OP_ML);

% Vẽ đồ thị so sánh một vài mẫu
figure;
plot(1:50, Max_Gain_True(1:50), 'b-o', 'LineWidth', 1.5); hold on;
plot(1:50, Max_Gain_ML_Choice(1:50), 'r-*', 'LineWidth', 1.5);
legend('Quét cạn (Toàn bộ N cổng)', 'Mạng Nơ-ron (Chỉ đo 5 cổng)');
xlabel('Chỉ số gói tin (Packet/Sample Index)');
ylabel('Cường độ kênh thu được (|g|^2)');
title('So sánh Cường độ kênh thu được: Optimal vs. Machine Learning');
grid on;