% stiefel_extrapolation_final_compact.m
% Stiefel流形外插实验 (紧凑PDF输出版)
% 功能：1. 执行相位对齐预测; 2. 打印详细数据表; 3. 生成两张紧凑PDF图片（无空白边距）

clc; clear; close all;

%% 1. 加载数据与重建变量
data_file = './LoS_5users_data.mat';
if ~isfile(data_file)
    error('未找到数据文件 %s，请先运行 run_quadriga_submit.m 生成数据。', data_file);
end
fprintf('加载数据: %s...\n', data_file);
load(data_file, 'H_seq', 'pos_seq', 'meta');

% 重建 F_dig (如果不存在)
if ~exist('F_dig', 'var')
    fprintf('正在从信道重建预编码矩阵 F_dig (SVD)...\n');
    [T_total, ~] = size(H_seq);
    U_total = length(H_seq{1});
    F_dig = cell(T_total, 1);
    for t = 1:T_total
        F_dig{t} = cell(U_total, 1);
        for u = 1:U_total
            [~,~,V] = svd(H_seq{t}{u}, 'econ');
            F_dig{t}{u} = V(:,1); 
        end
    end
end

% 实验参数
[T_total, U_count] = size(pos_seq(:, :, 1)); 
sigma2 = 1;
SNR_dB = 20; 
P_user = 10^(SNR_dB/10);

T_history = 10; 
T_future  = 5;

fprintf('配置: History=%d, Forecast=%d. 开始计算...\n', T_history, T_future);

%% 2. 执行外插实验 (带相位校准 Phase Alignment)
metrics = struct();
metrics.err_hold = zeros(T_future, U_count);
metrics.err_geo  = zeros(T_future, U_count);
metrics.rate_gt   = zeros(T_future, U_count);
metrics.rate_hold = zeros(T_future, U_count);
metrics.rate_geo  = zeros(T_future, U_count);

for u = 1:U_count
    Nt = size(F_dig{1}{u}, 1);
    F_u_gt = zeros(Nt, 1, T_total);
    for t = 1:T_total
        F_u_gt(:, 1, t) = F_dig{t}{u};
    end
    
    P_n_minus_1 = F_u_gt(:, 1, T_history-1);
    P_n         = F_u_gt(:, 1, T_history);
    
    % === 相位对齐 ===
    phase_diff = angle(P_n' * P_n_minus_1);
    P_n_minus_1_aligned = P_n_minus_1 * exp(-1j * phase_diff);
    
    % 计算流形速度
    V_back = stiefel_log_approx(P_n, P_n_minus_1_aligned);
    V_end = -1.0 * V_back; 
    
    P_pred_hold = repmat(P_n, 1, T_future);
    P_pred_geo = zeros(Nt, T_future);
    
    for k = 1:T_future
        delta_t = k;
        % 指数映射外插
        P_hat = stiefel_exp_approx(P_n, delta_t * V_end);
        P_pred_geo(:, k) = P_hat;
        
        % 获取真值
        P_true = F_u_gt(:, 1, T_history + k);
        H_true = H_seq{T_history + k}{u};
        
        % 记录指标
        metrics.err_hold(k, u) = stiefel_chordal(P_true, P_pred_hold(:, k));
        metrics.err_geo(k, u)  = stiefel_chordal(P_true, P_pred_geo(:, k));
        
        metrics.rate_gt(k, u)   = compute_rate(H_true, P_true, sigma2, P_user);
        metrics.rate_hold(k, u) = compute_rate(H_true, P_pred_hold(:, k), sigma2, P_user);
        metrics.rate_geo(k, u)  = compute_rate(H_true, P_pred_geo(:, k), sigma2, P_user);
    end
end

% 统计平均
avg_err_hold = mean(metrics.err_hold, 2);
avg_err_geo  = mean(metrics.err_geo, 2);
avg_rate_gt   = mean(metrics.rate_gt, 2);
avg_rate_hold = mean(metrics.rate_hold, 2);
avg_rate_geo  = mean(metrics.rate_geo, 2);

%% 3. 打印数据表 (Print Table)
fprintf('\n=== 优化后的外插结果 (Phase Aligned) ===\n');
fprintf('Step\t Hold误差\t Geo误差\t Rate提升\n');
for k = 1:T_future
    gain = (avg_rate_geo(k) - avg_rate_hold(k)) / avg_rate_hold(k) * 100;
    fprintf('%d\t %.4f\t %.4f\t %.2f%%\n', k, avg_err_hold(k), avg_err_geo(k), gain);
end
fprintf('======================================\n\n');

%% 4. 绘制并保存为紧凑 PDF (Compact PDF Output)

% --- 图1: 预测误差 (Chordal Error) ---
fig1 = figure('Color', 'w', 'Units', 'inches');
% 设置纸张大小为 5.5 x 4 英寸（紧凑尺寸）
fig1.PaperUnits = 'inches';
fig1.PaperSize = [5.5 4];
fig1.PaperPosition = [0 0 5.5 4];

ax1 = axes(fig1);
plot(ax1, 1:T_future, avg_err_hold, 'k--^', 'LineWidth', 1.5, 'MarkerSize', 8, 'DisplayName', 'Baseline (Hold)');
hold(ax1, 'on');
plot(ax1, 1:T_future, avg_err_geo, 'r-o', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Proposed (Geo)');
grid(ax1, 'on'); 
box(ax1, 'on');
ax1.FontSize = 11;
xlabel(ax1, 'Prediction Step $k$', 'Interpreter', 'latex', 'FontSize', 12);
ylabel(ax1, 'Chordal Distance Error', 'FontSize', 12);
title(ax1, 'Prediction Accuracy', 'FontSize', 13, 'FontWeight', 'bold');
legend(ax1, 'Location', 'NorthWest', 'FontSize', 10);
if max(avg_err_hold) > 1e-6
    ylim(ax1, [0, max(avg_err_hold)*1.15]);
end

% 保存为 PDF（无空白边距）
print(fig1, 'stiefel_result_error.pdf', '-dpdf', '-r300');
fprintf('已保存紧凑误差图 (PDF): stiefel_result_error.pdf\n');
close(fig1);

% --- 图2: 频谱效率 (Spectral Efficiency) ---
fig2 = figure('Color', 'w', 'Units', 'inches');
% 设置纸张大小为 5.5 x 4 英寸（紧凑尺寸）
fig2.PaperUnits = 'inches';
fig2.PaperSize = [5.5 4];
fig2.PaperPosition = [0 0 5.5 4];

ax2 = axes(fig2);
plot(ax2, 1:T_future, avg_rate_gt, 'g-s', 'LineWidth', 1.5, 'MarkerSize', 8, 'Color', [0 0.6 0], 'DisplayName', 'Ideal CSI');
hold(ax2, 'on');
plot(ax2, 1:T_future, avg_rate_geo, 'r-o', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Proposed');
plot(ax2, 1:T_future, avg_rate_hold, 'k--^', 'LineWidth', 1.5, 'MarkerSize', 8, 'DisplayName', 'Baseline');
grid(ax2, 'on'); 
box(ax2, 'on');
ax2.FontSize = 11;
xlabel(ax2, 'Prediction Step $k$', 'Interpreter', 'latex', 'FontSize', 12);
ylabel(ax2, 'Spectral Efficiency (bps/Hz)', 'FontSize', 12);
title(ax2, 'Rate Performance', 'FontSize', 13, 'FontWeight', 'bold');
legend(ax2, 'Location', 'SouthWest', 'FontSize', 10);

% 保存为 PDF（无空白边距）
print(fig2, 'stiefel_result_rate.pdf', '-dpdf', '-r300');
fprintf('已保存紧凑速率图 (PDF): stiefel_result_rate.pdf\n');
close(fig2);

fprintf('\n✓ 所有文件已生成完毕！\n');

%% 工具函数
function dist = stiefel_chordal(W1, W2)
    overlap = abs(W1' * W2);
    dist = sqrt(max(0, 1 - overlap^2));
end

function rate = compute_rate(H, w, sigma2, P)
    signal_power = P * norm(H * w)^2; 
    rate = log2(1 + signal_power / sigma2);
end

function U = stiefel_log_approx(X, Y)
    inner_prod = real(X' * Y);
    theta = acos(max(-1, min(1, inner_prod)));
    if abs(theta) < 1e-10, U = zeros(size(X)); else, U = (theta / sin(theta)) * (Y - inner_prod * X); end
end

function Y = stiefel_exp_approx(X, U)
    norm_U = norm(U);
    if norm_U < 1e-10, Y = X; else, Y = X * cos(norm_U) + (U / norm_U) * sin(norm_U); end
    Y = Y / norm(Y);
end