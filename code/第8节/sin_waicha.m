% stiefel_extrapolation_sinusoidal_single_user.m
% 测试曲线运动（正弦轨迹）下的外插性能 - 单用户版本
% 功能：生成Sin轨迹单个用户，执行相位对齐预测

clc; clear; close all;

%% 1. 生成正弦轨迹用户场景
fprintf('=== 生成正弦曲线轨迹场景（单用户）===\n');

% 实验参数
fc     = 3.5e9;         % 载频
lambda = 3e8 / fc;      % 波长
N_V_BS = 8; N_H_BS = 8;
N_t    = N_V_BS * N_H_BS;
h_BS   = 25;            % 基站高度
U      = 1;             % 用户数（改为1）
N_r    = 2;             % 接收天线数
h_UE   = 1.5;           % 用户高度

T_total   = 15;         % 总快照数
T_history = 10;         % 历史窗口
T_future  = 5;          % 预测窗口

% 轨迹参数
x_start = 100;          % 起始x坐标 (远离基站)
x_step  = 3.5;          % 每个时间步x方向前进距离

% 单用户的正弦参数
amplitude  = 10;        % 振幅（控制摆动幅度）
frequency  = 0.3;       % 频率（控制摆动频率）
phase      = 0;         % 初始相位

% 生成轨迹
pos_seq = zeros(T_total, 1, 3);
for t = 1:T_total
    x = x_start + (t-1) * x_step;  % x方向匀速前进
    y = amplitude * sin(frequency * (t-1) + phase); % y方向正弦摆动
    pos_seq(t, 1, :) = [x; y; h_UE];
end

fprintf('轨迹生成完成。用户初始位置约 %d m，按Sin曲线接近基站。\n', x_start);

%% 2. 生成信道矩阵（LoS模型）
fprintf('\n=== 生成LoS信道矩阵 ===\n');

H_seq = cell(T_total, 1);
for t = 1:T_total
    p_BS = [0; 0; h_BS];
    p_UE = squeeze(pos_seq(t, 1, :));
    d_vec = p_UE - p_BS;
    
    % 计算方位角和仰角
    azimuth   = atan2(d_vec(2), d_vec(1));
    elevation = atan2(d_vec(3), sqrt(d_vec(1)^2 + d_vec(2)^2));
    
    % 生成阵列响应向量
    a_t = get_UPA_response(N_V_BS, N_H_BS, azimuth, elevation);
    a_r = get_UPA_response(1, N_r, azimuth + pi, -elevation);
    
    H_seq{t} = a_r * a_t';
end

fprintf('信道生成完成。\n');

%% 3. 生成理想预编码矩阵（SVD分解）
fprintf('\n=== 生成理想预编码矩阵（Ground Truth）===\n');

F_dig = zeros(N_t, T_total);
for t = 1:T_total
    [~,~,V] = svd(H_seq{t}, 'econ');
    F_dig(:, t) = V(:, 1);
end

%% 4. 执行外插实验（相位对齐）
fprintf('\n=== 开始外插实验 ===\n');

sigma2 = 1;
SNR_dB = 20; 
P_user = 10^(SNR_dB/10);

% 初始化指标存储
err_hold = zeros(T_future, 1);
err_geo  = zeros(T_future, 1);
rate_gt   = zeros(T_future, 1);
rate_hold = zeros(T_future, 1);
rate_geo  = zeros(T_future, 1);

% 获取历史末端点
P_n_minus_1 = F_dig(:, T_history-1);
P_n         = F_dig(:, T_history);

% 相位对齐
phase_diff = angle(P_n' * P_n_minus_1);
P_n_minus_1_aligned = P_n_minus_1 * exp(-1j * phase_diff);

% 计算流形速度
V_back = stiefel_log_approx(P_n, P_n_minus_1_aligned);
V_end = -1.0 * V_back;

% Hold Baseline
P_pred_hold = repmat(P_n, 1, T_future);

% Geodesic Extrapolation
P_pred_geo = zeros(N_t, T_future);
for k = 1:T_future
    delta_t = k;
    P_hat = stiefel_exp_approx(P_n, delta_t * V_end);
    P_pred_geo(:, k) = P_hat;
    
    % 计算指标
    P_true = F_dig(:, T_history + k);
    H_true = H_seq{T_history + k};
    
    err_hold(k) = stiefel_chordal(P_true, P_pred_hold(:, k));
    err_geo(k)  = stiefel_chordal(P_true, P_pred_geo(:, k));
    
    rate_gt(k)   = compute_rate(H_true, P_true, sigma2, P_user);
    rate_hold(k) = compute_rate(H_true, P_pred_hold(:, k), sigma2, P_user);
    rate_geo(k)  = compute_rate(H_true, P_pred_geo(:, k), sigma2, P_user);
end

%% 5. 打印数据表
fprintf('\n=== 正弦轨迹外插结果 (Sinusoidal Trajectory - Single User) ===\n');
fprintf('Step\t Hold误差\t Geo误差\t Rate提升\n');
for k = 1:T_future
    gain = (rate_geo(k) - rate_hold(k)) / rate_hold(k) * 100;
    fprintf('%d\t %.4f\t %.4f\t %.2f%%\n', k, err_hold(k), err_geo(k), gain);
end
fprintf('==========================================\n\n');

%% 6. 可视化

% --- 图1: 用户轨迹 2D 投影图 ---
fig1 = figure('Color', 'w', 'Units', 'inches');
fig1.PaperUnits = 'inches';
fig1.PaperSize = [6 4.5];
fig1.PaperPosition = [0 0 6 4.5];

plot(0, 0, 'k^', 'MarkerSize', 12, 'MarkerFaceColor', 'y', 'DisplayName', 'Base Station');
hold on; grid on; box on;

traj = squeeze(pos_seq(:, 1, :));
plot(traj(:,1), traj(:,2), 'b-', 'LineWidth', 2.5, 'DisplayName', 'User 1');
% 标记起点
plot(traj(1,1), traj(1,2), 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 10, 'DisplayName', 'Start');
% 标记终点
plot(traj(end,1), traj(end,2), 'bs', 'MarkerFaceColor', 'b', 'MarkerSize', 10, 'DisplayName', 'End');

xlabel('x (m)', 'FontSize', 12);
ylabel('y (m)', 'FontSize', 12);
title('User Trajectory: Sinusoidal Path', 'FontSize', 13, 'FontWeight', 'bold');
legend('Location', 'NorthEast', 'FontSize', 10);
axis equal;
print(fig1, 'sinusoidal_trajectories.pdf', '-dpdf', '-r300');
fprintf('已保存轨迹图: sinusoidal_trajectories.pdf\n');
close(fig1);

% --- 图2: 预测误差 ---
fig2 = figure('Color', 'w', 'Units', 'inches');
fig2.PaperUnits = 'inches';
fig2.PaperSize = [5.5 4];
fig2.PaperPosition = [0 0 5.5 4];

ax2 = axes(fig2);
plot(ax2, 1:T_future, err_hold, 'k--^', 'LineWidth', 1.5, 'MarkerSize', 8, 'DisplayName', 'Baseline (Hold)');
hold(ax2, 'on');
plot(ax2, 1:T_future, err_geo, 'r-o', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Proposed (Geo)');
grid(ax2, 'on'); box(ax2, 'on');
ax2.FontSize = 11;
xlabel(ax2, 'Prediction Step $k$', 'Interpreter', 'latex', 'FontSize', 12);
ylabel(ax2, 'Chordal Distance Error', 'FontSize', 12);
title(ax2, 'Prediction Accuracy (Sinusoidal)', 'FontSize', 13, 'FontWeight', 'bold');
legend(ax2, 'Location', 'NorthWest', 'FontSize', 10);
if max(err_hold) > 1e-6, ylim(ax2, [0, max(err_hold)*1.15]); end
print(fig2, 'sinusoidal_result_error.pdf', '-dpdf', '-r300');
fprintf('已保存误差图: sinusoidal_result_error.pdf\n');
close(fig2);

% --- 图3: 频谱效率 ---
fig3 = figure('Color', 'w', 'Units', 'inches');
fig3.PaperUnits = 'inches';
fig3.PaperSize = [5.5 4];
fig3.PaperPosition = [0 0 5.5 4];

ax3 = axes(fig3);
plot(ax3, 1:T_future, rate_gt, 'g-s', 'LineWidth', 1.5, 'MarkerSize', 8, 'Color', [0 0.6 0], 'DisplayName', 'Ideal CSI');
hold(ax3, 'on');
plot(ax3, 1:T_future, rate_geo, 'r-o', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Proposed');
plot(ax3, 1:T_future, rate_hold, 'k--^', 'LineWidth', 1.5, 'MarkerSize', 8, 'DisplayName', 'Baseline');
grid(ax3, 'on'); box(ax3, 'on');
ax3.FontSize = 11;
xlabel(ax3, 'Prediction Step $k$', 'Interpreter', 'latex', 'FontSize', 12);
ylabel(ax3, 'Spectral Efficiency (bps/Hz)', 'FontSize', 12);
title(ax3, 'Rate Performance (Sinusoidal)', 'FontSize', 13, 'FontWeight', 'bold');
legend(ax3, 'Location', 'SouthWest', 'FontSize', 10);
print(fig3, 'sinusoidal_result_rate.pdf', '-dpdf', '-r300');
fprintf('已保存速率图: sinusoidal_result_rate.pdf\n');
close(fig3);

fprintf('\n✓ 所有文件已生成完毕！\n');

%% 诊断分析
fprintf('\n=== 诊断分析：误差与速率的关系 ===\n');
fprintf('Step\t Geo误差\t 理想速率\t Geo速率\t 速率损失(%%)\t 误差平方\n');

for k = 1:T_future
    rate_loss_pct = (rate_gt(k) - rate_geo(k)) / rate_gt(k) * 100;
    err_squared = err_geo(k)^2;
    fprintf('%d\t %.4f\t %.4f\t %.4f\t %.2f%%\t\t %.4f\n', ...
        k, err_geo(k), rate_gt(k), rate_geo(k), rate_loss_pct, err_squared);
end

% 计算相关性：误差平方 vs 速率损失
fprintf('\n速率损失与误差的关系：\n');
rate_losses = (rate_gt - rate_geo) ./ rate_gt * 100;
err_squared = err_geo.^2;
correlation = corrcoef(err_squared, rate_losses);
fprintf('相关系数（误差^2 vs 速率损失）: %.3f\n', correlation(1,2));

%% ==================== 工具函数 ====================

function a = get_UPA_response(Nv, Nh, az, el)
    % UPA阵列响应向量
    m_v = (0:Nv-1).'; 
    m_h = (0:Nh-1).';
    phase_v = exp(1j*pi*m_v*sin(el));
    phase_h = exp(1j*pi*m_h*cos(el).*sin(az));
    a = kron(phase_h, phase_v); 
    a = a/norm(a);
end

function dist = stiefel_chordal(W1, W2)
    % Chordal距离
    overlap = abs(W1' * W2);
    dist = sqrt(max(0, 1 - overlap^2));
end

function rate = compute_rate(H, w, sigma2, P)
    % 单用户速率
    signal_power = P * norm(H * w)^2; 
    rate = log2(1 + signal_power / sigma2);
end

function U = stiefel_log_approx(X, Y)
    % 对数映射
    inner_prod = real(X' * Y);
    theta = acos(max(-1, min(1, inner_prod)));
    if abs(theta) < 1e-10
        U = zeros(size(X));
    else
        U = (theta / sin(theta)) * (Y - inner_prod * X);
    end
end

function Y = stiefel_exp_approx(X, U)
    % 指数映射
    norm_U = norm(U);
    if norm_U < 1e-10
        Y = X;
    else
        Y = X * cos(norm_U) + (U / norm_U) * sin(norm_U);
    end
    Y = Y / norm(Y);
end