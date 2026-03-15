% stiefel_interpolation_experiment.m
% Stiefel流形上预编码数据的插值实验

clc; clear; close all;
addpath(genpath('.')); % 确保相关函数在路径中

%% ========== 加载原始数据 ==========
fprintf('加载原始数据...\n');
load('./LoS_5users_data.mat', 'F_dig', 'pos_seq', 'H_seq', 'meta');

T_snap = 15;  % 总快照数
U = 5;        % 用户数
N_t = size(F_dig{1}{1}, 1);  % 发送天线数

% 将F_dig转换为更方便的格式
F_dig_mat = zeros(N_t, U, T_snap);
for t = 1:T_snap
    for u = 1:U
        F_dig_mat(:, u, t) = F_dig{t}{u};
    end
end

%% ========== 选择插值策略 ==========
fprintf('\n选择插值策略:\n');
fprintf('1. 均匀采样插值\n');
fprintf('2. 自适应采样插值（基于Chordal距离）\n');
strategy = input('请选择策略 (1或2): ');

if strategy == 1
    % 策略1：均匀采样
    sample_ratio = 0.4;  % 采样率
    n_samples = max(3, round(T_snap * sample_ratio));  % 至少3个点
    sample_indices = round(linspace(1, T_snap, n_samples));
else
    % 策略2：基于Chordal距离的自适应采样
    chordal_dists = zeros(T_snap-1, U);
    for t = 1:(T_snap-1)
        for u = 1:U
            chordal_dists(t, u) = stiefel_chordal(F_dig{t}{u}, F_dig{t+1}{u});
        end
    end
    
    avg_dists = mean(chordal_dists, 2);
    [~, sorted_idx] = sort(avg_dists, 'descend');
    n_samples = min(6, T_snap);  % 最大采样点数
    sample_indices = sort([1, sorted_idx(1:n_samples-2)', T_snap]);
    sample_indices = unique(sample_indices);  % 去重
end

fprintf('采样点索引: %s\n', mat2str(sample_indices));
fprintf('采样点数: %d\n', length(sample_indices));

%% ========== Stiefel流形插值算法 ==========
fprintf('\n开始Stiefel流形插值...\n');

% 初始化插值结果
F_interp_mat = zeros(N_t, U, T_snap);

% 对每个用户分别进行插值
for u = 1:U
    fprintf('处理用户 %d/%d...\n', u, U);
    
    % 获取采样点数据
    X_samples = cell(length(sample_indices), 1);
    for i = 1:length(sample_indices)
        t = sample_indices(i);
        X_samples{i} = F_dig_mat(:, u, t);
    end
    
    % 为每个采样区间创建映射
    num_intervals = length(sample_indices) - 1;
    
    for interval_idx = 1:num_intervals
        % 获取区间端点
        t_start = sample_indices(interval_idx);
        t_end = sample_indices(interval_idx + 1);
        X_start = X_samples{interval_idx};
        X_end = X_samples{interval_idx + 1};
        
        % 计算区间内的测地线方向
        U_interval = stiefel_log_approx(X_start, X_end);
        
        % 计算区间内的所有插值点
        interval_points = t_start:t_end;
        num_points_in_interval = length(interval_points);
        
        for i = 1:num_points_in_interval
            t_curr = interval_points(i);
            
            if t_curr == t_start
                % 起点：直接使用采样点
                F_interp_mat(:, u, t_curr) = X_start;
            elseif t_curr == t_end
                % 终点：直接使用采样点
                F_interp_mat(:, u, t_curr) = X_end;
            else
                % 中间点：按比例插值
                t_param = (t_curr - t_start) / (t_end - t_start);
                F_interp_mat(:, u, t_curr) = stiefel_exp_approx(X_start, t_param * U_interval);
            end
        end
    end
end

%% ========== 验证插值曲线通过采样点 ==========
fprintf('\n验证插值曲线是否通过采样点...\n');
sampled_points = ismember(1:T_snap, sample_indices);

for idx = 1:length(sample_indices)
    t = sample_indices(idx);
    for u = 1:U
        % 检查插值结果是否与原始采样点相同
        F_orig = F_dig_mat(:, u, t);
        F_interp = F_interp_mat(:, u, t);
        error = norm(F_orig - F_interp);
        
        if error > 1e-10
            fprintf('警告: 用户%d在快照%d的插值误差为%g (应接近0)\n', u, t, error);
        end
    end
end
fprintf('验证完成。\n');

%% ========== 计算插值误差 ==========
fprintf('\n计算插值误差...\n');

% Chordal距离误差
chordal_errors = zeros(T_snap, U);
for t = 1:T_snap
    for u = 1:U
        if ~sampled_points(t)  % 只计算非采样点
            F_orig = F_dig_mat(:, u, t);
            F_interp = F_interp_mat(:, u, t);
            chordal_errors(t, u) = stiefel_chordal(F_orig, F_interp);
        end
    end
end

% 重构为矩阵形式（去除采样点）
valid_errors = chordal_errors(~sampled_points, :);
avg_error = mean(valid_errors(:));
max_error = max(valid_errors(:));
std_error = std(valid_errors(:));

fprintf('平均Chordal距离误差: %.4f\n', avg_error);
fprintf('最大Chordal距离误差: %.4f\n', max_error);
fprintf('误差标准差: %.4f\n', std_error);

%% ========== 性能评估：频谱效率 ==========
fprintf('\n评估插值预编码的频谱效率...\n');

SNR_dB = 20;  % 固定SNR
sigma2 = 1;
P_user = 10^(SNR_dB/10);

rates_orig = zeros(T_snap, U);
rates_interp = zeros(T_snap, U);

for t = 1:T_snap
    H_t = H_seq{t};
    
    for u = 1:U
        % 原始预编码的速率
        F_orig_u = F_dig_mat(:, u, t);
        rates_orig(t, u) = compute_rate(H_t{u}, F_orig_u, sigma2, P_user);
        
        % 插值预编码的速率
        F_interp_u = F_interp_mat(:, u, t);
        rates_interp(t, u) = compute_rate(H_t{u}, F_interp_u, sigma2, P_user);
    end
end

rate_loss = mean(abs(rates_orig - rates_interp), 2);
avg_rate_loss = mean(rate_loss(~sampled_points));

fprintf('平均速率损失: %.4f bps/Hz\n', avg_rate_loss);

%% ========== 可视化结果 ==========
fprintf('\n生成可视化结果...\n');

% 图1：采样点与插值点分布
figure('Position', [100, 100, 1200, 800]);
subplot(2, 2, 1);
hold on; grid on; box on;
colors = lines(U);

for u = 1:U
    % 原始轨迹
    traj = squeeze(pos_seq(:, u, 1:2));
    plot(traj(:,1), traj(:,2), '-', 'Color', colors(u,:)*0.7, ...
        'LineWidth', 1.5, 'DisplayName', sprintf('User %d Traj', u));
    
    % 采样点
    plot(traj(sample_indices,1), traj(sample_indices,2), 'o', ...
        'Color', colors(u,:), 'MarkerFaceColor', colors(u,:), ...
        'MarkerSize', 8, 'DisplayName', sprintf('User %d Samples', u));
end

xlabel('x (m)'); ylabel('y (m)');
title('采样点分布 (2D几何)');
legend('Location', 'bestoutside');
axis equal;

% 图2：Chordal距离误差（分用户）
subplot(2, 2, 2);
hold on; grid on; box on;
x_points = 1:T_snap;
sampled_points = ismember(1:T_snap, sample_indices);

for u = 1:U
    errors_u = chordal_errors(:, u);
    % 绘制插值点误差
    interp_idx = ~sampled_points;
    plot(x_points(interp_idx), errors_u(interp_idx), '-', ...
        'Color', colors(u,:), 'LineWidth', 1.5, ...
        'DisplayName', sprintf('User %d', u));
    % 标记采样点（误差为0）
    plot(x_points(sampled_points), zeros(sum(sampled_points),1), 'o', ...
        'Color', colors(u,:), 'MarkerFaceColor', colors(u,:), ...
        'MarkerSize', 6, 'HandleVisibility', 'off');
end

xlabel('快照索引'); ylabel('Chordal距离误差');
title('各用户插值误差');
legend('Location', 'best');
ylim([0, max(chordal_errors(:))*1.1]);


% 图3：频谱效率比较
subplot(2, 2, 3);
hold on; grid on; box on;

% 计算平均速率
avg_rate_orig = mean(rates_orig, 2);
avg_rate_interp = mean(rates_interp, 2);

% 绘制原始速率曲线
plot(1:T_snap, avg_rate_orig, 'k-^', 'LineWidth', 2, ...
    'MarkerSize', 6, 'DisplayName', '原始预编码');

% 绘制插值速率曲线
plot(1:T_snap, avg_rate_interp, 'r--o', 'LineWidth', 2, ...
    'MarkerSize', 6, 'DisplayName', '插值预编码');

% 高亮显示采样点（插值曲线与原始曲线重合的点）
plot(sample_indices, avg_rate_orig(sample_indices), 'ro', ...
    'MarkerSize', 10, 'MarkerFaceColor', 'r', ...
    'DisplayName', '采样点');

xlabel('快照索引'); ylabel('平均频谱效率 (bps/Hz)');
title(sprintf('频谱效率比较 @ SNR=%ddB', SNR_dB));
legend('Location', 'best');

% 图4：误差统计
subplot(2, 2, 4);
hold on; grid on; box on;

% 创建分组柱状图
x_pos = 1:U;
width = 0.35;

for u = 1:U
    % 当前用户的误差（排除采样点）
    user_errors = chordal_errors(~sampled_points, u);
    
    % 计算统计量
    avg_err = mean(user_errors);
    std_err = std(user_errors);
    
    % 绘制柱状图
    bar(x_pos(u) - width/2, avg_err, width, 'FaceColor', colors(u,:), ...
        'EdgeColor', colors(u,:)*0.7, 'LineWidth', 1.5);
    
    % 误差条
    errorbar(x_pos(u) - width/2, avg_err, std_err, 'k', 'LineWidth', 1.5);
end

xlabel('用户索引'); ylabel('Chordal距离误差');
title('各用户误差统计（均值±标准差）');
set(gca, 'XTick', 1:U, 'XTickLabel', arrayfun(@(x) sprintf('User %d', x), 1:U, 'UniformOutput', false));
xlim([0.5, U+0.5]);

% 添加整体统计信息
text(0.6, max(ylim)*0.9, sprintf('平均误差: %.4f', avg_error), ...
    'FontSize', 10, 'FontWeight', 'bold');
text(0.6, max(ylim)*0.85, sprintf('最大误差: %.4f', max_error), ...
    'FontSize', 10, 'FontWeight', 'bold');

sgtitle('Stiefel流形插值实验结果', 'FontSize', 14, 'FontWeight', 'bold');
saveas(gcf, './stiefel_interpolation_results.png');

%% ========== 保存结果 ==========
save('./interpolation_results.mat', 'F_interp_mat', 'sample_indices', ...
    'chordal_errors', 'rates_orig', 'rates_interp', 'avg_error', ...
    'max_error', 'std_error', 'avg_rate_loss');

fprintf('\n实验完成！结果已保存。\n');

%% ========== 工具函数 ==========

function d = stiefel_chordal(X, Y)
    % 计算Stiefel流形上的Chordal距离
    d = sqrt(max(0, 1 - abs(X'*Y)^2));
end

function rate = compute_rate(H, F, sigma2, P_user)
    % 计算单用户速率
    sig = P_user * norm(H*F)^2;
    rate = log2(1 + sig/sigma2);
end

function U = stiefel_log_approx(X, Y)
    % Stiefel流形对数映射（近似）
    % 对于St(n,1)（单位球面），使用球面对数映射
    
    n = size(X, 1);
    p = size(X, 2);
    
    if p == 1  % 单位球面情况
        cos_theta = real(X' * Y);
        cos_theta = max(-1, min(1, cos_theta));  % 数值稳定性
        theta = acos(cos_theta);
        
        if abs(theta) < 1e-10
            U = zeros(size(X));
        else
            U = (theta / sin(theta)) * (Y - cos_theta * X);
        end
    else
        % 一般Stiefel流形情况（简化处理）
        S = (eye(n) - X*X') * Y * pinv(X'*Y);
        U = X * real(logm(X'*Y)) + S;
    end
end

function Y = stiefel_exp_approx(X, U)
    % Stiefel流形指数映射（近似）
    % 对于St(n,1)（单位球面），使用球面指数映射
    
    n = size(X, 1);
    p = size(X, 2);
    
    if p == 1  % 单位球面情况
        norm_U = norm(U);
        if norm_U < 1e-10
            Y = X;
        else
            Y = cos(norm_U) * X + sin(norm_U) * (U / norm_U);
        end
    else
        % 一般Stiefel流形情况（使用QR分解）
        [Q, R] = qr(U, 0);
        Y = X + Q * R;
        [Y, ~] = qr(Y, 0);  % 重新正交化
    end
    
    % 确保在流形上（数值稳定性）
    if p == 1
        Y = Y / norm(Y);
    end
end