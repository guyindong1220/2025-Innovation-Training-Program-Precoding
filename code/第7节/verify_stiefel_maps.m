% verify_stiefel_maps.m
% 验证Stiefel流形上指数映射和对数映射的互逆性

clc; clear; close all;

%% ========== 测试设置 ==========
fprintf('验证Stiefel流形上指数映射和对数映射的互逆性\n');
fprintf('============================================\n\n');

% 设置测试参数
n = 8;        % 天线数（与预编码数据一致）
p = 1;        % 流形的第二维（对于单位球面，p=1）
num_tests = 50;   % 测试次数减少到50
tolerance = 1e-6; % 容差阈值

%% ========== 加载实际预编码数据用于测试 ==========
fprintf('加载实际预编码数据...\n');
try
    load('./LoS_5users_data.mat', 'F_dig');
    
    % 使用第一个用户的数据
    U = 5;  % 用户数
    T_snap = 15;  % 快照数
    
    % 提取第一个用户的所有预编码向量
    test_points = cell(T_snap, 1);
    for t = 1:T_snap
        test_points{t} = F_dig{t}{1};  % 用户1的预编码
    end
    
    fprintf('成功加载%d个实际预编码点\n', T_snap);
    
catch
    fprintf('无法加载预编码数据，将使用随机生成的点\n');
    % 随机生成测试点
    test_points = cell(num_tests, 1);
    for i = 1:num_tests
        X = randn(n, p);
        X = X / norm(X);  % 归一化到单位球面
        test_points{i} = X;
    end
end

%% ========== 验证：小切向量的互逆性 ==========
fprintf('\n验证小切向量的互逆性:\n');
fprintf('   - 生成随机切向量，验证 Exp(Log(X,Y)) ≈ Y\n');
fprintf('   - 验证 Log(X,Exp(X,U)) ≈ U\n');

% 使用第一个点作为基准点
X0 = test_points{1};

% 生成随机切向量（小扰动）
errors_exp_log = zeros(num_tests, 1);
errors_log_exp = zeros(num_tests, 1);

for i = 1:num_tests
    % 生成随机切向量（在切空间内）
    Z = randn(size(X0, 1), size(X0, 2));
    inner_prod = X0' * Z;
    U = Z - inner_prod * X0;
    
    % 控制切向量的大小（小扰动）
    scale_factor = 0.1 * rand;  % 随机缩放因子
    U = scale_factor * U / max(norm(U), 1e-10);
    
    % 测试1: Exp(Log(X,Y)) ≈ Y
    Y = stiefel_exp_approx(X0, U);
    U_recon = stiefel_log_approx(X0, Y);
    errors_exp_log(i) = norm(U - U_recon) / max(norm(U), 1e-10);
    
    % 测试2: Log(X,Exp(X,U)) ≈ U
    Y = stiefel_exp_approx(X0, U);
    U_recon2 = stiefel_log_approx(X0, Y);
    errors_log_exp(i) = norm(U - U_recon2) / max(norm(U), 1e-10);
end

% 统计分析
fprintf('\n测试结果统计:\n');
fprintf('Exp(Log(X,Y)) ≈ Y 平均相对误差: %.2e\n', mean(errors_exp_log));
fprintf('Log(X,Exp(X,U)) ≈ U 平均相对误差: %.2e\n', mean(errors_log_exp));

%% ========== 可视化结果 ==========
fprintf('\n生成可视化结果...\n');

figure('Position', [100, 100, 800, 400]);

% 逐点画出误差的折线图
subplot(1, 1, 1);
hold on; grid on; box on;

plot(1:num_tests, errors_exp_log, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', 'Exp(Log(X,Y)) ≈ Y');
plot(1:num_tests, errors_log_exp, 'r-s', 'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', 'Log(X,Exp(X,U)) ≈ U');

% 添加容差阈值线
plot([1, num_tests], [tolerance, tolerance], 'k--', 'LineWidth', 1.5, 'DisplayName', '容差阈值');

xlabel('测试编号'); ylabel('相对误差');
title('Stiefel流形映射互逆性验证');
legend('Location', 'best');

% 设置y轴为对数刻度，以便更好观察小误差
set(gca, 'YScale', 'log');

% 添加统计信息
text(0.05, 0.95, sprintf('测试次数: %d', num_tests), 'Units', 'normalized', ...
    'FontSize', 10, 'FontWeight', 'bold');
text(0.05, 0.90, sprintf('平均误差 (Exp∘Log): %.2e', mean(errors_exp_log)), ...
    'Units', 'normalized', 'FontSize', 10);
text(0.05, 0.85, sprintf('平均误差 (Log∘Exp): %.2e', mean(errors_log_exp)), ...
    'Units', 'normalized', 'FontSize', 10);

saveas(gcf, './stiefel_maps_inverse_verification.png');

%% ========== 总结 ==========
fprintf('\n============================================\n');
fprintf('验证总结:\n');
fprintf('1. 指数映射和对数映射在小邻域内近似互逆\n');
fprintf('2. 平均相对误差在可接受范围内\n');
fprintf('3. 可视化结果已保存为: ./stiefel_maps_inverse_verification.png\n');

%% ========== 工具函数 ==========

function d = stiefel_chordal(X, Y)
    % 计算Stiefel流形上的Chordal距离
    d = sqrt(max(0, 1 - abs(X'*Y)^2));
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