clear all; close all; clc;

%% 参数设置
rng(123); % 设置随机种子，确保结果可重现
K = 1;          % 单小区
I = 4;          % 4个用户
T = 8;          % 发射天线数
R = 2;          % 每个用户的接收天线数
d = 1;          % 每个用户的数据流数
SNR_dB = 10;    % 信噪比 (dB)
num_channels = 10; % 信道实现次数
max_iter = 1000; % 最大迭代次数
epsilon = 1e-4; % 收敛容差

% 权重设置 (所有用户平等)
alpha = ones(I, 1);

% 噪声功率 (所有用户相同)
sigma2 = 1;

% 发射功率
P = sigma2*10^(SNR_dB/10); 

%% 存储结果
sum_rate_wmmse = zeros(num_channels, 1);
sum_rate_mmse = zeros(num_channels, 1);

%% 主循环 - 对不同信道实现
for ch_idx = 1:num_channels
    fprintf('处理信道实现 %d/%d\n', ch_idx, num_channels);
    
    % 生成随机信道 (CN(0,1))
    H = cell(I, 1);
    for i = 1:I
        if mod(i,2) == 1
            % 好用户：信道增益大
            H{i} = 2.0 * (randn(R,T) + 1j*randn(R,T))/sqrt(2);
        else
            % 差用户：信道增益小
            H{i} = 0.3 * (randn(R,T) + 1j*randn(R,T))/sqrt(2);
        end
    end
    
    %% WMMSE 算法
    % 初始化
    V_wmmse = cell(I, 1);
    for i = 1:I
        V_wmmse{i} = sqrt(P/(I*d)) * (randn(T, d) + 1j * randn(T, d)) / sqrt(2);
    end
    
    prev_obj = -inf;
    for iter = 1:max_iter
        % 步骤1: 更新接收波束成形器 U (MMSE接收机)
        U = cell(I, 1);
        J = cell(I, 1);
        for i = 1:I
            % 计算干扰加噪声协方差矩阵
            J{i} = zeros(R, R);
            for j = 1:I
                J{i} = J{i} + H{i} * V_wmmse{j} * V_wmmse{j}' * H{i}';
            end
            J{i} = J{i} + sigma2 * eye(R);
            
            % MMSE接收机
            U{i} = J{i} \ (H{i} * V_wmmse{i});
        end
        
        % 步骤2: 更新权重矩阵 W
        W = cell(I, 1);
        E = cell(I, 1);
        for i = 1:I
            % 计算MSE矩阵
            E{i} = eye(d) - U{i}' * H{i} * V_wmmse{i};
            E{i} = E{i} * E{i}';
            for j = 1:I
                if j ~= i
                    E{i} = E{i} + U{i}' * H{i} * V_wmmse{j} * V_wmmse{j}' * H{i}' * U{i};
                end
            end
            E{i} = E{i} + sigma2 * U{i}' * U{i};
            
            % 最优权重矩阵
            W{i} = inv(E{i});
        end
        
        % 步骤3: 更新发射波束成形器 V
        % 使用二分法求解最优mu
        mu_min = 0;
        mu_max = 1e6;
        mu_tol = 1e-6;
        
        for bisect_iter = 1:100
            mu = (mu_min + mu_max) / 2;
            
            % 计算公共矩阵 A
            A = zeros(T, T);
            for i = 1:I
                A = A + alpha(i) * H{i}' * U{i} * W{i} * U{i}' * H{i};
            end
            A = A + mu * eye(T);
            
            % 计算新的V
            V_new = cell(I, 1);
            total_power = 0;
            for i = 1:I
                V_new{i} = A \ (alpha(i) * H{i}' * U{i} * W{i});
                total_power = total_power + trace(V_new{i} * V_new{i}');
            end
            
            % 更新二分法区间
            if total_power > P
                mu_min = mu;
            else
                mu_max = mu;
            end
            
            if (mu_max - mu_min) < mu_tol
                break;
            end
        end
        V_wmmse = V_new;
        
        % 计算当前目标函数值 (加权和速率)
        current_obj = 0;
        for i = 1:I
            % 计算用户i的干扰加噪声协方差矩阵
            J_i = zeros(R, R);
            for j = 1:I
                if j ~= i
                    J_i = J_i + H{i} * V_wmmse{j} * V_wmmse{j}' * H{i}';
                end
            end
            J_i = J_i + sigma2 * eye(R);
            
            % 计算用户i的速率
            signal_cov = H{i} * V_wmmse{i} * V_wmmse{i}' * H{i}';
            R_i = real(log2(det(eye(R) + signal_cov / J_i)));
            current_obj = current_obj + alpha(i) * R_i;
        end
        
        % 检查收敛
        if iter > 1 && abs(current_obj - prev_obj) < epsilon
            fprintf("WMMSE在第%d步收敛\n",iter);
            break;
        end
        prev_obj = current_obj;
    end
    
    % 记录WMMSE的和速率
    sum_rate_wmmse(ch_idx) = current_obj;
    
    %% 传统MMSE算法 (最小化Sum-MSE)
    % 初始化
    V_mmse = cell(I, 1);
    for i = 1:I
        V_mmse{i} = sqrt(P/(I*d)) * (randn(T, d) + 1j * randn(T, d)) / sqrt(2);
    end
    
    prev_obj_mmse = inf;
    for iter = 1:max_iter
        % 更新接收波束成形器 U (MMSE接收机)
        U_mmse = cell(I, 1);
        for i = 1:I
            J_mmse = zeros(R, R);
            for j = 1:I
                J_mmse = J_mmse + H{i} * V_mmse{j} * V_mmse{j}' * H{i}';
            end
            J_mmse = J_mmse + sigma2 * eye(R);
            U_mmse{i} = J_mmse \ (H{i} * V_mmse{i});
        end
        
        % 更新发射波束成形器 V (最小化Sum-MSE)
        mu_min_mmse = 0;
        mu_max_mmse = 1e6;
        
        for bisect_iter = 1:100
            mu_mmse = (mu_min_mmse + mu_max_mmse) / 2;
            
            % 计算公共矩阵 (对于Sum-MSE最小化)
            A_mmse = zeros(T, T);
            for i = 1:I
                A_mmse = A_mmse + H{i}' * U_mmse{i} * U_mmse{i}' * H{i};
            end
            A_mmse = A_mmse + mu_mmse * eye(T);
            
            % 计算新的V
            V_new_mmse = cell(I, 1);
            total_power_mmse = 0;
            for i = 1:I
                V_new_mmse{i} = A_mmse \ (H{i}' * U_mmse{i});
                total_power_mmse = total_power_mmse + trace(V_new_mmse{i} * V_new_mmse{i}');
            end
            
            if total_power_mmse > P
                mu_min_mmse = mu_mmse;
            else
                mu_max_mmse = mu_mmse;
            end
            
            if (mu_max_mmse - mu_min_mmse) < 1e-6
                break;
            end
        end
        V_mmse = V_new_mmse;
        
        % 计算当前Sum-MSE
        current_obj_mmse = 0;
        for i = 1:I
            E_mmse = eye(d) - U_mmse{i}' * H{i} * V_mmse{i};
            E_mmse = E_mmse * E_mmse';
            for j = 1:I
                if j ~= i
                    E_mmse = E_mmse + U_mmse{i}' * H{i} * V_mmse{j} * V_mmse{j}' * H{i}' * U_mmse{i};
                end
            end
            E_mmse = E_mmse + sigma2 * U_mmse{i}' * U_mmse{i};
            current_obj_mmse = current_obj_mmse + real(trace(E_mmse));
        end
        
        % 检查收敛
        if iter > 1 && abs(current_obj_mmse - prev_obj_mmse) < epsilon
            fprintf("MMSE在第%d步收敛\n",iter);
            break;
        end
        prev_obj_mmse = current_obj_mmse;
    end
    
    % 计算MMSE算法的和速率
    sum_rate_mmse_temp = 0;
    for i = 1:I
        J_i_mmse = zeros(R, R);
        for j = 1:I
            if j ~= i
                J_i_mmse = J_i_mmse + H{i} * V_mmse{j} * V_mmse{j}' * H{i}';
            end
        end
        J_i_mmse = J_i_mmse + sigma2 * eye(R);
        signal_cov_mmse = H{i} * V_mmse{i} * V_mmse{i}' * H{i}';
        R_i_mmse = real(log2(det(eye(R) + signal_cov_mmse / J_i_mmse)));
        sum_rate_mmse_temp = sum_rate_mmse_temp + R_i_mmse;
    end
    sum_rate_mmse(ch_idx) = sum_rate_mmse_temp;
end

%% 结果显示
fprintf('\n=== 性能对比结果 (单小区, %d个用户) ===\n', I);
fprintf('信道实现次数: %d\n', num_channels);
fprintf('SNR: %d dB\n', SNR_dB);
fprintf('发射天线: %d, 接收天线: %d\n', T, R);
fprintf('\n');

fprintf('WMMSE算法平均和速率: %.4f bps/Hz\n', mean(sum_rate_wmmse));
fprintf('MMSE算法平均和速率: %.4f bps/Hz\n', mean(sum_rate_mmse));
fprintf('性能增益: %.4f bps/Hz (%.2f%%)\n', ...
    mean(sum_rate_wmmse) - mean(sum_rate_mmse), ...
    (mean(sum_rate_wmmse) - mean(sum_rate_mmse)) / mean(sum_rate_mmse) * 100);

%% 绘制结果
figure;
plot(1:num_channels, sum_rate_wmmse, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 4);
hold on;
plot(1:num_channels, sum_rate_mmse, 'r-s', 'LineWidth', 1.5, 'MarkerSize', 4);
xlabel('信道实现索引');
ylabel('和速率 (bps/Hz)');
title(['单小区多用户MIMO系统和速率对比 (K=1, I=' num2str(I) ', SNR=' num2str(SNR_dB) 'dB)']);
legend('WMMSE算法', 'MMSE算法', 'Location', 'best');
grid on;

figure;
boxplot([sum_rate_wmmse, sum_rate_mmse], 'Labels', {'WMMSE', 'MMSE'});
ylabel('和速率 (bps/Hz)');
title('和速率分布对比');
grid on;