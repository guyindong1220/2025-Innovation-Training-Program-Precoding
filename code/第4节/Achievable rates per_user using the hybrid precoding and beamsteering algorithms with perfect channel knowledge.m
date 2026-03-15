clear all; close all; clc;

%% 系统参数设置
rng(123); % 设置随机数种子以保证结果可重现

% 基站参数
N_BS = 64; % 8x8 UPA = 64天线
Nx_BS = 8; Ny_BS = 8;
N_RF = 4;  % RF链数量

% 移动台参数
N_MS = 16; % 4x4 UPA = 16天线
Nx_MS = 4; Ny_MS = 4;
U = 4;     % 用户数量

% 仿真参数
SNR_dB = -10:5:20; % SNR范围 (dB)
SNR_lin = 10.^(SNR_dB/10); % 线性SNR
num_channels = 1000; % 信道实现次数

% 存储速率结果
rate_hybrid = zeros(length(SNR_dB), U);
rate_single_user = zeros(length(SNR_dB), U);
rate_beamsteering = zeros(length(SNR_dB), U);

%% 主仿真循环
for chan_idx = 1:num_channels
    if mod(chan_idx, 100) == 0
        fprintf('处理第 %d/%d 个信道...\n', chan_idx, num_channels);
    end
    
    %% 生成单路径信道
    H = cell(U, 1);
    alpha = zeros(U, 1);
    phi_az = zeros(U, 1); phi_el = zeros(U, 1); % BS AoD
    theta_az = zeros(U, 1); theta_el = zeros(U, 1); % MS AoA
    
    for u = 1:U
        % 随机生成信道参数
        alpha(u) = (randn(1) + 1j*randn(1))/sqrt(2); % 复高斯增益
        phi_az(u) = 2*pi*rand(1);   % BS方位角AoD [0, 2π]
        phi_el(u) = pi*rand(1) - pi/2; % BS俯仰角AoD [-π/2, π/2]
        theta_az(u) = 2*pi*rand(1); % MS方位角AoA [0, 2π]
        theta_el(u) = pi*rand(1) - pi/2; % MS俯仰角AoA [-π/2, π/2]
        
        % 生成阵列响应向量
        a_BS = get_UPA_response(phi_az(u), phi_el(u), Nx_BS, Ny_BS);
        a_MS = get_UPA_response(theta_az(u), theta_el(u), Nx_MS, Ny_MS);
        
        % 构建信道矩阵 (公式(4)，单路径情况)
        H{u} = sqrt(N_BS*N_MS) * alpha(u) * (a_MS * a_BS');
    end
    
    %% 第一阶段：模拟波束成形设计
    F_RF = zeros(N_BS, U);
    W = zeros(N_MS, U);
    
    for u = 1:U
        % 最优模拟波束成形和合并向量 (公式(7))
        a_BS = get_UPA_response(phi_az(u), phi_el(u), Nx_BS, Ny_BS);
        a_MS = get_UPA_response(theta_az(u), theta_el(u), Nx_MS, Ny_MS);
        
        F_RF(:, u) = a_BS;
        W(:, u) = a_MS;
    end
    
    %% 第二阶段：数字预编码设计
    % 计算有效信道 (公式(8))
    H_eff = zeros(U, U);
    for u = 1:U
        H_eff(u, :) = W(:, u)' * H{u} * F_RF;
    end
    
    % Zero-Forcing 数字预编码 (公式(10))
    F_BB = H_eff' / (H_eff * H_eff');
    
    % 功率归一化
    for u = 1:U
        F_BB(:, u) = F_BB(:, u) / norm(F_RF * F_BB(:, u), 'fro');
    end
    
    %% 计算各种方案的速率
    for snr_idx = 1:length(SNR_lin)
        P = SNR_lin(snr_idx); % 总发射功率
        sigma2 = 1; % 噪声功率归一化为1
        
        for u = 1:U
            %% 1. 混合预编码速率 (公式(5))
            signal_power = P/U * abs(W(:, u)' * H{u} * F_RF * F_BB(:, u))^2;
            interference_power = 0;
            for n = 1:U
                if n ~= u
                    interference_power = interference_power + ...
                        P/U * abs(W(:, u)' * H{u} * F_RF * F_BB(:, n))^2;
                end
            end
            SINR_hybrid = signal_power / (interference_power + sigma2);
            rate_hybrid(snr_idx, u) = rate_hybrid(snr_idx, u) + ...
                log2(1 + SINR_hybrid);
            
            %% 2. 单用户速率 (无干扰情况)
            signal_power_single = P * abs(W(:, u)' * H{u} * F_RF(:, u))^2;
            SINR_single = signal_power_single / sigma2;
            rate_single_user(snr_idx, u) = rate_single_user(snr_idx, u) + ...
                log2(1 + SINR_single);
            
            %% 3. 纯模拟波束成形速率 (有干扰)
            signal_power_bs = P/U * abs(W(:, u)' * H{u} * F_RF(:, u))^2;
            interference_power_bs = 0;
            for n = 1:U
                if n ~= u
                    interference_power_bs = interference_power_bs + ...
                        P/U * abs(W(:, u)' * H{u} * F_RF(:, n))^2;
                end
            end
            SINR_bs = signal_power_bs / (interference_power_bs + sigma2);
            rate_beamsteering(snr_idx, u) = rate_beamsteering(snr_idx, u) + ...
                log2(1 + SINR_bs);
        end
    end
end

%% 平均速率计算
rate_hybrid_avg = mean(rate_hybrid, 2) / num_channels;
rate_single_user_avg = mean(rate_single_user, 2) / num_channels;
rate_beamsteering_avg = mean(rate_beamsteering, 2) / num_channels;

%% 绘图 (复现Fig.3)
figure;
plot(SNR_dB, rate_single_user_avg, 'b-o', 'LineWidth', 2, 'MarkerSize', 8);
hold on;
plot(SNR_dB, rate_hybrid_avg, 'r-s', 'LineWidth', 2, 'MarkerSize', 8);
plot(SNR_dB, rate_beamsteering_avg, 'k-^', 'LineWidth', 2, 'MarkerSize', 8);
grid on;

xlabel('SNR (dB)', 'FontSize', 12);
ylabel('Per-user Achievable Rate (bps/Hz)', 'FontSize', 12);
title('Achievable Rates with Perfect Channel Knowledge', 'FontSize', 14);
legend('Single-user', 'Hybrid Precoding', 'Beamsteering', 'Location', 'northwest');
set(gca, 'FontSize', 11);

%% UPA阵列响应向量生成函数
function a = get_UPA_response(azimuth, elevation, Nx, Ny)
    % 生成UPA阵列响应向量
    % azimuth: 方位角 [0, 2π]
    % elevation: 俯仰角 [-π/2, π/2]
    % Nx, Ny: x和y方向的天线数量
    
    % x方向相位
    kx = 0:Nx-1;
    phase_x = exp(1j * pi * kx * sin(azimuth) * cos(elevation));
    
    % y方向相位
    ky = 0:Ny-1;
    phase_y = exp(1j * pi * ky * sin(elevation));
    
    % Kronecker积得到完整的阵列响应
    a = kron(phase_x(:), phase_y(:));
    a = a / norm(a); % 归一化
end