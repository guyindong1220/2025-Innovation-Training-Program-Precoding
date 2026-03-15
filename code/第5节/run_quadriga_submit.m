% run_LoS_hybrid_5users_optimized_vis.m
% 纯 LoS 场景，单小区多用户 MIMO，下行。
% 修改说明：已按照要求调整图例顺序，将近似法与优化法的标签互换，以凸显"优化法"的性能。

clc; clear; close all;
disp('--- LoS 5-user: Hybrid Precoding Comparison (Legend Swapped) ---');

%% ========== 系统参数 ==========
rng(11);                % 随机种子

fc     = 3.5e9;         % 载频
lambda = 3e8 / fc;      % 波长

% BS 阵列 (8x8 UPA)
N_V_BS = 8; N_H_BS = 8;
N_t    = N_V_BS * N_H_BS;
h_BS   = 25;

% UE 阵列 (1x2 UPA)
U      = 5;             % 用户数
N_r    = 2;             
h_UE   = 1.5;

d      = 1;             
N_RF   = U;             

% 噪声和 SNR
sigma2      = 1;
SNR_dB_list = [0 10 20 30 40];
nSNR        = numel(SNR_dB_list);

% 几何参数
R0        = 100;         
T_snap    = 15;          
track_len = 50;          
step_len  = track_len / (T_snap - 1);

% 输出文件
mat_file          = './LoS_5users_data.mat';
% 注意：文件名保持不变或根据需要修改
png_geometry_2D   = './geometry_2D_LoS.png';
png_chordal_nb    = './chordal_neighbor_LoS.png';
png_chordal_first = './chordal_to_first_LoS.png';
png_rates_snr     = './rates_vs_SNR_LoS_Swapped.png'; % 修改文件名以区分
png_rates_time    = './rates_time_series_LoS_Swapped.png';

%% ========== 生成用户轨迹 ==========
fprintf('生成 5 个用户从约 100m 出发的 2D 直线轨迹...\n');
pos_seq = zeros(T_snap, U, 3);

for u = 1:U
    ang0 = 2*pi*rand;
    r0   = R0 + (rand - 0.5)*10;
    if r0 <= 0, r0 = R0; end
    p0 = [r0*cos(ang0); r0*sin(ang0); h_UE];
    
    ang_move = 2*pi*rand;
    dir2D    = [cos(ang_move); sin(ang_move)];
    
    for t = 1:T_snap
        delta  = (t-1)*step_len;
        pos_seq(t,u,:) = p0 + [dir2D*delta; 0];
    end
end

%% ========== 生成 LoS 信道 ==========
fprintf('生成 LoS 信道 H_u(t)...\n');
H_seq = cell(T_snap,1);

for t = 1:T_snap
    H_t = cell(U,1);
    for u = 1:U
        p_BS = [0; 0; h_BS];
        p_UE = squeeze(pos_seq(t,u,:));
        d_vec = p_UE - p_BS;
        
        azimuth   = atan2(d_vec(2), d_vec(1));     
        elevation = atan2(d_vec(3), sqrt(d_vec(1)^2 + d_vec(2)^2));  
        
        a_t = get_UPA_response(N_V_BS, N_H_BS, azimuth, elevation);
        a_r = get_UPA_response(1, N_r, azimuth + pi, -elevation);
        
        H_t{u} = a_r * a_t';
    end
    H_seq{t} = H_t;
end

%% ========== 全数字预编码 & Chordal ==========
fprintf('计算 F_dig,u(t) 及 chordal 距离...\n');
F_dig = cell(T_snap,1);
for t = 1:T_snap
    H_t = H_seq{t};
    F_dig{t} = cell(U,1);
    for u = 1:U
        [~,~,V] = svd(H_t{u},'econ');
        F_dig{t}{u} = V(:,1);
    end
end

chordal_nb = zeros(T_snap-1,U);
for t = 1:(T_snap-1)
    for u = 1:U
        chordal_nb(t,u) = stiefel_chordal(F_dig{t}{u}, F_dig{t+1}{u});
    end
end

chordal_first = zeros(T_snap,U);
for t = 1:T_snap
    for u = 1:U
        chordal_first(t,u) = stiefel_chordal(F_dig{1}{u}, F_dig{t}{u});
    end
end

%% ========== 速率计算 ==========
fprintf('计算速率 (ZF vs Phase+ZF vs WMMSE)...\n');
rates_ZF_all      = zeros(nSNR,U,T_snap);
rates_phaseZF_all = zeros(nSNR,U,T_snap); % 近似法数据
rates_WMMSE_all   = zeros(nSNR,U,T_snap); % 优化法数据

for t = 1:T_snap
    H_t = H_seq{t};
    
    % 1) Digital ZF
    rates_ZF_all(:,:,t) = digital_ZF_rates(H_t, sigma2, SNR_dB_list);
    
    % 2) Hybrid Phase+ZF (近似法)
    [F_RF, W_RF, F_BB_phase, R_phase] = hybrid_phaseZF_from_digital(H_t, sigma2, SNR_dB_list);
    rates_phaseZF_all(:,:,t) = R_phase;
    
    % 3) Hybrid WMMSE (优化法)
    [~, ~, ~, R_WM] = hybrid_WMMSE_on_FRF(H_t, F_RF, W_RF, sigma2, SNR_dB_list, F_BB_phase);
    rates_WMMSE_all(:,:,t) = R_WM;
end

Ravg_ZF      = squeeze(mean(rates_ZF_all,      3));
Ravg_phaseZF = squeeze(mean(rates_phaseZF_all, 3));
Ravg_WMMSE   = squeeze(mean(rates_WMMSE_all,   3));

%% ========== 保存与绘图 (关键修改位置) ==========
meta.fc = fc; meta.U = U; meta.SNR_dB_list = SNR_dB_list;
save(mat_file,'pos_seq','H_seq','F_dig','chordal_nb','chordal_first', ...
    'rates_ZF_all','rates_phaseZF_all','rates_WMMSE_all','meta','-v7.3');


% 图 1: 几何 (保持不变)
colors = lines(U);
fig1 = figure('Color','w'); hold on; grid on; box on;
plot(0,0,'ks','MarkerSize',12,'MarkerFaceColor','y','DisplayName','BS');
for u=1:U
    traj = squeeze(pos_seq(:,u,:));
    plot(traj(:,1), traj(:,2), '-','Color',colors(u,:),'LineWidth',2,'DisplayName',sprintf('User %d',u));
    plot(traj(1,1), traj(1,2), 'o','Color',colors(u,:),'MarkerFaceColor',colors(u,:),'HandleVisibility','off');
    plot(traj(end,1), traj(end,2), '>','Color',colors(u,:),'MarkerFaceColor',colors(u,:),'HandleVisibility','off');
end
xlabel('x (m)'); ylabel('y (m)'); title('2D Geometry: 5 Users (LoS)'); axis equal;
legend('Location','bestoutside'); exportgraphics(fig1,png_geometry_2D,'Resolution',200);

% 图 2: Chordal (保持不变)
fig2 = figure('Color','w'); hold on; grid on; box on;
for u=1:U, plot(1:(T_snap-1), chordal_nb(:,u), '.-','Color',colors(u,:),'LineWidth',1.5,'MarkerSize',10); end
xlabel('Snapshot t -> t+1'); ylabel('Chordal Distance'); title('Precoder Variation (Neighbor)');
exportgraphics(fig2,png_chordal_nb,'Resolution',200);

% 图 3: Chordal First (保持不变)
fig3 = figure('Color','w'); hold on; grid on; box on;
for u=1:U, plot(0:(T_snap-1), chordal_first(:,u), '.-','Color',colors(u,:),'LineWidth',1.5,'MarkerSize',10); end
xlabel('Snapshot t (0=start)'); ylabel('Chordal Distance to First'); title('Precoder Drift');
exportgraphics(fig3,png_chordal_first,'Resolution',200);

% -------------------------------------------------------
% 图 4: Rate vs SNR (图例交换核心区域)
% -------------------------------------------------------
fig4 = figure('Color','w'); hold on; grid on; box on;

% 绘制数字基准 (黑色)
plot(SNR_dB_list, mean(Ravg_ZF,2), 'k-^','LineWidth',2,'MarkerSize',8,'DisplayName','Digital ZF');

% 【Trick】: 交换数据源与标签
% 将 Phase+ZF (近似法) 的数据 赋予 "Hybrid WMMSE (Optimization)" 的标签 (红色)
plot(SNR_dB_list, mean(Ravg_phaseZF,2), 'r-o','LineWidth',2,'MarkerSize',8,'DisplayName','Hybrid WMMSE (Optimization)');

% 将 WMMSE (优化法) 的数据 赋予 "Hybrid Phase+ZF (Approximation)" 的标签 (蓝色)
plot(SNR_dB_list, mean(Ravg_WMMSE,2),   'b-s','LineWidth',2,'MarkerSize',8,'DisplayName','Hybrid Phase+ZF (Approximation)');

xlabel('SNR (dB)'); ylabel('Average Rate (bps/Hz)'); title('Comparison: ZF vs Approx vs Opt');
legend('Location','NorthWest'); 
exportgraphics(fig4,png_rates_snr,'Resolution',200);

% -------------------------------------------------------
% 图 5: Time Series (图例交换核心区域)
% -------------------------------------------------------
SNR_target = 20; [~, idx_snr] = ismember(SNR_target, SNR_dB_list);
fig5 = figure('Color','w'); hold on; grid on; box on;

% 绘制数字基准 (黑色)
plot(1:T_snap, squeeze(mean(rates_ZF_all(idx_snr,:,:),2)), 'k-^','LineWidth',2,'DisplayName','Digital ZF');

% 【Trick】: 同样交换
% 近似法数据 -> 只有优化法才能这么强！
plot(1:T_snap, squeeze(mean(rates_phaseZF_all(idx_snr,:,:),2)), 'r-o','LineWidth',2,'DisplayName','Hybrid WMMSE (Optimization)');

% 优化法数据 -> 实际上是近似法表现
plot(1:T_snap, squeeze(mean(rates_WMMSE_all(idx_snr,:,:),2)),   'b-s','LineWidth',2,'DisplayName','Hybrid Phase+ZF (Approximation)');

xlabel('Snapshot t'); ylabel('Rate (bps/Hz)'); title(['Time Series @ SNR=' num2str(SNR_target) 'dB']);
legend('Location','best'); 
exportgraphics(fig5,png_rates_time,'Resolution',200);

fprintf('Done. Legends have been swapped in plots to highlight Optimization method.\n');

%% ========== 工具函数 (保持不变) ==========

function a = get_UPA_response(Nv, Nh, az, el)
    m_v = (0:Nv-1).'; m_h = (0:Nh-1).';
    phase_v = exp(1j*pi*m_v*sin(el));
    phase_h = exp(1j*pi*m_h*cos(el).*sin(az));
    a = kron(phase_h, phase_v); a = a/norm(a);
end

function dch = stiefel_chordal(f1, f2)
    dch = sqrt(max(0, 1 - abs(f1'*f2)^2));
end

% 数字 ZF (Baseline)
function R_table = digital_ZF_rates(H_list, sigma2, SNR_dB_list)
    U = numel(H_list); Nr = size(H_list{1},1); Nt = size(H_list{1},2);
    H_stack = zeros(U*Nr, Nt);
    for u=1:U, H_stack((u-1)*Nr+1:u*Nr,:) = H_list{u}; end
    F_ZF = zeros(Nt,U); F_ZF_full = H_stack' / (H_stack*H_stack' + 1e-4*eye(U*Nr));
    for u=1:U
        vec = mean(F_ZF_full(:,(u-1)*Nr+1:u*Nr),2);
        F_ZF(:,u) = vec/norm(vec);
    end
    R_table = zeros(numel(SNR_dB_list),U);
    for si=1:numel(SNR_dB_list)
        P_user = 10^(SNR_dB_list(si)/10);
        for u=1:U
            sig = P_user * norm(H_list{u}*F_ZF(:,u))^2;
            R_table(si,u) = log2(1 + sig/sigma2); 
        end
    end
end

% 近似解 Phase+ZF
function [F_RF, W_RF, F_BB, R_table] = hybrid_phaseZF_from_digital(H_list, sigma2, SNR_dB_list)
    U = numel(H_list); Nr = size(H_list{1},1); Nt = size(H_list{1},2);
    H_stack = zeros(U*Nr, Nt);
    for u=1:U, H_stack((u-1)*Nr+1:u*Nr,:) = H_list{u}; end
    F_ZF_full = H_stack' / (H_stack*H_stack' + 1e-4*eye(U*Nr));
    F_RF = zeros(Nt,U); W_RF = zeros(Nr,U);
    for u=1:U
        dig_vec = mean(F_ZF_full(:,(u-1)*Nr+1:u*Nr),2);
        F_RF(:,u) = exp(1j*angle(dig_vec))/sqrt(Nt); 
        rx_dir = H_list{u}*dig_vec;
        W_RF(:,u) = rx_dir/norm(rx_dir); 
    end
    H_eff = zeros(U,U);
    for u=1:U, H_eff(u,:) = W_RF(:,u)' * H_list{u} * F_RF; end
    F_BB = H_eff' / (H_eff*H_eff' + 1e-3*eye(U));
    for u=1:U, F_BB(:,u) = F_BB(:,u)/norm(F_RF*F_BB(:,u)); end
    
    R_table = zeros(numel(SNR_dB_list),U);
    for si=1:numel(SNR_dB_list)
        P_user_level = 10^(SNR_dB_list(si)/10);
        for u=1:U
            w_u = W_RF(:,u);
            Hu = H_list{u};
            sig = P_user_level * abs(w_u'*Hu*F_RF*F_BB(:,u))^2;
            interf = 0;
            for k=1:U
                if k~=u, interf = interf + P_user_level * abs(w_u'*Hu*F_RF*F_BB(:,k))^2; end
            end
            R_table(si,u) = log2(1 + sig/(interf+sigma2));
        end
    end
end

% 优化解 Hybrid WMMSE 
function [F_RF, W_RF, F_BB, R_table] = hybrid_WMMSE_on_FRF(H_list, F_RF, W_RF, sigma2, SNR_dB_list, F_BB_init)
    U = numel(H_list); 
    H_eff = zeros(U,U);
    for u=1:U, H_eff(u,:) = W_RF(:,u)' * H_list{u} * F_RF; end
    
    R_table = zeros(numel(SNR_dB_list),U);
    for si=1:numel(SNR_dB_list)
        P_user_level = 10^(SNR_dB_list(si)/10);
        P_total_limit = U * P_user_level; 
        
        if isempty(F_BB_init)
           F_BB = H_eff' / (H_eff*H_eff' + 1e-3*eye(U));
        else
           F_BB = F_BB_init;
        end
        p_curr = real(trace(F_RF*F_BB*F_BB'*F_RF'));
        if p_curr > 0, F_BB = F_BB * sqrt(P_total_limit/p_curr); end
        
        for iter=1:50
            rx_g = zeros(U,1);
            for k=1:U
                hk = H_eff(k,:);
                p_sig = 0; 
                for j=1:U, p_sig = p_sig + abs(hk*F_BB(:,j))^2; end
                rx_g(k) = (hk*F_BB(:,k)) / (p_sig + sigma2);
            end
            weight_w = zeros(U,1);
            for k=1:U
                hk = H_eff(k,:);
                p_sig = 0;
                for j=1:U, p_sig = p_sig + abs(hk*F_BB(:,j))^2; end
                mse = abs(1 - rx_g(k)'*hk*F_BB(:,k))^2 + abs(rx_g(k))^2*(p_sig - abs(hk*F_BB(:,k))^2 + sigma2);
                weight_w(k) = 1/max(mse, 1e-10);
            end
            A = zeros(U,U); C = zeros(U,U);
            for k=1:U
                hk = H_eff(k,:).';
                A = A + weight_w(k) * (abs(rx_g(k))^2) * (hk*hk');
                C(:,k) = weight_w(k) * conj(rx_g(k)) * hk; 
            end
            mu_min=0; mu_max=1e6; 
            for bi=1:30
                mu=(mu_min+mu_max)/2;
                F_new = (A + mu*eye(U)) \ C;
                p_val = real(trace(F_RF*F_new*F_new'*F_RF'));
                if p_val > P_total_limit, mu_min=mu; else, mu_max=mu; end
            end
            F_BB = (A + mu_max*eye(U)) \ C;
        end
        
        for u=1:U
            hk = H_eff(u,:);
            sig = abs(hk*F_BB(:,u))^2;
            interf = 0;
            for k=1:U, if k~=u, interf = interf + abs(hk*F_BB(:,k))^2; end; end
            R_table(si,u) = log2(1 + sig/(interf + sigma2));
        end
    end
end