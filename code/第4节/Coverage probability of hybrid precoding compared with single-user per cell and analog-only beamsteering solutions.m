% reproduce_fig5_fixed_full_v2.m
% 改进版：single-user per cell 采用随机信道模型，不再是水平线
clear; close all; rng(1);

%% 参数设置
area_side = 1000;
BS_density = 1/10000;
lambda_bs = BS_density;
num_realizations = 300;
U_to_test = [2,3,4,5];
rate_thresholds = linspace(0,3.5,40);
SNR_dB = 10; SNR = 10^(SNR_dB/10);

% 阵列参数
Mbs = 8; Nbs_cols = 8; Nbs = Mbs * Nbs_cols;
Nms_rows = 4; Nms_cols = 4; Nms_total = Nms_rows * Nms_cols;

% 阻塞模型
pLoS_at_1m = 0.9; decay_dist = 200;

steerUPA = @(M,N,az,el) kroneckerUPA_gen(M,N,az,el);

%% 累积变量
cov_hybrid = zeros(length(U_to_test), length(rate_thresholds));
cov_analog = zeros(length(U_to_test), length(rate_thresholds));
cov_single = zeros(1, length(rate_thresholds));

%% Monte Carlo 主循环
for real = 1:num_realizations
    % 基站 PPP
    Nbs_real = poissrnd(lambda_bs * area_side^2);
    if Nbs_real == 0, continue; end
    bs_pos = area_side * rand(Nbs_real,2);
    % 用户
    Nusers_real = max(50, 30 * Nbs_real);
    user_pos = area_side * rand(Nusers_real,2);
    dists = pdist2(user_pos, bs_pos);
    [minD, bs_idx] = min(dists,[],2);
    
    %% ---- Single-user per cell ----
    rates_single_local = [];
    for b = 1:Nbs_real
        assoc = find(bs_idx==b);
        if isempty(assoc), continue; end
        % 随机选一个用户
        uidx = assoc(randi(length(assoc)));
        d = minD(uidx);
        pLoS = pLoS_at_1m * exp(-d/decay_dist);
        isLoS = (rand < pLoS);
        alpha = (randn+1j*randn)/sqrt(2) * (isLoS*1.0 + (~isLoS)*0.2);
        az_bs = 2*pi*rand; az_ms = 2*pi*rand;
        a_bs = steerUPA(Mbs, Nbs_cols, az_bs, pi/2);
        a_ms = steerUPA(Nms_rows, Nms_cols, az_ms, pi/2);
        H = alpha * (a_ms * a_bs'); % 单用户信道
        % 全功率波束成形增益
        gain = abs(a_ms' * H * a_bs)^2;
        rate_single = log2(1 + SNR * gain);
        rates_single_local(end+1) = rate_single;
    end
    if ~isempty(rates_single_local)
        for t = 1:length(rate_thresholds)
            cov_single(t) = cov_single(t) + mean(rates_single_local >= rate_thresholds(t));
        end
    end
    
    %% ---- 多用户场景 ----
    for ni = 1:length(U_to_test)
        nserve = U_to_test(ni);
        rates_hybrid_local = [];
        rates_analog_local = [];
        
        for b = 1:Nbs_real
            assoc = find(bs_idx == b);
            if length(assoc) < nserve, continue; end
            pick = randsample(assoc, nserve);
            Hs = zeros(Nms_total, Nbs, nserve);
            Ws = zeros(Nms_total, nserve);
            FRF = zeros(Nbs, nserve);
            for uu = 1:nserve
                uidx = pick(uu);
                d = minD(uidx);
                pLoS = pLoS_at_1m * exp(-d/decay_dist);
                isLoS = (rand < pLoS);
                alpha = (randn+1j*randn)/sqrt(2) * (isLoS*1.0 + (~isLoS)*0.2);
                az_bs = 2*pi*rand; az_ms = 2*pi*rand;
                a_bs = steerUPA(Mbs, Nbs_cols, az_bs, pi/2);
                a_ms = steerUPA(Nms_rows, Nms_cols, az_ms, pi/2);
                Hs(:,:,uu) = alpha * (a_ms * a_bs');
                FRF(:,uu) = a_bs;
                Ws(:,uu) = a_ms;
            end
            
            % 有效信道
            H_eff = zeros(nserve, nserve);
            for uu=1:nserve
                H_eff(uu,:) = (Ws(:,uu)') * Hs(:,:,uu) * FRF;
            end
            % ZF 预编码
            if rank(H_eff) < nserve
                reg = 1e-3;
                Fbb = H_eff' * ((H_eff*H_eff' + reg*eye(nserve)) \ eye(nserve));
            else
                Fbb = H_eff' * ((H_eff*H_eff') \ eye(nserve));
            end
            for uu=1:nserve
                f = Fbb(:,uu); nf = sqrt(sum(abs(FRF * f).^2));
                if nf>0, Fbb(:,uu)=f/nf; end
            end
            
            % 计算速率
            for uu=1:nserve
                wu = Ws(:,uu);
                num_h = abs(wu' * Hs(:,:,uu) * FRF * Fbb(:,uu))^2 * (SNR/nserve);
                interf_h = 0;
                for nn=1:nserve
                    if nn~=uu
                        interf_h = interf_h + abs(wu' * Hs(:,:,uu) * FRF * Fbb(:,nn))^2 * (SNR/nserve);
                    end
                end
                rate_h = log2(1 + num_h/(interf_h + 1));
                
                f_analog = FRF(:,uu)/norm(FRF(:,uu));
                num_a = abs(wu' * Hs(:,:,uu) * f_analog)^2 * (SNR/nserve);
                interf_a = 0;
                for nn=1:nserve
                    if nn~=uu
                        f_n = FRF(:,nn)/norm(FRF(:,nn));
                        interf_a = interf_a + abs(wu' * Hs(:,:,uu) * f_n)^2 * (SNR/nserve);
                    end
                end
                rate_a = log2(1 + num_a/(interf_a + 1));
                
                rates_hybrid_local(end+1) = rate_h;
                rates_analog_local(end+1) = rate_a;
            end
        end
        
        if ~isempty(rates_hybrid_local)
            for t = 1:length(rate_thresholds)
                cov_hybrid(ni,t) = cov_hybrid(ni,t) + mean(rates_hybrid_local >= rate_thresholds(t));
                cov_analog(ni,t) = cov_analog(ni,t) + mean(rates_analog_local >= rate_thresholds(t));
            end
        end
    end
end

%% 平均
cov_single = cov_single / num_realizations;
cov_hybrid = cov_hybrid / num_realizations;
cov_analog = cov_analog / num_realizations;

%% 绘图
figure; hold on; box on;
plot(rate_thresholds, cov_single, '-k','LineWidth',2.4,'DisplayName','Single-user per cell');
colors = lines(length(U_to_test));
for i=1:length(U_to_test)
    plot(rate_thresholds, cov_hybrid(i,:), '-', 'LineWidth',1.6,'Color',colors(i,:),...
        'DisplayName',sprintf('Hybrid Precoding - %d users',U_to_test(i)));
end
for i=1:length(U_to_test)
    plot(rate_thresholds, cov_analog(i,:), '--', 'LineWidth',1.2,'Color',colors(i,:),...
        'DisplayName',sprintf('Analog-only Beamsteering - %d users',U_to_test(i)));
end
xlabel('Rate Threshold (bps/Hz)');
ylabel('Coverage Probability P(R\ge\eta)');
legend('Location','SouthWest');
grid on;
ylim([0 1]);
yticks(0:0.1:1);
xlim([min(rate_thresholds) max(rate_thresholds)]);

%% 辅助函数
function a = kroneckerUPA_gen(M,N,az,el)
    [mx,my] = meshgrid(0:M-1, 0:N-1);
    mx = mx(:); my = my(:);
    kx = pi * sin(el) * cos(az);
    ky = pi * sin(el) * sin(az);
    a = exp(1j*(mx*kx + my*ky));
    a = a / norm(a);
end
