% 初始化
F_RF = []; % 模拟预编码矩阵
F_Res = F_Opt; % 残差矩阵

for r = 1:NumRF
    % Step 4: 计算路径在最优预编码器上的投影
    Psi = At'*F_Res;
    % Step 5: 选择投影最大的路径
    [~,k] = max(diag(Psi*Psi'));
    % Step 6: 添加选中的向量到RF预编码器
    F_RF = [F_RF At(:,k)];
    % Step 7: 计算基带预编码器（最小二乘解）
    F_BB = (F_RF'*F_RF)\(F_RF'*F_Opt);
    % Step 8: 更新残差
    F_Res = (F_Opt-F_RF*F_BB)/norm(F_Opt-F_RF*F_BB,'fro');
end
% Step 10: 功率归一化
F_BB = sqrt(Ns)*(F_BB/norm(F_RF*F_BB,'fro'));

% 计算接收信号协方差矩阵
CovRx = (SNR/Ns)*H*F_RF*F_BB*F_BB'*F_RF'*H'+eye(Nr);
W_MMSE = ((1/sqrt(SNR))*(F_BB'*F_RF'*H'*H*F_RF*F_BB+(Ns/SNR)*eye(Ns))\(F_BB'*F_RF'*H'))';

W_RF = [];
W_Res = W_MMSE;

for r = 1:NumRF
    Psi = Ar'*CovRx*W_Res;
    [~,k] = max(diag(Psi*Psi'));
    W_RF = [W_RF Ar(:,k)];
    W_BB = (W_RF'*CovRx*W_RF)\(W_RF'*CovRx*W_MMSE);
    W_Res = (W_MMSE-W_RF*W_BB)/norm(W_MMSE-W_RF*W_BB,'fro');
end