% stiefel_utils.m
% Stiefel流形工具函数

function [U, info] = stiefel_log_numerical(X, Y, max_iter, tol)
    % Stiefel流形对数映射的数值计算方法
    % 使用梯度下降法
    
    if nargin < 3
        max_iter = 100;
    end
    if nargin < 4
        tol = 1e-6;
    end
    
    n = size(X, 1);
    p = size(X, 2);
    
    % 初始化：投影到切空间
    U = Y - X * (X' * Y);
    [U, ~] = qr(U, 0);  # 正交化
    
    % 迭代优化
    for iter = 1:max_iter
        % 计算指数映射
        Y_est = stiefel_exp_approx(X, U);
        
        % 计算残差
        R = Y_est - Y;
        res_norm = norm(R, 'fro');
        
        % 检查收敛
        if res_norm < tol
            break;
        end
        
        % 计算梯度（近似）
        JtR = stiefel_jacobian_transpose(X, U, R);
        G = stiefel_projection(X, JtR);
        
        % 线搜索
        alpha = 0.5;
        for ls_iter = 1:10
            U_new = U - alpha * G;
            Y_new = stiefel_exp_approx(X, U_new);
            res_new = norm(Y_new - Y, 'fro');
            
            if res_new < res_norm
                U = U_new;
                break;
            end
            alpha = alpha * 0.5;
        end
    end
    
    info.iter = iter;
    info.residual = res_norm;
end

function P = stiefel_projection(X, Z)
    % 投影到切空间 T_X St(n,p)
    P = Z - X * ((X' * Z + Z' * X) / 2);
end

function JtR = stiefel_jacobian_transpose(X, U, R)
    % 指数映射Jacobian的转置作用（近似）
    n = size(X, 1);
    p = size(X, 2);
    
    if p == 1
        % 球面情况的简化
        norm_U = norm(U);
        if norm_U < 1e-10
            JtR = R;
        else
            u_unit = U / norm_U;
            x_dot = -sin(norm_U) * norm_U * (X' * R) * u_unit + ...
                    cos(norm_U) * (R - (u_unit' * R) * u_unit);
            JtR = x_dot;
        end
    else
        % 一般情况的近似
        [Q, R_qr] = qr(U, 0);
        JtR = R * R_qr' + X * (X' * R * U');
    end
end

function d = riemannian_distance(X, Y)
    % 计算Stiefel流形上的黎曼距离
    % 对于球面情况，使用测地线距离
    
    p = size(X, 2);
    
    if p == 1
        % 球面情况
        cos_theta = real(X' * Y);
        cos_theta = max(-1, min(1, cos_theta));
        d = acos(cos_theta);
    else
        % 一般Stiefel流形情况（近似）
        [U, info] = stiefel_log_numerical(X, Y);
        d = norm(U, 'fro');
    end
end

function [curves, metrics] = evaluate_interpolation_quality(F_orig, F_interp, sample_indices)
    % 评估插值质量
    
    T = size(F_orig, 3);
    U = size(F_orig, 2);
    
    metrics = struct();
    
    % 1. Chordal距离
    chordal_errors = zeros(T, U);
    for t = 1:T
        for u = 1:U
            if ~ismember(t, sample_indices)
                chordal_errors(t, u) = stiefel_chordal(...
                    F_orig(:, u, t), F_interp(:, u, t));
            end
        end
    end
    metrics.chordal_mean = mean(chordal_errors(:));
    metrics.chordal_std = std(chordal_errors(:));
    metrics.chordal_max = max(chordal_errors(:));
    
    % 2. 黎曼距离
    riemann_errors = zeros(T, U);
    for t = 1:T
        for u = 1:U
            if ~ismember(t, sample_indices)
                riemann_errors(t, u) = riemannian_distance(...
                    F_orig(:, u, t), F_interp(:, u, t));
            end
        end
    end
    metrics.riemann_mean = mean(riemann_errors(:));
    metrics.riemann_std = std(riemann_errors(:));
    
    % 3. 角度误差（对于球面）
    if size(F_orig, 2) == 1
        angle_errors = zeros(T, 1);
        for t = 1:T
            if ~ismember(t, sample_indices)
                cos_err = abs(F_orig(:, 1, t)' * F_interp(:, 1, t));
                angle_errors(t) = acosd(min(1, cos_err));
            end
        end
        metrics.angle_mean = mean(angle_errors);
        metrics.angle_max = max(angle_errors);
    end
    
    % 4. 插值曲线的平滑度
    curves.smoothness = compute_curve_smoothness(F_interp);
    curves.variation = compute_total_variation(F_interp);
    
    fprintf('插值质量评估:\n');
    fprintf('  Chordal距离 - 均值: %.4f, 标准差: %.4f, 最大值: %.4f\n', ...
        metrics.chordal_mean, metrics.chordal_std, metrics.chordal_max);
    fprintf('  黎曼距离 - 均值: %.4f, 标准差: %.4f\n', ...
        metrics.riemann_mean, metrics.riemann_std);
    
    if isfield(metrics, 'angle_mean')
        fprintf('  角度误差 - 均值: %.2f°, 最大值: %.2f°\n', ...
            metrics.angle_mean, metrics.angle_max);
    end
end

function smoothness = compute_curve_smoothness(F)
    % 计算曲线的平滑度（二阶差分）
    T = size(F, 3);
    
    if T < 3
        smoothness = 0;
        return;
    end
    
    total_norm = 0;
    for t = 2:T-1
        % 近似二阶导数
        diff = F(:, :, t+1) - 2*F(:, :, t) + F(:, :, t-1);
        total_norm = total_norm + norm(diff, 'fro')^2;
    end
    
    smoothness = total_norm / (T-2);
end

function tv = compute_total_variation(F)
    % 计算曲线的总变差
    T = size(F, 3);
    
    if T < 2
        tv = 0;
        return;
    end
    
    total_variation = 0;
    for t = 1:T-1
        for u = 1:size(F, 2)
            d = stiefel_chordal(F(:, u, t), F(:, u, t+1));
            total_variation = total_variation + d;
        end
    end
    
    tv = total_variation;
end