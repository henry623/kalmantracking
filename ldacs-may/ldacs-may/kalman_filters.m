function [kf_output] = kalman_filters(kf_type, kf, z, dt)
    switch kf_type
        case 'standard'
            kf_output = standard_kf(kf, z, dt);
        case 'extended'
            kf_output = extended_kf(kf, z, dt);
        case 'unscented'
            kf_output = unscented_kf(kf, z, dt);
        otherwise
            error('Unknown Kalman filter type');
    end
end

function [kf_output] = standard_kf(kf, z, dt)
    % 预测步骤
    kf.F = [1 dt 0 0; 0 1 0 0; 0 0 1 dt; 0 0 0 1];
    kf.x = kf.F * kf.x;
    kf.P = kf.F * kf.P * kf.F' + kf.Q;
    
    % 更新步骤
    y = z - kf.H * kf.x;
    S = kf.H * kf.P * kf.H' + kf.R;
    K = kf.P * kf.H' / S;
    kf.x = kf.x + K * y;
    kf.P = (eye(4) - K * kf.H) * kf.P;
    
    kf_output = kf;
end

function [kf_output] = extended_kf(kf, z, dt)
    % 预测步骤
    kf.F = [1 dt 0 0; 0 1 0 0; 0 0 1 dt; 0 0 0 1];
    kf.x = kf.F * kf.x;
    kf.P = kf.F * kf.P * kf.F' + kf.Q;
    
    % 更新步骤
    h = @(x) [x(1); x(3)];  % 非线性测量函数
    H = [1 0 0 0; 0 0 1 0];  % 测量函数的雅可比矩阵
    y = z - h(kf.x);
    S = H * kf.P * H' + kf.R;
    K = kf.P * H' / S;
    kf.x = kf.x + K * y;
    kf.P = (eye(4) - K * H) * kf.P;
    
    kf_output = kf;
end

function [kf_output] = unscented_kf(kf, z, dt)
    n = 4;  % 状态维度
    alpha = 1e-3;
    kappa = 0;
    beta = 2;
    lambda = alpha^2 * (n + kappa) - n;
    
    % 计算权重
    Wm = [lambda / (n + lambda), repmat(1 / (2*(n + lambda)), 1, 2*n)];
    Wc = Wm;
    Wc(1) = Wc(1) + (1 - alpha^2 + beta);
    
    % 预测步骤
    kf.F = [1 dt 0 0; 0 1 0 0; 0 0 1 dt; 0 0 0 1];
    X = sqrtm((n + lambda) * kf.P);
    X = [zeros(n, 1), X, -X];
    X = kf.x + X;
    
    % 传播 sigma 点
    X = kf.F * X;
    
    % 计算预测均值和协方差
    kf.x = X * Wm';
    kf.P = zeros(n);
    for i = 1:2*n+1
        kf.P = kf.P + Wc(i) * (X(:,i) - kf.x) * (X(:,i) - kf.x)';
    end
    kf.P = kf.P + kf.Q;
    
    % 更新步骤
    Y = kf.H * X;
    y_pred = Y * Wm';
    Pyy = zeros(2);
    Pxy = zeros(n, 2);
    for i = 1:2*n+1
        Pyy = Pyy + Wc(i) * (Y(:,i) - y_pred) * (Y(:,i) - y_pred)';
        Pxy = Pxy + Wc(i) * (X(:,i) - kf.x) * (Y(:,i) - y_pred)';
    end
    Pyy = Pyy + kf.R;
    
    K = Pxy / Pyy;
    kf.x = kf.x + K * (z - y_pred);
    kf.P = kf.P - K * Pyy * K';
    
    kf_output = kf;
end
