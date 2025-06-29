function [kf_output] = kalman_filters2(kf_type, kf, z, dt, ideal_values)
    % 修改后的卡尔曼滤波器函数，状态变量为：
    % x = [码相位误差; 码频率误差; 载波相位误差; 载波频率误差]
    % z = [测量的码相位误差; 测量的载波相位误差]
    % ideal_values = [理想码相位; 理想码频率; 理想载波相位; 理想载波频率]
    
    switch kf_type
        case 'standard'
            kf_output = standard_kf(kf, z, dt);
        case 'extended'
            kf_output = extended_kf(kf, z, dt);
        case 'unscented'
            kf_output = unscented_kf(kf, z, dt);
        otherwise
            error('未知的卡尔曼滤波器类型');
    end
    
    % 将卡尔曼滤波器的输出（误差）与理想值相加，得到修正后的估计值
    kf_output.corrected_values = ideal_values + kf_output.x;
end

function [kf_output] = standard_kf(kf, z, dt)
    % 标准卡尔曼滤波器
    % 状态变量：[码相位误差; 码频率误差; 载波相位误差; 载波频率误差]
    
    % 预测步骤
    % 状态转移矩阵 F
    kf.F = [1 dt 0 0;   % 码相位误差 = 上一时刻码相位误差 + dt * 码频率误差
            0 1  0 0;   % 码频率误差 = 上一时刻码频率误差
            0 0  1 dt;  % 载波相位误差 = 上一时刻载波相位误差 + dt * 载波频率误差
            0 0  0 1];  % 载波频率误差 = 上一时刻载波频率误差
    
    % 预测状态
    kf.x = kf.F * kf.x;
    
    % 预测协方差
    kf.P = kf.F * kf.P * kf.F' + kf.Q;
    
    % 更新步骤
    % 计算卡尔曼增益
    y = z - kf.H * kf.x;  % 测量残差
    S = kf.H * kf.P * kf.H' + kf.R;  % 残差协方差
    K = kf.P * kf.H' / S;  % 卡尔曼增益
    
    % 更新状态和协方差
    kf.x = kf.x + K * y;
    kf.P = (eye(4) - K * kf.H) * kf.P;
    
    kf_output = kf;
end

function [kf_output] = extended_kf(kf, z, dt)
    % 扩展卡尔曼滤波器
    % 状态变量：[码相位误差; 码频率误差; 载波相位误差; 载波频率误差]
    
    % 预测步骤
    % 非线性状态转移函数
    f = @(x) [x(1) + dt * x(2);  % 码相位误差
              x(2);               % 码频率误差
              x(3) + dt * x(4);  % 载波相位误差
              x(4)];              % 载波频率误差
    
    % 状态转移函数的雅可比矩阵
    F = [1 dt 0 0;
         0 1  0 0;
         0 0  1 dt;
         0 0  0 1];
    
    % 预测状态
    kf.x = f(kf.x);
    
    % 预测协方差
    kf.P = F * kf.P * F' + kf.Q;
    
    % 更新步骤
    % 非线性测量函数
    h = @(x) [x(1); x(3)];  % 测量码相位误差和载波相位误差
    
    % 测量函数的雅可比矩阵
    H = [1 0 0 0;
         0 0 1 0];
    
    % 计算卡尔曼增益
    y = z - h(kf.x);  % 测量残差
    S = H * kf.P * H' + kf.R;  % 残差协方差
    K = kf.P * H' / S;  % 卡尔曼增益
    
    % 更新状态和协方差
    kf.x = kf.x + K * y;
    kf.P = (eye(4) - K * H) * kf.P;
    
    kf_output = kf;
end

function [kf_output] = unscented_kf(kf, z, dt)
    % 无迹卡尔曼滤波器
    % 状态变量：[码相位误差; 码频率误差; 载波相位误差; 载波频率误差]
    
    n = 4;  % 状态维度
    alpha = 1e-3;  % 控制 sigma 点的分布
    kappa = 0;  % 次要缩放参数
    beta = 2;  % 用于合并高阶项的参数（高斯分布的最优值为2）
    lambda = alpha^2 * (n + kappa) - n;
    
    % 计算权重
    Wm = [lambda / (n + lambda), repmat(1 / (2*(n + lambda)), 1, 2*n)];  % 均值权重
    Wc = Wm;  % 协方差权重
    Wc(1) = Wc(1) + (1 - alpha^2 + beta);  % 调整第一个协方差权重
    
    % 预测步骤
    % 非线性状态转移函数
    f = @(x) [x(1) + dt * x(2);  % 码相位误差
              x(2);               % 码频率误差
              x(3) + dt * x(4);  % 载波相位误差
              x(4)];              % 载波频率误差
    
    % 生成 sigma 点
    P_sqrt = sqrtm((n + lambda) * kf.P);  % 矩阵平方根
    X = [zeros(n, 1), P_sqrt, -P_sqrt];  % 2n+1 个 sigma 点
    X = kf.x + X;  % 中心化 sigma 点
    
    % 传播 sigma 点
    for i = 1:2*n+1
        X(:,i) = f(X(:,i));
    end
    
    % 计算预测均值
    kf.x = zeros(n, 1);
    for i = 1:2*n+1
        kf.x = kf.x + Wm(i) * X(:,i);
    end
    
    % 计算预测协方差
    kf.P = zeros(n);
    for i = 1:2*n+1
        kf.P = kf.P + Wc(i) * (X(:,i) - kf.x) * (X(:,i) - kf.x)';
    end
    kf.P = kf.P + kf.Q;
    
    % 更新步骤
    % 非线性测量函数
    h = @(x) [x(1); x(3)];  % 测量码相位误差和载波相位误差
    
    % 计算测量预测
    Y = zeros(2, 2*n+1);
    for i = 1:2*n+1
        Y(:,i) = h(X(:,i));
    end
    
    % 计算测量预测均值
    y_pred = zeros(2, 1);
    for i = 1:2*n+1
        y_pred = y_pred + Wm(i) * Y(:,i);
    end
    
    % 计算测量预测协方差和交叉协方差
    Pyy = zeros(2);
    Pxy = zeros(n, 2);
    for i = 1:2*n+1
        Pyy = Pyy + Wc(i) * (Y(:,i) - y_pred) * (Y(:,i) - y_pred)';
        Pxy = Pxy + Wc(i) * (X(:,i) - kf.x) * (Y(:,i) - y_pred)';
    end
    Pyy = Pyy + kf.R;
    
    % 计算卡尔曼增益
    K = Pxy / Pyy;
    
    % 更新状态和协方差
    kf.x = kf.x + K * (z - y_pred);
    kf.P = kf.P - K * Pyy * K';
    
    kf_output = kf;
end
