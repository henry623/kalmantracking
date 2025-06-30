function [kf_output] = kalman_filters(kf_type, kf, z, dt)
    % kalman_filters - 实现不同类型的卡尔曼滤波器
    %
    % 输入:
    %   kf_type - 卡尔曼滤波器类型: 'standard', 'extended', 或 'unscented'
    %   kf      - 卡尔曼滤波器状态结构体
    %   z       - 测量值 [码相位; 载波相位]
    %   dt      - 时间步长
    %
    % 输出:
    %   kf_output - 更新后的卡尔曼滤波器状态
    %
    % 状态向量 x 定义:
    %   x(1) - 载波频率误差 - 载波频率的估计误差
    %   x(2) - 载波相位误差 - 载波相位的估计误差
    %   x(3) - 码频率误差 - 码频率的估计误差
    %   x(4) - 码相位误差 - 码相位的估计误差
    %
    % 卡尔曼滤波器的作用是通过结合自己的预测对跟踪环路的输出结果进行纠正
    % 它通过估计上述四个参数的误差，并利用这些误差估计来校正跟踪环路的输出，
    % 从而提高跟踪精度和稳定性，特别是在低信噪比环境下。
    
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
    % standard_kf - 标准卡尔曼滤波器实现
    %
    % 状态向量 x = [载波频率误差; 载波相位误差; 码频率误差; 码相位误差]
    % 通过预测和更新步骤，估计跟踪环路中的误差，用于后续校正
    
    % 预测步骤 - 状态传播
    % 状态转移矩阵 F 将当前状态映射到下一个状态
    % 载波相位误差 += 载波频率误差 * dt
    % 码相位误差 += 码频率误差 * dt
    kf.F = [1  0  0  0;   % 载波频率误差保持不变
            dt 1  0  0;   % 载波相位误差 = 载波相位误差 + 载波频率误差*dt
            0  0  1  0;   % 码频率误差保持不变
            0  0  dt 1];  % 码相位误差 = 码相位误差 + 码频率误差*dt
            
    % 预测状态
    kf.x = kf.F * kf.x;
    
    % 预测协方差
    kf.P = kf.F * kf.P * kf.F' + kf.Q;
    
    % 更新步骤 - 使用测量值校正状态
    % 计算测量残差 (测量值 - 预测值)
    y = z - kf.H * kf.x;
    
    % 计算残差协方差
    S = kf.H * kf.P * kf.H' + kf.R;
    
    % 计算卡尔曼增益
    K = kf.P * kf.H' / S;
    
    % 更新状态估计
    kf.x = kf.x + K * y;
    
    % 更新状态协方差
    kf.P = (eye(size(kf.P)) - K * kf.H) * kf.P;
    
    kf_output = kf;
end

function [kf_output] = extended_kf(kf, z, dt)
    % extended_kf - 扩展卡尔曼滤波器实现
    %
    % 状态向量 x = [载波频率误差; 载波相位误差; 码频率误差; 码相位误差]
    % 处理非线性测量模型，更适合处理相位和频率之间的非线性关系
    
    % 预测步骤 - 状态传播
    kf.F = [1  0  0  0;   % 载波频率误差保持不变
            dt 1  0  0;   % 载波相位误差 = 载波相位误差 + 载波频率误差*dt
            0  0  1  0;   % 码频率误差保持不变
            0  0  dt 1];  % 码相位误差 = 码相位误差 + 码频率误差*dt
            
    % 预测状态 - 可以使用非线性状态转移函数
    % 这里使用线性状态转移作为简化
    kf.x = kf.F * kf.x;
    
    % 预测协方差
    kf.P = kf.F * kf.P * kf.F' + kf.Q;
    
    % 更新步骤 - 使用非线性测量模型
    % 定义非线性测量函数和其雅可比矩阵
    % 这里的非线性函数可以捕捉相位和频率之间的复杂关系
    h = @(x) [x(2) + 0.1*sin(x(1)); x(4) + 0.1*sin(x(3))];  % 非线性测量函数示例
    
    % 计算测量函数在当前状态估计处的雅可比矩阵
    H = [0.1*cos(kf.x(1)) 1 0 0;
         0 0 0.1*cos(kf.x(3)) 1];
    
    % 计算测量残差
    y = z - h(kf.x);
    
    % 计算残差协方差
    S = H * kf.P * H' + kf.R;
    
    % 计算卡尔曼增益
    K = kf.P * H' / S;
    
    % 更新状态估计
    kf.x = kf.x + K * y;
    
    % 更新状态协方差
    kf.P = (eye(size(kf.P)) - K * H) * kf.P;
    
    kf_output = kf;
end

function [kf_output] = unscented_kf(kf, z, dt)
    % unscented_kf - 无迹卡尔曼滤波器实现
    %
    % 状态向量 x = [载波频率误差; 载波相位误差; 码频率误差; 码相位误差]
    % 使用sigma点来处理非线性，无需计算雅可比矩阵，更适合处理强非线性系统
    
    n = size(kf.x, 1);  % 状态维度
    alpha = 1e-3;  % 控制sigma点分布的参数
    kappa = 0;     % 次要缩放参数
    beta = 2;      % 包含先验分布信息的参数
    lambda = alpha^2 * (n + kappa) - n;
    
    % 计算权重
    Wm = [lambda / (n + lambda), repmat(1 / (2*(n + lambda)), 1, 2*n)];  % 均值权重
    Wc = Wm;  % 协方差权重
    Wc(1) = Wc(1) + (1 - alpha^2 + beta);  % 调整协方差权重
    
    % 预测步骤 - 使用sigma点
    % 状态转移矩阵
    kf.F = [1  0  0  0;   % 载波频率误差保持不变
            dt 1  0  0;   % 载波相位误差 = 载波相位误差 + 载波频率误差*dt
            0  0  1  0;   % 码频率误差保持不变
            0  0  dt 1];  % 码相位误差 = 码相位误差 + 码频率误差*dt
    
    % 生成sigma点
    try
        X_sqrt = chol((n + lambda) * kf.P, 'lower');  % 尝试使用Cholesky分解
    catch
        % 如果Cholesky分解失败，使用SVD分解
        [U, S, ~] = svd((n + lambda) * kf.P);
        X_sqrt = U * sqrt(S);
    end
    
    X = [zeros(n, 1), X_sqrt, -X_sqrt];  % 2n+1个sigma点
    X = repmat(kf.x, 1, 2*n+1) + X;  % 中心化sigma点
    
    % 定义非线性状态转移函数
    % 这里使用线性状态转移作为简化
    f = @(x) kf.F * x;
    
    % 传播sigma点
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
    kf.P = kf.P + kf.Q;  % 添加过程噪声
    
    % 定义非线性测量函数
    % 这里使用一个简单的非线性函数作为示例
    h = @(x) [x(2) + 0.1*sin(x(1)); x(4) + 0.1*sin(x(3))];
    
    % 将sigma点映射到测量空间
    Y = zeros(size(z, 1), 2*n+1);
    for i = 1:2*n+1
        Y(:,i) = h(X(:,i));
    end
    
    % 计算预测测量
    y_pred = zeros(size(z));
    for i = 1:2*n+1
        y_pred = y_pred + Wm(i) * Y(:,i);
    end
    
    % 计算测量协方差和交叉协方差
    Pyy = zeros(size(z, 1));  % 测量协方差
    Pxy = zeros(n, size(z, 1));  % 状态-测量交叉协方差
    for i = 1:2*n+1
        Pyy = Pyy + Wc(i) * (Y(:,i) - y_pred) * (Y(:,i) - y_pred)';
        Pxy = Pxy + Wc(i) * (X(:,i) - kf.x) * (Y(:,i) - y_pred)';
    end
    Pyy = Pyy + kf.R;  % 添加测量噪声
    
    % 计算卡尔曼增益
    K = Pxy / Pyy;
    
    % 更新状态估计
    kf.x = kf.x + K * (z - y_pred);
    
    % 更新状态协方差
    kf.P = kf.P - K * Pyy * K';
    
    kf_output = kf;
end
