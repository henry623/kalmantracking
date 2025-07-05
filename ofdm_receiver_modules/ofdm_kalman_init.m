function kf_params = ofdm_kalman_init(simSettings)
% OFDM卡尔曼滤波器初始化模块
%
% 功能描述:
%   为OFDM接收机初始化卡尔曼滤波器参数
%   用于跟踪载波频率偏移、相位偏移和码延迟等状态
%
% 输入参数:
%   simSettings - 仿真设置参数结构体
%
% 输出参数:
%   kf_params - 卡尔曼滤波器参数结构体，包含：
%              * A: 状态转移矩阵
%              * H: 观测矩阵
%              * Q: 过程噪声协方差矩阵
%              * R: 观测噪声协方差矩阵
%              * P: 状态协方差矩阵
%              * x: 初始状态向量
%              * state_names: 状态变量名称
%
% 状态向量定义:
%   x = [code_delay; code_rate; carrier_freq; carrier_phase]
%   - code_delay: 码延迟 (samples)
%   - code_rate: 码速率偏差 (samples/sample)
%   - carrier_freq: 载波频率偏移 (Hz)
%   - carrier_phase: 载波相位偏移 (rad)
%
% 作者: OFDM接收机开发
% 日期: 2025年7月
% 版本: 1.0

%% 输入参数验证
if nargin < 1
    error('ofdm_kalman_init: 必须提供仿真设置参数');
end

% 检查必要字段
required_fields = {'kalman_Q', 'kalman_R', 'kalman_P0', 'fs'};
for i = 1:length(required_fields)
    if ~isfield(simSettings, required_fields{i})
        error('ofdm_kalman_init: simSettings缺少必要字段 %s', required_fields{i});
    end
end

fprintf('  - 初始化卡尔曼滤波器参数...\n');

%% 状态向量维度定义
n_states = 4;  % 状态数量
n_obs = 2;     % 观测数量 (DLL和PLL鉴别器输出)

%% 状态变量名称
kf_params.state_names = {'code_delay', 'code_rate', 'carrier_freq', 'carrier_phase'};
kf_params.obs_names = {'dll_discriminator', 'pll_discriminator'};

fprintf('    * 状态向量维度: %d\n', n_states);
fprintf('    * 观测向量维度: %d\n', n_obs);

%% 状态转移矩阵 A (4x4)
% 基于离散时间模型，采样间隔为 Ts
Ts = 1 / simSettings.fs;  % 采样间隔

% 状态转移矩阵
% x(k+1) = A * x(k) + w(k)
kf_params.A = [
    1,  Ts,  0,   0;    % code_delay(k+1) = code_delay(k) + Ts*code_rate(k)
    0,  1,   0,   0;    % code_rate(k+1) = code_rate(k)
    0,  0,   1,   Ts;   % carrier_freq(k+1) = carrier_freq(k) + Ts*carrier_phase_rate(k)
    0,  0,   0,   1     % carrier_phase(k+1) = carrier_phase(k)
];

fprintf('    * 状态转移矩阵: 4x4\n');

%% 观测矩阵 H (2x4)
% 观测方程: z(k) = H * x(k) + v(k)
% DLL鉴别器主要观测码延迟
% PLL鉴别器主要观测载波相位
kf_params.H = [
    1,  0,  0,  0;      % DLL鉴别器输出 ≈ code_delay
    0,  0,  0,  1       % PLL鉴别器输出 ≈ carrier_phase
];

fprintf('    * 观测矩阵: 2x4\n');

%% 过程噪声协方差矩阵 Q (4x4)
% 根据OFDM系统特性设置过程噪声
base_Q = simSettings.kalman_Q;

% 码延迟过程噪声 - OFDM符号边界相对稳定
code_delay_var = base_Q * 0.1;
code_rate_var = base_Q * 0.01;  % 码速率变化很慢

% 载波频率过程噪声 - 考虑振荡器稳定性
carrier_freq_var = base_Q * 10;  % 频率偏移变化
carrier_phase_var = base_Q * 100; % 相位噪声较大

kf_params.Q = diag([code_delay_var, code_rate_var, carrier_freq_var, carrier_phase_var]);

fprintf('    * 过程噪声协方差: 对角矩阵\n');
fprintf('      - 码延迟方差: %.2e\n', code_delay_var);
fprintf('      - 码速率方差: %.2e\n', code_rate_var);
fprintf('      - 载波频率方差: %.2e\n', carrier_freq_var);
fprintf('      - 载波相位方差: %.2e\n', carrier_phase_var);

%% 观测噪声协方差矩阵 R (2x2)
% 根据鉴别器性能设置观测噪声
base_R = simSettings.kalman_R;

% DLL鉴别器噪声 - 与信噪比相关
dll_noise_var = base_R;

% PLL鉴别器噪声 - 通常比DLL噪声大
pll_noise_var = base_R * 5;

kf_params.R = diag([dll_noise_var, pll_noise_var]);

fprintf('    * 观测噪声协方差: 对角矩阵\n');
fprintf('      - DLL鉴别器噪声方差: %.2e\n', dll_noise_var);
fprintf('      - PLL鉴别器噪声方差: %.2e\n', pll_noise_var);

%% 初始状态协方差矩阵 P (4x4)
base_P = simSettings.kalman_P0;

% 初始状态不确定性
code_delay_init_var = base_P * 100;    % 码延迟初始不确定性较大
code_rate_init_var = base_P * 0.1;     % 码速率初始不确定性较小
carrier_freq_init_var = base_P * 1000; % 载波频率初始不确定性大
carrier_phase_init_var = base_P * 10;  % 载波相位初始不确定性中等

kf_params.P = diag([code_delay_init_var, code_rate_init_var, ...
                   carrier_freq_init_var, carrier_phase_init_var]);

fprintf('    * 初始状态协方差: 对角矩阵\n');

%% 初始状态向量 x (4x1)
% 初始状态估计
kf_params.x = [
    0;      % 初始码延迟 (samples)
    0;      % 初始码速率偏差
    0;      % 初始载波频率偏移 (Hz)
    0       % 初始载波相位偏移 (rad)
];

fprintf('    * 初始状态向量: 零初始化\n');

%% 滤波器配置参数
kf_params.config.n_states = n_states;
kf_params.config.n_obs = n_obs;
kf_params.config.sampling_rate = simSettings.fs;
kf_params.config.update_rate = 1;  % 每个样本更新一次

%% 自适应参数
% 根据信噪比调整滤波器参数
if isfield(simSettings, 'SNR')
    SNR_dB = simSettings.SNR;
    
    if SNR_dB < -15
        % 低信噪比：增加过程噪声，减少观测权重
        kf_params.Q = kf_params.Q * 2;
        kf_params.R = kf_params.R * 0.5;
        fprintf('    * 低信噪比调整: 增加过程噪声，减少观测噪声\n');
        
    elseif SNR_dB > 5
        % 高信噪比：减少过程噪声，增加观测权重
        kf_params.Q = kf_params.Q * 0.5;
        kf_params.R = kf_params.R * 2;
        fprintf('    * 高信噪比调整: 减少过程噪声，增加观测噪声\n');
    end
end

%% 滤波器状态记录
kf_params.history.enabled = false;  % 默认不记录历史
kf_params.history.max_length = 1000;
kf_params.history.states = [];
kf_params.history.covariances = [];
kf_params.history.innovations = [];

%% 性能监控参数
kf_params.monitor.innovation_threshold = 3.0;  % 新息检验门限
kf_params.monitor.divergence_threshold = 10.0; % 发散检测门限
kf_params.monitor.reset_threshold = 100.0;     % 重置门限

%% 验证矩阵维度
validate_matrices(kf_params);

fprintf('    * 卡尔曼滤波器初始化完成\n');

end

function validate_matrices(kf_params)
% 验证卡尔曼滤波器矩阵维度的一致性

n_states = kf_params.config.n_states;
n_obs = kf_params.config.n_obs;

% 检查状态转移矩阵
if ~isequal(size(kf_params.A), [n_states, n_states])
    error('ofdm_kalman_init: 状态转移矩阵A维度错误');
end

% 检查观测矩阵
if ~isequal(size(kf_params.H), [n_obs, n_states])
    error('ofdm_kalman_init: 观测矩阵H维度错误');
end

% 检查过程噪声协方差矩阵
if ~isequal(size(kf_params.Q), [n_states, n_states])
    error('ofdm_kalman_init: 过程噪声协方差矩阵Q维度错误');
end

% 检查观测噪声协方差矩阵
if ~isequal(size(kf_params.R), [n_obs, n_obs])
    error('ofdm_kalman_init: 观测噪声协方差矩阵R维度错误');
end

% 检查状态协方差矩阵
if ~isequal(size(kf_params.P), [n_states, n_states])
    error('ofdm_kalman_init: 状态协方差矩阵P维度错误');
end

% 检查初始状态向量
if ~isequal(size(kf_params.x), [n_states, 1])
    error('ofdm_kalman_init: 初始状态向量x维度错误');
end

% 检查矩阵正定性
if any(eig(kf_params.Q) <= 0)
    warning('ofdm_kalman_init: 过程噪声协方差矩阵Q不是正定的');
end

if any(eig(kf_params.R) <= 0)
    warning('ofdm_kalman_init: 观测噪声协方差矩阵R不是正定的');
end

if any(eig(kf_params.P) <= 0)
    warning('ofdm_kalman_init: 初始状态协方差矩阵P不是正定的');
end

end