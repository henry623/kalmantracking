function simSettings = init(SNR_dB)
% OFDM接收机仿真参数初始化函数
%
% 功能描述:
%   初始化OFDM接收机仿真所需的所有参数
%   包括信号参数、OFDM参数、跟踪参数等
%
% 输入参数:
%   SNR_dB - 信噪比 (dB)，可选参数，默认为 -10 dB
%
% 输出参数:
%   simSettings - 仿真设置参数结构体
%
% 参数修改总结:
%   1. Nu: 512 -> 400 (留出保护子载波，避免频谱泄漏)
%   2. nSymbol: 10 -> 20 (提高处理增益，改善低信噪比性能)
%   3. CP_length: 128 -> 256 (增强抗多径能力)
%   4. famp: 4 -> 2 (降低计算复杂度，2倍插值已足够)
%   5. extra: 10 -> 20 (增强边界效应处理)
%   6. search_range: 1000 -> 2048 (覆盖完整信号周期)
%   7. search_step: 1 -> 2 (平衡精度和速度)
%
% 作者: OFDM接收机开发
% 日期: 2025年7月
% 版本: 1.0

%% 输入参数处理
if nargin < 1
    SNR_dB = -10;  % 默认信噪比
end

%% 基本仿真参数
simSettings.SNR = SNR_dB;                    % 信噪比 (dB)
simSettings.fs = 2e6;                        % 采样频率 (Hz)
simSettings.fc = 1e9;                        % 载波频率 (Hz)
simSettings.T = 1e-3;                        % 符号周期 (s)

%% OFDM参数
simSettings.NFFT = 1024;                     % FFT点数
simSettings.Nu = 400;                        % 每个频带的有效子载波数
                                            % 修改记录: 初始值=512, 修改值=400
                                            % 修改原因: 留出保护子载波，避免频谱泄漏和邻道干扰
simSettings.numBand = 2;                     % 频带数量
simSettings.nSymbol = 20;                    % OFDM符号数
                                            % 修改记录: 初始值=10, 修改值=20
                                            % 修改原因: 增加符号数提高处理增益，改善低信噪比下的性能
simSettings.CP_length = 256;                 % 循环前缀长度
                                            % 修改记录: 初始值=128, 修改值=256
                                            % 修改原因: 设为1/4 NFFT，增强抗多径能力和符号间干扰抑制

%% 信号处理参数
simSettings.famp = 2;                        % 插值倍数
                                            % 修改记录: 初始值=4, 修改值=2
                                            % 修改原因: 降低计算复杂度，对于OFDM系统2倍插值已足够
simSettings.extra = 20;                      % 保护间隔长度
                                            % 修改记录: 初始值=10, 修改值=20
                                            % 修改原因: 增加保护间隔以更好地处理边界效应和时延扩展
simSettings.codeLength = 10230;              % 码长度

%% 跟踪环路参数
simSettings.dll_bandwidth = 1.0;             % DLL噪声带宽 (Hz)
simSettings.pll_bandwidth = 10.0;            % PLL噪声带宽 (Hz)
simSettings.dll_damping = 0.707;             % DLL阻尼比
simSettings.pll_damping = 0.707;             % PLL阻尼比

%% 卡尔曼滤波器参数
simSettings.kalman_Q = 1e-6;                 % 过程噪声协方差
simSettings.kalman_R = 1e-3;                 % 观测噪声协方差
simSettings.kalman_P0 = 1e-2;                % 初始状态协方差

%% 搜索参数
simSettings.search_range = 2048;             % 搜索范围 (samples)
                                            % 修改记录: 初始值=1000, 修改值=2048
                                            % 修改原因: 设为约2个OFDM符号长度，确保能覆盖完整的信号周期
simSettings.search_step = 2;                 % 搜索步长 (samples)
                                            % 修改记录: 初始值=1, 修改值=2
                                            % 修改原因: 平衡搜索精度和计算速度，减少50%的计算量

%% 信道参数
simSettings.doppler_freq = 100;              % 多普勒频率 (Hz)
simSettings.multipath_delay = [0, 1e-6, 2e-6]; % 多径延迟 (s)
simSettings.multipath_gain = [1, 0.5, 0.3]; % 多径增益

%% 噪声参数
simSettings.noise_type = 'awgn';             % 噪声类型
simSettings.noise_seed = 12345;              % 噪声种子

%% 接收机参数
simSettings.acquisition_threshold = 0.3;     % 捕获门限
simSettings.tracking_threshold = 0.1;        % 跟踪门限
simSettings.lock_detector_threshold = 0.8;   % 锁定检测门限

%% 数据记录参数
simSettings.record_data = true;              % 是否记录数据
simSettings.plot_results = true;             % 是否绘制结果
simSettings.save_results = false;            % 是否保存结果

%% 调试参数
simSettings.debug_mode = false;              % 调试模式
simSettings.verbose = true;                  % 详细输出

%% 根据信噪比调整参数
if SNR_dB < -20
    % 极低信噪比环境
    simSettings.dll_bandwidth = 0.5;
    simSettings.pll_bandwidth = 5.0;
    simSettings.acquisition_threshold = 0.2;
    simSettings.tracking_threshold = 0.05;
    simSettings.kalman_Q = 1e-7;
    simSettings.kalman_R = 1e-2;
    fprintf('初始化: 极低信噪比模式 (SNR = %.1f dB)\n', SNR_dB);
    
elseif SNR_dB < -10
    % 低信噪比环境
    simSettings.dll_bandwidth = 0.8;
    simSettings.pll_bandwidth = 8.0;
    simSettings.acquisition_threshold = 0.25;
    simSettings.tracking_threshold = 0.08;
    simSettings.kalman_Q = 1e-6;
    simSettings.kalman_R = 5e-3;
    fprintf('初始化: 低信噪比模式 (SNR = %.1f dB)\n', SNR_dB);
    
elseif SNR_dB < 0
    % 中等信噪比环境
    simSettings.dll_bandwidth = 1.0;
    simSettings.pll_bandwidth = 10.0;
    simSettings.acquisition_threshold = 0.3;
    simSettings.tracking_threshold = 0.1;
    fprintf('初始化: 中等信噪比模式 (SNR = %.1f dB)\n', SNR_dB);
    
elseif SNR_dB < 10
    % 较高信噪比环境
    simSettings.dll_bandwidth = 1.5;
    simSettings.pll_bandwidth = 15.0;
    simSettings.acquisition_threshold = 0.4;
    simSettings.tracking_threshold = 0.15;
    simSettings.kalman_Q = 1e-5;
    simSettings.kalman_R = 1e-4;
    fprintf('初始化: 较高信噪比模式 (SNR = %.1f dB)\n', SNR_dB);
    
else
    % 高信噪比环境
    simSettings.dll_bandwidth = 2.0;
    simSettings.pll_bandwidth = 20.0;
    simSettings.acquisition_threshold = 0.5;
    simSettings.tracking_threshold = 0.2;
    simSettings.kalman_Q = 1e-4;
    simSettings.kalman_R = 1e-5;
    fprintf('初始化: 高信噪比模式 (SNR = %.1f dB)\n', SNR_dB);
end

%% 计算派生参数
% 采样间隔
simSettings.Ts = 1 / simSettings.fs;

% 子载波间隔
simSettings.subcarrier_spacing = simSettings.fs / simSettings.NFFT;

% OFDM符号持续时间
simSettings.symbol_duration = simSettings.NFFT / simSettings.fs;

% 循环前缀持续时间
simSettings.CP_duration = simSettings.CP_length / simSettings.fs;

% 总符号持续时间（包括CP）
simSettings.total_symbol_duration = simSettings.symbol_duration + simSettings.CP_duration;

% 信号总长度
simSettings.signal_length = simSettings.nSymbol * (simSettings.NFFT + simSettings.CP_length);

% 插值后信号长度
simSettings.interpolated_length = simSettings.signal_length * simSettings.famp;

% 最终信号长度（包括保护间隔）
simSettings.final_length = simSettings.interpolated_length + 2 * simSettings.extra * simSettings.famp;

%% 频率相关参数
% 多普勒频率对应的相位变化率
simSettings.doppler_phase_rate = 2 * pi * simSettings.doppler_freq;

% 载波频率对应的相位变化率
simSettings.carrier_phase_rate = 2 * pi * simSettings.fc;

%% 时间向量
simSettings.time_vector = (0:simSettings.final_length-1) * simSettings.Ts / simSettings.famp;

%% 验证参数合理性
validate_parameters(simSettings);

%% 显示初始化信息
if simSettings.verbose
    display_initialization_info(simSettings);
end

end

function validate_parameters(simSettings)
% 验证参数的合理性

% 检查OFDM参数
if simSettings.Nu > simSettings.NFFT
    error('init: 有效子载波数不能超过FFT点数');
end

if simSettings.numBand <= 0
    error('init: 频带数量必须为正数');
end

if simSettings.nSymbol <= 0
    error('init: OFDM符号数必须为正数');
end

% 检查插值参数
if simSettings.famp <= 0 || mod(simSettings.famp, 1) ~= 0
    error('init: 插值倍数必须为正整数');
end

% 检查跟踪参数
if simSettings.dll_bandwidth <= 0 || simSettings.pll_bandwidth <= 0
    error('init: 环路带宽必须为正数');
end

if simSettings.dll_damping <= 0 || simSettings.pll_damping <= 0
    error('init: 阻尼比必须为正数');
end

% 检查卡尔曼滤波器参数
if simSettings.kalman_Q <= 0 || simSettings.kalman_R <= 0 || simSettings.kalman_P0 <= 0
    error('init: 卡尔曼滤波器协方差参数必须为正数');
end

% 检查门限参数
if simSettings.acquisition_threshold <= 0 || simSettings.acquisition_threshold > 1
    warning('init: 捕获门限建议设置在0到1之间');
end

if simSettings.tracking_threshold <= 0 || simSettings.tracking_threshold > 1
    warning('init: 跟踪门限建议设置在0到1之间');
end

end

function display_initialization_info(simSettings)
% 显示初始化信息

fprintf('\n=== OFDM接收机参数初始化完成 ===\n');
fprintf('基本参数:\n');
fprintf('  - 信噪比: %.1f dB\n', simSettings.SNR);
fprintf('  - 采样频率: %.2f MHz\n', simSettings.fs / 1e6);
fprintf('  - 载波频率: %.2f GHz\n', simSettings.fc / 1e9);

fprintf('\nOFDM参数:\n');
fprintf('  - FFT点数: %d\n', simSettings.NFFT);
fprintf('  - 有效子载波数: %d\n', simSettings.Nu);
fprintf('  - 频带数量: %d\n', simSettings.numBand);
fprintf('  - OFDM符号数: %d\n', simSettings.nSymbol);
fprintf('  - 循环前缀长度: %d\n', simSettings.CP_length);

fprintf('\n信号处理参数:\n');
fprintf('  - 插值倍数: %d\n', simSettings.famp);
fprintf('  - 保护间隔: %d\n', simSettings.extra);
fprintf('  - 最终信号长度: %d samples\n', simSettings.final_length);

fprintf('\n跟踪环路参数:\n');
fprintf('  - DLL带宽: %.2f Hz\n', simSettings.dll_bandwidth);
fprintf('  - PLL带宽: %.2f Hz\n', simSettings.pll_bandwidth);
fprintf('  - DLL阻尼比: %.3f\n', simSettings.dll_damping);
fprintf('  - PLL阻尼比: %.3f\n', simSettings.pll_damping);

fprintf('\n门限参数:\n');
fprintf('  - 捕获门限: %.3f\n', simSettings.acquisition_threshold);
fprintf('  - 跟踪门限: %.3f\n', simSettings.tracking_threshold);
fprintf('  - 锁定检测门限: %.3f\n', simSettings.lock_detector_threshold);

fprintf('===============================\n\n');

end