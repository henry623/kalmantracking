function [tracking_params] = ofdm_tracking_init(simSettings)
% OFDM跟踪环路初始化模块
%
% 功能描述:
%   初始化OFDM接收机的跟踪环路参数，包括DLL（码跟踪环路）和PLL（载波跟踪环路）
%   的各种参数设置，为后续的信号跟踪处理做准备
%
% 输入参数:
%   simSettings - 仿真设置参数结构体，包含：
%                * fp: 码频率 (Hz)
%                * fs: 采样频率 (Hz)
%                * fi: 中频频率 (Hz)
%                * NFFT: FFT点数
%                * nSymbol: OFDM符号数
%                * dt: 时间间隔 (s)
%                * t_total: 总仿真时间 (s)
%
% 输出参数:
%   tracking_params - 跟踪参数结构体，包含：
%                    * DLL相关参数（码跟踪环路）
%                    * PLL相关参数（载波跟踪环路）
%                    * 环路滤波器参数
%                    * 初始状态变量
%                    * 性能监控参数
%
% 处理流程:
%   1. 输入参数验证
%   2. DLL参数计算和设置
%   3. PLL参数计算和设置
%   4. 环路滤波器系数计算
%   5. 初始状态变量设置
%   6. 性能监控参数设置
%   7. 参数合理性检查
%
% 技术细节:
%   - DLL使用早晚相关器结构进行码相位跟踪
%   - PLL使用Costas环路进行载波相位跟踪
%   - 支持FLL辅助的PLL结构提高动态性能
%   - 环路参数根据信号特性和噪声环境自适应调整
%
% 作者: OFDM接收机开发团队
% 日期: 2025年1月
% 版本: 1.0

%% 1. 输入参数验证
if nargin < 1
    error('ofdm_tracking_init: 必须提供仿真设置参数');
end

% 检查必要字段
required_fields = {'fp', 'fs', 'fi', 'NFFT', 'nSymbol', 'dt', 't_total'};
for i = 1:length(required_fields)
    if ~isfield(simSettings, required_fields{i})
        error('ofdm_tracking_init: simSettings缺少必要字段 %s', required_fields{i});
    end
end

fprintf('  - 输入参数验证通过\n');

% 提取基本参数
fp = simSettings.fp;           % 码频率
fs = simSettings.fs;           % 采样频率
fi = simSettings.fi;           % 中频频率
NFFT = simSettings.NFFT;       % FFT点数
nSymbol = simSettings.nSymbol; % OFDM符号数
dt = simSettings.dt;           % 时间间隔
t_total = simSettings.t_total; % 总仿真时间

% 计算基本时间参数
lenOFDM = NFFT * nSymbol;      % OFDM符号长度
PDI = lenOFDM / fs;            % 积分时间 (Predetection Integration time)

fprintf('    * 码频率: %.2f kHz\n', fp/1e3);
fprintf('    * 采样频率: %.2f MHz\n', fs/1e6);
fprintf('    * 积分时间: %.4f ms\n', PDI*1e3);

%% 2. DLL (Delay Lock Loop) 参数设置
fprintf('  - 设置DLL（码跟踪环路）参数...\n');

% DLL基本参数
dll_params = struct();

% 早晚相关器间距 (chips)
dll_params.earlyLateSpc = 0.5;  % 标准设置为0.5个码片

% DLL噪声带宽 (Hz) - 根据动态环境调整
if isfield(simSettings, 'SNR') && simSettings.SNR < -10
    % 低信噪比环境，使用较窄的噪声带宽
    dll_params.noiseBandwidth = 0.5;
elseif isfield(simSettings, 'SNR') && simSettings.SNR > 10
    % 高信噪比环境，可以使用较宽的噪声带宽
    dll_params.noiseBandwidth = 2.0;
else
    % 标准环境
    dll_params.noiseBandwidth = 1.0;
end

% DLL阻尼比 - 影响环路的稳定性和响应速度
dll_params.dampingRatio = 0.707;  % 临界阻尼，平衡稳定性和响应速度

% DLL环路增益 - 影响环路的跟踪精度
dll_params.loopGain = 1.0;

% 积分时间
dll_params.PDI = PDI;

% 计算DLL环路滤波器系数
fprintf('    * 计算DLL环路滤波器系数...\n');
[dll_params.tau1, dll_params.tau2] = calLoopCoef(dll_params.noiseBandwidth, ...
                                                  dll_params.dampingRatio, ...
                                                  dll_params.loopGain);

% DLL性能参数
dll_params.wn = dll_params.noiseBandwidth / 0.7845;  % 自然频率
dll_params.omegaOf = dll_params.noiseBandwidth / 0.53;  % 开环频率

fprintf('    * DLL噪声带宽: %.2f Hz\n', dll_params.noiseBandwidth);
fprintf('    * DLL阻尼比: %.3f\n', dll_params.dampingRatio);
fprintf('    * DLL tau1: %.6f, tau2: %.6f\n', dll_params.tau1, dll_params.tau2);

%% 3. PLL (Phase Lock Loop) 参数设置
fprintf('  - 设置PLL（载波跟踪环路）参数...\n');

% PLL基本参数
pll_params = struct();

% PLL噪声带宽 (Hz) - 通常比DLL宽一些
if isfield(simSettings, 'SNR') && simSettings.SNR < -10
    % 低信噪比环境
    pll_params.noiseBandwidth = 2.0;
elseif isfield(simSettings, 'SNR') && simSettings.SNR > 10
    % 高信噪比环境
    pll_params.noiseBandwidth = 8.0;
else
    % 标准环境
    pll_params.noiseBandwidth = 4.0;
end

% PLL阻尼比
pll_params.dampingRatio = 0.707;

% PLL环路增益 - Costas环路的增益通常较小
pll_params.loopGain = 0.25;

% 积分时间
pll_params.PDI = PDI;

% 计算PLL环路滤波器系数
fprintf('    * 计算PLL环路滤波器系数...\n');
[pll_params.tau1, pll_params.tau2] = calLoopCoef(pll_params.noiseBandwidth, ...
                                                  pll_params.dampingRatio, ...
                                                  pll_params.loopGain);

% PLL性能参数
pll_params.wn = pll_params.noiseBandwidth / 0.7845;  % 自然频率
pll_params.omegaOf = pll_params.noiseBandwidth / 0.53;  % 开环频率

fprintf('    * PLL噪声带宽: %.2f Hz\n', pll_params.noiseBandwidth);
fprintf('    * PLL阻尼比: %.3f\n', pll_params.dampingRatio);
fprintf('    * PLL tau1: %.6f, tau2: %.6f\n', pll_params.tau1, pll_params.tau2);

%% 4. FLL (Frequency Lock Loop) 参数设置
fprintf('  - 设置FLL（频率跟踪环路）参数...\n');

% FLL用于辅助PLL，特别是在初始捕获阶段
fll_params = struct();

% FLL噪声带宽 - 通常比PLL更宽
fll_params.noiseBandwidth = pll_params.noiseBandwidth * 2;

% FLL阻尼比
fll_params.dampingRatio = 0.707;

% FLL环路增益
fll_params.loopGain = 1.0;

% FLL性能参数
fll_params.wn = fll_params.noiseBandwidth / 0.7845;
fll_params.omegaOf = fll_params.noiseBandwidth / 0.53;
fll_params.a2 = 1.414;  % FLL特有参数

% FLL到PLL切换阈值
fll_params.switch_threshold = 10;  % 在第10个循环后切换到PLL

fprintf('    * FLL噪声带宽: %.2f Hz\n', fll_params.noiseBandwidth);
fprintf('    * FLL切换阈值: %d 个循环\n', fll_params.switch_threshold);

%% 5. 初始状态变量设置
fprintf('  - 设置初始状态变量...\n');

% 相位和频率初始值
initial_states = struct();

% 载波相关初始状态
initial_states.carrPhase = 0.0;      % 载波相位
initial_states.carrFreq = fi;        % 载波频率
initial_states.remCarrPhase = 0.0;   % 剩余载波相位

% 码相关初始状态
initial_states.codePhase = 0.0;      % 码相位
initial_states.codeFreq = fp;        % 码频率
initial_states.remCodePhase = 0.0;   % 剩余码相位

% 环路滤波器状态变量
initial_states.oldCodeNco = 0.0;     % 码NCO历史值
initial_states.oldCodeError = 0.0;   % 码误差历史值
initial_states.oldCarrNco = 0.0;     % 载波NCO历史值
initial_states.oldCarrError = 0.0;   % 载波误差历史值

% PLL状态变量
initial_states.oldPllX = 0.0;        % PLL X状态
initial_states.oldPllY = 0.0;        % PLL Y状态
initial_states.old_IP_SUM = 0.0;     % I通道累积
initial_states.old_QP_SUM = 0.0;     % Q通道累积

% FLL状态变量
initial_states.oldfllx = 0.0;        % FLL X状态
initial_states.oldflly = 0.0;        % FLL Y状态
initial_states.oldCarrFreq = 0.0;    % 载波频率历史值

% 读取索引
initial_states.readIndex = 0;        % 信号读取索引

fprintf('    * 初始载波频率: %.2f Hz\n', initial_states.carrFreq);
fprintf('    * 初始码频率: %.2f Hz\n', initial_states.codeFreq);

%% 6. 性能监控参数设置
fprintf('  - 设置性能监控参数...\n');

% 性能监控参数
performance_params = struct();

% 收敛检测参数
performance_params.convergence_window = 50;      % 收敛检测窗口长度
performance_params.convergence_threshold = 0.1;  % 收敛阈值

% 失锁检测参数
performance_params.lock_detector_threshold = 0.5;  % 失锁检测阈值
performance_params.lock_detector_window = 20;      % 失锁检测窗口

% 信噪比估计参数
performance_params.snr_estimation_window = 100;    % 信噪比估计窗口
performance_params.snr_update_interval = 10;       % 信噪比更新间隔

% 自适应参数调整
performance_params.adaptive_bandwidth = true;      % 是否启用自适应带宽
performance_params.bandwidth_scale_factor = 1.5;   % 带宽缩放因子

fprintf('    * 收敛检测窗口: %d\n', performance_params.convergence_window);
fprintf('    * 失锁检测阈值: %.2f\n', performance_params.lock_detector_threshold);

%% 7. 鉴相器参数设置
fprintf('  - 设置鉴相器参数...\n');

% 鉴相器参数
discriminator_params = struct();

% 码鉴相器类型和参数
discriminator_params.code_discriminator_type = 'early_late';  % 早晚门鉴相器
discriminator_params.code_discriminator_gain = 1.0;          % 码鉴相器增益

% 载波鉴相器类型和参数
discriminator_params.carrier_discriminator_type = 'costas';  % Costas鉴相器
discriminator_params.carrier_discriminator_gain = 1.0;       % 载波鉴相器增益

% 频率鉴相器参数
discriminator_params.frequency_discriminator_type = 'cross_product';  % 叉积鉴频器
discriminator_params.frequency_discriminator_gain = 1.0;              % 频率鉴相器增益

fprintf('    * 码鉴相器: %s\n', discriminator_params.code_discriminator_type);
fprintf('    * 载波鉴相器: %s\n', discriminator_params.carrier_discriminator_type);

%% 8. 构建输出结构体
fprintf('  - 构建跟踪参数结构体...\n');

tracking_params = struct();

% 基本参数
tracking_params.basic = struct();
tracking_params.basic.fp = fp;
tracking_params.basic.fs = fs;
tracking_params.basic.fi = fi;
tracking_params.basic.lenOFDM = lenOFDM;
tracking_params.basic.PDI = PDI;
tracking_params.basic.dt = dt;

% DLL参数
tracking_params.dll = dll_params;

% PLL参数
tracking_params.pll = pll_params;

% FLL参数
tracking_params.fll = fll_params;

% 初始状态
tracking_params.initial_states = initial_states;

% 性能监控参数
tracking_params.performance = performance_params;

% 鉴相器参数
tracking_params.discriminator = discriminator_params;

% 兼容性参数（与原有代码兼容）
tracking_params.earlyLateSpc = dll_params.earlyLateSpc;
tracking_params.dllNoiseBandwidth = dll_params.noiseBandwidth;
tracking_params.dllDampingRatio = dll_params.dampingRatio;
tracking_params.pllNoiseBandwidth = pll_params.noiseBandwidth;
tracking_params.pllDampingRatio = pll_params.dampingRatio;

fprintf('  - 跟踪参数结构体构建完成\n');

%% 9. 参数合理性检查
fprintf('  - 参数合理性检查...\n');

% 检查采样率与码率的关系
if fs < 2 * fp
    warning('ofdm_tracking_init: 采样频率可能过低，建议 fs >= 2*fp');
end

% 检查环路带宽的合理性
if dll_params.noiseBandwidth > fp / 10
    warning('ofdm_tracking_init: DLL噪声带宽可能过宽，建议 < fp/10');
end

if pll_params.noiseBandwidth > fs / 100
    warning('ofdm_tracking_init: PLL噪声带宽可能过宽，建议 < fs/100');
end

% 检查积分时间的合理性
if PDI < 1e-3
    warning('ofdm_tracking_init: 积分时间过短，可能影响跟踪性能');
elseif PDI > 0.1
    warning('ofdm_tracking_init: 积分时间过长，可能影响动态响应');
end

% 检查环路滤波器系数的稳定性
if dll_params.tau1 <= 0 || dll_params.tau2 <= 0
    error('ofdm_tracking_init: DLL环路滤波器系数无效');
end

if pll_params.tau1 <= 0 || pll_params.tau2 <= 0
    error('ofdm_tracking_init: PLL环路滤波器系数无效');
end

fprintf('    * 参数合理性检查通过\n');

%% 10. 生成参数报告
fprintf('  - 生成参数配置报告...\n');

% 创建参数报告
report = struct();
report.timestamp = datestr(now);
report.dll_bandwidth_hz = dll_params.noiseBandwidth;
report.pll_bandwidth_hz = pll_params.noiseBandwidth;
report.integration_time_ms = PDI * 1000;
report.early_late_spacing = dll_params.earlyLateSpc;
report.damping_ratio = dll_params.dampingRatio;

% 估算跟踪精度
code_tracking_error_std = sqrt(dll_params.noiseBandwidth * PDI / 2);  % 码跟踪误差标准差
carrier_tracking_error_std = sqrt(pll_params.noiseBandwidth * PDI / 2);  % 载波跟踪误差标准差

report.estimated_code_tracking_error = code_tracking_error_std;
report.estimated_carrier_tracking_error = carrier_tracking_error_std;

% 添加报告到输出结构体
tracking_params.configuration_report = report;

fprintf('    * 估算码跟踪误差标准差: %.6f\n', code_tracking_error_std);
fprintf('    * 估算载波跟踪误差标准差: %.6f\n', carrier_tracking_error_std);

fprintf('  - 跟踪环路初始化完成\n');

end