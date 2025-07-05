function [processed_signal, init_params] = ofdm_receiver_init(simSettings, received_signal)
% OFDM接收机初始化模块
%
% 功能描述:
%   负责接收机的参数初始化和输入信号的预处理工作
%   包括信号格式转换、长度检查、参数提取等基础处理
%
% 输入参数:
%   simSettings    - 仿真设置参数结构体，包含：
%                    * fp: 码频率 (Hz)
%                    * fs: 采样频率 (Hz) 
%                    * fi: 中频频率 (Hz)
%                    * fc: 载波频率 (Hz)
%                    * numBand: 频带数量
%                    * SNR: 信噪比 (dB)
%                    * dt: 时间间隔 (s)
%                    * t_total: 总仿真时间 (s)
%                    * NFFT: FFT点数
%                    * nSymbol: OFDM符号数
%   received_signal - 接收到的信号，支持两种格式：
%                    * [I; Q] 格式：2×N矩阵，第一行为I通道，第二行为Q通道
%                    * 复数格式：1×N复数向量
%
% 输出参数:
%   processed_signal - 预处理后的信号结构体，包含：
%                     * I: I通道信号 (实部)
%                     * Q: Q通道信号 (虚部)
%                     * length: 信号长度
%                     * format: 原始信号格式标识
%   init_params     - 初始化参数结构体，包含：
%                     * fp, fs, fi, fc: 各种频率参数
%                     * numBand: 频带数量
%                     * SNR: 信噪比
%                     * dt: 时间间隔
%                     * t_sim: 仿真时间
%                     * lenOFDM: OFDM符号长度
%                     * num_samples: 预期样本数
%
% 处理流程:
%   1. 输入参数验证
%   2. 信号格式识别和转换
%   3. 信号长度检查和调整
%   4. 参数提取和计算
%   5. 输出结构体构建
%
% 异常处理:
%   - 不支持的信号格式会抛出错误
%   - 信号长度不足时会进行零填充并给出警告
%   - 参数缺失时会使用默认值并给出提示
%
% 作者: OFDM接收机开发团队
% 日期: 2025年1月
% 版本: 1.0

%% 1. 输入参数验证
if nargin < 2
    error('ofdm_receiver_init: 必须提供仿真设置和接收信号');
end

% 检查simSettings结构体的必要字段
required_fields = {'fp', 'fs', 'fi', 'fc', 'numBand', 'SNR', 'dt', 't_total', 'NFFT', 'nSymbol'};
for i = 1:length(required_fields)
    if ~isfield(simSettings, required_fields{i})
        error('ofdm_receiver_init: simSettings缺少必要字段 %s', required_fields{i});
    end
end

fprintf('  - 输入参数验证通过\n');

%% 2. 信号格式识别和转换
fprintf('  - 识别输入信号格式...\n');

% 获取信号维度信息
[rows, cols] = size(received_signal);

if rows == 2 && cols > 1
    % I/Q分离格式：2×N矩阵
    rawSignalI = received_signal(1, :);
    rawSignalQ = received_signal(2, :);
    signal_format = 'IQ_separate';
    fprintf('    * 检测到I/Q分离格式信号\n');
    
elseif rows == 1 && cols > 1 && ~isreal(received_signal)
    % 复数格式：1×N复数向量
    rawSignalI = real(received_signal);
    rawSignalQ = imag(received_signal);
    signal_format = 'complex';
    fprintf('    * 检测到复数格式信号\n');
    
elseif cols == 1 && rows > 1 && ~isreal(received_signal)
    % 复数格式：N×1复数向量（转置）
    received_signal = received_signal.';
    rawSignalI = real(received_signal);
    rawSignalQ = imag(received_signal);
    signal_format = 'complex_transposed';
    fprintf('    * 检测到复数格式信号（已转置）\n');
    
elseif rows == 1 && cols > 1 && isreal(received_signal)
    % 实数信号，假设为I通道，Q通道置零
    rawSignalI = received_signal;
    rawSignalQ = zeros(size(received_signal));
    signal_format = 'real_only';
    fprintf('    * 检测到实数信号，Q通道置零\n');
    warning('ofdm_receiver_init: 输入为实数信号，Q通道已置零，可能影响性能');
    
else
    error('ofdm_receiver_init: 不支持的信号格式，维度为 %dx%d', rows, cols);
end

% 确保信号为行向量
if size(rawSignalI, 1) > 1
    rawSignalI = rawSignalI.';
    rawSignalQ = rawSignalQ.';
end

original_length = length(rawSignalI);
fprintf('    * 原始信号长度: %d 采样点\n', original_length);

%% 3. 信号长度检查和调整
fprintf('  - 检查信号长度...\n');

% 计算预期信号长度
t0 = 0:simSettings.dt:simSettings.t_total;
num_samples = length(t0);
lenOFDM = simSettings.NFFT * simSettings.nSymbol;

% 估算每个处理周期需要的采样点数
samples_per_cycle = ceil(lenOFDM * simSettings.fs / simSettings.fp);
expected_length = num_samples * samples_per_cycle;

fprintf('    * 预期信号长度: %d 采样点\n', expected_length);
fprintf('    * 处理周期数: %d\n', num_samples);
fprintf('    * 每周期采样点: %d\n', samples_per_cycle);

if original_length < expected_length
    % 信号长度不足，进行零填充
    padding_length = expected_length - original_length;
    rawSignalI = [rawSignalI, zeros(1, padding_length)];
    rawSignalQ = [rawSignalQ, zeros(1, padding_length)];
    
    fprintf('    * 警告: 信号长度不足，已零填充 %d 个采样点\n', padding_length);
    warning('ofdm_receiver_init: 信号长度不足，已进行零填充，可能影响处理结果');
    
elseif original_length > expected_length
    % 信号长度过长，截取所需长度
    rawSignalI = rawSignalI(1:expected_length);
    rawSignalQ = rawSignalQ(1:expected_length);
    
    fprintf('    * 信息: 信号长度过长，已截取至 %d 个采样点\n', expected_length);
    
else
    fprintf('    * 信号长度匹配，无需调整\n');
end

final_length = length(rawSignalI);

%% 4. 信号质量初步检查
fprintf('  - 信号质量初步检查...\n');

% 检查信号动态范围
I_max = max(abs(rawSignalI));
Q_max = max(abs(rawSignalQ));
I_rms = rms(rawSignalI);
Q_rms = rms(rawSignalQ);

fprintf('    * I通道: 最大值=%.4f, RMS=%.4f\n', I_max, I_rms);
fprintf('    * Q通道: 最大值=%.4f, RMS=%.4f\n', Q_max, Q_rms);

% 检查信号是否存在饱和
saturation_threshold = 0.95; % 饱和阈值
if I_max > saturation_threshold || Q_max > saturation_threshold
    warning('ofdm_receiver_init: 信号可能存在饱和，最大幅度接近1.0');
end

% 检查信号是否过小
min_signal_threshold = 1e-6;
if I_rms < min_signal_threshold && Q_rms < min_signal_threshold
    warning('ofdm_receiver_init: 信号幅度过小，可能影响处理性能');
end

% 估算初始信噪比
signal_power = I_rms^2 + Q_rms^2;
estimated_snr = 10 * log10(signal_power / (2 * 1e-6)); % 假设噪声功率
fprintf('    * 估算信噪比: %.2f dB\n', estimated_snr);

%% 5. 参数提取和计算
fprintf('  - 提取和计算系统参数...\n');

% 基本频率参数
fp = simSettings.fp;           % 码频率
fs = simSettings.fs;           % 采样频率
fi = simSettings.fi;           % 中频频率
fc = simSettings.fc;           % 载波频率

% 系统参数
numBand = simSettings.numBand; % 频带数量
SNR = simSettings.SNR;         % 信噪比
dt = simSettings.dt;           % 时间间隔
t_sim = simSettings.t_total;   % 仿真时间

% OFDM参数
NFFT = simSettings.NFFT;       % FFT点数
nSymbol = simSettings.nSymbol; % 符号数
lenOFDM = NFFT * nSymbol;      % OFDM符号长度

% 计算采样率相关参数
oversampling_ratio = fs / fp;  % 过采样率
samples_per_symbol = NFFT * oversampling_ratio / numBand; % 每符号采样数

fprintf('    * 码频率: %.2f kHz\n', fp/1e3);
fprintf('    * 采样频率: %.2f MHz\n', fs/1e6);
fprintf('    * 过采样率: %.2f\n', oversampling_ratio);
fprintf('    * OFDM符号长度: %d\n', lenOFDM);
fprintf('    * 每符号采样数: %.1f\n', samples_per_symbol);

%% 6. 构建输出结构体
fprintf('  - 构建输出结构体...\n');

% 预处理后的信号结构体
processed_signal = struct();
processed_signal.I = rawSignalI;
processed_signal.Q = rawSignalQ;
processed_signal.length = final_length;
processed_signal.format = signal_format;
processed_signal.original_length = original_length;
processed_signal.padding_applied = (final_length > original_length);
processed_signal.truncation_applied = (original_length > expected_length);

% 信号统计信息
processed_signal.statistics = struct();
processed_signal.statistics.I_max = I_max;
processed_signal.statistics.Q_max = Q_max;
processed_signal.statistics.I_rms = I_rms;
processed_signal.statistics.Q_rms = Q_rms;
processed_signal.statistics.estimated_snr = estimated_snr;
processed_signal.statistics.signal_power = signal_power;

% 初始化参数结构体
init_params = struct();

% 频率参数
init_params.fp = fp;
init_params.fs = fs;
init_params.fi = fi;
init_params.fc = fc;

% 系统参数
init_params.numBand = numBand;
init_params.SNR = SNR;
init_params.dt = dt;
init_params.t_sim = t_sim;

% OFDM参数
init_params.NFFT = NFFT;
init_params.nSymbol = nSymbol;
init_params.lenOFDM = lenOFDM;

% 计算参数
init_params.num_samples = num_samples;
init_params.samples_per_cycle = samples_per_cycle;
init_params.expected_length = expected_length;
init_params.oversampling_ratio = oversampling_ratio;
init_params.samples_per_symbol = samples_per_symbol;

% 处理信息
init_params.processing_info = struct();
init_params.processing_info.signal_format = signal_format;
init_params.processing_info.length_adjustment = struct();
init_params.processing_info.length_adjustment.original = original_length;
init_params.processing_info.length_adjustment.final = final_length;
init_params.processing_info.length_adjustment.expected = expected_length;
init_params.processing_info.length_adjustment.padding = (final_length > original_length);
init_params.processing_info.length_adjustment.truncation = (original_length > expected_length);

fprintf('  - 初始化完成\n');

%% 7. 最终验证
fprintf('  - 最终验证...\n');

% 验证输出信号的完整性
if length(processed_signal.I) ~= length(processed_signal.Q)
    error('ofdm_receiver_init: I/Q通道长度不匹配');
end

if any(~isfinite(processed_signal.I)) || any(~isfinite(processed_signal.Q))
    error('ofdm_receiver_init: 信号包含无效值（NaN或Inf）');
end

% 验证参数的合理性
if init_params.fs <= init_params.fp
    warning('ofdm_receiver_init: 采样频率可能过低，建议fs > 2*fp');
end

if init_params.num_samples < 10
    warning('ofdm_receiver_init: 处理样本数过少，可能影响跟踪性能');
end

fprintf('    * 验证通过，初始化模块完成\n');

end