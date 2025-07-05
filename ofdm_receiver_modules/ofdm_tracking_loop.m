function tracking_output = ofdm_tracking_loop(processed_signal, local_signals, tracking_params, kf_params, simSettings)
% OFDM信号跟踪循环处理模块
%
% 功能描述:
%   执行OFDM信号的跟踪处理，包括码跟踪和载波跟踪
%   支持传统环路滤波器和卡尔曼滤波器两种跟踪方式
%
% 输入参数:
%   processed_signal - 预处理后的接收信号
%   local_signals   - 本地参考信号
%   tracking_params - 跟踪参数
%   kf_params      - 卡尔曼滤波器参数 (可选)
%   simSettings    - 仿真设置参数
%
% 输出参数:
%   tracking_output - 跟踪结果结构体
%
% 作者: OFDM接收机开发
% 日期: 2025年7月
% 版本: 1.0

%% 输入参数验证
if nargin < 4
    error('ofdm_tracking_loop: 至少需要提供前4个输入参数');
end

if nargin < 5
    simSettings = struct();
end

use_kalman = (nargin >= 4) && ~isempty(kf_params);

fprintf('  - 开始OFDM信号跟踪处理...\n');
fprintf('    * 使用卡尔曼滤波器: %s\n', char(string(use_kalman)));

%% 提取信号和参数
I_signal = processed_signal.I;
Q_signal = processed_signal.Q;
signal_length = length(I_signal);

yr0 = local_signals.yr0;
yi0 = local_signals.yi0;

% 提取跟踪参数
dll_params = tracking_params.dll;
pll_params = tracking_params.pll;
basic_params = tracking_params.basic;

fprintf('    * 信号长度: %d 采样点\n', signal_length);
fprintf('    * 本地信号大小: %dx%d\n', size(yr0));

%% 初始化跟踪变量
% 时间向量
time_samples = 1:signal_length;
num_loops = floor(signal_length / basic_params.lenOFDM);

fprintf('    * 跟踪循环数: %d\n', num_loops);

% 初始化输出数组
OutCarrFreq = zeros(1, num_loops);
OutCodeFreq = zeros(1, num_loops);
OutCarrPhase = zeros(1, num_loops);
OutCodePhase = zeros(1, num_loops);
OutSNR = zeros(1, num_loops);
OutLockIndicator = zeros(1, num_loops);

% 鉴别器输出
DLL_discriminator = zeros(1, num_loops);
PLL_discriminator = zeros(1, num_loops);

% 环路滤波器状态
dll_filter_state = [0, 0];  % [积分器, 比例器]
pll_filter_state = [0, 0];

% 当前估计值
current_code_freq = basic_params.fp;
current_carr_freq = basic_params.fi;
current_code_phase = 0;
current_carr_phase = 0;

%% 卡尔曼滤波器初始化
if use_kalman
    kf_state = kf_params.x;
    kf_covariance = kf_params.P;
    fprintf('    * 卡尔曼滤波器已初始化\n');
end

%% 主跟踪循环
fprintf('  - 执行主跟踪循环...\n');

for loop_idx = 1:num_loops
    if mod(loop_idx, 100) == 0
        fprintf('    * 处理循环 %d/%d (%.1f%%)\n', loop_idx, num_loops, 100*loop_idx/num_loops);
    end
    
    % 计算当前循环的样本范围
    start_sample = (loop_idx - 1) * basic_params.lenOFDM + 1;
    end_sample = min(start_sample + basic_params.lenOFDM - 1, signal_length);
    
    if end_sample <= start_sample
        break;
    end
    
    % 提取当前循环的信号段
    current_I = I_signal(start_sample:end_sample);
    current_Q = Q_signal(start_sample:end_sample);
    
    %% 生成本地载波
    sample_indices = start_sample:end_sample;
    local_carr_I = cos(2*pi*current_carr_freq/basic_params.fs * sample_indices + current_carr_phase);
    local_carr_Q = sin(2*pi*current_carr_freq/basic_params.fs * sample_indices + current_carr_phase);
    
    %% 载波剥离
    baseband_I = current_I .* local_carr_I + current_Q .* local_carr_Q;
    baseband_Q = current_Q .* local_carr_I - current_I .* local_carr_Q;
    
    %% 码相关处理
    % 计算本地码的起始位置
    code_samples_per_chip = basic_params.fs / current_code_freq;
    code_start_sample = mod(current_code_phase, size(yr0, 2));
    
    % 生成早、准时、晚本地码
    [early_code_I, prompt_code_I, late_code_I] = generate_local_codes(yr0, yi0, ...
        code_start_sample, length(current_I), dll_params.earlyLateSpc);
    [early_code_Q, prompt_code_Q, late_code_Q] = generate_local_codes(yi0, yr0, ...
        code_start_sample, length(current_I), dll_params.earlyLateSpc);
    
    %% 相关器输出
    % 早期相关器
    IE = sum(baseband_I .* early_code_I);
    QE = sum(baseband_Q .* early_code_Q);
    
    % 准时相关器
    IP = sum(baseband_I .* prompt_code_I);
    QP = sum(baseband_Q .* prompt_code_Q);
    
    % 晚期相关器
    IL = sum(baseband_I .* late_code_I);
    QL = sum(baseband_Q .* late_code_Q);
    
    %% 鉴别器计算
    % DLL鉴别器 (Early-Late Power)
    early_power = IE^2 + QE^2;
    late_power = IL^2 + QL^2;
    dll_disc = (early_power - late_power) / (early_power + late_power + eps);
    
    % PLL鉴别器 (Costas Loop)
    pll_disc = atan2(QP, IP);
    
    % 存储鉴别器输出
    DLL_discriminator(loop_idx) = dll_disc;
    PLL_discriminator(loop_idx) = pll_disc;
    
    %% 环路滤波器更新
    if use_kalman
        % 卡尔曼滤波器更新
        [kf_state, kf_covariance] = kalman_filter_update(kf_state, kf_covariance, ...
            [dll_disc; pll_disc], kf_params);
        
        % 从卡尔曼滤波器状态提取估计值
        code_freq_correction = kf_state(2);
        carr_freq_correction = kf_state(3);
        code_phase_correction = kf_state(1);
        carr_phase_correction = kf_state(4);
    else
        % 传统环路滤波器
        % DLL环路滤波器
        dll_filter_state(1) = dll_filter_state(1) + dll_params.tau1 * dll_disc;
        dll_filter_output = dll_filter_state(1) + dll_params.tau2 * dll_disc;
        
        % PLL环路滤波器
        pll_filter_state(1) = pll_filter_state(1) + pll_params.tau1 * pll_disc;
        pll_filter_output = pll_filter_state(1) + pll_params.tau2 * pll_disc;
        
        % 频率和相位修正
        code_freq_correction = dll_filter_output;
        carr_freq_correction = pll_filter_output;
        code_phase_correction = dll_filter_output * basic_params.PDI;
        carr_phase_correction = pll_filter_output * basic_params.PDI;
    end
    
    %% 更新载波和码的频率/相位
    current_code_freq = basic_params.fp + code_freq_correction;
    current_carr_freq = basic_params.fi + carr_freq_correction;
    current_code_phase = current_code_phase + code_phase_correction;
    current_carr_phase = current_carr_phase + carr_phase_correction;
    
    % 相位归一化
    current_carr_phase = mod(current_carr_phase, 2*pi);
    current_code_phase = mod(current_code_phase, size(yr0, 2));
    
    %% 性能指标计算
    % 信噪比估计
    signal_power = IP^2 + QP^2;
    noise_power = var([IE, QE, IL, QL]) / 2;
    snr_linear = signal_power / (noise_power + eps);
    snr_db = 10 * log10(snr_linear + eps);
    
    % 锁定指示器
    lock_indicator = signal_power / (signal_power + noise_power + eps);
    
    %% 存储结果
    OutCarrFreq(loop_idx) = current_carr_freq;
    OutCodeFreq(loop_idx) = current_code_freq;
    OutCarrPhase(loop_idx) = current_carr_phase;
    OutCodePhase(loop_idx) = current_code_phase;
    OutSNR(loop_idx) = snr_db;
    OutLockIndicator(loop_idx) = lock_indicator;
end

%% 计算平均性能指标
valid_indices = 1:loop_idx;
Average_SNR = mean(OutSNR(valid_indices));
Average_Lock = mean(OutLockIndicator(valid_indices));

fprintf('    * 跟踪循环完成\n');
fprintf('    * 平均信噪比: %.2f dB\n', Average_SNR);
fprintf('    * 平均锁定指示器: %.4f\n', Average_Lock);

%% 构建输出结构体
tracking_output = struct();

% 主要跟踪结果
tracking_output.OutCarrFreq = OutCarrFreq(valid_indices);
tracking_output.OutCodeFreq = OutCodeFreq(valid_indices);
tracking_output.OutCarrPhase = OutCarrPhase(valid_indices);
tracking_output.OutCodePhase = OutCodePhase(valid_indices);
tracking_output.OutSNR = OutSNR(valid_indices);
tracking_output.OutLockIndicator = OutLockIndicator(valid_indices);

% 鉴别器输出
tracking_output.DLL_discriminator = DLL_discriminator(valid_indices);
tracking_output.PLL_discriminator = PLL_discriminator(valid_indices);

% 性能指标
tracking_output.Average_SNR = Average_SNR;
tracking_output.Average_Lock = Average_Lock;
tracking_output.num_loops = length(valid_indices);

% 卡尔曼滤波器结果
if use_kalman
    tracking_output.kalman_states = kf_state;
    tracking_output.kalman_covariance = kf_covariance;
end

% 处理信息
tracking_output.processing_info = struct();
tracking_output.processing_info.use_kalman = use_kalman;
tracking_output.processing_info.signal_length = signal_length;
tracking_output.processing_info.samples_per_loop = basic_params.lenOFDM;

fprintf('  - OFDM信号跟踪处理完成\n');

end

function [early_code, prompt_code, late_code] = generate_local_codes(code_I, code_Q, start_sample, length_needed, early_late_spacing)
% 生成早、准时、晚本地码
%
% 输入参数:
%   code_I, code_Q - 本地码的I和Q分量
%   start_sample   - 起始采样点
%   length_needed  - 需要的长度
%   early_late_spacing - 早晚间距

code_length = size(code_I, 2);

% 计算早、准时、晚的起始位置
prompt_start = mod(start_sample - 1, code_length) + 1;
early_start = mod(prompt_start - early_late_spacing - 1, code_length) + 1;
late_start = mod(prompt_start + early_late_spacing - 1, code_length) + 1;

% 生成码序列
prompt_code = generate_code_sequence(code_I, prompt_start, length_needed);
early_code = generate_code_sequence(code_I, early_start, length_needed);
late_code = generate_code_sequence(code_I, late_start, length_needed);

end

function code_sequence = generate_code_sequence(code_matrix, start_pos, length_needed)
% 从码矩阵生成指定长度的码序列

if isempty(code_matrix)
    code_sequence = zeros(1, length_needed);
    return;
end

[num_bands, code_length] = size(code_matrix);

% 使用第一个频带的码
code_vector = code_matrix(1, :);

% 生成重复的码序列
num_repeats = ceil(length_needed / code_length);
extended_code = repmat(code_vector, 1, num_repeats);

% 提取所需长度的序列
end_pos = start_pos + length_needed - 1;
if end_pos <= length(extended_code)
    code_sequence = extended_code(start_pos:end_pos);
else
    % 处理边界情况
    code_sequence = zeros(1, length_needed);
    available_length = min(length_needed, length(extended_code) - start_pos + 1);
    code_sequence(1:available_length) = extended_code(start_pos:start_pos+available_length-1);
end

end

function [updated_state, updated_covariance] = kalman_filter_update(state, covariance, observation, kf_params)
% 卡尔曼滤波器更新步骤

% 预测步骤
predicted_state = kf_params.A * state;
predicted_covariance = kf_params.A * covariance * kf_params.A' + kf_params.Q;

% 更新步骤
innovation = observation - kf_params.H * predicted_state;
innovation_covariance = kf_params.H * predicted_covariance * kf_params.H' + kf_params.R;
kalman_gain = predicted_covariance * kf_params.H' / innovation_covariance;

% 状态和协方差更新
updated_state = predicted_state + kalman_gain * innovation;
updated_covariance = (eye(size(covariance)) - kalman_gain * kf_params.H) * predicted_covariance;

end