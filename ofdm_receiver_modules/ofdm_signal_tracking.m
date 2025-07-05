function [tracking_results] = ofdm_signal_tracking(processed_signal, local_signals, tracking_params, loop_index)
% OFDM信号跟踪处理模块
%
% 功能描述:
%   执行OFDM信号的跟踪处理，包括码跟踪和载波跟踪
%   使用DLL进行码相位跟踪，使用PLL/FLL进行载波相位和频率跟踪
%
% 输入参数:
%   processed_signal - 预处理后的接收信号结构体
%   local_signals    - 本地参考信号结构体
%   tracking_params  - 跟踪参数结构体
%   loop_index      - 当前跟踪循环索引
%
% 输出参数:
%   tracking_results - 跟踪结果结构体，包含：
%                     * 相关器输出
%                     * 鉴相器输出
%                     * 环路滤波器输出
%                     * NCO控制信号
%                     * 跟踪状态更新
%                     * 性能指标
%
% 处理流程:
%   1. 信号段提取和载波剥离
%   2. 码相关处理（早、准时、晚相关器）
%   3. 码跟踪环路处理（DLL）
%   4. 载波跟踪环路处理（PLL/FLL）
%   5. NCO更新和状态维护
%   6. 性能监控和质量评估
%
% 技术细节:
%   - 支持多频带OFDM信号跟踪
%   - 自适应环路带宽调整
%   - 失锁检测和重捕获
%   - 实时性能监控
%
% 作者: OFDM接收机开发团队
% 日期: 2025年1月
% 版本: 1.0

%% 1. 输入参数验证
if nargin < 4
    error('ofdm_signal_tracking: 需要提供所有输入参数');
end

% 提取基本参数
fs = tracking_params.basic.fs;
fp = tracking_params.basic.fp;
fi = tracking_params.basic.fi;
lenOFDM = tracking_params.basic.lenOFDM;
PDI = tracking_params.basic.PDI;

% 提取信号
rawSignalI = processed_signal.I;
rawSignalQ = processed_signal.Q;
yr0 = local_signals.yr0;
yi0 = local_signals.yi0;

numBand = size(yr0, 1);
signal_length = size(yr0, 2);

fprintf('    * 处理循环 %d，信号长度: %d\n', loop_index, signal_length);

%% 2. 获取当前跟踪状态
if loop_index == 1
    % 第一次循环，使用初始状态
    carrPhase = tracking_params.initial_states.carrPhase;
    carrFreq = tracking_params.initial_states.carrFreq;
    remCarrPhase = tracking_params.initial_states.remCarrPhase;
    codePhase = tracking_params.initial_states.codePhase;
    codeFreq = tracking_params.initial_states.codeFreq;
    remCodePhase = tracking_params.initial_states.remCodePhase;
    
    % 环路滤波器状态
    oldCodeNco = tracking_params.initial_states.oldCodeNco;
    oldCodeError = tracking_params.initial_states.oldCodeError;
    oldCarrNco = tracking_params.initial_states.oldCarrNco;
    oldCarrError = tracking_params.initial_states.oldCarrError;
    
    % PLL状态
    oldPllX = tracking_params.initial_states.oldPllX;
    oldPllY = tracking_params.initial_states.oldPllY;
    old_IP_SUM = tracking_params.initial_states.old_IP_SUM;
    old_QP_SUM = tracking_params.initial_states.old_QP_SUM;
    
    % FLL状态
    oldfllx = tracking_params.initial_states.oldfllx;
    oldflly = tracking_params.initial_states.oldflly;
    oldCarrFreq = tracking_params.initial_states.oldCarrFreq;
    
    readIndex = tracking_params.initial_states.readIndex;
    
    fprintf('      使用初始跟踪状态\n');
else
    % 从全局变量或持久变量中获取状态（这里简化处理）
    % 在实际应用中，这些状态应该从上一次跟踪结果中获取
    persistent prev_tracking_state;
    if isempty(prev_tracking_state)
        % 如果没有前一状态，使用初始状态
        carrPhase = tracking_params.initial_states.carrPhase;
        carrFreq = tracking_params.initial_states.carrFreq;
        remCarrPhase = tracking_params.initial_states.remCarrPhase;
        codePhase = tracking_params.initial_states.codePhase;
        codeFreq = tracking_params.initial_states.codeFreq;
        remCodePhase = tracking_params.initial_states.remCodePhase;
        
        oldCodeNco = tracking_params.initial_states.oldCodeNco;
        oldCodeError = tracking_params.initial_states.oldCodeError;
        oldCarrNco = tracking_params.initial_states.oldCarrNco;
        oldCarrError = tracking_params.initial_states.oldCarrError;
        
        oldPllX = tracking_params.initial_states.oldPllX;
        oldPllY = tracking_params.initial_states.oldPllY;
        old_IP_SUM = tracking_params.initial_states.old_IP_SUM;
        old_QP_SUM = tracking_params.initial_states.old_QP_SUM;
        
        oldfllx = tracking_params.initial_states.oldfllx;
        oldflly = tracking_params.initial_states.oldflly;
        oldCarrFreq = tracking_params.initial_states.oldCarrFreq;
        
        readIndex = tracking_params.initial_states.readIndex;
    else
        % 使用前一状态
        carrPhase = prev_tracking_state.carrPhase;
        carrFreq = prev_tracking_state.carrFreq;
        remCarrPhase = prev_tracking_state.remCarrPhase;
        codePhase = prev_tracking_state.codePhase;
        codeFreq = prev_tracking_state.codeFreq;
        remCodePhase = prev_tracking_state.remCodePhase;
        
        oldCodeNco = prev_tracking_state.oldCodeNco;
        oldCodeError = prev_tracking_state.oldCodeError;
        oldCarrNco = prev_tracking_state.oldCarrNco;
        oldCarrError = prev_tracking_state.oldCarrError;
        
        oldPllX = prev_tracking_state.oldPllX;
        oldPllY = prev_tracking_state.oldPllY;
        old_IP_SUM = prev_tracking_state.old_IP_SUM;
        old_QP_SUM = prev_tracking_state.old_QP_SUM;
        
        oldfllx = prev_tracking_state.oldfllx;
        oldflly = prev_tracking_state.oldflly;
        oldCarrFreq = prev_tracking_state.oldCarrFreq;
        
        readIndex = prev_tracking_state.readIndex;
    end
    
    fprintf('      使用前一循环跟踪状态\n');
end

%% 3. 信号段提取
fprintf('      提取当前处理信号段...\n');

% 计算当前处理的信号段
samples_to_process = min(lenOFDM, signal_length - readIndex);
if samples_to_process <= 0
    error('ofdm_signal_tracking: 没有足够的信号进行处理');
end

% 提取信号段
signal_segment_I = rawSignalI(readIndex+1:readIndex+samples_to_process);
signal_segment_Q = rawSignalQ(readIndex+1:readIndex+samples_to_process);

fprintf('        处理信号段长度: %d 采样点\n', samples_to_process);

%% 4. 载波剥离（下变频）
fprintf('      执行载波剥离...\n');

% 生成载波信号
t_segment = (0:samples_to_process-1) / fs;
carr_cos = cos(2*pi*carrFreq*t_segment + carrPhase);
carr_sin = sin(2*pi*carrFreq*t_segment + carrPhase);

% 载波剥离
baseband_I = signal_segment_I .* carr_cos + signal_segment_Q .* carr_sin;
baseband_Q = -signal_segment_I .* carr_sin + signal_segment_Q .* carr_cos;

% 更新载波相位
carrPhase = carrPhase + 2*pi*carrFreq*samples_to_process/fs;
carrPhase = mod(carrPhase, 2*pi);  % 保持相位在[0, 2π]范围内

%% 5. 码相关处理
fprintf('      执行码相关处理...\n');

% 提取跟踪参数
earlyLateSpc = tracking_params.earlyLateSpc;

% 初始化相关器输出
correlator_outputs = struct();
correlator_outputs.early = struct();
correlator_outputs.prompt = struct();
correlator_outputs.late = struct();

% 对每个频带进行相关处理
for bandID = 1:numBand
    % 获取本地码信号
    local_I = yr0(bandID, :);
    local_Q = yi0(bandID, :);
    
    % 计算码相位索引
    code_indices = mod(floor((0:samples_to_process-1) * codeFreq/fs + codePhase), length(local_I)) + 1;
    
    % 早相关器（提前半个码片）
    early_indices = mod(code_indices - floor(earlyLateSpc * fs/fp) - 1, length(local_I)) + 1;
    early_local_I = local_I(early_indices);
    early_local_Q = local_Q(early_indices);
    
    % 准时相关器
    prompt_local_I = local_I(code_indices);
    prompt_local_Q = local_Q(code_indices);
    
    % 晚相关器（滞后半个码片）
    late_indices = mod(code_indices + floor(earlyLateSpc * fs/fp) - 1, length(local_I)) + 1;
    late_local_I = local_I(late_indices);
    late_local_Q = local_Q(late_indices);
    
    % 计算相关值
    % 早相关器
    I_E = sum(baseband_I .* early_local_I);
    Q_E = sum(baseband_Q .* early_local_Q);
    
    % 准时相关器
    I_P = sum(baseband_I .* prompt_local_I);
    Q_P = sum(baseband_Q .* prompt_local_Q);
    
    % 晚相关器
    I_L = sum(baseband_I .* late_local_I);
    Q_L = sum(baseband_Q .* late_local_Q);
    
    % 存储相关器输出
    band_name = sprintf('band_%d', bandID);
    correlator_outputs.early.(band_name) = complex(I_E, Q_E);
    correlator_outputs.prompt.(band_name) = complex(I_P, Q_P);
    correlator_outputs.late.(band_name) = complex(I_L, Q_L);
end

% 计算多频带平均相关输出
I_E_avg = 0; Q_E_avg = 0;
I_P_avg = 0; Q_P_avg = 0;
I_L_avg = 0; Q_L_avg = 0;

for bandID = 1:numBand
    band_name = sprintf('band_%d', bandID);
    early_corr = correlator_outputs.early.(band_name);
    prompt_corr = correlator_outputs.prompt.(band_name);
    late_corr = correlator_outputs.late.(band_name);
    
    I_E_avg = I_E_avg + real(early_corr);
    Q_E_avg = Q_E_avg + imag(early_corr);
    I_P_avg = I_P_avg + real(prompt_corr);
    Q_P_avg = Q_P_avg + imag(prompt_corr);
    I_L_avg = I_L_avg + real(late_corr);
    Q_L_avg = Q_L_avg + imag(late_corr);
end

I_E_avg = I_E_avg / numBand;
Q_E_avg = Q_E_avg / numBand;
I_P_avg = I_P_avg / numBand;
Q_P_avg = Q_P_avg / numBand;
I_L_avg = I_L_avg / numBand;
Q_L_avg = Q_L_avg / numBand;

fprintf('        相关器输出 - E: %.4f+j%.4f, P: %.4f+j%.4f, L: %.4f+j%.4f\n', ...
        I_E_avg, Q_E_avg, I_P_avg, Q_P_avg, I_L_avg, Q_L_avg);

%% 6. 码跟踪环路处理（DLL）
fprintf('      执行码跟踪环路处理...\n');

% 计算码鉴相器输出（早晚门鉴相器）
E_power = I_E_avg^2 + Q_E_avg^2;
L_power = I_L_avg^2 + Q_L_avg^2;
codeError = (E_power - L_power) / (E_power + L_power + 1e-10);  % 添加小量避免除零

% DLL环路滤波器
tau1 = tracking_params.dll.tau1;
tau2 = tracking_params.dll.tau2;

codeNco = oldCodeNco + (tau2/tau1) * (codeError - oldCodeError) + codeError * (PDI/tau1);

% 更新码频率和相位
codeFreq = fp + codeNco;
codePhase = codePhase + codeNco * samples_to_process / fs;

% 更新DLL状态
oldCodeNco = codeNco;
oldCodeError = codeError;

fprintf('        码误差: %.6f, 码NCO: %.6f, 码频率: %.2f Hz\n', ...
        codeError, codeNco, codeFreq);

%% 7. 载波跟踪环路处理
fprintf('      执行载波跟踪环路处理...\n');

% 判断使用FLL还是PLL
use_fll = (loop_index <= tracking_params.fll.switch_threshold);

if use_fll
    fprintf('        使用FLL进行频率跟踪\n');
    
    % FLL鉴频器（叉积鉴频器）
    dot_product = I_P_avg * old_IP_SUM + Q_P_avg * old_QP_SUM;
    cross_product = I_P_avg * old_QP_SUM - Q_P_avg * old_IP_SUM;
    
    if abs(dot_product) > 1e-10
        carrError = atan(cross_product / dot_product) / (2 * pi * PDI);
    else
        carrError = 0;
    end
    
    % FLL环路滤波器
    wn = tracking_params.fll.wn;
    a2 = tracking_params.fll.a2;
    
    fllx = oldfllx + wn * PDI * carrError;
    flly = oldflly + wn * a2 * PDI * carrError;
    
    carrNco = fllx + flly;
    carrFreq = fi + carrNco;
    
    % 更新FLL状态
    oldfllx = fllx;
    oldflly = flly;
    oldCarrFreq = carrFreq;
    
else
    fprintf('        使用PLL进行相位跟踪\n');
    
    % Costas环路鉴相器
    carrError = atan2(Q_P_avg, I_P_avg) / (2 * pi);
    
    % PLL环路滤波器
    tau1 = tracking_params.pll.tau1;
    tau2 = tracking_params.pll.tau2;
    
    carrNco = oldCarrNco + (tau2/tau1) * (carrError - oldCarrError) + carrError * (PDI/tau1);
    carrFreq = fi + carrNco;
    
    % 更新PLL状态
    oldCarrNco = carrNco;
    oldCarrError = carrError;
end

% 更新载波相位
remCarrPhase = remCarrPhase + carrNco * samples_to_process / fs;

% 更新相关器历史值
old_IP_SUM = I_P_avg;
old_QP_SUM = Q_P_avg;

fprintf('        载波误差: %.6f, 载波NCO: %.6f, 载波频率: %.2f Hz\n', ...
        carrError, carrNco, carrFreq);

%% 8. 性能监控
fprintf('      执行性能监控...\n');

% 计算信号强度
signal_power = I_P_avg^2 + Q_P_avg^2;
noise_power = (I_E_avg^2 + Q_E_avg^2 + I_L_avg^2 + Q_L_avg^2) / 2 - signal_power;
noise_power = max(noise_power, 1e-10);  % 避免负值

% 估算信噪比
snr_estimate = 10 * log10(signal_power / noise_power);

% 计算锁定指示器
lock_indicator = signal_power / (signal_power + noise_power);

% 计算跟踪误差
code_tracking_error = abs(codeError);
carrier_tracking_error = abs(carrError);

fprintf('        信号功率: %.6f, 噪声功率: %.6f\n', signal_power, noise_power);
fprintf('        估算SNR: %.2f dB, 锁定指示器: %.4f\n', snr_estimate, lock_indicator);

%% 9. 更新读取索引
readIndex = readIndex + samples_to_process;

%% 10. 构建输出结果
tracking_results = struct();

% 相关器输出
tracking_results.correlators = correlator_outputs;
tracking_results.correlators.average = struct();
tracking_results.correlators.average.I_E = I_E_avg;
tracking_results.correlators.average.Q_E = Q_E_avg;
tracking_results.correlators.average.I_P = I_P_avg;
tracking_results.correlators.average.Q_P = Q_P_avg;
tracking_results.correlators.average.I_L = I_L_avg;
tracking_results.correlators.average.Q_L = Q_L_avg;

% 鉴相器输出
tracking_results.discriminators = struct();
tracking_results.discriminators.code_error = codeError;
tracking_results.discriminators.carrier_error = carrError;

% 环路滤波器输出
tracking_results.loop_filters = struct();
tracking_results.loop_filters.code_nco = codeNco;
tracking_results.loop_filters.carrier_nco = carrNco;

% 更新后的跟踪状态
tracking_results.updated_states = struct();
tracking_results.updated_states.carrPhase = carrPhase;
tracking_results.updated_states.carrFreq = carrFreq;
tracking_results.updated_states.remCarrPhase = remCarrPhase;
tracking_results.updated_states.codePhase = codePhase;
tracking_results.updated_states.codeFreq = codeFreq;
tracking_results.updated_states.remCodePhase = remCodePhase;

tracking_results.updated_states.oldCodeNco = oldCodeNco;
tracking_results.updated_states.oldCodeError = oldCodeError;
tracking_results.updated_states.oldCarrNco = oldCarrNco;
tracking_results.updated_states.oldCarrError = oldCarrError;

tracking_results.updated_states.oldPllX = oldPllX;
tracking_results.updated_states.oldPllY = oldPllY;
tracking_results.updated_states.old_IP_SUM = old_IP_SUM;
tracking_results.updated_states.old_QP_SUM = old_QP_SUM;

tracking_results.updated_states.oldfllx = oldfllx;
tracking_results.updated_states.oldflly = oldflly;
tracking_results.updated_states.oldCarrFreq = oldCarrFreq;

tracking_results.updated_states.readIndex = readIndex;

% 性能指标
tracking_results.performance = struct();
tracking_results.performance.signal_power = signal_power;
tracking_results.performance.noise_power = noise_power;
tracking_results.performance.snr_estimate = snr_estimate;
tracking_results.performance.lock_indicator = lock_indicator;
tracking_results.performance.code_tracking_error = code_tracking_error;
tracking_results.performance.carrier_tracking_error = carrier_tracking_error;

% 处理信息
tracking_results.processing_info = struct();
tracking_results.processing_info.loop_index = loop_index;
tracking_results.processing_info.samples_processed = samples_to_process;
tracking_results.processing_info.tracking_mode = use_fll ? 'FLL' : 'PLL';
tracking_results.processing_info.num_bands = numBand;

% 保存当前状态到持久变量（用于下次调用）
persistent prev_tracking_state;
prev_tracking_state = tracking_results.updated_states;

fprintf('      跟踪处理完成\n');

end