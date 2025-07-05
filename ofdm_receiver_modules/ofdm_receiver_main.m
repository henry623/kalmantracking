function [output] = ofdm_receiver_main(simSettings, received_signal, enable_kalman)
% OFDM接收机主函数 - 模块化设计
% 
% 功能描述:
%   这是OFDM接收机的主控制函数，协调各个模块完成完整的接收处理流程
%   采用模块化设计，便于维护和扩展
%
% 输入参数:
%   simSettings    - 仿真设置参数结构体
%   received_signal - 接收到的信号 [I; Q] 或复数信号
%   enable_kalman  - 是否启用卡尔曼滤波器，默认为true
%
% 输出参数:
%   output - 包含所有处理结果的结构体，包括：
%            * 跟踪环路输出
%            * 卡尔曼滤波器结果
%            * OFDM解调数据
%            * 性能评估指标
%
% 处理流程:
%   1. 参数初始化和信号预处理
%   2. 本地参考信号生成
%   3. 跟踪环路初始化
%   4. 卡尔曼滤波器初始化（可选）
%   5. 主跟踪循环处理
%   6. OFDM解调
%   7. 性能评估
%   8. 结果可视化
%
% 作者: OFDM接收机开发
% 日期: 2025年7月
% 版本: 1.0

%% 输入参数检查和默认值设置
if nargin < 3
    enable_kalman = true; % 默认启用卡尔曼滤波器
end

if nargin < 2
    error('必须提供仿真设置和接收信号');
end

fprintf('=== OFDM接收机启动 ===\n');
fprintf('模块化OFDM接收机 v1.0\n');
fprintf('处理开始时间: %s\n', datestr(now));

%% 模块1: 参数初始化和信号预处理
fprintf('\n[模块1] 参数初始化和信号预处理...\n');
[processed_signal, init_params] = ofdm_receiver_init(simSettings, received_signal);

fprintf('- 信号长度: %d 采样点\n', length(processed_signal.I));
fprintf('- 采样频率: %.2f MHz\n', init_params.fs/1e6);
fprintf('- 码频率: %.2f kHz\n', init_params.fp/1e3);
fprintf('- 频带数量: %d\n', init_params.numBand);

%% 模块2: 本地参考信号生成
fprintf('\n[模块2] 本地参考信号生成...\n');
local_signals = ofdm_generate_local_signals(simSettings);

fprintf('- 本地码长度: %d\n', length(simSettings.code));
fprintf('- OFDM符号数: %d\n', simSettings.nSymbol);
fprintf('- 子载波数: %d\n', simSettings.Nu);

%% 模块3: 跟踪环路初始化
fprintf('\n[模块3] 跟踪环路初始化...\n');
tracking_params = ofdm_tracking_init(simSettings);

fprintf('- DLL噪声带宽: %.2f Hz\n', tracking_params.dllNoiseBandwidth);
fprintf('- PLL噪声带宽: %.2f Hz\n', tracking_params.pllNoiseBandwidth);
fprintf('- 早晚相关器间距: %.2f\n', tracking_params.earlyLateSpc);

%% 模块4: 卡尔曼滤波器初始化（可选）
kf_params = [];
if enable_kalman
    fprintf('\n[模块4] 卡尔曼滤波器初始化...\n');
    kf_params = ofdm_kalman_init(simSettings);
    fprintf('- 标准卡尔曼滤波器: 已初始化\n');
    fprintf('- 扩展卡尔曼滤波器: 已初始化\n');
    fprintf('- 无迹卡尔曼滤波器: 已初始化\n');
else
    fprintf('\n[模块4] 跳过卡尔曼滤波器初始化\n');
end

%% 模块5: 主跟踪循环处理
fprintf('\n[模块5] 主跟踪循环处理...\n');
tracking_output = ofdm_tracking_loop(processed_signal, local_signals, ...
                                   tracking_params, kf_params, simSettings);

fprintf('- 处理样本数: %d\n', length(tracking_output.OutCarrFreq));
fprintf('- 平均信噪比: %.2f dB\n', tracking_output.Average_SNR);

%% 模块6: OFDM解调
fprintf('\n[模块6] OFDM解调处理...\n');
demod_output = ofdm_demodulation(tracking_output, simSettings);

if ~isempty(demod_output.DemodulatedData)
    fprintf('- 解调符号数: %d\n', length(demod_output.DemodulatedData));
    fprintf('- 资源网格大小: %dx%dx%d\n', size(demod_output.ResourceGrid));
else
    fprintf('- OFDM解调失败\n');
end

%% 模块7: 性能评估
fprintf('\n[模块7] 性能评估...\n');
performance_metrics = ofdm_performance_evaluation(tracking_output, demod_output, simSettings);

fprintf('- 载波频率稳定性: %.4f Hz\n', performance_metrics.CarrFreqStd);
fprintf('- 码频率稳定性: %.4f Hz\n', performance_metrics.CodeFreqStd);
fprintf('- 跟踪收敛时间: %.2f s\n', performance_metrics.ConvergenceTime);

%% 模块8: 结果整合和输出
fprintf('\n[模块8] 结果整合...\n');

% 整合所有模块的输出结果
output = struct();

% 跟踪结果
output.Tracking = tracking_output;

% 解调结果
output.Demodulation = demod_output;

% 性能指标
output.Performance = performance_metrics;

% 处理参数
output.Parameters = struct();
output.Parameters.Init = init_params;
output.Parameters.Tracking = tracking_params;
output.Parameters.KalmanFilter = kf_params;
output.Parameters.SimSettings = simSettings;

% 处理时间戳
output.ProcessingInfo = struct();
output.ProcessingInfo.StartTime = datestr(now);
output.ProcessingInfo.EnableKalman = enable_kalman;
output.ProcessingInfo.SignalLength = length(processed_signal.I);

fprintf('- 结果整合完成\n');

%% 模块9: 结果可视化
fprintf('\n[模块9] 结果可视化...\n');
ofdm_visualization(output, enable_kalman);

fprintf('- 生成图表数量: %d\n', 5 + 2*enable_kalman);

%% 处理完成
fprintf('\n=== OFDM接收机处理完成 ===\n');
fprintf('总处理时间: %.2f 秒\n', toc);
fprintf('处理结束时间: %s\n', datestr(now));

% 显示关键性能指标摘要
fprintf('\n=== 性能指标摘要 ===\n');
fprintf('平均信噪比: %.2f dB\n', output.Tracking.Average_SNR);
fprintf('载波频率标准差: %.4f Hz\n', output.Performance.CarrFreqStd);
fprintf('码频率标准差: %.4f Hz\n', output.Performance.CodeFreqStd);
fprintf('跟踪收敛时间: %.2f s\n', output.Performance.ConvergenceTime);

if enable_kalman
    fprintf('\n=== 卡尔曼滤波器性能 ===\n');
    fprintf('标准KF改善度: %.2f%%\n', output.Performance.KF_Standard_Improvement);
    fprintf('扩展KF改善度: %.2f%%\n', output.Performance.KF_Extended_Improvement);
    fprintf('无迹KF改善度: %.2f%%\n', output.Performance.KF_Unscented_Improvement);
end

end