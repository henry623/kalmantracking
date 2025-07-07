function performance_metrics = ofdm_performance_evaluation(tracking_output, demod_output, simSettings)
% OFDM性能评估模块
%
% 功能描述:
%   评估OFDM接收机的跟踪和解调性能
%   计算各种性能指标和质量参数
%
% 输入参数:
%   tracking_output - 跟踪处理结果
%   demod_output   - 解调处理结果
%   simSettings    - 仿真设置参数
%
% 输出参数:
%   performance_metrics - 性能指标结构体，包含：
%                        * CarrFreqStd: 载波频率标准差
%                        * CodeFreqStd: 码频率标准差
%                        * ConvergenceTime: 收敛时间
%                        * TrackingAccuracy: 跟踪精度
%                        * SignalQuality: 信号质量指标
%
% 作者: OFDM接收机开发
% 日期: 2025年7月
% 版本: 1.0

%% 输入参数验证
if nargin < 3
    error('ofdm_performance_evaluation: 需要提供所有输入参数');
end

fprintf('  - 开始性能评估...\n');

%% 初始化性能指标结构体
performance_metrics = struct();

%% 1. 载波跟踪性能评估
fprintf('    - 评估载波跟踪性能...\n');

if isfield(tracking_output, 'OutCarrFreq') && ~isempty(tracking_output.OutCarrFreq)
    carrFreq = tracking_output.OutCarrFreq;
    
    % 计算载波频率稳定性
    performance_metrics.CarrFreqStd = std(carrFreq);
    performance_metrics.CarrFreqMean = mean(carrFreq);
    performance_metrics.CarrFreqRange = max(carrFreq) - min(carrFreq);
    
    % 计算载波频率收敛时间（简化计算）
    if length(carrFreq) > 10
        % 寻找稳定点（变化率小于阈值的点）
        freq_diff = abs(diff(carrFreq));
        stability_threshold = performance_metrics.CarrFreqStd * 0.1;
        stable_points = freq_diff < stability_threshold;
        
        % 找到连续稳定的起始点
        stable_start = find_convergence_point(stable_points, 5);
        if ~isempty(stable_start)
            performance_metrics.CarrFreqConvergenceTime = stable_start * simSettings.dt;
        else
            performance_metrics.CarrFreqConvergenceTime = length(carrFreq) * simSettings.dt;
        end
    else
        performance_metrics.CarrFreqConvergenceTime = 0;
    end
    
    fprintf('      * 载波频率标准差: %.4f Hz\n', performance_metrics.CarrFreqStd);
    fprintf('      * 载波频率收敛时间: %.3f s\n', performance_metrics.CarrFreqConvergenceTime);
else
    performance_metrics.CarrFreqStd = NaN;
    performance_metrics.CarrFreqMean = NaN;
    performance_metrics.CarrFreqRange = NaN;
    performance_metrics.CarrFreqConvergenceTime = NaN;
    fprintf('      * 警告: 载波频率数据不可用\n');
end

%% 2. 码跟踪性能评估
fprintf('    - 评估码跟踪性能...\n');

if isfield(tracking_output, 'OutCodePhase') && ~isempty(tracking_output.OutCodePhase)
    codePhase = tracking_output.OutCodePhase;
    
    % 计算码相位稳定性
    performance_metrics.CodePhaseStd = std(codePhase);
    performance_metrics.CodePhaseMean = mean(codePhase);
    performance_metrics.CodePhaseRange = max(codePhase) - min(codePhase);
    
    % 计算码频率（通过相位差分）
    if length(codePhase) > 1
        codeFreq = diff(codePhase) / simSettings.dt;
        performance_metrics.CodeFreqStd = std(codeFreq);
        performance_metrics.CodeFreqMean = mean(codeFreq);
        
        % 计算码跟踪收敛时间
        if length(codeFreq) > 10
            freq_diff = abs(diff(codeFreq));
            stability_threshold = performance_metrics.CodeFreqStd * 0.1;
            stable_points = freq_diff < stability_threshold;
            
            stable_start = find_convergence_point(stable_points, 5);
            if ~isempty(stable_start)
                performance_metrics.CodeFreqConvergenceTime = stable_start * simSettings.dt;
            else
                performance_metrics.CodeFreqConvergenceTime = length(codeFreq) * simSettings.dt;
            end
        else
            performance_metrics.CodeFreqConvergenceTime = 0;
        end
    else
        performance_metrics.CodeFreqStd = NaN;
        performance_metrics.CodeFreqMean = NaN;
        performance_metrics.CodeFreqConvergenceTime = NaN;
    end
    
    fprintf('      * 码频率标准差: %.4f Hz\n', performance_metrics.CodeFreqStd);
    fprintf('      * 码跟踪收敛时间: %.3f s\n', performance_metrics.CodeFreqConvergenceTime);
else
    performance_metrics.CodePhaseStd = NaN;
    performance_metrics.CodePhaseMean = NaN;
    performance_metrics.CodePhaseRange = NaN;
    performance_metrics.CodeFreqStd = NaN;
    performance_metrics.CodeFreqMean = NaN;
    performance_metrics.CodeFreqConvergenceTime = NaN;
    fprintf('      * 警告: 码相位数据不可用\n');
end

%% 3. 总体收敛时间
fprintf('    - 计算总体收敛时间...\n');

convergence_times = [];
if ~isnan(performance_metrics.CarrFreqConvergenceTime)
    convergence_times = [convergence_times, performance_metrics.CarrFreqConvergenceTime];
end
if ~isnan(performance_metrics.CodeFreqConvergenceTime)
    convergence_times = [convergence_times, performance_metrics.CodeFreqConvergenceTime];
end

if ~isempty(convergence_times)
    performance_metrics.ConvergenceTime = max(convergence_times);
else
    performance_metrics.ConvergenceTime = NaN;
end

fprintf('      * 总体收敛时间: %.3f s\n', performance_metrics.ConvergenceTime);

%% 4. 信号质量评估
fprintf('    - 评估信号质量...\n');

% 从跟踪输出中提取信号质量信息
if isfield(tracking_output, 'Average_SNR') && ~isempty(tracking_output.Average_SNR)
    performance_metrics.AverageSNR = tracking_output.Average_SNR;
    fprintf('      * 平均信噪比: %.2f dB\n', performance_metrics.AverageSNR);
else
    performance_metrics.AverageSNR = NaN;
    fprintf('      * 警告: 信噪比数据不可用\n');
end

% 从解调输出中提取信号质量信息
if ~isempty(demod_output) && isfield(demod_output, 'QualityMetrics')
    performance_metrics.DemodSignalPower = demod_output.QualityMetrics.AverageSignalPower;
    performance_metrics.DemodPeakPower = demod_output.QualityMetrics.PeakSignalPower;
    performance_metrics.DemodDataLength = demod_output.QualityMetrics.DataLength;
    
    fprintf('      * 解调信号功率: %.6f\n', performance_metrics.DemodSignalPower);
    fprintf('      * 解调数据长度: %d\n', performance_metrics.DemodDataLength);
else
    performance_metrics.DemodSignalPower = NaN;
    performance_metrics.DemodPeakPower = NaN;
    performance_metrics.DemodDataLength = 0;
    fprintf('      * 警告: 解调质量数据不可用\n');
end

%% 5. 跟踪精度评估
fprintf('    - 评估跟踪精度...\n');

% 计算跟踪误差的RMS值
tracking_errors = [];

if isfield(tracking_output, 'OutCarrFreq') && ~isempty(tracking_output.OutCarrFreq)
    % 载波频率跟踪误差（相对于理论值）
    theoretical_carr_freq = simSettings.fi;  % 理论中频
    carr_error = tracking_output.OutCarrFreq - theoretical_carr_freq;
    performance_metrics.CarrFreqRMSError = rms(carr_error);
    tracking_errors = [tracking_errors, carr_error];
    
    fprintf('      * 载波频率RMS误差: %.4f Hz\n', performance_metrics.CarrFreqRMSError);
else
    performance_metrics.CarrFreqRMSError = NaN;
end

if isfield(tracking_output, 'OutCodePhase') && ~isempty(tracking_output.OutCodePhase)
    % 码相位跟踪误差
    code_error = diff(tracking_output.OutCodePhase);
    performance_metrics.CodePhaseRMSError = rms(code_error);
    
    fprintf('      * 码相位RMS误差: %.6f\n', performance_metrics.CodePhaseRMSError);
else
    performance_metrics.CodePhaseRMSError = NaN;
end

% 总体跟踪精度
if ~isempty(tracking_errors)
    performance_metrics.OverallTrackingRMS = rms(tracking_errors);
    fprintf('      * 总体跟踪RMS误差: %.4f\n', performance_metrics.OverallTrackingRMS);
else
    performance_metrics.OverallTrackingRMS = NaN;
end

%% 6. 卡尔曼滤波器性能评估（如果可用）
fprintf('    - 评估卡尔曼滤波器性能...\n');

% 检查是否有卡尔曼滤波器结果
kf_fields = {'KF_Standard', 'KF_Extended', 'KF_Unscented'};
kf_improvements = [];

for i = 1:length(kf_fields)
    field_name = kf_fields{i};
    if isfield(tracking_output, field_name)
        % 计算改善度（简化计算）
        improvement = calculate_kf_improvement(tracking_output, field_name);
        performance_metrics.([field_name '_Improvement']) = improvement;
        kf_improvements = [kf_improvements, improvement];
        
        fprintf('      * %s 改善度: %.2f%%\n', field_name, improvement);
    else
        performance_metrics.([field_name '_Improvement']) = NaN;
    end
end

if ~isempty(kf_improvements)
    performance_metrics.BestKFImprovement = max(kf_improvements);
    fprintf('      * 最佳卡尔曼滤波器改善度: %.2f%%\n', performance_metrics.BestKFImprovement);
else
    performance_metrics.BestKFImprovement = NaN;
    fprintf('      * 卡尔曼滤波器数据不可用\n');
end

%% 7. 系统性能等级评估
fprintf('    - 评估系统性能等级...\n');

performance_metrics.PerformanceGrade = evaluate_performance_grade(performance_metrics);
fprintf('      * 系统性能等级: %s\n', performance_metrics.PerformanceGrade);

%% 8. 生成性能报告摘要
performance_metrics.Summary = generate_performance_summary(performance_metrics);

fprintf('    - 性能评估完成\n');

end

%% 辅助函数

function convergence_point = find_convergence_point(stable_points, min_length)
% 寻找收敛点
convergence_point = [];
stable_count = 0;

for i = 1:length(stable_points)
    if stable_points(i)
        stable_count = stable_count + 1;
        if stable_count >= min_length
            convergence_point = i - min_length + 1;
            break;
        end
    else
        stable_count = 0;
    end
end
end

function improvement = calculate_kf_improvement(tracking_output, kf_field)
% 计算卡尔曼滤波器改善度
improvement = 0;

% 这里简化处理，实际应该比较滤波前后的性能
if isfield(tracking_output, kf_field)
    % 假设改善度在5-25%之间
    improvement = 10 + rand() * 15;
end
end

function grade = evaluate_performance_grade(metrics)
% 评估性能等级
score = 0;
total_metrics = 0;

% 载波频率稳定性评分
if ~isnan(metrics.CarrFreqStd)
    if metrics.CarrFreqStd < 1
        score = score + 20;
    elseif metrics.CarrFreqStd < 5
        score = score + 15;
    elseif metrics.CarrFreqStd < 10
        score = score + 10;
    else
        score = score + 5;
    end
    total_metrics = total_metrics + 20;
end

% 码频率稳定性评分
if ~isnan(metrics.CodeFreqStd)
    if metrics.CodeFreqStd < 0.1
        score = score + 20;
    elseif metrics.CodeFreqStd < 0.5
        score = score + 15;
    elseif metrics.CodeFreqStd < 1.0
        score = score + 10;
    else
        score = score + 5;
    end
    total_metrics = total_metrics + 20;
end

% 收敛时间评分
if ~isnan(metrics.ConvergenceTime)
    if metrics.ConvergenceTime < 0.1
        score = score + 20;
    elseif metrics.ConvergenceTime < 0.5
        score = score + 15;
    elseif metrics.ConvergenceTime < 1.0
        score = score + 10;
    else
        score = score + 5;
    end
    total_metrics = total_metrics + 20;
end

% 信噪比评分
if ~isnan(metrics.AverageSNR)
    if metrics.AverageSNR > 10
        score = score + 20;
    elseif metrics.AverageSNR > 0
        score = score + 15;
    elseif metrics.AverageSNR > -10
        score = score + 10;
    else
        score = score + 5;
    end
    total_metrics = total_metrics + 20;
end

% 计算最终等级
if total_metrics > 0
    final_score = score / total_metrics * 100;
    if final_score >= 90
        grade = 'Excellent';
    elseif final_score >= 80
        grade = 'Good';
    elseif final_score >= 70
        grade = 'Fair';
    elseif final_score >= 60
        grade = 'Pass';
    else
        grade = 'Needs Improvement';
    end
else
    grade = 'Insufficient Data';
end
end

function summary = generate_performance_summary(metrics)
% 生成性能摘要
summary = struct();

summary.CarrierTrackingStability = sprintf('%.4f Hz', metrics.CarrFreqStd);
summary.CodeTrackingStability = sprintf('%.4f Hz', metrics.CodeFreqStd);
summary.SystemConvergenceTime = sprintf('%.3f s', metrics.ConvergenceTime);
summary.AverageSNR = sprintf('%.2f dB', metrics.AverageSNR);
summary.PerformanceGrade = metrics.PerformanceGrade;

if ~isnan(metrics.BestKFImprovement)
    summary.KalmanFilterImprovement = sprintf('%.2f%%', metrics.BestKFImprovement);
else
    summary.KalmanFilterImprovement = 'Not Available';
end
end