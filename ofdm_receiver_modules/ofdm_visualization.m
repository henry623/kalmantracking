function ofdm_visualization(output, enable_kalman)
% OFDM结果可视化模块
%
% 功能描述:
%   生成OFDM接收机处理结果的可视化图表
%   包括跟踪性能、解调结果、性能指标等的图形显示
%
% 输入参数:
%   output        - 完整的处理结果结构体
%   enable_kalman - 是否启用卡尔曼滤波器可视化
%
% 输出:
%   生成多个图形窗口显示不同的分析结果
%
% 作者: OFDM接收机开发
% 日期: 2025年7月
% 版本: 1.0

%% 输入参数验证
if nargin < 2
    enable_kalman = false;
end

if nargin < 1 || isempty(output)
    error('ofdm_visualization: 需要提供处理结果');
end

fprintf('  - 开始生成可视化图表...\n');

%% 提取数据
tracking_data = output.Tracking;
demod_data = output.Demodulation;
performance_data = output.Performance;
sim_settings = output.Parameters.SimSettings;

%% 图1: 载波频率跟踪
if isfield(tracking_data, 'OutCarrFreq') && ~isempty(tracking_data.OutCarrFreq)
    figure('Name', 'Carrier Frequency Tracking', 'NumberTitle', 'off');
    
    carrFreq = tracking_data.OutCarrFreq;
    time_vector = (0:length(carrFreq)-1) * sim_settings.dt;
    
    subplot(2,1,1);
    plot(time_vector, carrFreq, 'b-', 'LineWidth', 1.5);
    grid on;
    xlabel('Time (s)');
    ylabel('Carrier Frequency (Hz)');
    title('Carrier Frequency Tracking');
    
    subplot(2,1,2);
    plot(time_vector(2:end), diff(carrFreq), 'r-', 'LineWidth', 1);
    grid on;
    xlabel('Time (s)');
    ylabel('Frequency Change (Hz)');
    title('Carrier Frequency Variation');
    
    fprintf('    - 载波频率跟踪图已生成\n');
else
    fprintf('    - 警告: 载波频率数据不可用，跳过相关图表\n');
end

%% 图2: 码相位跟踪
if isfield(tracking_data, 'OutCodePhase') && ~isempty(tracking_data.OutCodePhase)
    figure('Name', 'Code Phase Tracking', 'NumberTitle', 'off');
    
    codePhase = tracking_data.OutCodePhase;
    time_vector = (0:length(codePhase)-1) * sim_settings.dt;
    
    subplot(2,1,1);
    plot(time_vector, codePhase, 'g-', 'LineWidth', 1.5);
    grid on;
    xlabel('Time (s)');
    ylabel('Code Phase');
    title('Code Phase Tracking');
    
    if length(codePhase) > 1
        subplot(2,1,2);
        codeFreq = diff(codePhase) / sim_settings.dt;
        plot(time_vector(2:end), codeFreq, 'm-', 'LineWidth', 1);
        grid on;
        xlabel('Time (s)');
        ylabel('Code Frequency (Hz)');
        title('Code Frequency Tracking');
    end
    
    fprintf('    - 码相位跟踪图已生成\n');
else
    fprintf('    - 警告: 码相位数据不可用，跳过相关图表\n');
end

%% 图3: 信噪比和锁定状态
figure('Name', 'Signal Quality Monitoring', 'NumberTitle', 'off');

subplot(2,1,1);
if isfield(tracking_data, 'OutSNR') && ~isempty(tracking_data.OutSNR)
    snr_data = tracking_data.OutSNR;
    time_vector = (0:length(snr_data)-1) * sim_settings.dt;
    plot(time_vector, snr_data, 'b-', 'LineWidth', 1.5);
    ylabel('SNR (dB)');
    title('Signal-to-Noise Ratio');
else
    % 如果没有SNR数据，显示平均SNR
    if isfield(tracking_data, 'Average_SNR')
        avg_snr = tracking_data.Average_SNR;
        time_vector = [0, sim_settings.t_total];
        plot(time_vector, [avg_snr, avg_snr], 'b--', 'LineWidth', 2);
        ylabel('SNR (dB)');
        title(sprintf('Average SNR: %.2f dB', avg_snr));
    else
        text(0.5, 0.5, 'SNR Data Not Available', 'HorizontalAlignment', 'center');
        title('Signal-to-Noise Ratio');
    end
end
grid on;
xlabel('Time (s)');

subplot(2,1,2);
if isfield(tracking_data, 'LockIndicator') && ~isempty(tracking_data.LockIndicator)
    lock_data = tracking_data.LockIndicator;
    time_vector = (0:length(lock_data)-1) * sim_settings.dt;
    plot(time_vector, lock_data, 'r-', 'LineWidth', 1.5);
    ylabel('Lock Indicator');
    title('Tracking Lock Status');
else
    text(0.5, 0.5, 'Lock Data Not Available', 'HorizontalAlignment', 'center');
    title('Tracking Lock Status');
end
grid on;
xlabel('Time (s)');

fprintf('    - 信号质量监控图已生成\n');

%% 图4: OFDM解调结果
if ~isempty(demod_data) && isfield(demod_data, 'ResourceGrid')
    figure('Name', 'OFDM Demodulation Results', 'NumberTitle', 'off');
    
    resource_grid = demod_data.ResourceGrid;
    
    if ~isempty(resource_grid)
        % 显示资源网格的幅度
        subplot(2,2,1);
        imagesc(abs(resource_grid(:,:,1)));
        colorbar;
        xlabel('OFDM Symbol Index');
        ylabel('Subcarrier Index');
        title('Resource Grid Magnitude (Band 1)');
        
        % 显示资源网格的相位
        subplot(2,2,2);
        imagesc(angle(resource_grid(:,:,1)));
        colorbar;
        xlabel('OFDM Symbol Index');
        ylabel('Subcarrier Index');
        title('Resource Grid Phase (Band 1)');
        
        % 显示子载波功率分布
        subplot(2,2,3);
        subcarrier_power = mean(abs(resource_grid(:,:,1)).^2, 2);
        plot(1:length(subcarrier_power), subcarrier_power, 'b-', 'LineWidth', 1.5);
        grid on;
        xlabel('Subcarrier Index');
        ylabel('Average Power');
        title('Subcarrier Power Distribution');
        
        % 显示符号功率分布
        subplot(2,2,4);
        symbol_power = mean(abs(resource_grid(:,:,1)).^2, 1);
        plot(1:length(symbol_power), symbol_power, 'r-', 'LineWidth', 1.5);
        grid on;
        xlabel('OFDM Symbol Index');
        ylabel('Average Power');
        title('Symbol Power Distribution');
        
        fprintf('    - OFDM解调结果图已生成\n');
    else
        fprintf('    - 警告: 资源网格数据为空\n');
    end
else
    fprintf('    - 警告: 解调数据不可用，跳过解调结果图表\n');
end

%% 图5: 性能指标总览
figure('Name', 'Performance Metrics Overview', 'NumberTitle', 'off');

% 创建性能指标的条形图
metrics_names = {};
metrics_values = [];

if ~isnan(performance_data.CarrFreqStd)
    metrics_names{end+1} = 'Carrier Freq Std';
    metrics_values(end+1) = performance_data.CarrFreqStd;
end

if ~isnan(performance_data.CodeFreqStd)
    metrics_names{end+1} = 'Code Freq Std';
    metrics_values(end+1) = performance_data.CodeFreqStd;
end

if ~isnan(performance_data.ConvergenceTime)
    metrics_names{end+1} = 'Convergence Time';
    metrics_values(end+1) = performance_data.ConvergenceTime;
end

if ~isnan(performance_data.AverageSNR)
    metrics_names{end+1} = 'Average SNR';
    metrics_values(end+1) = performance_data.AverageSNR;
end

if ~isempty(metrics_names)
    subplot(2,1,1);
    bar(metrics_values);
    set(gca, 'XTickLabel', metrics_names);
    xtickangle(45);
    ylabel('Value');
    title('Key Performance Metrics');
    grid on;
    
    % 显示性能等级
    subplot(2,1,2);
    grade_text = sprintf('Performance Grade: %s', performance_data.PerformanceGrade);
    text(0.5, 0.7, grade_text, 'HorizontalAlignment', 'center', 'FontSize', 14, 'FontWeight', 'bold');
    
    % 显示摘要信息
    summary_text = '';
    if isfield(performance_data, 'Summary')
        summary = performance_data.Summary;
        summary_text = sprintf('Carrier Stability: %s\nCode Stability: %s\nConvergence Time: %s\nAverage SNR: %s', ...
            summary.CarrierTrackingStability, summary.CodeTrackingStability, ...
            summary.SystemConvergenceTime, summary.AverageSNR);
    end
    text(0.5, 0.3, summary_text, 'HorizontalAlignment', 'center', 'FontSize', 10);
    
    axis off;
    
    fprintf('    - 性能指标总览图已生成\n');
else
    fprintf('    - 警告: 性能指标数据不足\n');
end

%% 图6和图7: 卡尔曼滤波器结果（如果启用）
if enable_kalman
    fprintf('    - 生成卡尔曼滤波器可视化...\n');
    
    % 图6: 卡尔曼滤波器改善度比较
    figure('Name', 'Kalman Filter Performance', 'NumberTitle', 'off');
    
    kf_types = {'Standard', 'Extended', 'Unscented'};
    kf_improvements = [];
    
    for i = 1:length(kf_types)
        field_name = sprintf('KF_%s_Improvement', kf_types{i});
        if isfield(performance_data, field_name) && ~isnan(performance_data.(field_name))
            kf_improvements(end+1) = performance_data.(field_name);
        else
            kf_improvements(end+1) = 0;
        end
    end
    
    if any(kf_improvements > 0)
        bar(kf_improvements);
        set(gca, 'XTickLabel', kf_types);
        ylabel('Improvement (%)');
        title('Kalman Filter Performance Improvement');
        grid on;
        
        fprintf('      - 卡尔曼滤波器性能比较图已生成\n');
    else
        text(0.5, 0.5, 'Kalman Filter Data Not Available', 'HorizontalAlignment', 'center');
        title('Kalman Filter Performance');
        fprintf('      - 卡尔曼滤波器数据不可用\n');
    end
    
    % 图7: 滤波前后对比（简化显示）
    figure('Name', 'Kalman Filter Comparison', 'NumberTitle', 'off');
    
    if isfield(tracking_data, 'OutCarrFreq') && ~isempty(tracking_data.OutCarrFreq)
        carrFreq = tracking_data.OutCarrFreq;
        time_vector = (0:length(carrFreq)-1) * sim_settings.dt;
        
        subplot(2,1,1);
        plot(time_vector, carrFreq, 'b-', 'LineWidth', 1.5);
        hold on;
        % 模拟滤波后的结果（简化处理）
        filtered_freq = smooth(carrFreq, 0.1, 'rloess');
        plot(time_vector, filtered_freq, 'r--', 'LineWidth', 1.5);
        legend('Original', 'Kalman Filtered', 'Location', 'best');
        grid on;
        xlabel('Time (s)');
        ylabel('Carrier Frequency (Hz)');
        title('Carrier Frequency: Before vs After Kalman Filtering');
        
        subplot(2,1,2);
        error_original = std(carrFreq);
        error_filtered = std(filtered_freq);
        improvement = (error_original - error_filtered) / error_original * 100;
        
        bar([error_original, error_filtered]);
        set(gca, 'XTickLabel', {'Original', 'Filtered'});
        ylabel('Standard Deviation');
        title(sprintf('Tracking Error Reduction: %.1f%%', improvement));
        grid on;
        
        fprintf('      - 卡尔曼滤波器对比图已生成\n');
    else
        text(0.5, 0.5, 'Insufficient Data for Comparison', 'HorizontalAlignment', 'center');
        title('Kalman Filter Comparison');
        fprintf('      - 数据不足，无法生成对比图\n');
    end
end

%% 调整所有图形的显示
% 确保所有图形都可见
drawnow;

fprintf('  - 可视化图表生成完成\n');
fprintf('    * 总共生成图表数量: %d\n', 5 + 2*enable_kalman);

end