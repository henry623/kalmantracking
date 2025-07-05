function ofdm_results_output(all_tracking_results, simSettings, output_options)
% OFDM接收机结果输出和可视化模块
%
% 功能描述:
%   处理和可视化OFDM接收机的跟踪结果，生成性能分析图表和数据报告
%   支持多种输出格式和可视化选项
%
% 输入参数:
%   all_tracking_results - 所有跟踪循环的结果数组
%   simSettings         - 仿真设置参数结构体
%   output_options      - 输出选项结构体，包含：
%                        * save_figures: 是否保存图形文件
%                        * show_plots: 是否显示图形
%                        * save_data: 是否保存数据文件
%                        * output_dir: 输出目录
%                        * file_prefix: 文件名前缀
%
% 输出:
%   - 跟踪性能图表（相关器输出、误差曲线、频率跟踪等）
%   - 统计分析报告
%   - 数据文件（可选）
%
% 处理流程:
%   1. 数据提取和预处理
%   2. 跟踪性能分析
%   3. 可视化图表生成
%   4. 统计报告生成
%   5. 文件保存（可选）
%
% 作者: OFDM接收机开发团队
% 日期: 2025年1月
% 版本: 1.0

%% 1. 输入参数验证和默认设置
if nargin < 2
    error('ofdm_results_output: 至少需要提供跟踪结果和仿真设置');
end

if nargin < 3 || isempty(output_options)
    % 使用默认输出选项
    output_options = struct();
    output_options.save_figures = false;
    output_options.show_plots = true;
    output_options.save_data = false;
    output_options.output_dir = './results';
    output_options.file_prefix = 'ofdm_tracking';
end

% 检查必要字段
if ~isfield(output_options, 'save_figures'), output_options.save_figures = false; end
if ~isfield(output_options, 'show_plots'), output_options.show_plots = true; end
if ~isfield(output_options, 'save_data'), output_options.save_data = false; end
if ~isfield(output_options, 'output_dir'), output_options.output_dir = './results'; end
if ~isfield(output_options, 'file_prefix'), output_options.file_prefix = 'ofdm_tracking'; end

fprintf('开始生成OFDM跟踪结果输出...\n');

% 创建输出目录
if output_options.save_figures || output_options.save_data
    if ~exist(output_options.output_dir, 'dir')
        mkdir(output_options.output_dir);
        fprintf('  - 创建输出目录: %s\n', output_options.output_dir);
    end
end

%% 2. 数据提取和预处理
fprintf('  - 提取和预处理跟踪数据...\n');

num_loops = length(all_tracking_results);
if num_loops == 0
    error('ofdm_results_output: 没有跟踪结果数据');
end

% 初始化数据数组
time_vector = zeros(1, num_loops);
code_errors = zeros(1, num_loops);
carrier_errors = zeros(1, num_loops);
code_ncos = zeros(1, num_loops);
carrier_ncos = zeros(1, num_loops);
code_freqs = zeros(1, num_loops);
carrier_freqs = zeros(1, num_loops);
snr_estimates = zeros(1, num_loops);
lock_indicators = zeros(1, num_loops);
signal_powers = zeros(1, num_loops);
noise_powers = zeros(1, num_loops);

% 相关器输出
I_E_values = zeros(1, num_loops);
Q_E_values = zeros(1, num_loops);
I_P_values = zeros(1, num_loops);
Q_P_values = zeros(1, num_loops);
I_L_values = zeros(1, num_loops);
Q_L_values = zeros(1, num_loops);

% 跟踪模式
tracking_modes = cell(1, num_loops);

% 提取数据
dt = simSettings.dt;
for i = 1:num_loops
    result = all_tracking_results{i};
    
    time_vector(i) = (i-1) * dt;
    
    % 鉴相器输出
    code_errors(i) = result.discriminators.code_error;
    carrier_errors(i) = result.discriminators.carrier_error;
    
    % 环路滤波器输出
    code_ncos(i) = result.loop_filters.code_nco;
    carrier_ncos(i) = result.loop_filters.carrier_nco;
    
    % 频率跟踪
    code_freqs(i) = result.updated_states.codeFreq;
    carrier_freqs(i) = result.updated_states.carrFreq;
    
    % 性能指标
    snr_estimates(i) = result.performance.snr_estimate;
    lock_indicators(i) = result.performance.lock_indicator;
    signal_powers(i) = result.performance.signal_power;
    noise_powers(i) = result.performance.noise_power;
    
    % 相关器输出
    I_E_values(i) = result.correlators.average.I_E;
    Q_E_values(i) = result.correlators.average.Q_E;
    I_P_values(i) = result.correlators.average.I_P;
    Q_P_values(i) = result.correlators.average.Q_P;
    I_L_values(i) = result.correlators.average.I_L;
    Q_L_values(i) = result.correlators.average.Q_L;
    
    % 跟踪模式
    tracking_modes{i} = result.processing_info.tracking_mode;
end

fprintf('    * 提取了 %d 个跟踪循环的数据\n', num_loops);

%% 3. 跟踪性能分析
fprintf('  - 分析跟踪性能...\n');

% 计算统计量
stats = struct();

% 码跟踪统计
stats.code_error_mean = mean(code_errors);
stats.code_error_std = std(code_errors);
stats.code_error_rms = rms(code_errors);

% 载波跟踪统计
stats.carrier_error_mean = mean(carrier_errors);
stats.carrier_error_std = std(carrier_errors);
stats.carrier_error_rms = rms(carrier_errors);

% 频率跟踪统计
stats.code_freq_mean = mean(code_freqs);
stats.code_freq_std = std(code_freqs);
stats.carrier_freq_mean = mean(carrier_freqs);
stats.carrier_freq_std = std(carrier_freqs);

% 信噪比统计
stats.snr_mean = mean(snr_estimates);
stats.snr_std = std(snr_estimates);
stats.snr_min = min(snr_estimates);
stats.snr_max = max(snr_estimates);

% 锁定性能
stats.lock_mean = mean(lock_indicators);
stats.lock_min = min(lock_indicators);
stats.lock_percentage = sum(lock_indicators > 0.5) / num_loops * 100;

fprintf('    * 码跟踪误差RMS: %.6f\n', stats.code_error_rms);
fprintf('    * 载波跟踪误差RMS: %.6f\n', stats.carrier_error_rms);
fprintf('    * 平均信噪比: %.2f dB\n', stats.snr_mean);
fprintf('    * 锁定百分比: %.1f%%\n', stats.lock_percentage);

%% 4. 生成可视化图表
fprintf('  - 生成可视化图表...\n');

% 设置图形属性
if output_options.show_plots
    figure_visibility = 'on';
else
    figure_visibility = 'off';
end

%% 4.1 跟踪误差时间序列图
fig1 = figure('Name', 'OFDM跟踪误差时间序列', 'Position', [100, 100, 1200, 800], 'Visible', figure_visibility);

subplot(2, 1, 1);
plot(time_vector, code_errors, 'b-', 'LineWidth', 1.5);
title('码跟踪误差');
xlabel('时间 (s)');
ylabel('码误差');
grid on;
legend('DLL误差', 'Location', 'best');

subplot(2, 1, 2);
plot(time_vector, carrier_errors, 'r-', 'LineWidth', 1.5);
title('载波跟踪误差');
xlabel('时间 (s)');
ylabel('载波误差');
grid on;
legend('PLL/FLL误差', 'Location', 'best');

sgtitle('OFDM信号跟踪误差时间序列');

%% 4.2 频率跟踪图
fig2 = figure('Name', 'OFDM频率跟踪', 'Position', [150, 150, 1200, 800], 'Visible', figure_visibility);

subplot(2, 1, 1);
plot(time_vector, (code_freqs - simSettings.fp)/1e3, 'b-', 'LineWidth', 1.5);
title('码频率跟踪');
xlabel('时间 (s)');
ylabel('码频率偏差 (kHz)');
grid on;
legend('码频率偏差', 'Location', 'best');

subplot(2, 1, 2);
plot(time_vector, (carrier_freqs - simSettings.fi)/1e3, 'r-', 'LineWidth', 1.5);
title('载波频率跟踪');
xlabel('时间 (s)');
ylabel('载波频率偏差 (kHz)');
grid on;
legend('载波频率偏差', 'Location', 'best');

sgtitle('OFDM信号频率跟踪');

%% 4.3 相关器输出图
fig3 = figure('Name', 'OFDM相关器输出', 'Position', [200, 200, 1200, 800], 'Visible', figure_visibility);

subplot(2, 2, 1);
plot(time_vector, I_E_values, 'g-', time_vector, I_P_values, 'b-', time_vector, I_L_values, 'r-', 'LineWidth', 1.5);
title('I通道相关器输出');
xlabel('时间 (s)');
ylabel('相关值');
legend('早相关', '准时相关', '晚相关', 'Location', 'best');
grid on;

subplot(2, 2, 2);
plot(time_vector, Q_E_values, 'g-', time_vector, Q_P_values, 'b-', time_vector, Q_L_values, 'r-', 'LineWidth', 1.5);
title('Q通道相关器输出');
xlabel('时间 (s)');
ylabel('相关值');
legend('早相关', '准时相关', '晚相关', 'Location', 'best');
grid on;

subplot(2, 2, 3);
plot(time_vector, sqrt(I_P_values.^2 + Q_P_values.^2), 'b-', 'LineWidth', 1.5);
title('信号幅度');
xlabel('时间 (s)');
ylabel('幅度');
grid on;

subplot(2, 2, 4);
plot(time_vector, atan2(Q_P_values, I_P_values)*180/pi, 'r-', 'LineWidth', 1.5);
title('信号相位');
xlabel('时间 (s)');
ylabel('相位 (度)');
grid on;

sgtitle('OFDM相关器输出分析');

%% 4.4 性能监控图
fig4 = figure('Name', 'OFDM性能监控', 'Position', [250, 250, 1200, 800], 'Visible', figure_visibility);

subplot(2, 2, 1);
plot(time_vector, snr_estimates, 'b-', 'LineWidth', 1.5);
title('信噪比估计');
xlabel('时间 (s)');
ylabel('SNR (dB)');
grid on;
legend(sprintf('平均SNR: %.2f dB', stats.snr_mean), 'Location', 'best');

subplot(2, 2, 2);
plot(time_vector, lock_indicators, 'r-', 'LineWidth', 1.5);
hold on;
plot([time_vector(1), time_vector(end)], [0.5, 0.5], 'k--', 'LineWidth', 1);
title('锁定指示器');
xlabel('时间 (s)');
ylabel('锁定指示器');
legend('锁定指示器', '锁定阈值', 'Location', 'best');
grid on;

subplot(2, 2, 3);
semilogy(time_vector, signal_powers, 'b-', time_vector, noise_powers, 'r-', 'LineWidth', 1.5);
title('信号和噪声功率');
xlabel('时间 (s)');
ylabel('功率');
legend('信号功率', '噪声功率', 'Location', 'best');
grid on;

subplot(2, 2, 4);
plot(time_vector, code_ncos, 'b-', time_vector, carrier_ncos, 'r-', 'LineWidth', 1.5);
title('NCO控制信号');
xlabel('时间 (s)');
ylabel('NCO输出');
legend('码NCO', '载波NCO', 'Location', 'best');
grid on;

sgtitle('OFDM接收机性能监控');

%% 4.5 跟踪模式切换图
fig5 = figure('Name', 'OFDM跟踪模式', 'Position', [300, 300, 1200, 400], 'Visible', figure_visibility);

% 将跟踪模式转换为数值
mode_values = zeros(1, num_loops);
for i = 1:num_loops
    if strcmp(tracking_modes{i}, 'FLL')
        mode_values(i) = 1;
    else
        mode_values(i) = 2;
    end
end

plot(time_vector, mode_values, 'ko-', 'LineWidth', 2, 'MarkerSize', 6);
title('跟踪模式切换');
xlabel('时间 (s)');
ylabel('跟踪模式');
ylim([0.5, 2.5]);
yticks([1, 2]);
yticklabels({'FLL', 'PLL'});
grid on;

%% 4.6 统计分析图
fig6 = figure('Name', 'OFDM统计分析', 'Position', [350, 350, 1200, 800], 'Visible', figure_visibility);

subplot(2, 2, 1);
histogram(code_errors, 30, 'Normalization', 'probability');
title('码跟踪误差分布');
xlabel('码误差');
ylabel('概率');
grid on;

subplot(2, 2, 2);
histogram(carrier_errors, 30, 'Normalization', 'probability');
title('载波跟踪误差分布');
xlabel('载波误差');
ylabel('概率');
grid on;

subplot(2, 2, 3);
histogram(snr_estimates, 30, 'Normalization', 'probability');
title('信噪比分布');
xlabel('SNR (dB)');
ylabel('概率');
grid on;

subplot(2, 2, 4);
% 绘制跟踪误差的功率谱密度
[psd_code, f_code] = pwelch(code_errors, [], [], [], 1/dt);
[psd_carrier, f_carrier] = pwelch(carrier_errors, [], [], [], 1/dt);
semilogy(f_code, psd_code, 'b-', f_carrier, psd_carrier, 'r-', 'LineWidth', 1.5);
title('跟踪误差功率谱密度');
xlabel('频率 (Hz)');
ylabel('PSD');
legend('码误差PSD', '载波误差PSD', 'Location', 'best');
grid on;

sgtitle('OFDM跟踪统计分析');

%% 5. 生成文本报告
fprintf('  - 生成统计报告...\n');

report_text = sprintf(['OFDM接收机跟踪性能报告\n' ...
                      '========================\n\n' ...
                      '仿真参数:\n' ...
                      '  码频率: %.2f kHz\n' ...
                      '  采样频率: %.2f MHz\n' ...
                      '  中频频率: %.2f MHz\n' ...
                      '  信噪比: %.2f dB\n' ...
                      '  跟踪循环数: %d\n' ...
                      '  总仿真时间: %.3f s\n\n' ...
                      '码跟踪性能:\n' ...
                      '  平均误差: %.6f\n' ...
                      '  误差标准差: %.6f\n' ...
                      '  误差RMS: %.6f\n' ...
                      '  频率标准差: %.3f Hz\n\n' ...
                      '载波跟踪性能:\n' ...
                      '  平均误差: %.6f\n' ...
                      '  误差标准差: %.6f\n' ...
                      '  误差RMS: %.6f\n' ...
                      '  频率标准差: %.3f Hz\n\n' ...
                      '信号质量:\n' ...
                      '  平均信噪比: %.2f dB\n' ...
                      '  信噪比标准差: %.2f dB\n' ...
                      '  信噪比范围: %.2f ~ %.2f dB\n' ...
                      '  平均锁定指示器: %.4f\n' ...
                      '  锁定百分比: %.1f%%\n\n' ...
                      '报告生成时间: %s\n'], ...
                     simSettings.fp/1e3, simSettings.fs/1e6, simSettings.fi/1e6, ...
                     simSettings.SNR, num_loops, time_vector(end), ...
                     stats.code_error_mean, stats.code_error_std, stats.code_error_rms, ...
                     stats.code_freq_std, ...
                     stats.carrier_error_mean, stats.carrier_error_std, stats.carrier_error_rms, ...
                     stats.carrier_freq_std, ...
                     stats.snr_mean, stats.snr_std, stats.snr_min, stats.snr_max, ...
                     stats.lock_mean, stats.lock_percentage, ...
                     datestr(now));

% 显示报告
fprintf('\n%s\n', report_text);

%% 6. 保存文件
if output_options.save_figures
    fprintf('  - 保存图形文件...\n');
    
    % 保存所有图形
    figures = [fig1, fig2, fig3, fig4, fig5, fig6];
    figure_names = {'tracking_errors', 'frequency_tracking', 'correlator_outputs', ...
                   'performance_monitoring', 'tracking_modes', 'statistical_analysis'};
    
    for i = 1:length(figures)
        filename = fullfile(output_options.output_dir, ...
                           sprintf('%s_%s.png', output_options.file_prefix, figure_names{i}));
        saveas(figures(i), filename);
        fprintf('    * 保存图形: %s\n', filename);
    end
end

if output_options.save_data
    fprintf('  - 保存数据文件...\n');
    
    % 保存跟踪数据
    tracking_data = struct();
    tracking_data.time_vector = time_vector;
    tracking_data.code_errors = code_errors;
    tracking_data.carrier_errors = carrier_errors;
    tracking_data.code_freqs = code_freqs;
    tracking_data.carrier_freqs = carrier_freqs;
    tracking_data.snr_estimates = snr_estimates;
    tracking_data.lock_indicators = lock_indicators;
    tracking_data.correlators = struct();
    tracking_data.correlators.I_E = I_E_values;
    tracking_data.correlators.Q_E = Q_E_values;
    tracking_data.correlators.I_P = I_P_values;
    tracking_data.correlators.Q_P = Q_P_values;
    tracking_data.correlators.I_L = I_L_values;
    tracking_data.correlators.Q_L = Q_L_values;
    tracking_data.statistics = stats;
    tracking_data.simSettings = simSettings;
    
    data_filename = fullfile(output_options.output_dir, ...
                            sprintf('%s_data.mat', output_options.file_prefix));
    save(data_filename, 'tracking_data');
    fprintf('    * 保存数据文件: %s\n', data_filename);
    
    % 保存文本报告
    report_filename = fullfile(output_options.output_dir, ...
                              sprintf('%s_report.txt', output_options.file_prefix));
    fid = fopen(report_filename, 'w');
    if fid > 0
        fprintf(fid, '%s', report_text);
        fclose(fid);
        fprintf('    * 保存报告文件: %s\n', report_filename);
    end
end

fprintf('OFDM跟踪结果输出完成\n');

end