% 测试OFDM接收机模块化系统
% 
% 本脚本演示如何使用新创建的模块化OFDM接收机系统
% 包括初始化、本地信号生成、跟踪参数设置、信号跟踪和结果输出
%
% 作者: OFDM接收机开发团队
% 日期: 2025年1月
% 版本: 1.0

clear; close all; clc;

fprintf('=== OFDM接收机模块化系统测试 ===\n\n');

%% 1. 初始化仿真参数
fprintf('1. 初始化仿真参数...\n');

% 使用现有的初始化函数
simSettings = init(-15);  % SNR = -15 dB

% 添加OFDM接收机需要的额外参数
if ~isfield(simSettings, 'NFFT')
    simSettings.NFFT = 1024;  % FFT点数
end
if ~isfield(simSettings, 'nSymbol')
    simSettings.nSymbol = 10;  % OFDM符号数
end
if ~isfield(simSettings, 't_total')
    simSettings.t_total = 1.0;  % 总仿真时间 (秒)
end

% 加载本地码
load('weil10230_signed.mat');
simSettings.code = weil10230_signed(1, :);  % 使用第一个PRN码

fprintf('   - 仿真参数初始化完成\n');
fprintf('   - 码频率: %.2f kHz\n', simSettings.fp/1e3);
fprintf('   - 采样频率: %.2f MHz\n', simSettings.fs/1e6);
fprintf('   - 信噪比: %.1f dB\n', simSettings.SNR);
fprintf('   - 本地码长度: %d\n', length(simSettings.code));

%% 2. 生成测试信号
fprintf('\n2. 生成测试信号...\n');

try
    % 生成OFDM调制信号作为接收信号
    [resource_grid, test_signal] = generateCrossOFDM(simSettings);
    
    % 添加噪声
    signal_power = mean(abs(test_signal).^2, 'all');
    noise_power = signal_power / (10^(simSettings.SNR/10));
    noise = sqrt(noise_power/2) * (randn(size(test_signal)) + 1j*randn(size(test_signal)));
    
    % 接收信号 = 信号 + 噪声
    received_signal = test_signal + noise;
    
    fprintf('   - 测试信号生成完成\n');
    fprintf('   - 信号大小: %dx%d\n', size(received_signal));
    fprintf('   - 信号功率: %.6f\n', signal_power);
    fprintf('   - 噪声功率: %.6f\n', noise_power);
    
catch ME
    fprintf('   - 错误: 信号生成失败 - %s\n', ME.message);
    return;
end

%% 3. 测试接收机初始化模块
fprintf('\n3. 测试接收机初始化模块...\n');

try
    [processed_signal, init_params] = ofdm_receiver_init(simSettings, received_signal);
    
    fprintf('   - 接收机初始化成功\n');
    fprintf('   - 处理后信号长度: %d\n', processed_signal.length);
    fprintf('   - 信号格式: %s\n', processed_signal.format);
    
catch ME
    fprintf('   - 错误: 接收机初始化失败 - %s\n', ME.message);
    return;
end

%% 4. 测试本地信号生成模块
fprintf('\n4. 测试本地信号生成模块...\n');

try
    local_signals = ofdm_generate_local_signals(simSettings);
    
    fprintf('   - 本地信号生成成功\n');
    fprintf('   - 本地信号大小: %dx%d\n', size(local_signals.yr0));
    fprintf('   - 码类型: %s\n', local_signals.code_info.type);
    
catch ME
    fprintf('   - 错误: 本地信号生成失败 - %s\n', ME.message);
    return;
end

%% 5. 测试跟踪参数初始化模块
fprintf('\n5. 测试跟踪参数初始化模块...\n');

try
    tracking_params = ofdm_tracking_init(simSettings);
    
    fprintf('   - 跟踪参数初始化成功\n');
    fprintf('   - DLL噪声带宽: %.2f Hz\n', tracking_params.dll.noiseBandwidth);
    fprintf('   - PLL噪声带宽: %.2f Hz\n', tracking_params.pll.noiseBandwidth);
    
catch ME
    fprintf('   - 错误: 跟踪参数初始化失败 - %s\n', ME.message);
    return;
end

%% 6. 测试信号跟踪处理模块
fprintf('\n6. 测试信号跟踪处理模块...\n');

% 计算跟踪循环数
total_samples = processed_signal.length;
samples_per_loop = tracking_params.basic.lenOFDM;
num_tracking_loops = floor(total_samples / samples_per_loop);

fprintf('   - 总样本数: %d\n', total_samples);
fprintf('   - 每循环样本数: %d\n', samples_per_loop);
fprintf('   - 跟踪循环数: %d\n', num_tracking_loops);

% 存储所有跟踪结果
all_tracking_results = cell(1, num_tracking_loops);

try
    for loop_idx = 1:num_tracking_loops
        fprintf('   - 执行跟踪循环 %d/%d\n', loop_idx, num_tracking_loops);
        
        tracking_result = ofdm_signal_tracking(processed_signal, local_signals, ...
                                              tracking_params, loop_idx);
        
        all_tracking_results{loop_idx} = tracking_result;
        
        % 更新跟踪参数（在实际应用中，这里应该更新tracking_params的状态）
        % 这里简化处理，实际的状态更新在ofdm_signal_tracking函数内部通过持久变量处理
    end
    
    fprintf('   - 信号跟踪处理完成\n');
    
catch ME
    fprintf('   - 错误: 信号跟踪处理失败 - %s\n', ME.message);
    return;
end

%% 7. 测试结果输出模块
fprintf('\n7. 测试结果输出模块...\n');

% 设置输出选项
output_options = struct();
output_options.save_figures = true;   % 保存图形
output_options.show_plots = true;     % 显示图形
output_options.save_data = true;      % 保存数据
output_options.output_dir = './test_results';  % 输出目录
output_options.file_prefix = 'ofdm_test';      % 文件前缀

try
    ofdm_results_output(all_tracking_results, simSettings, output_options);
    
    fprintf('   - 结果输出完成\n');
    
catch ME
    fprintf('   - 错误: 结果输出失败 - %s\n', ME.message);
    return;
end

%% 8. 性能总结
fprintf('\n8. 性能总结...\n');

% 提取最后几个循环的性能指标
last_n = min(10, num_tracking_loops);
last_results = all_tracking_results(end-last_n+1:end);

% 计算平均性能
avg_snr = 0;
avg_lock = 0;
avg_code_error = 0;
avg_carrier_error = 0;

for i = 1:last_n
    result = last_results{i};
    avg_snr = avg_snr + result.performance.snr_estimate;
    avg_lock = avg_lock + result.performance.lock_indicator;
    avg_code_error = avg_code_error + abs(result.discriminators.code_error);
    avg_carrier_error = avg_carrier_error + abs(result.discriminators.carrier_error);
end

avg_snr = avg_snr / last_n;
avg_lock = avg_lock / last_n;
avg_code_error = avg_code_error / last_n;
avg_carrier_error = avg_carrier_error / last_n;

fprintf('   - 最后%d个循环的平均性能:\n', last_n);
fprintf('     * 平均信噪比: %.2f dB\n', avg_snr);
fprintf('     * 平均锁定指示器: %.4f\n', avg_lock);
fprintf('     * 平均码跟踪误差: %.6f\n', avg_code_error);
fprintf('     * 平均载波跟踪误差: %.6f\n', avg_carrier_error);

% 判断跟踪性能
if avg_lock > 0.7
    fprintf('     * 跟踪状态: 良好锁定\n');
elseif avg_lock > 0.5
    fprintf('     * 跟踪状态: 基本锁定\n');
else
    fprintf('     * 跟踪状态: 锁定困难\n');
end

%% 9. 与原始系统比较（可选）
fprintf('\n9. 与原始系统比较...\n');

try
    % 调用原始的跟踪函数进行比较
    fprintf('   - 运行原始跟踪系统进行比较...\n');
    
    % 这里可以调用原始的track函数进行比较
    % [trackResults, channel] = track(acqResults, simSettings);
    
    fprintf('   - 注意: 原始系统比较需要完整的捕获结果，此处跳过\n');
    
catch ME
    fprintf('   - 原始系统比较失败: %s\n', ME.message);
end

%% 10. 测试完成
fprintf('\n=== OFDM接收机模块化系统测试完成 ===\n');
fprintf('测试结果已保存到: %s\n', output_options.output_dir);
fprintf('请查看生成的图形和报告文件\n\n');

% 显示模块化系统的优势
fprintf('模块化系统优势:\n');
fprintf('1. 代码结构清晰，易于理解和维护\n');
fprintf('2. 模块间接口明确，便于测试和调试\n');
fprintf('3. 支持独立的性能分析和可视化\n');
fprintf('4. 易于扩展和修改特定功能\n');
fprintf('5. 提供详细的处理过程信息\n\n');

fprintf('建议后续工作:\n');
fprintf('1. 根据实际需求调整跟踪参数\n');
fprintf('2. 添加更多的性能监控指标\n');
fprintf('3. 实现自适应参数调整功能\n');
fprintf('4. 优化计算效率和内存使用\n');
fprintf('5. 添加更多的错误处理和恢复机制\n');