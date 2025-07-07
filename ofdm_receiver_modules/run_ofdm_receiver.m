% 运行OFDM接收机主函数的脚本
%
% 这个脚本设置必要的参数并调用ofdm_receiver_main函数
%
% 作者: OFDM接收机开发
% 日期: 2025年7月
% 版本: 1.0

clear; close all; clc;

fprintf('=== OFDM接收机运行脚本 ===\n\n');

%% 1. 初始化仿真参数
fprintf('1. 初始化仿真参数...\n');

% 使用现有的初始化函数
simSettings = init(-20);  % SNR = -15 dB

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
load('../weil10230_signed.mat');
simSettings.code = weil10230_signed(1, :);  % 使用第一个PRN码

fprintf('   ✓ 仿真参数初始化完成\n');
fprintf('   - 信噪比: %.1f dB\n', simSettings.SNR);
fprintf('   - 采样频率: %.2f MHz\n', simSettings.fs/1e6);
fprintf('   - 码频率: %.2f kHz\n', simSettings.fp/1e3);

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
    
    fprintf('   ✓ 测试信号生成完成\n');
    fprintf('   - 信号维度: %s\n', mat2str(size(received_signal)));
    fprintf('   - 信号功率: %.6f\n', signal_power);
    fprintf('   - 噪声功率: %.6f\n', noise_power);
    
catch ME
    fprintf('   ✗ 错误: 信号生成失败 - %s\n', ME.message);
    return;
end

%% 3. 运行OFDM接收机主函数
fprintf('\n3. 运行OFDM接收机主函数...\n');

try
    % 调用主函数
    % 参数: (仿真设置, 接收信号, 是否启用卡尔曼滤波)
    enable_kalman = true;  % 暂时禁用卡尔曼滤波器，先测试基本功能
    
    fprintf('   调用 ofdm_receiver_main...\n');
    output = ofdm_receiver_main(simSettings, received_signal, enable_kalman);
    
    fprintf('   ✅ OFDM接收机处理完成！\n');
    
catch ME
    fprintf('   ✗ 错误: OFDM接收机处理失败 - %s\n', ME.message);
    fprintf('   错误堆栈:\n');
    for i = 1:length(ME.stack)
        fprintf('     %s (行 %d): %s\n', ME.stack(i).file, ME.stack(i).line, ME.stack(i).name);
    end
    return;
end

%% 4. 显示处理结果
fprintf('\n4. 处理结果总结:\n');

if exist('output', 'var') && isstruct(output)
    % 显示主要结果
    fields = fieldnames(output);
    fprintf('   输出结构包含以下字段:\n');
    for i = 1:length(fields)
        fprintf('   - %s\n', fields{i});
    end
    
    % 如果有跟踪结果，显示一些关键指标
    if isfield(output, 'tracking_results')
        fprintf('\n   跟踪处理结果:\n');
        if isfield(output.tracking_results, 'snr_estimate')
            fprintf('   - 估算信噪比: %.2f dB\n', output.tracking_results.snr_estimate(end));
        end
        if isfield(output.tracking_results, 'lock_indicator')
            fprintf('   - 锁定指示器: %.4f\n', output.tracking_results.lock_indicator(end));
        end
    end
    
    % 如果有性能评估结果
    if isfield(output, 'performance')
        fprintf('\n   性能评估结果:\n');
        perf_fields = fieldnames(output.performance);
        for i = 1:length(perf_fields)
            fprintf('   - %s: 已计算\n', perf_fields{i});
        end
    end
    
else
    fprintf('   ⚠ 输出结果格式异常\n');
end

fprintf('\n=== OFDM接收机运行完成 ===\n');
fprintf('总处理时间: %.2f 秒\n', toc);

%% 5. 保存结果（可选）
save_results = false;  % 设置为true以保存结果

if save_results && exist('output', 'var')
    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    filename = sprintf('ofdm_results_%s.mat', timestamp);
    save(filename, 'output', 'simSettings');
    fprintf('结果已保存到: %s\n', filename);
end