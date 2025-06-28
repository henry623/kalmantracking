% 主仿真程序：高动态性航空信道中的信号跟踪

% 添加所有模块文件所在的路径
addpath(pwd);

% 1. 设置仿真参数
fs = 1e6;  % 采样频率 (Hz)
duration = 1;  % 仿真持续时间 (s)
t = 0:1/fs:duration-1/fs;  % 时间向量
fc = 1.5e9;  % 载波频率 (Hz)
weil_length = 10230;  % Weil码长度
use_custom_code = true;  % 使用自定义Weil码
v = 1000;  % 初始速度 (m/s)
a = 100;  % 加速度 (m/s^2)
snr_db = -10;  % 信噪比 (dB)

% PLL参数
pll_params.bw = 50;  % 环路带宽 (Hz)
pll_params.damping = 0.707;  % 阻尼系数

% DLL参数
dll_params.bw = 5;  % 环路带宽 (Hz)
dll_params.damping = 0.707;  % 阻尼系数
dll_params.early_late_spacing = 0.5;  % 早晚码间隔 (码片)

% 2. 生成信号
[signal, t] = signal_generation(fs, duration, fc, weil_length, use_custom_code);

% 3. 生成信道响应
[channel_response, doppler_shift] = channel_model(t, fc, v, a);

% 4. 信号通过信道传输
received_signal = signal_transmission(signal, channel_response, snr_db);

% 5. 信号跟踪
[tracked_signal, estimated_phase, estimated_freq] = signal_tracking(received_signal, t, fs, fc, pll_params, dll_params);

% 6. 计算真实相位和频率
c = 3e8;  % 光速 (m/s)
true_phase = 2 * pi * (fc * v * t + 0.5 * fc * a * t.^2) / c;
true_freq = fc * (v + a * t) / c;

% 7. 可视化结果
result_visualization(t, signal, received_signal, tracked_signal, estimated_phase, estimated_freq, true_phase, true_freq);

% 8. 计算和显示性能指标
phase_error = wrapToPi(estimated_phase - true_phase);
freq_error = estimated_freq - true_freq;

phase_rmse = sqrt(mean(phase_error.^2));
freq_rmse = sqrt(mean(freq_error.^2));

fprintf('相位估计RMSE: %.4f rad\n', phase_rmse);
fprintf('频率估计RMSE: %.4f Hz\n', freq_rmse);

% 9. 绘制多普勒频移
figure;
plot(t, doppler_shift, 'b', t, estimated_freq, 'r');
title('多普勒频移估计');
xlabel('时间 (s)');
ylabel('频率 (Hz)');
legend('真实多普勒频移', '估计频率');
grid on;

% 10. 计算并显示多普勒频移估计误差
doppler_error = estimated_freq - doppler_shift;
doppler_rmse = sqrt(mean(doppler_error.^2));
fprintf('多普勒频移估计RMSE: %.4f Hz\n', doppler_rmse);
