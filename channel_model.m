function [channel_response, doppler_shift] = channel_model(t, fc, v, a)
% 高动态性航空信道模型
%
% 输入参数：
%   t: 时间向量 (s)
%   fc: 载波频率 (Hz)
%   v: 初始速度 (m/s)
%   a: 加速度 (m/s^2)
%
% 输出参数：
%   channel_response: 信道响应
%   doppler_shift: 多普勒频移

% 光速 (m/s)
c = 3e8;

% 计算时变多普勒频移
doppler_shift = fc * (v + a * t) / c;

% 计算相位变化
phase_shift = 2 * pi * (fc * v * t + 0.5 * fc * a * t.^2) / c;

% 生成信道响应
channel_response = exp(1i * phase_shift);

% 添加多径效应（简化模型）
num_paths = 3;
delays = [0, 1e-6, 2e-6];  % 延迟时间
amplitudes = [1, 0.8, 0.6];  % 振幅衰减

for i = 2:num_paths
    delayed_phase = 2 * pi * (fc * v * (t - delays(i)) + 0.5 * fc * a * (t - delays(i)).^2) / c;
    channel_response = channel_response + amplitudes(i) * exp(1i * delayed_phase);
end

% 归一化信道响应
channel_response = channel_response / sqrt(sum(abs(amplitudes).^2));

end