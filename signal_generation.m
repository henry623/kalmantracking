function [signal, t] = signal_generation(fs, duration, fc, weil_length)
% 信号生成模块：生成基于Weil码的定位信号
%
% 输入参数：
%   fs: 采样频率 (Hz)
%   duration: 信号持续时间 (s)
%   fc: 载波频率 (Hz)
%   weil_length: Weil码长度
%
% 输出参数：
%   signal: 生成的信号
%   t: 时间向量

% 生成时间向量
t = 0:1/fs:duration-1/fs;

% 生成Weil码
weil_code = generate_weil_code(weil_length);

% 将Weil码重复以覆盖整个信号持续时间
code_samples = round(fs * duration / weil_length);
repeated_code = repmat(weil_code, 1, ceil(code_samples / weil_length));
repeated_code = repeated_code(1:code_samples);

% 生成载波信号
carrier = cos(2*pi*fc*t);

% 调制信号（BPSK调制）
signal = carrier .* repeated_code;

end

function weil_code = generate_weil_code(length)
% 生成Weil码
%
% 输入参数：
%   length: Weil码长度（必须是质数）
%
% 输出参数：
%   weil_code: 生成的Weil码序列

if ~isprime(length)
    error('Weil码长度必须是质数');
end

weil_code = zeros(1, length);
for m = 0:length-1
    weil_code(m+1) = exp(1i * pi * mod(m^2, length) / length);
end

% 将复数序列转换为实数序列（1和-1）
weil_code = real(weil_code) > 0;
weil_code = 2 * weil_code - 1;

end