function [signal, t] = signal_generation(fs, duration, fc, weil_length, use_custom_code)
% 信号生成模块：生成基于Weil码的定位信号
%
% 输入参数：
%   fs: 采样频率 (Hz)
%   duration: 信号持续时间 (s)
%   fc: 载波频率 (Hz)
%   weil_length: Weil码长度
%   use_custom_code: 是否使用自定义Weil码 (布尔值)
%
% 输出参数：
%   signal: 生成的信号
%   t: 时间向量

% 生成时间向量
t = 0:1/fs:duration-1/fs;

% 生成或加载Weil码
if use_custom_code
    weil_code = load_custom_weil_code(weil_length);
else
    weil_code = generate_weil_code(weil_length);
end

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
%   length: Weil码长度
%
% 输出参数：
%   weil_code: 生成的Weil码序列

% 如果长度不是质数，找到最接近的较小质数
if ~isprime(length)
    length = max(primes(length));
    warning('Weil码长度已调整为最接近的较小质数: %d', length);
end

weil_code = zeros(1, length);
for m = 0:length-1
    weil_code(m+1) = exp(1i * pi * mod(m^2, length) / length);
end

% 将复数序列转换为实数序列（1和-1）
weil_code = real(weil_code) > 0;
weil_code = 2 * weil_code - 1;

end

function weil_code = load_custom_weil_code(length)
% 加载自定义Weil码
%
% 输入参数：
%   length: 期望的Weil码长度
%
% 输出参数：
%   weil_code: 加载的Weil码序列

% 加载.mat文件
data = load('custom_weil_code.mat');

% 假设.mat文件中的Weil码变量名为'weil_code'
if ~isfield(data, 'weil_code')
    error('custom_weil_code.mat 文件中没有找到 weil_code 变量');
end

custom_code = data.weil_code;

% 检查加载的Weil码长度是否符合要求
if length(custom_code) ~= length
    warning('加载的Weil码长度 (%d) 与要求的长度 (%d) 不符。将进行截断或重复。', length(custom_code), length);
    if length(custom_code) > length
        weil_code = custom_code(1:length);
    else
        weil_code = repmat(custom_code, 1, ceil(length / length(custom_code)));
        weil_code = weil_code(1:length);
    end
else
    weil_code = custom_code;
end

end
