function interpolated_signal = interpo(input_signal, interpolation_factor)
% 信号插值函数
%
% 功能描述:
%   对输入信号进行插值处理，提高信号的采样率
%   支持多行信号矩阵的插值处理
%
% 输入参数:
%   input_signal         - 输入信号矩阵 (numBand × signal_length)
%   interpolation_factor - 插值倍数 (正整数)
%
% 输出参数:
%   interpolated_signal - 插值后的信号矩阵 (numBand × signal_length*interpolation_factor)
%
% 插值方法:
%   使用零填充插值法，在频域进行插值以保持信号特性
%   对于每一行信号，先进行FFT，然后零填充，最后IFFT
%
% 作者: OFDM接收机开发
% 日期: 2025年7月
% 版本: 1.0

%% 输入参数验证
if nargin < 2
    error('interpo: 需要提供输入信号和插值倍数');
end

if interpolation_factor <= 0 || mod(interpolation_factor, 1) ~= 0
    error('interpo: 插值倍数必须为正整数');
end

if interpolation_factor == 1
    % 插值倍数为1，直接返回原信号
    interpolated_signal = input_signal;
    return;
end

if isempty(input_signal)
    error('interpo: 输入信号不能为空');
end

%% 获取信号维度信息
[num_bands, signal_length] = size(input_signal);

if signal_length == 1
    % 如果是列向量，转置为行向量处理
    input_signal = input_signal.';
    [num_bands, signal_length] = size(input_signal);
    transpose_output = true;
else
    transpose_output = false;
end

%% 初始化输出信号
output_length = signal_length * interpolation_factor;
interpolated_signal = zeros(num_bands, output_length);

%% 对每个频带进行插值处理
for band_idx = 1:num_bands
    % 获取当前频带的信号
    current_signal = input_signal(band_idx, :);
    
    % 检查信号是否为复数
    is_complex = ~isreal(current_signal);
    
    if is_complex
        % 复数信号：分别对实部和虚部进行插值
        real_part = real(current_signal);
        imag_part = imag(current_signal);
        
        % 插值实部
        interp_real = perform_interpolation(real_part, interpolation_factor);
        
        % 插值虚部
        interp_imag = perform_interpolation(imag_part, interpolation_factor);
        
        % 合并复数信号
        interpolated_signal(band_idx, :) = interp_real + 1j * interp_imag;
    else
        % 实数信号：直接插值
        interpolated_signal(band_idx, :) = perform_interpolation(current_signal, interpolation_factor);
    end
end

%% 如果输入是列向量，转置输出
if transpose_output
    interpolated_signal = interpolated_signal.';
end

end

function interpolated_data = perform_interpolation(data, factor)
% 执行实际的插值操作
%
% 输入参数:
%   data   - 输入数据向量
%   factor - 插值倍数
%
% 输出参数:
%   interpolated_data - 插值后的数据向量

data_length = length(data);

% 选择插值方法
if data_length < 100
    % 对于短信号，使用简单的线性插值
    interpolated_data = linear_interpolation(data, factor);
else
    % 对于长信号，使用频域插值
    interpolated_data = frequency_domain_interpolation(data, factor);
end

end

function interpolated_data = linear_interpolation(data, factor)
% 线性插值方法
%
% 适用于短信号或实时处理场景

data_length = length(data);
output_length = data_length * factor;

% 创建原始采样点索引
original_indices = 1:data_length;

% 创建插值后的采样点索引
new_indices = linspace(1, data_length, output_length);

% 执行线性插值
interpolated_data = interp1(original_indices, data, new_indices, 'linear', 'extrap');

end

function interpolated_data = frequency_domain_interpolation(data, factor)
% 频域插值方法
%
% 通过零填充频域实现高质量插值

data_length = length(data);

% 计算FFT
data_fft = fft(data);

% 计算零填充后的FFT长度
padded_length = data_length * factor;

% 创建零填充的频域信号
padded_fft = zeros(1, padded_length);

% 计算频域数据的放置位置
half_length = floor(data_length / 2);

if mod(data_length, 2) == 0
    % 偶数长度
    % 放置低频部分
    padded_fft(1:half_length+1) = data_fft(1:half_length+1);
    % 放置高频部分
    padded_fft(end-half_length+1:end) = data_fft(half_length+2:end);
    
    % 处理奈奎斯特频率分量
    if half_length > 0
        padded_fft(half_length+1) = padded_fft(half_length+1) / 2;
        padded_fft(end-half_length+1) = padded_fft(half_length+1);
    end
else
    % 奇数长度
    % 放置低频部分
    padded_fft(1:half_length+1) = data_fft(1:half_length+1);
    % 放置高频部分
    padded_fft(end-half_length+1:end) = data_fft(half_length+2:end);
end

% 执行IFFT并调整幅度
interpolated_data = ifft(padded_fft) * factor;

% 确保输出为实数（如果输入为实数）
if isreal(data)
    interpolated_data = real(interpolated_data);
end

end