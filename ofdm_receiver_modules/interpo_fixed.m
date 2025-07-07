function interpolated_signal = interpo_fixed(input_signal, interpolation_factor)
% 修复版本的信号插值函数
%
% 功能描述:
%   对输入信号进行插值处理，提高信号的采样率
%   支持多行信号矩阵的插值处理，完全解决维度不匹配问题
%
% 输入参数:
%   input_signal         - 输入信号矩阵 (numBand × signal_length)
%   interpolation_factor - 插值倍数 (正整数)
%
% 输出参数:
%   interpolated_signal - 插值后的信号矩阵 (numBand × signal_length*interpolation_factor)
%
% 作者: OFDM接收机开发
% 日期: 2025年7月
% 版本: 2.0 (修复版)

%% 输入参数验证
if nargin < 2
    error('interpo_fixed: 需要提供输入信号和插值倍数');
end

if interpolation_factor <= 0 || mod(interpolation_factor, 1) ~= 0
    error('interpo_fixed: 插值倍数必须为正整数');
end

if interpolation_factor == 1
    % 插值倍数为1，直接返回原信号
    interpolated_signal = input_signal;
    return;
end

if isempty(input_signal)
    error('interpo_fixed: 输入信号不能为空');
end

%% 获取信号维度信息并标准化输入
original_input = input_signal;
transpose_output = false;

% 确保输入是矩阵格式 (行数 = 频带数, 列数 = 信号长度)
if isvector(input_signal)
    if size(input_signal, 1) > size(input_signal, 2)
        % 列向量，转置为行向量
        input_signal = input_signal.';
        transpose_output = true;
    end
    % 现在确保是行向量 (1 × N)
end

[num_bands, signal_length] = size(input_signal);

%% 计算输出维度
output_length = signal_length * interpolation_factor;
interpolated_signal = zeros(num_bands, output_length);

%% 对每个频带进行插值处理
for band_idx = 1:num_bands
    % 获取当前频带的信号
    current_signal = input_signal(band_idx, :);
    
    % 确保current_signal是行向量
    if size(current_signal, 1) > size(current_signal, 2)
        current_signal = current_signal.';
    end
    
    % 检查信号是否为复数
    if ~isreal(current_signal)
        % 复数信号：分别对实部和虚部进行插值
        % 数据验证和清理
        current_signal(~isfinite(current_signal)) = 0;
        real_part = real(current_signal);
        imag_part = imag(current_signal);
        
        % 插值实部和虚部
        interp_real = safe_interpolation(real_part, interpolation_factor, output_length);
        interp_imag = safe_interpolation(imag_part, interpolation_factor, output_length);
        
        % 合并复数信号
        interpolated_signal(band_idx, :) = interp_real + 1j * interp_imag;
    else
        % 实数信号：直接插值
        interpolated_signal(band_idx, :) = safe_interpolation(current_signal, interpolation_factor, output_length);
    end
end

%% 如果输入是列向量，转置输出
if transpose_output && num_bands == 1
    interpolated_signal = interpolated_signal.';
end

end

function interpolated_data = safe_interpolation(data, factor, target_length)
% 安全的插值函数，确保输出长度严格匹配目标长度
%
% 输入参数:
%   data          - 输入数据向量 (1 × N)
%   factor        - 插值倍数
%   target_length - 目标输出长度
%
% 输出参数:
%   interpolated_data - 插值后的数据向量 (1 × target_length)

data_length = length(data);

% 确保data是行向量
if size(data, 1) > size(data, 2)
    data = data.';
end

try
    if data_length <= 10
        % 对于很短的信号，使用简单的线性插值
        original_indices = 1:data_length;
        new_indices = linspace(1, data_length, target_length);
        interpolated_data = interp1(original_indices, data, new_indices, 'linear', 'extrap');
    else
        % 对于较长的信号，使用频域插值
        interpolated_data = frequency_domain_interpolation_safe(data, factor, target_length);
    end
    
    % 确保输出是行向量
    if size(interpolated_data, 1) > size(interpolated_data, 2)
        interpolated_data = interpolated_data.';
    end
    
    % 严格确保输出长度匹配
    current_length = length(interpolated_data);
    if current_length ~= target_length
        if current_length > target_length
            % 截取到目标长度
            interpolated_data = interpolated_data(1:target_length);
        else
            % 零填充到目标长度
            temp_data = zeros(1, target_length);
            temp_data(1:current_length) = interpolated_data;
            interpolated_data = temp_data;
        end
    end
    
catch ME
    % 如果插值失败，使用简单的重复填充方法
    warning('interpo_fixed:interpolationFailed', '插值失败，使用备用方法: %s', ME.message);
    interpolated_data = fallback_interpolation(data, target_length);
end

end

function interpolated_data = frequency_domain_interpolation_safe(data, factor, target_length)
% 安全的频域插值方法
%
% 输入参数:
%   data          - 输入数据向量
%   factor        - 插值倍数
%   target_length - 目标长度
%
% 输出参数:
%   interpolated_data - 插值后的数据向量

data_length = length(data);

% 计算FFT
data_fft = fft(data);

% 创建零填充的频域信号
padded_fft = zeros(1, target_length);

% 计算频域数据的放置位置
half_length = floor(data_length / 2);

if mod(data_length, 2) == 0
    % 偶数长度
    % 放置低频部分
    padded_fft(1:half_length+1) = data_fft(1:half_length+1);
    % 放置高频部分
    if half_length > 0 && (half_length+2) <= data_length
        end_start = max(1, target_length - half_length + 1);
        end_indices = end_start:target_length;
        source_indices = (half_length+2):data_length;
        
        % 确保索引不越界
        copy_length = min(length(end_indices), length(source_indices));
        if copy_length > 0
            padded_fft(end_indices(1:copy_length)) = data_fft(source_indices(1:copy_length));
        end
    end
    
    % 处理奈奎斯特频率分量
    if half_length > 0 && half_length+1 <= target_length
        padded_fft(half_length+1) = padded_fft(half_length+1) / 2;
        if target_length - half_length + 1 <= target_length
            padded_fft(target_length - half_length + 1) = padded_fft(half_length+1);
        end
    end
else
    % 奇数长度
    % 放置低频部分
    padded_fft(1:half_length+1) = data_fft(1:half_length+1);
    % 放置高频部分
    if half_length > 0 && (half_length+2) <= data_length
        end_start = max(1, target_length - half_length + 1);
        end_indices = end_start:target_length;
        source_indices = (half_length+2):data_length;
        
        % 确保索引不越界
        copy_length = min(length(end_indices), length(source_indices));
        if copy_length > 0
            padded_fft(end_indices(1:copy_length)) = data_fft(source_indices(1:copy_length));
        end
    end
end

% 执行IFFT并调整幅度
interpolated_data = ifft(padded_fft) * factor;

% 确保输出为实数（如果输入为实数）
if isreal(data)
    interpolated_data = real(interpolated_data);
end

end

function interpolated_data = fallback_interpolation(data, target_length)
% 备用插值方法：简单的重复和截取
%
% 输入参数:
%   data          - 输入数据向量
%   target_length - 目标长度
%
% 输出参数:
%   interpolated_data - 插值后的数据向量

data_length = length(data);

if target_length <= data_length
    % 目标长度较短，直接截取
    interpolated_data = data(1:target_length);
else
    % 目标长度较长，重复数据
    repeat_times = ceil(target_length / data_length);
    extended_data = repmat(data, 1, repeat_times);
    interpolated_data = extended_data(1:target_length);
end

end