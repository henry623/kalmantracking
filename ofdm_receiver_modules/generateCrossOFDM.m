function [resource_grid, rawSignal] = generateCrossOFDM(simSettings)
% 生成交叉OFDM调制信号
%
% 功能描述:
%   根据输入的码序列生成OFDM调制信号，支持多频带和多符号处理
%   将码序列映射到OFDM子载波上，并生成时域信号
%
% 输入参数:
%   simSettings - 仿真设置参数结构体，包含：
%                * code: 输入码序列
%                * Nu: 每个频带的子载波数量
%                * numBand: 频带数量
%                * NFFT: FFT点数
%                * nSymbol: OFDM符号数
%
% 输出参数:
%   resource_grid - OFDM资源网格 (nSymbol × NFFT × numBand)
%   rawSignal    - 时域OFDM信号 (numBand × nSymbol*NFFT)
%
% 处理流程:
%   1. 参数验证和提取
%   2. 码序列预处理和扩展
%   3. 子载波映射
%   4. IFFT变换生成时域信号
%   5. 多频带信号组合
%
% 作者: OFDM接收机开发
% 日期: 2025年7月
% 版本: 1.0

%% 1. 输入参数验证
if nargin < 1
    error('generateCrossOFDM: 必须提供仿真设置参数');
end

% 检查必要字段
required_fields = {'code', 'Nu', 'numBand', 'NFFT', 'nSymbol'};
for i = 1:length(required_fields)
    if ~isfield(simSettings, required_fields{i})
        error('generateCrossOFDM: simSettings缺少必要字段 %s', required_fields{i});
    end
end

% 提取参数
code = simSettings.code;
Nu = simSettings.Nu;
numBand = simSettings.numBand;
NFFT = simSettings.NFFT;
nSymbol = simSettings.nSymbol;

%% 2. 参数合理性检查
if length(code) == 0
    error('generateCrossOFDM: 码序列不能为空');
end

if Nu <= 0 || numBand <= 0 || NFFT <= 0 || nSymbol <= 0
    error('generateCrossOFDM: OFDM参数必须为正数');
end

if Nu > NFFT
    error('generateCrossOFDM: 子载波数不能超过FFT点数');
end

%% 3. 码序列预处理
fprintf('  - 码序列预处理...\n');

% 计算总的数据子载波数
total_subcarriers = Nu * numBand;
total_data_symbols = total_subcarriers * nSymbol;

% 扩展码序列以填充所有数据位置
code_length = length(code);
if total_data_symbols > code_length
    % 重复码序列
    repeat_times = ceil(total_data_symbols / code_length);
    extended_code = repmat(code, 1, repeat_times);
    data_symbols = extended_code(1:total_data_symbols);
    fprintf('    * 码序列重复 %d 次以填充 %d 个数据符号\n', repeat_times, total_data_symbols);
else
    % 截取码序列
    data_symbols = code(1:total_data_symbols);
    fprintf('    * 截取码序列前 %d 个符号\n', total_data_symbols);
end

% 将数据符号重新整形为 (nSymbol × total_subcarriers)
% 确保数据符号是行向量
if size(data_symbols, 1) > size(data_symbols, 2)
    data_symbols = data_symbols.';
end

% 安全的reshape操作
try
    data_matrix = reshape(data_symbols, total_subcarriers, nSymbol).';
catch ME
    % 如果reshape失败，尝试调整数据长度
    fprintf('    * 警告: reshape失败，调整数据长度\n');
    required_length = total_subcarriers * nSymbol;
    if length(data_symbols) < required_length
        % 数据不足，重复填充
        data_symbols = repmat(data_symbols, 1, ceil(required_length / length(data_symbols)));
        data_symbols = data_symbols(1:required_length);
    else
        % 数据过多，截取
        data_symbols = data_symbols(1:required_length);
    end
    data_matrix = reshape(data_symbols, total_subcarriers, nSymbol).';
end

%% 4. 初始化资源网格
fprintf('  - 初始化OFDM资源网格...\n');

% 资源网格: nSymbol × NFFT × numBand
resource_grid = zeros(nSymbol, NFFT, numBand);

%% 5. 子载波映射
fprintf('  - 执行子载波映射...\n');

for bandID = 1:numBand
    fprintf('    * 处理频带 %d/%d\n', bandID, numBand);
    
    for symbolID = 1:nSymbol
        % 获取当前符号的数据
        start_idx = (bandID - 1) * Nu + 1;
        end_idx = bandID * Nu;
        symbol_data = data_matrix(symbolID, start_idx:end_idx);
        
        % 子载波映射策略：将数据放在中心子载波上
        % 避免DC子载波和边缘子载波
        dc_index = NFFT/2 + 1;  % DC子载波索引
        
        % 计算数据子载波的起始和结束索引
        data_start = dc_index - floor(Nu/2);
        data_end = data_start + Nu - 1;
        
        % 确保索引在有效范围内
        if data_start < 1
            data_start = 2;  % 避开DC
            data_end = data_start + Nu - 1;
        end
        
        if data_end > NFFT
            data_end = NFFT;
            data_start = data_end - Nu + 1;
        end
        
        % 跳过DC子载波
        if data_start <= dc_index && data_end >= dc_index
            % 分两段放置数据，跳过DC
            pre_dc_length = dc_index - data_start;
            post_dc_length = Nu - pre_dc_length;
            
            if pre_dc_length > 0
                resource_grid(symbolID, data_start:dc_index-1, bandID) = symbol_data(1:pre_dc_length);
            end
            if post_dc_length > 0
                resource_grid(symbolID, dc_index+1:dc_index+post_dc_length, bandID) = symbol_data(pre_dc_length+1:end);
            end
        else
            % 直接放置数据
            actual_length = min(Nu, data_end - data_start + 1);
            resource_grid(symbolID, data_start:data_start+actual_length-1, bandID) = symbol_data(1:actual_length);
        end
    end
end

fprintf('    * 子载波映射完成\n');

%% 6. IFFT变换生成时域信号
fprintf('  - 执行IFFT变换...\n');

% 初始化时域信号
rawSignal = zeros(numBand, nSymbol * NFFT);

for bandID = 1:numBand
    fprintf('    * IFFT处理频带 %d/%d\n', bandID, numBand);
    
    band_time_signal = [];
    
    for symbolID = 1:nSymbol
        % 获取当前符号的频域数据
        freq_data = squeeze(resource_grid(symbolID, :, bandID));
        
        % IFFT变换
        time_data = ifft(freq_data, NFFT);
        
        % 确保time_data是行向量，然后添加到时域信号
        if size(time_data, 1) > size(time_data, 2)
            time_data = time_data.';  % 转置为行向量
        end
        band_time_signal = [band_time_signal, time_data];
    end
    
    % 存储频带信号
    rawSignal(bandID, :) = band_time_signal;
end

fprintf('    * IFFT变换完成\n');

%% 7. 信号质量检查
fprintf('  - 信号质量检查...\n');

for bandID = 1:numBand
    band_signal = rawSignal(bandID, :);
    
    % 计算信号统计量
    signal_power = mean(abs(band_signal).^2);
    signal_peak = max(abs(band_signal));
    signal_rms = rms(band_signal);
    
    % 计算峰均比 (PAPR)
    papr_db = 20 * log10(signal_peak / signal_rms);
    
    fprintf('    * 频带 %d: 功率=%.6f, PAPR=%.2f dB\n', bandID, signal_power, papr_db);
    
    % 检查异常值
    if signal_power < 1e-10
        warning('generateCrossOFDM: 频带 %d 信号功率过低', bandID);
    end
    
    if papr_db > 15
        warning('generateCrossOFDM: 频带 %d 峰均比过高 (%.2f dB)', bandID, papr_db);
    end
end

%% 8. 最终验证
fprintf('  - 最终验证...\n');

% 检查输出维度
expected_resource_grid_size = [nSymbol, NFFT, numBand];
actual_resource_grid_size = size(resource_grid);
if ~isequal(expected_resource_grid_size, actual_resource_grid_size)
    error('generateCrossOFDM: 资源网格大小不匹配，期望 %s，实际 %s', ...
          mat2str(expected_resource_grid_size), mat2str(actual_resource_grid_size));
end

expected_signal_size = [numBand, nSymbol * NFFT];
actual_signal_size = size(rawSignal);
if ~isequal(expected_signal_size, actual_signal_size)
    error('generateCrossOFDM: 时域信号大小不匹配，期望 %s，实际 %s', ...
          mat2str(expected_signal_size), mat2str(actual_signal_size));
end

% 检查信号中是否包含无效值
if any(~isfinite(resource_grid), 'all')
    error('generateCrossOFDM: 资源网格包含无效值（NaN或Inf）');
end

if any(~isfinite(rawSignal), 'all')
    error('generateCrossOFDM: 时域信号包含无效值（NaN或Inf）');
end

fprintf('    * 验证通过，OFDM信号生成完成\n');

end