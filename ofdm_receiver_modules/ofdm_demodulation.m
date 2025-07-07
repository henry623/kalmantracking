function demod_output = ofdm_demodulation(tracking_output, simSettings)
% OFDM解调模块
%
% 功能描述:
%   对跟踪后的OFDM信号进行解调处理
%   包括FFT变换、子载波解映射、资源网格提取等
%
% 输入参数:
%   tracking_output - 跟踪处理结果结构体
%   simSettings    - 仿真设置参数
%
% 输出参数:
%   demod_output   - 解调结果结构体，包含：
%                   * DemodulatedData: 解调后的数据
%                   * ResourceGrid: 资源网格
%                   * SubcarrierData: 子载波数据
%                   * ChannelEstimate: 信道估计
%
% 作者: OFDM接收机开发
% 日期: 2025年7月
% 版本: 1.0

%% 输入参数验证
if nargin < 2
    error('ofdm_demodulation: 需要提供跟踪输出和仿真设置');
end

fprintf('  - 开始OFDM解调处理...\n');

%% 提取参数
NFFT = simSettings.NFFT;
Nu = simSettings.Nu;
nSymbol = simSettings.nSymbol;
CP_length = simSettings.CP_length;
numBand = simSettings.numBand;

%% 初始化输出结构体
demod_output = struct();
demod_output.DemodulatedData = [];
demod_output.ResourceGrid = [];
demod_output.SubcarrierData = [];
demod_output.ChannelEstimate = [];

%% 检查跟踪输出是否有效
if ~isfield(tracking_output, 'OutCarrFreq') || isempty(tracking_output.OutCarrFreq)
    fprintf('    - 警告: 跟踪输出无效，跳过解调处理\n');
    return;
end

%% 提取跟踪后的信号
% 这里简化处理，假设跟踪输出包含载波频率和码相位信息
carrFreq = tracking_output.OutCarrFreq;
codePhase = tracking_output.OutCodePhase;

% 获取信号长度
signal_length = length(carrFreq);
fprintf('    - 信号长度: %d 采样点\n', signal_length);

%% 计算OFDM符号数量
samples_per_symbol = NFFT + CP_length;
available_symbols = floor(signal_length / samples_per_symbol);
actual_symbols = min(available_symbols, nSymbol);

fprintf('    - 可用符号数: %d, 处理符号数: %d\n', available_symbols, actual_symbols);

if actual_symbols == 0
    fprintf('    - 警告: 信号长度不足以解调一个OFDM符号\n');
    return;
end

%% 初始化资源网格
resource_grid = zeros(Nu, actual_symbols, numBand);
subcarrier_data = cell(numBand, 1);

%% 对每个频带进行解调
for bandID = 1:numBand
    fprintf('    - 处理频带 %d/%d\n', bandID, numBand);
    
    % 初始化当前频带的子载波数据
    band_subcarriers = zeros(Nu, actual_symbols);
    
    % 对每个OFDM符号进行处理
    for symIdx = 1:actual_symbols
        % 计算当前符号的起始位置
        symbol_start = (symIdx - 1) * samples_per_symbol + CP_length + 1;
        symbol_end = symbol_start + NFFT - 1;
        
        % 检查索引范围
        if symbol_end > signal_length
            fprintf('      - 警告: 符号 %d 超出信号范围，跳过\n', symIdx);
            break;
        end
        
        % 提取当前符号（去除循环前缀）
        symbol_samples = symbol_start:symbol_end;
        
        % 生成时域信号（简化处理，使用载波频率信息）
        if length(symbol_samples) == NFFT
            % 构造复数信号（简化处理）
            time_signal = exp(1j * 2 * pi * carrFreq(symbol_samples) / simSettings.fs);
            
            % 执行FFT变换
            freq_domain = fft(time_signal, NFFT);
            
            % 提取有效子载波
            % 假设有效子载波位于中心位置
            start_idx = (NFFT - Nu) / 2 + 1;
            end_idx = start_idx + Nu - 1;
            
            if start_idx >= 1 && end_idx <= NFFT
                valid_subcarriers = freq_domain(start_idx:end_idx);
                band_subcarriers(:, symIdx) = valid_subcarriers;
            else
                fprintf('      - 警告: 子载波索引超出范围\n');
            end
        end
    end
    
    % 存储当前频带的数据
    subcarrier_data{bandID} = band_subcarriers;
    resource_grid(:, :, bandID) = band_subcarriers;
end

%% 数据解调（简化处理）
fprintf('    - 执行数据解调...\n');

% 将所有频带的数据合并
all_data = [];
for bandID = 1:numBand
    band_data = subcarrier_data{bandID};
    % 数据验证和清理
    if ~isnumeric(band_data)
        band_data = zeros(size(band_data));
    end
    band_data(~isfinite(band_data)) = 0;
    
    % 简单的幅度检测
    demod_symbols = sign(real(band_data(:)));
    all_data = [all_data; demod_symbols];
end

%% 信道估计（简化处理）
fprintf('    - 执行信道估计...\n');

% 简单的信道估计：假设信道为平坦衰落
channel_estimate = zeros(Nu, actual_symbols, numBand);
for bandID = 1:numBand
    for symIdx = 1:actual_symbols
        % 使用导频符号进行信道估计（这里简化为单位增益）
        channel_estimate(:, symIdx, bandID) = ones(Nu, 1);
    end
end

%% 构建输出结构体
demod_output.DemodulatedData = all_data;
demod_output.ResourceGrid = resource_grid;
demod_output.SubcarrierData = subcarrier_data;
demod_output.ChannelEstimate = channel_estimate;

% 添加处理信息
demod_output.ProcessingInfo = struct();
demod_output.ProcessingInfo.ProcessedSymbols = actual_symbols;
demod_output.ProcessingInfo.ProcessedBands = numBand;
demod_output.ProcessingInfo.SubcarriersPerBand = Nu;
demod_output.ProcessingInfo.FFTSize = NFFT;
demod_output.ProcessingInfo.CPLength = CP_length;

% 添加质量指标
demod_output.QualityMetrics = struct();
if ~isempty(all_data)
    demod_output.QualityMetrics.DataLength = length(all_data);
    demod_output.QualityMetrics.AverageSignalPower = mean(abs(resource_grid(:)).^2);
    demod_output.QualityMetrics.PeakSignalPower = max(abs(resource_grid(:)).^2);
else
    demod_output.QualityMetrics.DataLength = 0;
    demod_output.QualityMetrics.AverageSignalPower = 0;
    demod_output.QualityMetrics.PeakSignalPower = 0;
end

fprintf('    - OFDM解调完成\n');
fprintf('      * 解调数据长度: %d\n', length(all_data));
fprintf('      * 资源网格大小: %dx%dx%d\n', size(resource_grid));
fprintf('      * 平均信号功率: %.6f\n', demod_output.QualityMetrics.AverageSignalPower);

end