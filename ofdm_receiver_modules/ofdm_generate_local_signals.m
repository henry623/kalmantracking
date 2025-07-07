function [local_signals] = ofdm_generate_local_signals(simSettings)
% OFDM本地参考信号生成模块
%
% 功能描述:
%   生成OFDM接收机所需的本地参考信号，包括本地码的OFDM调制
%   以及插值处理，为后续的相关运算提供参考信号
%
% 输入参数:
%   simSettings - 仿真设置参数结构体，包含：
%                * code: 本地伪随机码序列
%                * Nu: 每个频带的子载波数量
%                * numBand: 频带数量
%                * NFFT: FFT点数
%                * nSymbol: OFDM符号数
%                * famp: 插值倍数（采样率倍数）
%                * extra: 前后额外保护长度
%
% 输出参数:
%   local_signals - 本地信号结构体，包含：
%                  * yr0: 实部本地信号 (numBand × length)
%                  * yi0: 虚部本地信号 (numBand × length)
%                  * resource_grid: OFDM资源网格 (nSymbol × NFFT × numBand)
%                  * rawSignal: 原始OFDM时域信号 (numBand × nSymbol*NFFT)
%                  * code_info: 码信息结构体
%                  * processing_info: 处理信息结构体
%
% 处理流程:
%   1. 输入参数验证
%   2. 调用OFDM调制函数生成资源网格和时域信号
%   3. 对时域信号进行插值处理
%   4. 添加循环前缀和保护间隔
%   5. 信号质量检查和统计
%   6. 构建输出结构体
%
% 技术细节:
%   - 使用现有的generateCrossOFDM函数进行OFDM调制
%   - 使用interpo函数进行信号插值
%   - 添加循环保护间隔以处理边界效应
%   - 支持多频带OFDM信号生成
%
% 作者: OFDM接收机开发
% 日期: 2025年7月
% 版本: 1.0

%% 1. 输入参数验证
if nargin < 1
    error('ofdm_generate_local_signals: 必须提供仿真设置参数');
end

% 检查必要字段
required_fields = {'code', 'Nu', 'numBand', 'NFFT', 'nSymbol', 'famp', 'extra'};
for i = 1:length(required_fields)
    if ~isfield(simSettings, required_fields{i})
        error('ofdm_generate_local_signals: simSettings缺少必要字段 %s', required_fields{i});
    end
end

fprintf('  - 输入参数验证通过\n');

% 提取关键参数
code = simSettings.code;
Nu = simSettings.Nu;
numBand = simSettings.numBand;
NFFT = simSettings.NFFT;
nSymbol = simSettings.nSymbol;
famp = simSettings.famp;
extra = simSettings.extra;

fprintf('    * 本地码长度: %d\n', length(code));
fprintf('    * 子载波数: %d\n', Nu);
fprintf('    * 频带数: %d\n', numBand);
fprintf('    * FFT点数: %d\n', NFFT);
fprintf('    * OFDM符号数: %d\n', nSymbol);

%% 2. 码序列分析和验证
fprintf('  - 分析本地码序列...\n');

% 检查码序列的基本属性
code_length = length(code);
code_unique_values = unique(code);
code_mean = mean(code);
code_std = std(code);

fprintf('    * 码长度: %d\n', code_length);
fprintf('    * 唯一值数量: %d\n', length(code_unique_values));
fprintf('    * 码均值: %.4f\n', code_mean);
fprintf('    * 码标准差: %.4f\n', code_std);

% 检查码的取值范围
if all(ismember(code_unique_values, [-1, 1]))
    code_type = 'bipolar'; % 双极性码 (-1, +1)
    fprintf('    * 码类型: 双极性码 (-1, +1)\n');
elseif all(ismember(code_unique_values, [0, 1]))
    code_type = 'unipolar'; % 单极性码 (0, 1)
    fprintf('    * 码类型: 单极性码 (0, 1)\n');
else
    code_type = 'other';
    fprintf('    * 码类型: 其他类型，取值范围 [%.2f, %.2f]\n', min(code), max(code));
end

% 检查码长度与OFDM参数的匹配性
total_subcarriers = Nu * numBand;
required_symbols = ceil(code_length / total_subcarriers);

if required_symbols > nSymbol
    warning('ofdm_generate_local_signals: 码长度过长，需要 %d 个OFDM符号，但只有 %d 个', ...
            required_symbols, nSymbol);
elseif required_symbols < nSymbol
    fprintf('    * 信息: 码将在 %d 个OFDM符号中重复或填充\n', nSymbol);
end

%% 3. 生成OFDM调制信号
fprintf('  - 生成OFDM调制信号...\n');

try
    % 调用现有的OFDM调制函数
    [resource_grid, rawSignal] = generateCrossOFDM(simSettings);
    
    fprintf('    * OFDM调制成功\n');
    fprintf('    * 资源网格大小: %dx%dx%d\n', size(resource_grid));
    fprintf('    * 时域信号大小: %dx%d\n', size(rawSignal));
    
catch ME
    error('ofdm_generate_local_signals: OFDM调制失败 - %s', ME.message);
end

% 分离实部和虚部
% 数据验证和清理
if ~isnumeric(rawSignal)
    error('ofdm_generate_local_signals: rawSignal不是数值类型');
end

% 处理非有限值
rawSignal(~isfinite(rawSignal)) = 0;

yr_original = real(rawSignal);
yi_original = imag(rawSignal);

% 检查生成信号的质量
signal_power_real = mean(yr_original.^2, 'all');
signal_power_imag = mean(yi_original.^2, 'all');
total_signal_power = signal_power_real + signal_power_imag;

fprintf('    * 实部信号功率: %.6f\n', signal_power_real);
fprintf('    * 虚部信号功率: %.6f\n', signal_power_imag);
fprintf('    * 总信号功率: %.6f\n', total_signal_power);

%% 4. 信号插值处理
fprintf('  - 进行信号插值处理...\n');

if famp <= 1
    warning('ofdm_generate_local_signals: 插值倍数 <= 1，跳过插值处理');
    yr_interpolated = yr_original;
    yi_interpolated = yi_original;
else
    try
        % 对实部和虚部分别进行插值
        yr_interpolated = interpo_fixed(yr_original, famp);
        yi_interpolated = interpo_fixed(yi_original, famp);
        
        fprintf('    * 插值倍数: %d\n', famp);
        fprintf('    * 插值后信号大小: %dx%d\n', size(yr_interpolated));
        
        % 检查插值后的信号质量
        interp_power_real = mean(yr_interpolated.^2, 'all');
        interp_power_imag = mean(yi_interpolated.^2, 'all');
        
        fprintf('    * 插值后实部功率: %.6f\n', interp_power_real);
        fprintf('    * 插值后虚部功率: %.6f\n', interp_power_imag);
        
    catch ME
        error('ofdm_generate_local_signals: 信号插值失败 - %s', ME.message);
    end
end

%% 5. 添加循环保护间隔
fprintf('  - 添加循环保护间隔...\n');

if extra <= 0
    warning('ofdm_generate_local_signals: 保护间隔长度 <= 0，跳过保护间隔添加');
    yr0 = yr_interpolated;
    yi0 = yi_interpolated;
else
    try
        % 计算保护间隔长度
        guard_length = extra * famp;
        signal_length = size(yr_interpolated, 2);
        
        if guard_length >= signal_length
            warning('ofdm_generate_local_signals: 保护间隔长度过长，调整为信号长度的1/4');
            guard_length = floor(signal_length / 4);
        end
        
        % 添加循环前缀和后缀
        % 前缀：取信号末尾部分
        prefix_real = yr_interpolated(:, end-guard_length+1:end);
        prefix_imag = yi_interpolated(:, end-guard_length+1:end);
        
        % 后缀：取信号开头部分
        suffix_real = yr_interpolated(:, 1:guard_length);
        suffix_imag = yi_interpolated(:, 1:guard_length);
        
        % 组合最终信号
        yr0 = [prefix_real, yr_interpolated, suffix_real];
        yi0 = [prefix_imag, yi_interpolated, suffix_imag];
        
        fprintf('    * 保护间隔长度: %d 采样点\n', guard_length);
        fprintf('    * 最终信号长度: %d 采样点\n', size(yr0, 2));
        
    catch ME
        error('ofdm_generate_local_signals: 添加保护间隔失败 - %s', ME.message);
    end
end

%% 6. 信号质量最终检查
fprintf('  - 最终信号质量检查...\n');

% 检查每个频带的信号质量
for bandID = 1:numBand
    band_real = yr0(bandID, :);
    band_imag = yi0(bandID, :);
    
    % 计算统计量
    real_rms = rms(band_real);
    imag_rms = rms(band_imag);
    real_peak = max(abs(band_real));
    imag_peak = max(abs(band_imag));
    
    % 计算峰均比 (PAPR)
    real_papr = 20 * log10(real_peak / real_rms);
    imag_papr = 20 * log10(imag_peak / imag_rms);
    
    fprintf('    * 频带 %d: 实部RMS=%.4f, 虚部RMS=%.4f\n', bandID, real_rms, imag_rms);
    fprintf('    * 频带 %d: 实部PAPR=%.2fdB, 虚部PAPR=%.2fdB\n', bandID, real_papr, imag_papr);
    
    % 检查异常值
    if real_rms < 1e-10 || imag_rms < 1e-10
        warning('ofdm_generate_local_signals: 频带 %d 信号功率过低', bandID);
    end
    
    if real_papr > 15 || imag_papr > 15
        warning('ofdm_generate_local_signals: 频带 %d 峰均比过高', bandID);
    end
end

%% 7. 构建输出结构体
fprintf('  - 构建输出结构体...\n');

local_signals = struct();

% 主要信号输出
local_signals.yr0 = yr0;  % 实部本地信号
local_signals.yi0 = yi0;  % 虚部本地信号
local_signals.resource_grid = resource_grid;  % OFDM资源网格
local_signals.rawSignal = rawSignal;  % 原始OFDM时域信号

% 中间处理结果（用于调试和分析）
local_signals.intermediate = struct();
local_signals.intermediate.yr_original = yr_original;
local_signals.intermediate.yi_original = yi_original;
local_signals.intermediate.yr_interpolated = yr_interpolated;
local_signals.intermediate.yi_interpolated = yi_interpolated;

% 码信息
local_signals.code_info = struct();
local_signals.code_info.code = code;
local_signals.code_info.length = code_length;
local_signals.code_info.type = code_type;
local_signals.code_info.unique_values = code_unique_values;
local_signals.code_info.mean = code_mean;
local_signals.code_info.std = code_std;
local_signals.code_info.required_symbols = required_symbols;

% OFDM参数信息
local_signals.ofdm_info = struct();
local_signals.ofdm_info.Nu = Nu;
local_signals.ofdm_info.numBand = numBand;
local_signals.ofdm_info.NFFT = NFFT;
local_signals.ofdm_info.nSymbol = nSymbol;
local_signals.ofdm_info.total_subcarriers = total_subcarriers;

% 处理参数信息
local_signals.processing_info = struct();
local_signals.processing_info.famp = famp;
local_signals.processing_info.extra = extra;
local_signals.processing_info.guard_length = extra * famp;
local_signals.processing_info.final_length = size(yr0, 2);
local_signals.processing_info.interpolation_applied = (famp > 1);
local_signals.processing_info.guard_interval_applied = (extra > 0);

% 信号质量统计
local_signals.quality_metrics = struct();
for bandID = 1:numBand
    band_name = sprintf('band_%d', bandID);
    
    band_real = yr0(bandID, :);
    band_imag = yi0(bandID, :);
    
    local_signals.quality_metrics.(band_name) = struct();
    local_signals.quality_metrics.(band_name).real_rms = rms(band_real);
    local_signals.quality_metrics.(band_name).imag_rms = rms(band_imag);
    local_signals.quality_metrics.(band_name).real_peak = max(abs(band_real));
    local_signals.quality_metrics.(band_name).imag_peak = max(abs(band_imag));
    local_signals.quality_metrics.(band_name).real_papr = 20 * log10(max(abs(band_real)) / rms(band_real));
    local_signals.quality_metrics.(band_name).imag_papr = 20 * log10(max(abs(band_imag)) / rms(band_imag));
end

% 整体信号统计
local_signals.overall_stats = struct();
local_signals.overall_stats.total_power_real = mean(yr0.^2, 'all');
local_signals.overall_stats.total_power_imag = mean(yi0.^2, 'all');
local_signals.overall_stats.total_power = local_signals.overall_stats.total_power_real + ...
                                         local_signals.overall_stats.total_power_imag;
local_signals.overall_stats.snr_estimate = 10 * log10(local_signals.overall_stats.total_power / 1e-12);

fprintf('  - 本地信号生成完成\n');

%% 8. 最终验证
fprintf('  - 最终验证...\n');

% 验证输出信号的完整性
if size(yr0, 1) ~= numBand || size(yi0, 1) ~= numBand
    error('ofdm_generate_local_signals: 输出信号频带数不匹配');
end

if size(yr0, 2) ~= size(yi0, 2)
    error('ofdm_generate_local_signals: 实部和虚部信号长度不匹配');
end

% 检查信号中是否包含无效值
if any(~isfinite(yr0), 'all') || any(~isfinite(yi0), 'all')
    error('ofdm_generate_local_signals: 输出信号包含无效值（NaN或Inf）');
end

% 检查资源网格的完整性
expected_grid_size = [nSymbol, NFFT, numBand];
actual_grid_size = size(resource_grid);
if ~isequal(expected_grid_size, actual_grid_size)
    warning('ofdm_generate_local_signals: 资源网格大小不匹配，期望 %s，实际 %s', ...
            mat2str(expected_grid_size), mat2str(actual_grid_size));
end

fprintf('    * 验证通过，本地信号生成模块完成\n');

end