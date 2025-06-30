function verify_code_correlation(simSettings, prn_idx, compare_prns)
% 验证本地码的自相关特性
% 输入:
%   simSettings - 仿真设置参数
%   prn_idx - 要分析的PRN码索引，默认为1
%   compare_prns - 是否比较不同PRN码的自相关特性，默认为false
%
% 本函数验证码的自相关特性，即当码对齐时（零延迟）相关值很大，
% 而当有任何移位时相关值应该很小。这是跟踪系统中非常重要的特性。

if nargin < 1
    % 如果没有提供仿真设置，则使用默认设置
    simSettings = init(-20); % SNR=-20dB
end

if nargin < 2
    prn_idx = 1; % 默认使用第一个PRN码
end

if nargin < 3
    compare_prns = false; % 默认不比较不同PRN码
end

% 加载本地码
load('weil10230_signed.mat');

% 获取PRN码
code = weil10230_signed(prn_idx, :);
simSettings.code = code;

fprintf('分析PRN %d的自相关特性...\n', prn_idx);

% 生成OFDM调制后的信号
[resource_grid, rawSignal] = generateCrossOFDM(simSettings);
yr0 = real(rawSignal);
yi0 = imag(rawSignal);
yr0 = interpo(yr0, simSettings.famp);
yi0 = interpo(yi0, simSettings.famp);
yr0 = [yr0(:,end-simSettings.extra*simSettings.famp+1:end) yr0 yr0(:,1:1+simSettings.extra*simSettings.famp)];
yi0 = [yi0(:,end-simSettings.extra*simSettings.famp+1:end) yi0 yi0(:,1:1+simSettings.extra*simSettings.famp)];

% 计算自相关
max_lag = min(5000, floor(length(yr0)/2)); % 最大延迟
auto_corr = zeros(simSettings.numBand, 2*max_lag+1);

for bandID = 1:simSettings.numBand
    % 实部自相关
    signal_real = yr0(bandID, :);
    auto_corr_real = xcorr(signal_real, signal_real, max_lag);
    
    % 虚部自相关
    signal_imag = yi0(bandID, :);
    auto_corr_imag = xcorr(signal_imag, signal_imag, max_lag);
    
    % 合并实部和虚部的自相关
    auto_corr(bandID, :) = auto_corr_real + auto_corr_imag;
end

% 计算每个频带的平均自相关
avg_auto_corr = mean(auto_corr, 1);

% 归一化自相关
norm_auto_corr = avg_auto_corr / max(avg_auto_corr);

% 找到主峰值（零延迟）
[main_peak, peak_idx] = max(norm_auto_corr);
assert(peak_idx == max_lag + 1, '主峰值应该在零延迟处');

% 找到最大旁瓣（非零延迟）
temp_corr = norm_auto_corr;
exclude_range = max(1, peak_idx-10):min(length(temp_corr), peak_idx+10);
temp_corr(exclude_range) = 0;
[side_lobe, side_idx] = max(temp_corr);

% 计算主峰与旁瓣的比值
peak_to_sidelobe = main_peak / side_lobe;

% 可视化结果
figure('Color', 'w', 'Position', [100, 100, 1200, 600]);

% 绘制完整的自相关
subplot(2, 2, 1);
plot(-max_lag:max_lag, norm_auto_corr, 'LineWidth', 1.5);
title('本地码自相关特性');
xlabel('延迟');
ylabel('归一化相关值');
grid on;

% 放大显示主峰附近的区域
subplot(2, 2, 2);
zoom_range = 100;
plot(-zoom_range:zoom_range, norm_auto_corr(max_lag+1-zoom_range:max_lag+1+zoom_range), 'LineWidth', 1.5);
title('主峰附近放大');
xlabel('延迟');
ylabel('归一化相关值');
grid on;

% 绘制自相关的对数尺度
subplot(2, 2, 3);
semilogy(-max_lag:max_lag, abs(norm_auto_corr), 'LineWidth', 1.5);
title('自相关特性（对数尺度）');
xlabel('延迟');
ylabel('归一化相关值（对数）');
grid on;

% 绘制自相关的频谱
subplot(2, 2, 4);
fft_size = 2^nextpow2(length(norm_auto_corr));
auto_corr_fft = abs(fftshift(fft(norm_auto_corr, fft_size)));
auto_corr_fft = auto_corr_fft / max(auto_corr_fft);
freq_axis = (-fft_size/2:fft_size/2-1) / fft_size;
plot(freq_axis, auto_corr_fft, 'LineWidth', 1.5);
title('自相关频谱');
xlabel('归一化频率');
ylabel('幅度');
grid on;

% 分析不同频带的自相关特性
figure('Color', 'w', 'Position', [100, 100, 1200, 800]);
for bandID = 1:simSettings.numBand
    subplot(simSettings.numBand, 1, bandID);
    band_auto_corr = auto_corr(bandID, :);
    band_auto_corr = band_auto_corr / max(band_auto_corr);
    plot(-max_lag:max_lag, band_auto_corr, 'LineWidth', 1.5);
    title(sprintf('频带 %d 的自相关特性', bandID));
    xlabel('延迟');
    ylabel('归一化相关值');
    grid on;
end

% 分析不同移位的自相关特性
figure('Color', 'w', 'Position', [100, 100, 1200, 600]);

% 选择几个代表性的移位
shifts = [0, 1, 2, 5, 10, 20, 50, 100];
num_shifts = length(shifts);
colors = jet(num_shifts);

% 为每个移位计算相关值
subplot(1, 2, 1);
hold on;
for i = 1:num_shifts
    shift = shifts(i);
    shifted_idx = max_lag + 1 + shift;
    if shifted_idx <= length(norm_auto_corr)
        stem(shift, norm_auto_corr(shifted_idx), 'Color', colors(i,:), 'LineWidth', 2, 'Marker', 'o', 'MarkerSize', 8);
    end
end
title('不同移位的自相关值');
xlabel('移位');
ylabel('归一化相关值');
grid on;
hold off;

% 对数尺度
subplot(1, 2, 2);
hold on;
for i = 1:num_shifts
    shift = shifts(i);
    shifted_idx = max_lag + 1 + shift;
    if shifted_idx <= length(norm_auto_corr)
        stem(shift, abs(norm_auto_corr(shifted_idx)), 'Color', colors(i,:), 'LineWidth', 2, 'Marker', 'o', 'MarkerSize', 8);
    end
end
set(gca, 'YScale', 'log');
title('不同移位的自相关值（对数尺度）');
xlabel('移位');
ylabel('归一化相关值（对数）');
grid on;
hold off;

% 打印结果
fprintf('主峰值（零延迟）: %.4f\n', main_peak);
fprintf('最大旁瓣值（延迟 = %d）: %.4f\n', side_idx - max_lag - 1, side_lobe);
fprintf('主峰与旁瓣的比值: %.4f\n', peak_to_sidelobe);

% 分析不同移位的相关值
fprintf('\n不同移位的相关值:\n');
fprintf('移位\t相关值\t\t相对于主峰的比值\n');
for shift = [0, 1, 2, 5, 10, 20, 50, 100]
    shifted_idx = max_lag + 1 + shift;
    if shifted_idx <= length(norm_auto_corr)
        corr_val = norm_auto_corr(shifted_idx);
        ratio = corr_val / main_peak;
        fprintf('%d\t%.6f\t%.6f\n', shift, corr_val, ratio);
    end
end

% 分析码的周期性
figure('Color', 'w', 'Position', [100, 100, 1200, 400]);
code_length = length(code);
code_periods = floor(length(yr0(1,:)) / code_length);

% 计算多个码周期的自相关
if code_periods > 1
    extended_code = repmat(code, 1, code_periods);
    extended_code_corr = xcorr(extended_code, extended_code, code_length);
    extended_code_corr = extended_code_corr / max(extended_code_corr);
    
    plot(-code_length:code_length, extended_code_corr, 'LineWidth', 1.5);
    title('多个码周期的自相关');
    xlabel('延迟');
    ylabel('归一化相关值');
    grid on;
    
    % 分析码的周期性
    [peaks, locs] = findpeaks(extended_code_corr, 'MinPeakHeight', 0.5);
    if length(peaks) > 1
        avg_period = mean(diff(locs));
        fprintf('\n码的周期性分析:\n');
        fprintf('检测到的峰值数量: %d\n', length(peaks));
        fprintf('平均周期: %.2f\n', avg_period);
    end
end

% 分析OFDM调制对自相关特性的影响
figure('Color', 'w', 'Position', [100, 100, 1200, 800]);

% 原始码的自相关
subplot(2, 1, 1);
code_corr = xcorr(code, code, min(max_lag, length(code)-1));
code_corr = code_corr / max(code_corr);
plot(-min(max_lag, length(code)-1):min(max_lag, length(code)-1), code_corr, 'LineWidth', 1.5);
title('原始码的自相关');
xlabel('延迟');
ylabel('归一化相关值');
grid on;

% OFDM调制后的自相关
subplot(2, 1, 2);
plot(-max_lag:max_lag, norm_auto_corr, 'LineWidth', 1.5);
title('OFDM调制后的自相关');
xlabel('延迟');
ylabel('归一化相关值');
grid on;

% 比较原始码和OFDM调制后的自相关特性
[~, code_peak_idx] = max(code_corr);
code_temp = code_corr;
code_exclude_range = max(1, code_peak_idx-10):min(length(code_temp), code_peak_idx+10);
code_temp(code_exclude_range) = 0;
code_side_lobe = max(code_temp);
code_peak_to_sidelobe = 1 / code_side_lobe;

fprintf('\n自相关特性比较:\n');
fprintf('原始码的主峰与旁瓣比值: %.4f\n', code_peak_to_sidelobe);
fprintf('OFDM调制后的主峰与旁瓣比值: %.4f\n', peak_to_sidelobe);
fprintf('改善比例: %.4f\n', peak_to_sidelobe / code_peak_to_sidelobe);

% 分析不同子载波配置对自相关特性的影响
original_Nu = simSettings.Nu;
original_numBand = simSettings.numBand;

if simSettings.Nu > 10 && simSettings.numBand > 1
    fprintf('\n分析不同子载波配置对自相关特性的影响...\n');
    
    % 测试不同的子载波数量
    Nu_values = [original_Nu/2, original_Nu, original_Nu*1.5];
    Nu_values = round(Nu_values);
    Nu_values = Nu_values(Nu_values > 0);
    
    % 测试不同的频带数量
    numBand_values = [1, original_numBand, original_numBand*2];
    numBand_values = numBand_values(numBand_values > 0);
    
    % 初始化结果存储
    peak_to_sidelobe_values = zeros(length(Nu_values), length(numBand_values));
    
    figure('Color', 'w', 'Position', [100, 100, 1200, 800]);
    plot_idx = 1;
    
    for i = 1:length(Nu_values)
        for j = 1:length(numBand_values)
            % 修改配置
            test_settings = simSettings;
            test_settings.Nu = Nu_values(i);
            test_settings.numBand = numBand_values(j);
            
            % 重新计算nSymbol
            test_settings.nSymbol = ceil(length(code)/(test_settings.Nu*test_settings.numBand));
            
            % 生成OFDM调制后的信号
            [~, test_rawSignal] = generateCrossOFDM(test_settings);
            test_yr0 = real(test_rawSignal);
            test_yi0 = imag(test_rawSignal);
            test_yr0 = interpo(test_yr0, test_settings.famp);
            test_yi0 = interpo(test_yi0, test_settings.famp);
            test_yr0 = [test_yr0(:,end-test_settings.extra*test_settings.famp+1:end) test_yr0 test_yr0(:,1:1+test_settings.extra*test_settings.famp)];
            test_yi0 = [test_yi0(:,end-test_settings.extra*test_settings.famp+1:end) test_yi0 test_yi0(:,1:1+test_settings.extra*test_settings.famp)];
            
            % 计算自相关
            test_max_lag = min(1000, floor(length(test_yr0)/2));
            test_auto_corr = zeros(test_settings.numBand, 2*test_max_lag+1);
            
            for bandID = 1:test_settings.numBand
                % 实部自相关
                signal_real = test_yr0(bandID, :);
                auto_corr_real = xcorr(signal_real, signal_real, test_max_lag);
                
                % 虚部自相关
                signal_imag = test_yi0(bandID, :);
                auto_corr_imag = xcorr(signal_imag, signal_imag, test_max_lag);
                
                % 合并实部和虚部的自相关
                test_auto_corr(bandID, :) = auto_corr_real + auto_corr_imag;
            end
            
            % 计算每个频带的平均自相关
            test_avg_auto_corr = mean(test_auto_corr, 1);
            
            % 归一化自相关
            test_norm_auto_corr = test_avg_auto_corr / max(test_avg_auto_corr);
            
            % 找到主峰值（零延迟）
            [test_main_peak, test_peak_idx] = max(test_norm_auto_corr);
            
            % 找到最大旁瓣（非零延迟）
            test_temp_corr = test_norm_auto_corr;
            test_exclude_range = max(1, test_peak_idx-10):min(length(test_temp_corr), test_peak_idx+10);
            test_temp_corr(test_exclude_range) = 0;
            [test_side_lobe, ~] = max(test_temp_corr);
            
            % 计算主峰与旁瓣的比值
            test_peak_to_sidelobe = test_main_peak / test_side_lobe;
            peak_to_sidelobe_values(i, j) = test_peak_to_sidelobe;
            
            % 绘制自相关
            subplot(length(Nu_values), length(numBand_values), plot_idx);
            plot(-test_max_lag:test_max_lag, test_norm_auto_corr, 'LineWidth', 1.5);
            title(sprintf('Nu=%d, numBand=%d', Nu_values(i), numBand_values(j)));
            xlabel('延迟');
            ylabel('归一化相关值');
            grid on;
            
            fprintf('Nu=%d, numBand=%d: 主峰与旁瓣比值 = %.4f\n', Nu_values(i), numBand_values(j), test_peak_to_sidelobe);
            
            plot_idx = plot_idx + 1;
        end
    end
    
    % 恢复原始设置
    simSettings.Nu = original_Nu;
    simSettings.numBand = original_numBand;
    
    % 绘制热图
    figure('Color', 'w', 'Position', [100, 100, 600, 500]);
    imagesc(numBand_values, Nu_values, peak_to_sidelobe_values);
    colorbar;
    title('不同子载波配置的主峰与旁瓣比值');
    xlabel('频带数量');
    ylabel('每个频带的子载波数量');
    colormap('jet');
    
    % 找到最佳配置
    [max_val, max_idx] = max(peak_to_sidelobe_values(:));
    [max_i, max_j] = ind2sub(size(peak_to_sidelobe_values), max_idx);
    fprintf('\n最佳子载波配置: Nu=%d, numBand=%d, 主峰与旁瓣比值=%.4f\n', ...
        Nu_values(max_i), numBand_values(max_j), max_val);
end

% 如果需要比较不同PRN码的自相关特性
if compare_prns
    fprintf('\n比较不同PRN码的自相关特性...\n');
    
    % 选择要比较的PRN码
    prn_indices = [1, 2, 3, 4, 5];
    num_prns = length(prn_indices);
    
    % 初始化结果存储
    prn_peak_to_sidelobe = zeros(num_prns, 1);
    
    figure('Color', 'w', 'Position', [100, 100, 1200, 800]);
    
    for i = 1:num_prns
        prn_i = prn_indices(i);
        
        % 获取PRN码
        prn_code = weil10230_signed(prn_i, :);
        test_settings = simSettings;
        test_settings.code = prn_code;
        
        % 生成OFDM调制后的信号
        [~, prn_rawSignal] = generateCrossOFDM(test_settings);
        prn_yr0 = real(prn_rawSignal);
        prn_yi0 = imag(prn_rawSignal);
        prn_yr0 = interpo(prn_yr0, test_settings.famp);
        prn_yi0 = interpo(prn_yi0, test_settings.famp);
        prn_yr0 = [prn_yr0(:,end-test_settings.extra*test_settings.famp+1:end) prn_yr0 prn_yr0(:,1:1+test_settings.extra*test_settings.famp)];
        prn_yi0 = [prn_yi0(:,end-test_settings.extra*test_settings.famp+1:end) prn_yi0 prn_yi0(:,1:1+test_settings.extra*test_settings.famp)];
        
        % 计算自相关
        prn_max_lag = min(1000, floor(length(prn_yr0)/2));
        prn_auto_corr = zeros(test_settings.numBand, 2*prn_max_lag+1);
        
        for bandID = 1:test_settings.numBand
            % 实部自相关
            signal_real = prn_yr0(bandID, :);
            auto_corr_real = xcorr(signal_real, signal_real, prn_max_lag);
            
            % 虚部自相关
            signal_imag = prn_yi0(bandID, :);
            auto_corr_imag = xcorr(signal_imag, signal_imag, prn_max_lag);
            
            % 合并实部和虚部的自相关
            prn_auto_corr(bandID, :) = auto_corr_real + auto_corr_imag;
        end
        
        % 计算每个频带的平均自相关
        prn_avg_auto_corr = mean(prn_auto_corr, 1);
        
        % 归一化自相关
        prn_norm_auto_corr = prn_avg_auto_corr / max(prn_avg_auto_corr);
        
        % 找到主峰值（零延迟）
        [prn_main_peak, prn_peak_idx] = max(prn_norm_auto_corr);
        
        % 找到最大旁瓣（非零延迟）
        prn_temp_corr = prn_norm_auto_corr;
        prn_exclude_range = max(1, prn_peak_idx-10):min(length(prn_temp_corr), prn_peak_idx+10);
        prn_temp_corr(prn_exclude_range) = 0;
        [prn_side_lobe, ~] = max(prn_temp_corr);
        
        % 计算主峰与旁瓣的比值
        prn_peak_to_sidelobe(i) = prn_main_peak / prn_side_lobe;
        
        % 绘制自相关
        subplot(num_prns, 1, i);
        plot(-prn_max_lag:prn_max_lag, prn_norm_auto_corr, 'LineWidth', 1.5);
        title(sprintf('PRN %d 的自相关特性', prn_i));
        xlabel('延迟');
        ylabel('归一化相关值');
        grid on;
        
        fprintf('PRN %d: 主峰与旁瓣比值 = %.4f\n', prn_i, prn_peak_to_sidelobe(i));
    end
    
    % 找到最佳PRN码
    [max_prn_val, max_prn_idx] = max(prn_peak_to_sidelobe);
    fprintf('\n最佳PRN码: PRN %d, 主峰与旁瓣比值=%.4f\n', ...
        prn_indices(max_prn_idx), max_prn_val);
    
    % 计算互相关
    fprintf('\n计算PRN码之间的互相关...\n');
    
    % 初始化互相关矩阵
    cross_corr_matrix = zeros(num_prns, num_prns);
    
    for i = 1:num_prns
        for j = 1:num_prns
            if i ~= j
                % 获取PRN码
                prn_i_code = weil10230_signed(prn_indices(i), :);
                prn_j_code = weil10230_signed(prn_indices(j), :);
                
                % 计算互相关
                cross_corr = xcorr(prn_i_code, prn_j_code, 0);
                auto_corr_i = xcorr(prn_i_code, prn_i_code, 0);
                auto_corr_j = xcorr(prn_j_code, prn_j_code, 0);
                
                % 归一化互相关
                norm_cross_corr = cross_corr / sqrt(auto_corr_i * auto_corr_j);
                
                cross_corr_matrix(i, j) = abs(norm_cross_corr);
            end
        end
    end
    
    % 绘制互相关矩阵
    figure('Color', 'w', 'Position', [100, 100, 600, 500]);
    imagesc(prn_indices, prn_indices, cross_corr_matrix);
    colorbar;
    title('PRN码之间的互相关');
    xlabel('PRN索引');
    ylabel('PRN索引');
    colormap('jet');
    
    % 找到最小互相关
    cross_corr_matrix(cross_corr_matrix == 0) = 1; % 忽略对角线上的零
    [min_cross_val, min_cross_idx] = min(cross_corr_matrix(:));
    [min_i, min_j] = ind2sub(size(cross_corr_matrix), min_cross_idx);
    fprintf('最小互相关: PRN %d 与 PRN %d, 互相关值=%.4f\n', ...
        prn_indices(min_i), prn_indices(min_j), min_cross_val);
end

end
