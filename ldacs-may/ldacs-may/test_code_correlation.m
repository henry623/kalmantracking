% 测试本地码的自相关特性
clear all; close all;

% 加载仿真设置
simSettings = init(-20); % SNR=-20dB

% 选择测试模式
fprintf('请选择测试模式：\n');
fprintf('1. 基本自相关分析（单个PRN码）\n');
fprintf('2. 比较不同PRN码的自相关特性\n');
fprintf('3. 分析不同子载波配置的影响\n');
fprintf('4. 全面分析（包括所有上述测试）\n');

test_mode = input('请输入选项（1-4）: ');

% 选择PRN码
prn_idx = 1; % 默认使用第一个PRN码
if test_mode ~= 1
    prn_idx = input('请输入要分析的PRN码索引（默认为1）: ');
    if isempty(prn_idx)
        prn_idx = 1;
    end
end

% 根据测试模式调用验证函数
switch test_mode
    case 1
        % 基本自相关分析
        verify_code_correlation(simSettings, prn_idx, false);
        
    case 2
        % 比较不同PRN码的自相关特性
        verify_code_correlation(simSettings, prn_idx, true);
        
    case 3
        % 分析不同子载波配置的影响
        % 修改simSettings以确保可以测试不同配置
        original_Nu = simSettings.Nu;
        original_numBand = simSettings.numBand;
        
        % 确保Nu和numBand足够大以便进行比较
        if simSettings.Nu < 10
            simSettings.Nu = 50;
        end
        if simSettings.numBand < 2
            simSettings.numBand = 4;
        end
        
        verify_code_correlation(simSettings, prn_idx, false);
        
        % 恢复原始设置
        simSettings.Nu = original_Nu;
        simSettings.numBand = original_numBand;
        
    case 4
        % 全面分析
        % 确保Nu和numBand足够大以便进行比较
        original_Nu = simSettings.Nu;
        original_numBand = simSettings.numBand;
        
        if simSettings.Nu < 10
            simSettings.Nu = 50;
        end
        if simSettings.numBand < 2
            simSettings.numBand = 4;
        end
        
        verify_code_correlation(simSettings, prn_idx, true);
        
        % 恢复原始设置
        simSettings.Nu = original_Nu;
        simSettings.numBand = original_numBand;
        
    otherwise
        fprintf('无效的选项，执行基本自相关分析。\n');
        verify_code_correlation(simSettings, prn_idx, false);
end

% 输出提示信息
fprintf('\n自相关分析完成。请查看图表了解详细结果。\n');
fprintf('主要关注点：\n');
fprintf('1. 零延迟（完全对齐）时的相关峰值\n');
fprintf('2. 非零延迟时相关值的快速衰减\n');
fprintf('3. 主峰与最大旁瓣的比值\n');
fprintf('4. 原始码与OFDM调制后的自相关特性比较\n');

if test_mode >= 2
    fprintf('5. 不同PRN码的自相关特性比较\n');
    fprintf('6. PRN码之间的互相关分析\n');
end

if test_mode == 3 || test_mode == 4
    fprintf('7. 不同子载波配置对自相关特性的影响\n');
    fprintf('8. 最佳子载波配置建议\n');
end
