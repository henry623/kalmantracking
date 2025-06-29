% 5月对跟踪代码的修改
clear all; close all;

% 读取参数
SNR = -20; % 设置当前仿真的信噪比
simSettings = init(SNR); % 设置参数
flag = 1; % 要不要重新生成信号

% 生成信号
[grid, rawSignal] = generateCrossOFDM(simSettings); % 资源网格，行是时间，列是频率
yr0 = real(rawSignal); % 实部信号
yi0 = imag(rawSignal); % 虚部信号
yr0 = interpo(yr0, simSettings.famp); % 实部时域采样结果
yi0 = interpo(yi0, simSettings.famp); % 虚部时域采样结果
yr0 = [yr0(:,end-simSettings.extra*simSettings.famp+1:end) yr0 yr0(:,1:1+simSettings.extra*simSettings.famp)];
yi0 = [yi0(:,end-simSettings.extra*simSettings.famp+1:end) yi0 yi0(:,1:1+simSettings.extra*simSettings.famp)];

% 生成或加载信号
if flag == 1
    for bandID = 1:simSettings.numBand
        yr = yr0(bandID,:);
        yi = yi0(bandID,:);
        record_E(:,:,bandID) = generateSignal(simSettings, yr, yi);
    end
    filenm1 = ['output/signal_SNR'  num2str(SNR) '.mat'];
    save(filenm1, 'record_E');
else
    load(['output/signal_SNR'  num2str(SNR) '.mat']);
end

% 跟踪信号
output = track(simSettings, yr0, yi0);
savenm = ['output/tracked_SNR'  num2str(SNR) '.mat'];
save(savenm, 'output');

% 加载理想跟踪结果（SNR = inf 的情况）
load("output/tracked_SNRinf.mat");
ideal_codePhase = output.OutCodePhase;
ideal_carrPhase = output.OutCarrPhase;
ideal_codeFreq = output.OutCodeFreq;
ideal_carrFreq = output.OutCarrFreq;

% 加载当前 SNR 的跟踪结果
loadnm = ['output/tracked_SNR'  num2str(SNR) '.mat'];
load(loadnm);
codePhase = output.OutCodePhase;
carrPhase = output.OutCarrPhase;
codeFreq = output.OutCodeFreq;
carrFreq = output.OutCarrFreq;

% 提取相关数据
I_E = output.I_E;
Q_E = output.Q_E;
I_P = output.I_P;
Q_P = output.Q_P;
I_L = output.I_L;
Q_L = output.Q_L;

% 计算各通道的幅度
Amp_E = sqrt(I_E.^2 + Q_E.^2); % 早（Early）通道幅度
Amp_P = sqrt(I_P.^2 + Q_P.^2); % 准时（Prompt）通道幅度
Amp_L = sqrt(I_L.^2 + Q_L.^2); % 晚（Late）通道幅度

% 创建时间轴
num_samples = length(Amp_E);
t = 1:num_samples;

% 计算相位误差和频率误差
code_phase_error = codePhase - ideal_codePhase;
carr_phase_error = carrPhase - ideal_carrPhase;
code_freq_error = codeFreq - ideal_codeFreq;
carr_freq_error = carrFreq - ideal_carrFreq;

% 计算 RMSE（从第200个周期开始）
start_index = 200;
rmse_codePhase = sqrt(mean(code_phase_error(start_index:end).^2));
rmse_carrPhase = sqrt(mean(carr_phase_error(start_index:end).^2));
rmse_codeFreq = sqrt(mean(code_freq_error(start_index:end).^2));
rmse_carrFreq = sqrt(mean(carr_freq_error(start_index:end).^2));

rmse_codePhase_dis = simSettings.Ts * simSettings.c * rmse_codePhase;

% 绘制结果
plot_results(t, I_P, 'Bits of the navigation message', 't', 'I prompt', []);
plot_results(t, [Amp_E; Amp_P; Amp_L], 'Correlation', 'Sample Index', 'Amplitude', {'Early', 'Prompt', 'Late'});
plot_results(t, code_phase_error, 'Code Phase Error', 'Time (s)', 'Phase Error', []);
plot_results(t, carr_phase_error, 'Carrier Phase Error', 'Time (s)', 'Phase Error', []);
plot_results(t, code_freq_error, 'Code Frequency Error', 'Time (s)', 'Frequency Error (Hz)', []);
plot_results(t, carr_freq_error, 'Carrier Frequency Error', 'Time (s)', 'Frequency Error (Hz)', []);

% 比较原始跟踪结果和不同卡尔曼滤波器改进后的结果
figure('Color', 'w');
subplot(2,2,1);
plot(t, code_phase_error, 'b', ...
     t, output.KF_Standard_CodePhase - ideal_codePhase, 'r', ...
     t, output.KF_Extended_CodePhase - ideal_codePhase, 'g', ...
     t, output.KF_Unscented_CodePhase - ideal_codePhase, 'm');
title('码相位误差比较');
legend('原始跟踪', '标准卡尔曼滤波', '扩展卡尔曼滤波', 'UKF');
xlabel('时间 (s)');
ylabel('码相位误差');

subplot(2,2,2);
plot(t, carr_phase_error, 'b', ...
     t, output.KF_Standard_CarrPhase - ideal_carrPhase, 'r', ...
     t, output.KF_Extended_CarrPhase - ideal_carrPhase, 'g', ...
     t, output.KF_Unscented_CarrPhase - ideal_carrPhase, 'm');
title('载波相位误差比较');
legend('原始跟踪', '标准卡尔曼滤波', '扩展卡尔曼滤波', 'UKF');
xlabel('时间 (s)');
ylabel('载波相位误差');

subplot(2,2,3);
plot(t, code_freq_error, 'b', ...
     t, output.KF_Standard_CodeFreq - ideal_codeFreq, 'r', ...
     t, output.KF_Extended_CodeFreq - ideal_codeFreq, 'g', ...
     t, output.KF_Unscented_CodeFreq - ideal_codeFreq, 'm');
title('码频率误差比较');
legend('原始跟踪', '标准卡尔曼滤波', '扩展卡尔曼滤波', 'UKF');
xlabel('时间 (s)');
ylabel('码频率误差');

subplot(2,2,4);
plot(t, carr_freq_error, 'b', ...
     t, output.KF_Standard_CarrFreq - ideal_carrFreq, 'r', ...
     t, output.KF_Extended_CarrFreq - ideal_carrFreq, 'g', ...
     t, output.KF_Unscented_CarrFreq - ideal_carrFreq, 'm');
title('载波频率误差比较');
legend('原始跟踪', '标准卡尔曼滤波', '扩展卡尔曼滤波', 'UKF');
xlabel('时间 (s)');
ylabel('载波频率误差');

% 新增：分别展示原始跟踪和三种卡尔曼滤波器的跟踪误差结果
filter_types = {'Original', 'Standard', 'Extended', 'Unscented'};
for i = 1:length(filter_types)
    figure('Color', 'w');
    sgtitle([filter_types{i} ' 跟踪误差结果']);
    
    subplot(2,2,1);
    if i == 1
        plot(t, code_phase_error);
    else
        plot(t, output.(['KF_' filter_types{i} '_CodePhase']) - ideal_codePhase);
    end
    title('码相位误差');
    xlabel('时间 (s)');
    ylabel('误差');
    
    subplot(2,2,2);
    if i == 1
        plot(t, code_freq_error);
    else
        plot(t, output.(['KF_' filter_types{i} '_CodeFreq']) - ideal_codeFreq);
    end
    title('码频率误差');
    xlabel('时间 (s)');
    ylabel('误差 (Hz)');
    
    subplot(2,2,3);
    if i == 1
        plot(t, carr_phase_error);
    else
        plot(t, output.(['KF_' filter_types{i} '_CarrPhase']) - ideal_carrPhase);
    end
    title('载波相位误差');
    xlabel('时间 (s)');
    ylabel('误差');
    
    subplot(2,2,4);
    if i == 1
        plot(t, carr_freq_error);
    else
        plot(t, output.(['KF_' filter_types{i} '_CarrFreq']) - ideal_carrFreq);
    end
    title('载波频率误差');
    xlabel('时间 (s)');
    ylabel('误差 (Hz)');
end

% 新增：显示卡尔曼滤波器估计的误差
figure('Color', 'w');
sgtitle('卡尔曼滤波器估计的误差');

subplot(2,2,1);
plot(t, output.KF_Standard_Error_CodePhase, 'r', ...
     t, output.KF_Extended_Error_CodePhase, 'g', ...
     t, output.KF_Unscented_Error_CodePhase, 'm');
title('码相位误差估计');
legend('标准卡尔曼滤波', '扩展卡尔曼滤波', 'UKF');
xlabel('时间 (s)');
ylabel('误差估计');

subplot(2,2,2);
plot(t, output.KF_Standard_Error_CodeFreq, 'r', ...
     t, output.KF_Extended_Error_CodeFreq, 'g', ...
     t, output.KF_Unscented_Error_CodeFreq, 'm');
title('码频率误差估计');
legend('标准卡尔曼滤波', '扩展卡尔曼滤波', 'UKF');
xlabel('时间 (s)');
ylabel('误差估计 (Hz)');

subplot(2,2,3);
plot(t, output.KF_Standard_Error_CarrPhase, 'r', ...
     t, output.KF_Extended_Error_CarrPhase, 'g', ...
     t, output.KF_Unscented_Error_CarrPhase, 'm');
title('载波相位误差估计');
legend('标准卡尔曼滤波', '扩展卡尔曼滤波', 'UKF');
xlabel('时间 (s)');
ylabel('误差估计');

subplot(2,2,4);
plot(t, output.KF_Standard_Error_CarrFreq, 'r', ...
     t, output.KF_Extended_Error_CarrFreq, 'g', ...
     t, output.KF_Unscented_Error_CarrFreq, 'm');
title('载波频率误差估计');
legend('标准卡尔曼滤波', '扩展卡尔曼滤波', 'UKF');
xlabel('时间 (s)');
ylabel('误差估计 (Hz)');

% 打印 RMSE 结果
fprintf('原始跟踪码相位RMSE: %.4f\n', rmse_codePhase);
fprintf('原始跟踪载波相位RMSE: %.4f\n', rmse_carrPhase);
fprintf('原始跟踪码频率RMSE: %.4f Hz\n', rmse_codeFreq);
fprintf('原始跟踪载波频率RMSE: %.4f Hz\n', rmse_carrFreq);
fprintf('原始跟踪码相位距离RMSE: %.4f m\n', rmse_codePhase_dis);

% 计算并打印卡尔曼滤波器的 RMSE
filter_types = {'Standard', 'Extended', 'Unscented'};
for i = 1:length(filter_types)
    kf_type = filter_types{i};
    
    kf_codePhase_diff = output.(['KF_' kf_type '_CodePhase']) - ideal_codePhase;
    kf_carrPhase_diff = output.(['KF_' kf_type '_CarrPhase']) - ideal_carrPhase;
    kf_codeFreq_diff = output.(['KF_' kf_type '_CodeFreq']) - ideal_codeFreq;
    kf_carrFreq_diff = output.(['KF_' kf_type '_CarrFreq']) - ideal_carrFreq;
    
    rmse_kf_codePhase = sqrt(mean(kf_codePhase_diff(start_index:end).^2));
    rmse_kf_carrPhase = sqrt(mean(kf_carrPhase_diff(start_index:end).^2));
    rmse_kf_codeFreq = sqrt(mean(kf_codeFreq_diff(start_index:end).^2));
    rmse_kf_carrFreq = sqrt(mean(kf_carrFreq_diff(start_index:end).^2));
    
    rmse_kf_codePhase_dis = simSettings.Ts * simSettings.c * rmse_kf_codePhase;
    
    fprintf('%s卡尔曼滤波码相位RMSE: %.4f\n', kf_type, rmse_kf_codePhase);
    fprintf('%s卡尔曼滤波载波相位RMSE: %.4f\n', kf_type, rmse_kf_carrPhase);
    fprintf('%s卡尔曼滤波码频率RMSE: %.4f Hz\n', kf_type, rmse_kf_codeFreq);
    fprintf('%s卡尔曼滤波载波频率RMSE: %.4f Hz\n', kf_type, rmse_kf_carrFreq);
    fprintf('%s卡尔曼滤波码相位距离RMSE: %.4f m\n', kf_type, rmse_kf_codePhase_dis);
    
    % 计算改进百分比
    improvement_codePhase = (rmse_codePhase - rmse_kf_codePhase) / rmse_codePhase * 100;
    improvement_carrPhase = (rmse_carrPhase - rmse_kf_carrPhase) / rmse_carrPhase * 100;
    improvement_codeFreq = (rmse_codeFreq - rmse_kf_codeFreq) / rmse_codeFreq * 100;
    improvement_carrFreq = (rmse_carrFreq - rmse_kf_carrFreq) / rmse_carrFreq * 100;
    
    fprintf('%s卡尔曼滤波码相位改进: %.2f%%\n', kf_type, improvement_codePhase);
    fprintf('%s卡尔曼滤波载波相位改进: %.2f%%\n', kf_type, improvement_carrPhase);
    fprintf('%s卡尔曼滤波码频率改进: %.2f%%\n', kf_type, improvement_codeFreq);
    fprintf('%s卡尔曼滤波载波频率改进: %.2f%%\n', kf_type, improvement_carrFreq);
end

% 绘图函数定义
function plot_results(t, data, title_str, xlabel_str, ylabel_str, legend_str)
    figure('Color', 'w');
    plot(t, data, 'LineWidth', 1.5);
    title(title_str);
    xlabel(xlabel_str);
    ylabel(ylabel_str);
    if ~isempty(legend_str)
        legend(legend_str, 'Location', 'best');
    end
    grid on;
end
