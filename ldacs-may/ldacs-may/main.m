% 5月对跟踪代码的修改
clear all; close all;

% 读取参数
SNR = -20;
simSettings = init(SNR); % 设置参数
flag = 1; % 要不要重新生成

[grid, rawSignal] = generateCrossOFDM(simSettings); % 资源网格，行是时间，列是频率

yr0 = real(rawSignal); % 实部信号
yi0 = imag(rawSignal); % 虚部信号

yr0=interpo(yr0,simSettings.famp); % 实部时域采样结果
yi0=interpo(yi0,simSettings.famp); % 虚部时域采样结果

yr0 = [yr0(:,end-simSettings.extra*simSettings.famp+1:end) yr0 yr0(:,1:1+simSettings.extra*simSettings.famp)];
yi0 = [yi0(:,end-simSettings.extra*simSettings.famp+1:end) yi0 yi0(:,1:1+simSettings.extra*simSettings.famp)];

if flag==1
    for bandID=1:simSettings.numBand
        yr = yr0(bandID,:);
        yi = yi0(bandID,:);
        record_E(:,:,bandID) = generateSignal(simSettings, yr, yi);
    end
    filenm1 = ['output/signal_SNR'  num2str(SNR) '.mat' ];
    save(filenm1,'record_E');
else
    load(['output/signal_SNR'  num2str(SNR) '.mat' ]);
end

output = track(simSettings, yr0, yi0);
savenm=['output/tracked_SNR'  num2str(SNR) '.mat'];
save(savenm,'output');

load("output/tracked_SNRinf.mat");
codePhase0 = output.OutCodePhase;
carrPhase0 = output.OutCarrPhase;
codeFreq0 = output.OutCodeFreq;

loadnm = ['output/tracked_SNR'  num2str(SNR) '.mat' ];
load(loadnm);
codePhase = output.OutCodePhase;
carrPhase = output.OutCarrPhase;
codeFreq = output.OutCodeFreq;

I_E = output.I_E;
Q_E = output.Q_E;
I_P = output.I_P;
Q_P = output.Q_P;
I_L = output.I_L;
Q_L = output.Q_L;

% 计算各通道的幅度：sqrt(I² + Q²)
Amp_E = sqrt(I_E.^2 + Q_E.^2); % 早（Early）通道幅度
Amp_P = sqrt(I_P.^2 + Q_P.^2); % 准时（Prompt）通道幅度
Amp_L = sqrt(I_L.^2 + Q_L.^2); % 晚（Late）通道幅度

% 创建时间轴（假设数据为时间序列，若无时间信息可用样本索引代替）
num_samples = length(Amp_E);% 获取数据长度
t = 1:num_samples;             % 样本索引

% 创建一个新的图形窗口，背景颜色为白色
figure('Color', 'w');
plot(t,I_P);
xlabel('t');
ylabel('I prompt');
title('Bits of the navigation message');

% 创建一个新的图形窗口，背景颜色为白色
figure('Color', 'w');
% 绘制 I_P 和 Q_P 的散点图
plot(I_P, Q_P, '.');  % 使用 '.' 来绘制点的散点图
xlabel('I prompt');
ylabel('Q prompt');
title('Discrete-Time Scatter Plot');
axis equal;

% 创建新图形窗口
figure('Color', 'w');            % 白色背景窗口
hold on;                         % 保持绘图状态
% 绘制三条通道幅度曲线
plot(t , Amp_E, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Early');
plot(t , Amp_P, 'g-', 'LineWidth', 1.5, 'DisplayName', 'Prompt');
plot(t , Amp_L, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Late');
% 添加图形标注
grid on;                         % 显示网格
xlabel('Sample Index');          % 横轴标签
ylabel('Amplitude');             % 纵轴标签
title('Correlation');             % 图形标题
legend('show', 'Location', 'best'); % 显示图例（自动选择最佳位置）
hold off;                        % 释放绘图状态

figure();
c = simSettings.c;
Ts = simSettings.Ts;
diff = codePhase-codePhase0;
diff = diff(1:end);
plot(diff);
title('Code Phase Difference');
xlabel('Sample Index');
ylabel('Phase Difference');

figure();
plot(codeFreq-simSettings.fp,'r-', 'LineWidth', 1.5, 'DisplayName', '-10'); 
hold on; 
plot(codeFreq0-simSettings.fp,'b-', 'LineWidth', 1.5, 'DisplayName', 'inf');
legend('show', 'Location', 'best');
title('Code Frequency Deviation');
xlabel('Sample Index');
ylabel('Frequency Deviation (Hz)');

rmse_codePhase = sqrt(sum(diff.^2)/size(diff,2));
rmse_codePhase_dis = Ts*c*sqrt(sum(diff.^2)/size(diff,2));

% 比较原始跟踪结果和不同卡尔曼滤波器改进后的结果
figure;
subplot(2,1,1);
plot(t, output.OutCodePhase, 'b', t, output.KF_Standard_CodePhase, 'r', t, output.KF_Extended_CodePhase, 'g', t, output.KF_Unscented_CodePhase, 'm');
title('码相位比较');
legend('原始跟踪', '标准卡尔曼滤波', '扩展卡尔曼滤波', 'UKF');
xlabel('时间 (s)');
ylabel('码相位');

subplot(2,1,2);
plot(t, output.OutCarrPhase, 'b', t, output.KF_Standard_CarrPhase, 'r', t, output.KF_Extended_CarrPhase, 'g', t, output.KF_Unscented_CarrPhase, 'm');
title('载波相位比较');
legend('原始跟踪', '标准卡尔曼滤波', '扩展卡尔曼滤波', 'UKF');
xlabel('时间 (s)');
ylabel('载波相位');

% 计算不同卡尔曼滤波器改进后的RMSE
kalman_standard_diff = output.KF_Standard_CodePhase - codePhase0;
kalman_extended_diff = output.KF_Extended_CodePhase - codePhase0;
kalman_unscented_diff = output.KF_Unscented_CodePhase - codePhase0;

kalman_standard_diff = kalman_standard_diff(200:end);
kalman_extended_diff = kalman_extended_diff(200:end);
kalman_unscented_diff = kalman_unscented_diff(200:end);

rmse_kalman_standard = sqrt(sum(kalman_standard_diff.^2)/size(kalman_standard_diff,2));
rmse_kalman_extended = sqrt(sum(kalman_extended_diff.^2)/size(kalman_extended_diff,2));
rmse_kalman_unscented = sqrt(sum(kalman_unscented_diff.^2)/size(kalman_unscented_diff,2));

rmse_kalman_standard_dis = Ts*c*rmse_kalman_standard;
rmse_kalman_extended_dis = Ts*c*rmse_kalman_extended;
rmse_kalman_unscented_dis = Ts*c*rmse_kalman_unscented;

fprintf('原始跟踪码相位RMSE: %.4f\n', rmse_codePhase);
fprintf('标准卡尔曼滤波码相位RMSE: %.4f\n', rmse_kalman_standard);
fprintf('扩展卡尔曼滤波码相位RMSE: %.4f\n', rmse_kalman_extended);
fprintf('UKF码相位RMSE: %.4f\n', rmse_kalman_unscented);
fprintf('原始跟踪码相位距离RMSE: %.4f m\n', rmse_codePhase_dis);
fprintf('标准卡尔曼滤波码相位距离RMSE: %.4f m\n', rmse_kalman_standard_dis);
fprintf('扩展卡尔曼滤波码相位距离RMSE: %.4f m\n', rmse_kalman_extended_dis);
fprintf('UKF码相位距离RMSE: %.4f m\n', rmse_kalman_unscented_dis);

% 绘制I路bits of navigation message
figure;
plot(t, output.I_P, 'b', t, output.KF_Standard_CodePhase, 'r', t, output.KF_Extended_CodePhase, 'g', t, output.KF_Unscented_CodePhase, 'm');
title('I路bits of navigation message');
legend('原始跟踪', '标准卡尔曼滤波', '扩展卡尔曼滤波', 'UKF');
xlabel('时间 (s)');
ylabel('幅度/相位');
xlim([0, max(t)]); % 显示全部范围

% 绘制码频率偏差
figure;
plot(t, codeFreq - simSettings.fp, 'b-', 'LineWidth', 1.5, 'DisplayName', '原始跟踪');
hold on;
plot(t, output.KF_Standard_CodeFreq - simSettings.fp, 'r-', 'LineWidth', 1.5, 'DisplayName', '标准卡尔曼滤波');
plot(t, output.KF_Extended_CodeFreq - simSettings.fp, 'g-', 'LineWidth', 1.5, 'DisplayName', '扩展卡尔曼滤波');
plot(t, output.KF_Unscented_CodeFreq - simSettings.fp, 'm-', 'LineWidth', 1.5, 'DisplayName', 'UKF');
title('码频率偏差');
legend('show', 'Location', 'best');
xlabel('时间 (s)');
ylabel('频率偏差 (Hz)');
xlim([0, max(t)]); % 显示全部范围

% 绘制码相位误差
figure;
plot(t, codePhase - codePhase0, 'b-', 'LineWidth', 1.5, 'DisplayName', '原始跟踪');
hold on;
plot(t, output.KF_Standard_CodePhase - codePhase0, 'r-', 'LineWidth', 1.5, 'DisplayName', '标准卡尔曼滤波');
plot(t, output.KF_Extended_CodePhase - codePhase0, 'g-', 'LineWidth', 1.5, 'DisplayName', '扩展卡尔曼滤波');
plot(t, output.KF_Unscented_CodePhase - codePhase0, 'm-', 'LineWidth', 1.5, 'DisplayName', 'UKF');
title('码相位误差');
legend('show', 'Location', 'best');
xlabel('时间 (s)');
ylabel('相位误差');
xlim([0, max(t)]); % 显示全部范围
