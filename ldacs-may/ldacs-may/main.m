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
%disp(num_samples);%1785
t = 1:num_samples;             % 样本索引

% 创建一个新的图形窗口，背景颜色为白色
figure('Color', 'w');

plot(t,I_P);

% 设置横坐标标签
xlabel('t');

% 设置纵坐标标签
ylabel('I prompt');

% 设置图形标题
title('Bits of the navigation message');



% 若已知采样率，可生成实际时间轴，例如：
% fs = 1000;                   % 假设采样率1kHz
% t = (0:num_samples-1)/fs;    % 时间向量（单位：秒）
% 创建一个新的图形窗口，背景颜色为白色
figure('Color', 'w');

% 绘制 I_P 和 Q_P 的散点图
plot(I_P, Q_P, '.');  % 使用 '.' 来绘制点的散点图

% 设置横坐标标签
xlabel('I prompt');

% 设置纵坐标标签
ylabel('Q prompt');

% 设置图形标题
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
title('corelation');             % 图形标题（按需求保留拼写）
legend('show', 'Location', 'best'); % 显示图例（自动选择最佳位置）

hold off;                        % 释放绘图状态
figure();
c = simSettings.c;
Ts = simSettings.Ts;

diff = codePhase-codePhase0;
diff = diff(1:end);
plot(diff);
figure();
plot(codeFreq-simSettings.fp,'r-', 'LineWidth', 1.5, 'DisplayName', '-10'); 
hold on; 
plot(codeFreq0-simSettings.fp,'b-', 'LineWidth', 1.5, 'DisplayName', 'inf');
legend('show', 'Location', 'best'); % 显示图例（自动选择最佳位置）
rmse_codePhase = sqrt(sum(diff.^2)/size(diff,2))
rmse_codePhase_dis = Ts*c*sqrt(sum(diff.^2)/size(diff,2))

% 比较原始跟踪结果和卡尔曼滤波器改进后的结果
figure;
subplot(2,1,1);
plot(t0, output.OutCodePhase, 'b', t0, output.KF_CodePhase, 'r');
title('码相位比较');
legend('原始跟踪', '卡尔曼滤波');
xlabel('时间 (s)');
ylabel('码相位');

subplot(2,1,2);
plot(t0, output.OutCarrPhase, 'b', t0, output.KF_CarrPhase, 'r');
title('载波相位比较');
legend('原始跟踪', '卡尔曼滤波');
xlabel('时间 (s)');
ylabel('载波相位');

% 计算卡尔曼滤波器改进后的RMSE
kalman_diff = output.KF_CodePhase - codePhase0;
kalman_diff = kalman_diff(200:end);
rmse_kalman_codePhase = sqrt(sum(kalman_diff.^2)/size(kalman_diff,2))
rmse_kalman_codePhase_dis = Ts*c*sqrt(sum(kalman_diff.^2)/size(kalman_diff,2))

fprintf('原始跟踪码相位RMSE: %.4f\n', rmse_codePhase);
fprintf('卡尔曼滤波码相位RMSE: %.4f\n', rmse_kalman_codePhase);
fprintf('原始跟踪码相位距离RMSE: %.4f m\n', rmse_codePhase_dis);
fprintf('卡尔曼滤波码相位距离RMSE: %.4f m\n', rmse_kalman_codePhase_dis);

% 绘制I路bits of navigation message
figure;
plot(t0(1:800), output.I_P(1:800), 'b', t0(1:800), output.KF_CodePhase(1:800), 'r');
title('I路bits of navigation message (前800个周期)');
legend('原始跟踪', '卡尔曼滤波');
xlabel('时间 (s)');
ylabel('幅度/相位');
