function result_visualization(t, original_signal, received_signal, tracked_signal, estimated_phase, estimated_freq, true_phase, true_freq)
% 结果可视化模块：展示信号跟踪结果
%
% 输入参数：
%   t: 时间向量
%   original_signal: 原始发送信号
%   received_signal: 接收到的信号
%   tracked_signal: 跟踪后的信号
%   estimated_phase: 估计的相位
%   estimated_freq: 估计的频率
%   true_phase: 真实相位
%   true_freq: 真实频率

% 创建新图形窗口
figure('Name', '信号跟踪结果', 'NumberTitle', 'off', 'Position', [100, 100, 1200, 800]);

% 绘制信号波形
subplot(3, 2, 1);
plot(t, real(original_signal), 'b', t, real(received_signal), 'r');
title('原始信号与接收信号对比');
xlabel('时间 (s)');
ylabel('幅度');
legend('原始信号', '接收信号');

subplot(3, 2, 2);
plot(t, real(received_signal), 'r', t, real(tracked_signal), 'g');
title('接收信号与跟踪信号对比');
xlabel('时间 (s)');
ylabel('幅度');
legend('接收信号', '跟踪信号');

% 绘制相位估计结果
subplot(3, 2, 3);
plot(t, unwrap(true_phase), 'b', t, unwrap(estimated_phase), 'r');
title('相位估计结果');
xlabel('时间 (s)');
ylabel('相位 (rad)');
legend('真实相位', '估计相位');

% 绘制频率估计结果
subplot(3, 2, 4);
plot(t, true_freq, 'b', t, estimated_freq, 'r');
title('频率估计结果');
xlabel('时间 (s)');
ylabel('频率 (Hz)');
legend('真实频率', '估计频率');

% 绘制相位误差
subplot(3, 2, 5);
phase_error = wrapToPi(estimated_phase - true_phase);
plot(t, phase_error);
title('相位估计误差');
xlabel('时间 (s)');
ylabel('相位误差 (rad)');

% 绘制频率误差
subplot(3, 2, 6);
freq_error = estimated_freq - true_freq;
plot(t, freq_error);
title('频率估计误差');
xlabel('时间 (s)');
ylabel('频率误差 (Hz)');

% 调整子图间距
tight_subplot(3, 2, [0.1 0.05], [0.1 0.05], [0.05 0.02]);

% 计算和显示性能指标
phase_rmse = sqrt(mean(phase_error.^2));
freq_rmse = sqrt(mean(freq_error.^2));

fprintf('相位估计RMSE: %.4f rad\n', phase_rmse);
fprintf('频率估计RMSE: %.4f Hz\n', freq_rmse);

end

function h = tight_subplot(m, n, margin_h, margin_w, space_h, space_w)
% 创建紧凑的子图布局
%
% 输入参数：
%   m, n: 子图的行数和列数
%   margin_h, margin_w: 图形边缘的高度和宽度 [下 上] [左 右]
%   space_h, space_w: 子图之间的垂直和水平间距

if nargin < 3; margin_h = 0.05; end
if nargin < 4; margin_w = 0.05; end
if nargin < 5; space_h = 0.05; end
if nargin < 6; space_w = 0.05; end

% 计算每个子图的位置
ax_w = (1 - sum(margin_w) - (n-1)*space_w) / n;
ax_h = (1 - sum(margin_h) - (m-1)*space_h) / m;

py = 1 - margin_h(2) - ax_h;
h = zeros(m*n, 1);

for i = 1:m
    px = margin_w(1);
    for j = 1:n
        pos = [px py ax_w ax_h];
        h((i-1)*n+j) = subplot('Position', pos);
        px = px + ax_w + space_w;
    end
    py = py - ax_h - space_h;
end

end