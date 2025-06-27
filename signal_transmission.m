function received_signal = signal_transmission(signal, channel_response, snr_db)
% 信号传输模块：将信号通过信道并添加噪声
%
% 输入参数：
%   signal: 发送的信号
%   channel_response: 信道响应
%   snr_db: 信噪比 (dB)
%
% 输出参数：
%   received_signal: 接收到的信号

% 将信号通过信道
transmitted_signal = signal .* channel_response;

% 计算信号功率
signal_power = mean(abs(transmitted_signal).^2);

% 计算噪声功率
noise_power = signal_power / (10^(snr_db/10));

% 生成复高斯白噪声
noise = sqrt(noise_power/2) * (randn(size(signal)) + 1i*randn(size(signal)));

% 添加噪声到信号
received_signal = transmitted_signal + noise;

end