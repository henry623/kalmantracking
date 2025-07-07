function simSettings = init(SNR)

% 常规
simSettings.c = 3e8; % 光速
simSettings.fi = 0; % 中频频率
simSettings.fc = 1.5e9; % 载波频率

% OFDM调制参数
simSettings.Nu = 50; % 每个频段可用的子载波数量
simSettings.NFFT = 64; % FFT数量
simSettings.fsub = 9.765625e3; % 子载波间隔
simSettings.To = 1/simSettings.fsub; % 一个OFDM符号的时间

% 码参数
load('weil10230_signed.mat');
simSettings.code = weil10230_signed(1,:);
% simSettings.code = ones(size(weil10230_signed(1,:)));
% simSettings.code(simSettings.code==-1) = 0;
simSettings.Lc = length(simSettings.code); % 码的长度

% 仿真参数
simSettings.SNR = SNR; % 信噪比

simSettings.t_total = 20; % 仿真的时间长度

simSettings.p0 = [2500,2500,1000]; % 初始位置
simSettings.v0 = [100,0,0]*1000/3600; % 速度
simSettings.BS = simSettings.p0; % 基站位置
% simSettings.BS = [0,0,0]; % 基站位置

simSettings.famp = 3; % 采样率是码率的倍数（没有对非整数做处理）
simSettings.numBand = 4; % 跨的频段数量
simSettings.extra = 5; % 前后额外取的量

% 耦合参数
simSettings.nSymbol = ceil(simSettings.Lc/(simSettings.Nu*simSettings.numBand)); % 一个完整的码周期所需要的符号数
simSettings.Tp = simSettings.To*simSettings.nSymbol; % 在OFDM里传完一整个码周期所需的时间
simSettings.Ts = simSettings.To/simSettings.NFFT/simSettings.numBand; % 一个数据的时间
simSettings.fp = 1/simSettings.Ts; % 码频率
% simSettings.fs = simSettings.fp*simSettings.famp; % 采样频率
simSettings.fs = simSettings.fp*2; % 采样频率
simSettings.dt = simSettings.nSymbol*simSettings.To*2; % 仿真的时间间隔
% simSettings.dt = 0.1;