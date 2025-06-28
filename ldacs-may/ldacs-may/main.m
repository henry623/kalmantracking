% 5月对跟踪代码的修改
clear all; close all;

% 读取参数
SNR = inf;
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

c = simSettings.c;
Ts = simSettings.Ts;

diff = codePhase-codePhase0;
diff = diff(200:end);
plot(diff);
figure();
plot(codeFreq-simSettings.fp); hold on; plot(codeFreq0-simSettings.fp);
rmse_codePhase = sqrt(sum(diff.^2)/size(diff,2))
rmse_codePhase_dis = Ts*c*sqrt(sum(diff.^2)/size(diff,2))