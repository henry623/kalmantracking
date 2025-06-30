function [output]=track2(simSettings, yr, yi)
% 跟踪环路 - 修改版本，使用误差状态的卡尔曼滤波器

% 添加 kalman_filters2 函数的路径
addpath('ldacs-may/ldacs-may');

% 仿真参数
fp = simSettings.fp;
fs = simSettings.fs;
fi = simSettings.fi;
numBand = simSettings.numBand;
SNR = simSettings.SNR;
dt = simSettings.dt;
t_sim = simSettings.t_total;
t0 = 0:dt:t_sim;
lenOFDM = simSettings.NFFT*simSettings.nSymbol;

% 初始化卡尔曼滤波器
% 状态变量：[载波频率误差; 载波相位误差; 码频率误差; 码相位误差]
kf.x = [0; 0; 0; 0];  % 初始状态都设为0（无误差）
kf.P = diag([1e-4, 1e-2, 1e-4, 1e-2]);  % 初始协方差矩阵
kf.Q = diag([1e-8, 1e-6, 1e-8, 1e-6]);  % 过程噪声协方差
kf.R = diag([1e-2, 1e-2]);  % 测量噪声协方差
kf.H = [0 1 0 0; 0 0 0 1];  % 测量矩阵 - 只测量载波相位误差和码相位误差

% 初始化不同类型的卡尔曼滤波器
kf_standard = kf;
kf_extended = kf;
kf_unscented = kf;

% 读取生成信号文件
loadnm = ['output/signal_SNR'  num2str(SNR) '.mat' ];
record_E = load(loadnm);
record_E = record_E.record_E;

% 输出数据准备
output.OutCarrFreq = zeros(1,size(t0,2));
output.OutCodeFreq = zeros(1,size(t0,2));
output.OutCarrPhase = zeros(1,size(t0,2));
output.OutCodePhase = zeros(1,size(t0,2));
output.OutCarrNCO = zeros(1,size(t0,2));
output.OutCodeNCO = zeros(1,size(t0,2));
output.OutBlksize = zeros(1,size(t0,2));
output.OutRemCode = zeros(1,size(t0,2));

%% 环路数据准备
earlyLateSpc = 0.5;
dllNoiseBandwidth = 1;
dllDampingRatio = 0.707;
PDIcode = lenOFDM / fs;
[tau1code, tau2code] = calLoopCoef(dllNoiseBandwidth, dllDampingRatio, 1.0);
wn = dllNoiseBandwidth/0.7845;
omegaOf = dllNoiseBandwidth/0.53;

pllNoiseBandwidth = 4;
pllDampingRatio = 0.707;
PDIcarr = lenOFDM / fs;
[tau1carr, tau2carr] = calLoopCoef(pllNoiseBandwidth, pllDampingRatio, 0.25);
carrPhase = 0.0;
codePhase = 0;
carrFreq = simSettings.fi;
codeFreq = fp;

remCarrPhase = 0;
remCodePhase = 0;

oldCodeNco = 0.0;
oldCodeError = 0.0;
carrErrorFreq = 0.0;

oldPllX = 0.0;
oldPllY = 0.0;
old_IP_SUM = 0.0;
old_QP_SUM = 0.0;
a2 = 1.414;
oldfllx = 0.0;
oldflly = 0.0;
oldCarrFreq = 0;

oldCarrNco = 0.0;
oldCarrError = 0.0;

readIndex = 0;

trkCorrCnt = 1;

loopCntToProcess = floor(t_sim / simSettings.To); 

% 理想值（用于计算误差）
ideal_codePhase = 0;
ideal_carrPhase = 0;
ideal_codeFreq = fp;
ideal_carrFreq = fi;

%% 跟踪
for loopCnt = 1:size(t0,2)
    
    codePhaseStep = codeFreq / fs;
    
    blksize = ceil((lenOFDM-remCodePhase) / codePhaseStep);
    
    readIndex = readIndex(end)+1:readIndex(end)+1+blksize-1;

    I_E = 0;
    Q_E = 0;
    I_P = 0;
    Q_P = 0;
    I_L = 0;
    Q_L = 0;

    for bandID = 1:numBand

        yl_r = yr(bandID,:);
        yl_i = yi(bandID,:);

        rawSignal = record_E(:,readIndex,bandID);
        rawSignalI = rawSignal(1,:);
        rawSignalQ = rawSignal(2,:);
        
        % E
        tcode = (remCodePhase-earlyLateSpc) : ...
                codePhaseStep : ...
                ((blksize-1)*codePhaseStep+remCodePhase-earlyLateSpc);
        tcode2 = round(tcode*simSettings.famp+1)+simSettings.extra*simSettings.famp;
        earlyCode_r = yl_r(tcode2);
        earlyCode_i = yl_i(tcode2);
        earlyCode = earlyCode_r+earlyCode_i;
        
        % L
        tcode = (remCodePhase+earlyLateSpc) : ...
                codePhaseStep : ...
                ((blksize-1)*codePhaseStep+remCodePhase+earlyLateSpc);
        tcode2 = round(tcode*simSettings.famp+1)+simSettings.extra*simSettings.famp;
        lateCode_r = yl_r(tcode2);
        lateCode_i = yl_i(tcode2);
        lateCode = lateCode_r+lateCode_i;
        
        % P     
        tcode = remCodePhase : ...
                codePhaseStep : ...
                ((blksize-1)*codePhaseStep+remCodePhase);
        tcode2 = round(tcode*simSettings.famp+1)+simSettings.extra*simSettings.famp;
        promptCode_r = yl_r(tcode2);        
        promptCode_i = yl_i(tcode2);
        promptCode = promptCode_r+promptCode_i;
        
        time = (0:blksize) ./ fs;
        
        % Get the argument to sin/cos functions
        trigarg = ((carrFreq * 2.0 * pi) .* time) + 2.0 * pi * remCarrPhase;
        
        % Finally compute the signal to mix the collected data to bandband
        cosCarr = cos(trigarg(1:blksize));
        sinCarr = sin(trigarg(1:blksize));
    
        iBasebandSignal = rawSignalI .* cosCarr + rawSignalQ .* sinCarr;
        qBasebandSignal = rawSignalQ .* cosCarr - rawSignalI .* sinCarr;

        % Now get early, prompt, and late values for each
        I_E = I_E+sum(earlyCode .* iBasebandSignal);
        Q_E = Q_E+sum(earlyCode .* qBasebandSignal);
        I_P = I_P+sum(promptCode .* iBasebandSignal);
        Q_P = Q_P+sum(promptCode .* qBasebandSignal);
        I_L = I_L+sum(lateCode .* iBasebandSignal);
        Q_L = Q_L+sum(lateCode .* qBasebandSignal);
    end
    remCarrPhase = rem(trigarg(blksize+1), (2 * pi))/(2.0 * pi);
    carrPhase = carrPhase + ((carrFreq) * (blksize / fs));

    remCodePhase = tcode(blksize) + codePhaseStep - lenOFDM;
    codePhase = codePhase + (blksize * codePhaseStep);

    % Implement carrier loop discriminator (phase detector)
    carrError = atan(Q_P / I_P) / (2.0 * pi);

    if loopCnt == 1
        pDot = 1;
    else
        pDot = old_IP_SUM * I_P + old_QP_SUM * Q_P;
    end

    pCross = old_IP_SUM * Q_P - old_QP_SUM * I_P;
    old_IP_SUM = I_P;
    old_QP_SUM = Q_P;
    carrErrorFreq = atan(pCross/pDot) / PDIcarr/2/pi;

    if(loopCnt <= 0)
        carrErrorFreq = carrErrorFreq * 1;
        fllx = carrErrorFreq*(omegaOf*omegaOf)*PDIcarr + oldfllx;
        flly = 0.5*(fllx + oldfllx) + a2*omegaOf*carrErrorFreq;
        carrFreq = oldCarrFreq+(flly-oldflly)*PDIcarr;
        oldfllx = fllx;
        oldflly = flly;
        oldCarrFreq = carrFreq;            
    else              
        carrError = carrError * 1;
        carrErrorFreq = carrErrorFreq * 1;                
        pllX = oldPllX + carrErrorFreq*(omegaOf*omegaOf)*PDIcarr + carrError*(wn*wn*wn)*PDIcarr;
        pllY = oldPllY + carrErrorFreq*1.414*omegaOf*PDIcarr + 0.5*(pllX + oldPllX)*PDIcarr+1.1*carrError* wn*wn*PDIcarr;
        carrNco = carrError* 2.4*wn+0.5*(pllY+oldPllY);
        oldPllX = pllX;
        oldPllY = pllY; 
        carrFreq = fi + carrNco;
    end

    %% Find DLL error and update code NCO -------------------------------------
    codeError = (sqrt(I_E * I_E + Q_E * Q_E) - sqrt(I_L * I_L + Q_L * Q_L)) / ...
        (sqrt(I_E * I_E + Q_E * Q_E) + sqrt(I_L * I_L + Q_L * Q_L));

    % Implement code loop filter and generate NCO command
    codeNco = oldCodeNco + (tau2code/tau1code) * ...
        (codeError - oldCodeError) + codeError * (PDIcode/tau1code);
    oldCodeNco = codeNco;
    oldCodeError = codeError;

    % Modify code freq based on NCO command
    codeFreq = fp - codeNco;

    % 计算当前测量值与理想值之间的误差
    carr_phase_error = carrPhase - ideal_carrPhase;
    code_phase_error = codePhase - ideal_codePhase;
    
    % 使用不同类型的卡尔曼滤波器
    z = [carr_phase_error; code_phase_error];  % 测量值：载波相位误差和码相位误差
    ideal_values = [ideal_carrFreq; ideal_carrPhase; ideal_codeFreq; ideal_codePhase];  % 理想值
    
    % 使用标准卡尔曼滤波器
    kf_standard = kalman_filters2('standard', kf_standard, z, dt, ideal_values);
    
    % 使用扩展卡尔曼滤波器
    kf_extended = kalman_filters2('extended', kf_extended, z, dt, ideal_values);
    
    % 使用无迹卡尔曼滤波器
    kf_unscented = kalman_filters2('unscented', kf_unscented, z, dt, ideal_values);
    
    % 保存不同卡尔曼滤波器的结果
    output.KF_Standard_CarrFreq(loopCnt) = kf_standard.corrected_values(1);
    output.KF_Standard_CarrPhase(loopCnt) = kf_standard.corrected_values(2);
    output.KF_Standard_CodeFreq(loopCnt) = kf_standard.corrected_values(3);
    output.KF_Standard_CodePhase(loopCnt) = kf_standard.corrected_values(4);
    
    output.KF_Extended_CarrFreq(loopCnt) = kf_extended.corrected_values(1);
    output.KF_Extended_CarrPhase(loopCnt) = kf_extended.corrected_values(2);
    output.KF_Extended_CodeFreq(loopCnt) = kf_extended.corrected_values(3);
    output.KF_Extended_CodePhase(loopCnt) = kf_extended.corrected_values(4);
    
    output.KF_Unscented_CarrFreq(loopCnt) = kf_unscented.corrected_values(1);
    output.KF_Unscented_CarrPhase(loopCnt) = kf_unscented.corrected_values(2);
    output.KF_Unscented_CodeFreq(loopCnt) = kf_unscented.corrected_values(3);
    output.KF_Unscented_CodePhase(loopCnt) = kf_unscented.corrected_values(4);

    % 保存原始跟踪结果
    output.OutCarrFreq(loopCnt) = carrFreq;
    output.OutCodeFreq(loopCnt) = codeFreq;
    output.OutCarrPhase(loopCnt) = carrPhase;
    output.OutCodePhase(loopCnt) = codePhase;
    output.OutCarrNCO(loopCnt) = carrError;
    output.OutCodeNCO(loopCnt) = codeError;
    output.OutBlksize(loopCnt) = blksize;
    output.OutRemCode(loopCnt) = remCodePhase;
    output.I_E(loopCnt) = I_E;
    output.Q_E(loopCnt) = Q_E;
    output.I_P(loopCnt) = I_P;
    output.Q_P(loopCnt) = Q_P;
    output.I_L(loopCnt) = I_L;
    output.Q_L(loopCnt) = Q_L;
    
    % 更新理想值（随时间累积）
    ideal_codePhase = ideal_codePhase + ideal_codeFreq * dt;
    ideal_carrPhase = ideal_carrPhase + ideal_carrFreq * dt;
end

% 添加卡尔曼滤波器的最终状态到输出
output.KF_Standard_FinalState = kf_standard.x;
output.KF_Standard_FinalCovariance = kf_standard.P;
output.KF_Extended_FinalState = kf_extended.x;
output.KF_Extended_FinalCovariance = kf_extended.P;
output.KF_Unscented_FinalState = kf_unscented.x;
output.KF_Unscented_FinalCovariance = kf_unscented.P;

end
