 function [output]=track(simSettings, yr, yi)
% 跟踪环路

% 添加 kalman_filters 函数的路径
addpath('ldacs-may/ldacs-may');

% 仿真参数
fp = simSettings.fp;
fs = simSettings.fs;
fi = simSettings.fi;
numBand = simSettings.numBand;
SNR = simSettings.SNR;
dt = simSettings.dt;
t_sim = simSettings.t_total-1;
t0 = 0:dt:t_sim;
lenOFDM = simSettings.NFFT*simSettings.nSymbol;
% lenOFDM = 10230;

% 初始化卡尔曼滤波器
kf.x = [0; fp; 0; fi];  % 初始状态 [码相位; 码频率; 载波相位; 载波频率]
kf.P = diag([1e-2, 1e-4, 1e-2, 1e-4]);  % 初始协方差矩阵
kf.Q = diag([1e-6, 1e-8, 1e-6, 1e-8]);  % 过程噪声协方差
kf.R = diag([1e-2, 1e-2]);  % 测量噪声协方差
kf.H = [1 0 0 0; 0 0 1 0];  % 测量矩阵

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
dllNoiseBandwidth       = 1;
dllDampingRatio         = 0.707;
PDIcode = lenOFDM / fs;
[tau1code, tau2code] = calLoopCoef(dllNoiseBandwidth, ...
                                    dllDampingRatio, ...
                                    1.0);
wn = dllNoiseBandwidth/0.7845;
omegaOf = dllNoiseBandwidth/0.53;

pllNoiseBandwidth       = 4;
pllDampingRatio         = 0.707;
PDIcarr = lenOFDM / fs;
[tau1carr, tau2carr] = calLoopCoef(pllNoiseBandwidth, ...
                                    pllDampingRatio, ...
                                    0.25);
carrPhase = 0.0;
codePhase = 0;
carrFreq = simSettings.fi;
codeFreq = fp;

remCarrPhase=0;
remCodePhase=0;

oldCodeNco   = 0.0;
oldCodeError = 0.0;
carrErrorFreq = 0.0;

oldPllX =0.0;
oldPllY =0.0;
old_IP_SUM = 0.0;
old_QP_SUM = 0.0;
a2=1.414;
oldfllx = 0.0;
oldflly = 0.0;
oldCarrFreq = 0;

oldCarrNco   = 0.0;
oldCarrError = 0.0;

readIndex=0;

trkCorrCnt=1;

loopCntToProcess = floor(t_sim / simSettings.To); 

%% 跟踪
for loopCnt=1:size(t0,2)
    
    codePhaseStep = codeFreq / fs;
    
    blksize = ceil((lenOFDM-remCodePhase) / codePhaseStep);
    
    readIndex = readIndex(end)+1:readIndex(end)+1+blksize-1;

    I_E = 0;
    Q_E = 0;
    I_P = 0;
    Q_P = 0;
    I_L = 0;
    Q_L = 0;

    for bandID=1:numBand

        yl_r = yr(bandID,:);
        yl_i = yi(bandID,:);

        rawSignal=record_E(:,readIndex,bandID);
        rawSignalI = rawSignal(1,:);
        rawSignalQ = rawSignal(2,:);
        
        % E
        tcode       = (remCodePhase-earlyLateSpc) : ...
                      codePhaseStep : ...
                      ((blksize-1)*codePhaseStep+remCodePhase-earlyLateSpc);
        tcode2 = round(tcode*simSettings.famp+1)+simSettings.extra*simSettings.famp;
        earlyCode_r   = yl_r(tcode2);
        earlyCode_i   = yl_i(tcode2);
        earlyCode = earlyCode_r+earlyCode_i;
        
        % L
        tcode       = (remCodePhase+earlyLateSpc) : ...
                      codePhaseStep : ...
                      ((blksize-1)*codePhaseStep+remCodePhase+earlyLateSpc);
        tcode2 = round(tcode*simSettings.famp+1)+simSettings.extra*simSettings.famp;
        lateCode_r   = yl_r(tcode2);
        lateCode_i   = yl_i(tcode2);
        lateCode = lateCode_r+lateCode_i;
        
        % P     
        tcode       = remCodePhase : ...
                      codePhaseStep : ...
                      ((blksize-1)*codePhaseStep+remCodePhase);
        tcode2 = round(tcode*simSettings.famp+1)+simSettings.extra*simSettings.famp;
        promptCode_r  = yl_r(tcode2);        
        promptCode_i  = yl_i(tcode2);
        promptCode = promptCode_r+promptCode_i;
        
        time    = (0:blksize) ./ fs;
        
        % Get the argument to sin/cos functions
        trigarg = ((carrFreq * 2.0 * pi) .* time) + 2.0 * pi * remCarrPhase;
        
        % Finally compute the signal to mix the collected data to bandband
        cosCarr = cos(trigarg(1:blksize));
        sinCarr = sin(trigarg(1:blksize));
    
        iBasebandSignal = rawSignalI .* cosCarr + rawSignalQ .* sinCarr;
        qBasebandSignal = rawSignalQ .* cosCarr - rawSignalI .* sinCarr;

        % Now get early, prompt, and late values for each
        I_E = I_E+sum(earlyCode  .* iBasebandSignal);
        Q_E = Q_E+sum(earlyCode  .* qBasebandSignal);
        I_P = I_P+sum(promptCode .* iBasebandSignal);
        Q_P = Q_P+sum(promptCode .* qBasebandSignal);
        I_L = I_L+sum(lateCode   .* iBasebandSignal);
        Q_L = Q_L+sum(lateCode   .* qBasebandSignal);
    end
    remCarrPhase = rem(trigarg(blksize+1), (2 * pi))/(2.0 * pi);
    carrPhase = carrPhase + ((carrFreq) * (blksize / fs));

    remCodePhase = tcode(blksize) + codePhaseStep - lenOFDM;
    codePhase = codePhase + (blksize * codePhaseStep);

    % Implement carrier loop discriminator (phase detector)
    carrError = atan(Q_P / I_P) / (2.0 * pi);

    if loopCnt ==1
        pDot=1;
    else
        pDot   = old_IP_SUM * I_P + old_QP_SUM * Q_P;
    end

    pCross = old_IP_SUM * Q_P - old_QP_SUM * I_P;
    old_IP_SUM = I_P;
    old_QP_SUM = Q_P;
    carrErrorFreq = atan(pCross/pDot)  / PDIcarr/2/pi;%
    % carrErrorFreq=pCross*sign(pDot)/(IP_SUM^2+QP_SUM^2)/PDIcarr;

    if(loopCnt<=0)
        carrErrorFreq  = carrErrorFreq * 1;
                        fllx = carrErrorFreq*(omegaOf*omegaOf)*PDIcarr + oldfllx;
                        flly = 0.5*(fllx + oldfllx) + a2*omegaOf*carrErrorFreq;
                        carrFreq=oldCarrFreq+(flly-oldflly)*PDIcarr;
                        oldfllx = fllx;
                        oldflly=flly;
                        oldCarrFreq=carrFreq;            
    else              

    % carrNco = carrNco ;
        carrError = carrError * 1;
        carrErrorFreq  = carrErrorFreq * 1;                
        pllX = oldPllX + carrErrorFreq*(omegaOf*omegaOf)*PDIcarr + carrError*(wn*wn*wn)*PDIcarr;
        pllY = oldPllY + carrErrorFreq*1.414*omegaOf*PDIcarr + 0.5*(pllX + oldPllX)*PDIcarr+1.1*carrError* wn*wn*PDIcarr;
        carrNco =carrError* 2.4*wn+0.5*(pllY+oldPllY);
        oldPllX = pllX;
        oldPllY = pllY; 
        carrFreq = fi + carrNco;
        % carrFreq = fi;
    end
%         carrNco = oldCarrNco + (tau2carr/tau1carr) * ...
%             (carrError - oldCarrError) + carrError * (PDIcarr/tau1carr);
% 
% %         oldy = y;
%         oldCarrNco   = carrNco;
%         oldCarrError = carrError;
% 
%         % Modify carrier freq based on NCO command
%         carrFreq = fi + carrNco;

%% Find DLL error and update code NCO -------------------------------------
    codeError = (sqrt(I_E * I_E + Q_E * Q_E) - sqrt(I_L * I_L + Q_L * Q_L)) / ...
        (sqrt(I_E * I_E + Q_E * Q_E) + sqrt(I_L * I_L + Q_L * Q_L));

    % Implement code loop filter and generate NCO command
    
%         a = olda + (PDIcarr * wn^3)/2 * (codeError - oldCodeError) + (1.1 * wn^2) * (codeError - oldCodeError);
%         codeNco = oldCodeNco + (PDIcarr/2) * (a + olda) + (2.4 * wn) * (codeError - oldCodeError);
    
    codeNco = oldCodeNco + (tau2code/tau1code) * ...
        (codeError - oldCodeError) + codeError * (PDIcode/tau1code);
%         olda = a;
    oldCodeNco   = codeNco;
    oldCodeError = codeError;

    % Modify code freq based on NCO command
    codeFreq = fp - codeNco;
    % codeFreq = fp;

    % 使用不同类型的卡尔曼滤波器
    z = [codePhase; carrPhase];
    
    kf_standard = kalman_filters('standard', kf_standard, z, dt);
    kf_extended = kalman_filters('extended', kf_extended, z, dt);
    kf_unscented = kalman_filters('unscented', kf_unscented, z, dt);
    
    % 使用标准卡尔曼滤波器的估计结果
    kf_codePhase = kf_standard.x(1);
    kf_codeFreq = kf_standard.x(2);
    kf_carrPhase = kf_standard.x(3);
    kf_carrFreq = kf_standard.x(4);
    
    % 保存不同卡尔曼滤波器的结果
    output.KF_Standard_CodePhase(loopCnt) = kf_standard.x(1);
    output.KF_Standard_CodeFreq(loopCnt) = kf_standard.x(2);
    output.KF_Standard_CarrPhase(loopCnt) = kf_standard.x(3);
    output.KF_Standard_CarrFreq(loopCnt) = kf_standard.x(4);
    
    output.KF_Extended_CodePhase(loopCnt) = kf_extended.x(1);
    output.KF_Extended_CodeFreq(loopCnt) = kf_extended.x(2);
    output.KF_Extended_CarrPhase(loopCnt) = kf_extended.x(3);
    output.KF_Extended_CarrFreq(loopCnt) = kf_extended.x(4);
    
    output.KF_Unscented_CodePhase(loopCnt) = kf_unscented.x(1);
    output.KF_Unscented_CodeFreq(loopCnt) = kf_unscented.x(2);
    output.KF_Unscented_CarrPhase(loopCnt) = kf_unscented.x(3);
    output.KF_Unscented_CarrFreq(loopCnt) = kf_unscented.x(4);

    output.OutCarrFreq(loopCnt)=carrFreq;
    output.OutCodeFreq(loopCnt)=codeFreq;
    output.OutCarrPhase(loopCnt)=carrPhase;
    output.OutCodePhase(loopCnt)=codePhase;
    output.OutCarrNCO(loopCnt)=carrError;
    output.OutCodeNCO(loopCnt)=codeError;
    output.OutBlksize(loopCnt)=blksize;
    output.OutRemCode(loopCnt)=remCodePhase;
    output.I_E(loopCnt)=I_E;
    output.Q_E(loopCnt)=Q_E;
    output.I_P(loopCnt)=I_P;
    output.Q_P(loopCnt)=Q_P;
    output.I_L(loopCnt)=I_L;
    output.Q_L(loopCnt)=Q_L;
    
    % 添加卡尔曼滤波器的结果到输出
    output.KF_CodePhase(loopCnt) = kf_codePhase;
    output.KF_CodeFreq(loopCnt) = kf_codeFreq;
    output.KF_CarrPhase(loopCnt) = kf_carrPhase;
    output.KF_CarrFreq(loopCnt) = kf_carrFreq;
end

% 添加卡尔曼滤波器的最终状态到输出
output.KF_FinalState = kf.x;
output.KF_FinalCovariance = kf.P;
