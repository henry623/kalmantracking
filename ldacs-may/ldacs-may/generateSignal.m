function [record_E]=generateSignal(simSettings,yr,yi)

%% 仿真设置
SNR = simSettings.SNR;
c = simSettings.c;
dt = simSettings.dt;
t_total = simSettings.t_total;

t0 = 0:dt:t_total;
p0 = simSettings.p0;
v0 = simSettings.v0;

BS = simSettings.BS;

p = p0+v0.*t0'; % 实际位置
v = v0+0.*t0';
code = simSettings.code;
lenOFDM = simSettings.NFFT*simSettings.nSymbol;
% lenOFDM = 10230;


fp = simSettings.fp;
fs = simSettings.fs;
fi = simSettings.fi;
fc = simSettings.fc;
remCodePhase = 0;
remCarrPhase = 0;

%% generate
samplesPerMs = ceil(fs*dt)-1;

record_E0=[];
for t=1:size(t0,2)

    rawSignal_ms = [];
    pos = p(t,:);
    vel = v(t,:);

    for bs=1:size(BS,1)


        %% 信源发送的信号
        % while(size(rawSignal_ms,2) < samplesPerMs) 
            codeFreq = fp;
            carrFreq = fi;
    
            codePhaseStep = codeFreq / fs;
    
            blksize = ceil((lenOFDM-remCodePhase) / codePhaseStep);
    
            rawSignal=1;
    
            tcode       = remCodePhase : ...
                          codePhaseStep : ...
                          ((blksize-1)*codePhaseStep+remCodePhase);
            tcode2 = round(tcode*simSettings.famp+1)+simSettings.extra*simSettings.famp;

            remCodePhase = tcode(blksize) + codePhaseStep - lenOFDM;
    
            rawSignalr   = rawSignal * yr(tcode2);
            rawSignali   = rawSignal * yi(tcode2);
    
            time    = (0:blksize) ./ fs;
            trigarg = ((carrFreq * 2.0 * pi) .* time) + 2.0 * pi * remCarrPhase;
            remCarrPhase = rem(trigarg(blksize+1), (2 * pi))/(2.0 * pi);
    
            rawSignalr = rawSignalr .* cos(trigarg(1:blksize));
            rawSignali = rawSignali .* sin(trigarg(1:blksize));
    
            rawSignal = [rawSignalr;rawSignali];
    
            % rawSignal = rawSignalr + rawSignali;
    
            rawSignal_ms = [rawSignal_ms rawSignal];
    
        % end

        % rawSignal_Ebs=zeros(1,blksize);
        % x = rawSignal_ms;
        % samplingN = size(x,2);
        % rawSignal_ms = rawSignal_ms(samplesPerMs+1:end);
    
        % %% 接收机收到的信号
        % fsc=(-samplingN/2:samplingN/2-1).*fs/samplingN;
        % % 信道矩阵
        % k = 2*pi.*fsc/c;
        % 
        % rho=norm(pos-BS(bs,:));
        % 
        % alpha = 10.^((32.45+20*log10(rho/1e3)+20*log10(fc/1e6))/10);
        % tau = rho/c;
        % H = alpha.*exp(-1i.*k.*tau*c);
        % 
        % X = fftshift(fft(x,samplingN)./samplingN);
        % 
        % Y = H.*X;
        % y = ifft(ifftshift(Y).*samplingN);
        % 
        % rawSignal_E0=real(y(1:blksize));
        % 
        % rawSignal_Ebs=rawSignal_Ebs+rawSignal_E0;
        % rho = norm(pos-BS(bs,:));
        % alpha = 10.^((32.45+20*log10(rho/1e3)+20*log10(fc/1e6))/10);
        rawSignal_Ebs=rawSignal_ms;
    end

    rawSignal_E0=awgn(rawSignal_Ebs,SNR,'measured');

%     比特量化
    % amp = max(abs(rawSignal_E(:)));
    % q = amp/(2^(2-1)); 
    % for i=1:samplesPerMs
    %     for j=1:2
    %         if(rawSignal_E(i)>=0) 
    %             rawSignal_E(j,i) =  min(ceil(rawSignal_E(i)/q), 2^(3-1));
    %         else
    %             rawSignal_E(j,i) = -min(ceil(-rawSignal_E(i)/q), 2^(3-1));
    %         end
    %     end
    % end
    
    record_E0=[record_E0,rawSignal_E0];

end

rho=norm(p(1,:)-BS(1,:));
tau = rho/c;
tau_index = ceil(tau*fs)+1;
record_E = record_E0(:,tau_index:end);
filenm1 = ['output/signal_SNR'  num2str(SNR) '.mat' ];