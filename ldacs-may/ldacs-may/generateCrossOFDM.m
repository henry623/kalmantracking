% 生成OFDM网格与时域信号，跨频段
function [grid, rawSignal]=generateCrossOFDM(simSettings)

nSymbol = simSettings.nSymbol;
% code = repmat(simSettings.code,1,3);
code = simSettings.code;
Nu = simSettings.Nu;
numBand = simSettings.numBand;
NFFT = simSettings.NFFT;

grid = zeros(nSymbol,NFFT,numBand);
rawSignal = zeros(1,nSymbol*NFFT);

for i=1:nSymbol
    if i==nSymbol
        temp = code(1+(i-1)*Nu*numBand:end);
        for j=1:ceil(size(temp,2)/Nu)
            grid(i,8:8+size(temp,2)-1,j) = temp;
        end
    else
        temp = code(1+(i-1)*Nu*numBand:i*Nu*numBand);
        for j=1:numBand
            grid(i,8:8+Nu-1,j) = temp(1+(j-1)*Nu:1+j*Nu-1);
        end
    end

    % temp = code(1+(i-1)*Nu*numBand:i*Nu*numBand);
    % for j=1:numBand
    %     grid(i,8:8+Nu-1,j) = temp(1+(j-1)*Nu:1+j*Nu-1);
    % end

    for j=1:numBand
        rawSignal(j,(i-1)*NFFT+1:i*NFFT) = ifft(grid(i,:,j),NFFT); % 对资源网格进行IFFT
    end
end