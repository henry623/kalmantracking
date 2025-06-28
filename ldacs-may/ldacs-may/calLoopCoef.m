function [taul1,taul2]=calLoopCoef(noiseBandwidth,Zeta,K)
% Bl---noiseBandwidth
%zeta--dampingRatio
%K--loopGain


Wn=noiseBandwidth*8*Zeta/(1+4*Zeta*Zeta);
taul1=K/(Wn*Wn);
taul2=Zeta*2/Wn;