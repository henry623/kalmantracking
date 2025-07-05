function [tau1, tau2] = calLoopCoef(noiseBandwidth, dampingRatio, loopGain)
% 计算环路滤波器系数
%
% 功能描述:
%   根据环路噪声带宽、阻尼比和环路增益计算二阶环路滤波器的系数
%   用于DLL和PLL的环路滤波器设计
%
% 输入参数:
%   noiseBandwidth - 环路噪声带宽 (Hz)
%   dampingRatio   - 阻尼比 (无量纲)
%   loopGain      - 环路增益 (无量纲)
%
% 输出参数:
%   tau1 - 环路滤波器系数1
%   tau2 - 环路滤波器系数2
%
% 数学原理:
%   基于二阶环路的传递函数设计
%   Wn = 8*dampingRatio*noiseBandwidth / (4*dampingRatio^2 + 1)
%   tau1 = loopGain / (Wn * Wn)
%   tau2 = 2.0 * dampingRatio / Wn
%
% 作者: OFDM接收机开发
% 日期: 2025年7月
% 版本: 1.0

%% 输入参数验证
if nargin < 3
    error('calLoopCoef: 需要提供所有三个输入参数');
end

if noiseBandwidth <= 0
    error('calLoopCoef: 噪声带宽必须为正数');
end

if dampingRatio <= 0
    error('calLoopCoef: 阻尼比必须为正数');
end

if loopGain <= 0
    error('calLoopCoef: 环路增益必须为正数');
end

%% 计算自然频率
% 根据噪声带宽和阻尼比计算自然频率
% 这是基于二阶环路系统的标准公式
Wn = 8 * dampingRatio * noiseBandwidth / (4 * dampingRatio^2 + 1);

%% 计算环路滤波器系数
% tau1 对应积分器系数
tau1 = loopGain / (Wn * Wn);

% tau2 对应比例系数
tau2 = 2.0 * dampingRatio / Wn;

%% 参数合理性检查
if tau1 <= 0 || tau2 <= 0
    warning('calLoopCoef: 计算得到的系数可能不合理，请检查输入参数');
end

if tau1 > 1000 || tau2 > 1000
    warning('calLoopCoef: 计算得到的系数过大，可能导致环路不稳定');
end

end