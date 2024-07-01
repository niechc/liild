function [deformRate] = deform_rate(T_0, rT)
%SUPPORT_RATE 此处显示有关此函数的摘要
%   此处显示详细说明
deformRate = -95.4 .* 0.42 .^ rT .* T_0 ...
             - 141.6 .* 0.51 .^ rT + 100;
end
