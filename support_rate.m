function [supportRate] = support_rate(T_0, rT)
%SUPPORT_RATE 此处显示有关此函数的摘要
%   此处显示详细说明
supportRate = -238.7 .* 0.57 .^ rT .* (T_0 - 1.15) .^ 2 ...
             - 36.2 .* 0.71 .^ rT + 100;
end

