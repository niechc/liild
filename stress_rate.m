function [stressRate] = stress_rate(T_0, rT)
%STRESS_RATE 此处显示有关此函数的摘要
%   此处显示详细说明
stressRate = -219.8 .* 0.49 .^ rT .* (T_0 - 1.15) .^ 2 ...
             - 44.7 .* 0.58 .^ rT + 100;
end

