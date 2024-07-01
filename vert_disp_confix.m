function [rate] = vert_disp_confix(ke, k1h,xy)
%VERT_DISP_CONFIX 此处显示有关此函数的摘要
%   此处显示详细说明
    rate = 1;
    rate = max(rate, 0.9179 + 0.0956 * ke / k1h + 0.0263 * xy);
    rate = rate / 0.97;
end

