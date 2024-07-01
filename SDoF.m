function [u,v] = SDoF(ag, zeta, omg, dt)
%SDOF ag zeta
%   此处显示详细说明
omg_d = omg * sqrt(1.0 - zeta ^ 2);
n = length(ag);
u = zeros(1,n);
v = zeros(1,n);

B1 = exp(-zeta * omg * dt) * cos(omg_d * dt);
B2 = exp(-zeta * omg * dt) * sin(omg_d * dt);

omg_2 = 1.0 / omg / omg;
omg_3 = 1.0 / omg / omg / omg;
for i = 1:n - 1
    u_i = u(i);
    v_i = v(i);
    p_i = -ag(i);
    alpha_i = (-ag(i+1) + ag(i)) / dt;

    A0 = p_i * omg_2 - 2.0 * zeta * alpha_i * omg_3;
    A1 = alpha_i * omg_2;
    A2 = u_i - A0;
    A3 = (v_i + zeta * omg * A2 - A1) / omg_d;

    u(i+1) = A0 + A1 * dt + A2 * B1 + A3 * B2;
    v(i+1) = A1 + (omg_d * A3 - zeta * omg * A2) * B1...
            -(omg_d * A2 + zeta * omg * A3 ) * B2;
end

end

