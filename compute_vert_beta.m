function [vbeta] = compute_vert_beta(kb, zeta_b, ks, zeta_s,  ms, mb, Tg, ag)
    m = ms + mb;
    omega_s = sqrt(ks / ms);
    omega_b = sqrt(kb / m);

    ep = (omega_b / omega_s) ^ 2;
    ga = ms / m;

    Delta = sqrt(4 * ep * ga - 2 * ep + ep ^ 2 + 1);
    C = -(1 + ep - Delta) / 2 / (ga - 1);
    D = -(1 + ep + Delta) / 2 / (ga - 1);

    omega_1 = sqrt(C * omega_s ^ 2);
    omega_2 = sqrt(D * omega_s ^ 2);

    x1 = (ep - C) / ga / C;
    x2 = D / (1 - D);
%     zeta_1 = (sqrt(ep) * zeta_b + ga * x1 ^ 2 * zeta_s) / ...
%         sqrt(C) / (1 + 2 * ga * x1 + ga * x1 ^2);
%     zeta_2 = (sqrt(ep) * zeta_b + ga * x2 ^ 2 * zeta_s) / ...
%         sqrt(D) / (1 + 2 * ga * x2 + ga * x2 ^2);
    zeta_1 = zeta_b ;
    zeta_2 = zeta_s ;
    L1 = (1 + ga * x1) / (1 + 2 * ga * x1 + ga * x1 ^2);
    L2 = (1 + ga * x2) / (1 + 2 * ga * x2 + ga * x2 ^2);

    vSA_1 = rsa(Tg, zeta_1, 2 * pi / omega_1, ag) * 9806;
    vSA_2 = rsa(Tg, zeta_2, 2 * pi / omega_2, ag) * 9806;
    q1 = L1 * vSA_1 / (omega_1 ^ 2);
    q2 = L2 * vSA_2 / (omega_2 ^ 2);

    vSA_s = rsa(Tg, zeta_s, 2 * pi / omega_s, ag) * 9806;
    Ri = sqrt((x1 * q1) ^ 2 + (x2 * q2) ^2) * ks;
    disp(Ri)
    Rf = ms * vSA_s;
    disp(Rf)
    vbeta = Ri / Rf;
end

