function [kb, zeta_b, beta] = compute_hori_beta(ke, kp, xy, ks, zeta_s, zeta_h, ms, mb, Tg, ag)
    m = ms + mb;
    xm = 100;
    times = 0;
    while 1
        kb = 0;
        zeta_b = 0;
        for i=1:length(ke)
            if xm <= xy(i)
                kb = kb + ke(i);
            else
                kbi = (ke(i) * xy(i) + kp(i) * (xm - xy(i))) / xm;
                zeta_bi = 2 * xy(i) * (ke(i) - kp(i)) * ...
                    (xm - xy(i)) / pi / kbi / xm / xm;
                kb = kb + kbi;
                zeta_b = zeta_bi * kb;
            end
        end
        zeta_b = zeta_b / kb + zeta_h;
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
        
%         zeta_1 = (sqrt(ep) * zeta_b + ga * x1 ^ 2 * zeta_s) / ...
%             sqrt(C) / (1 + 2 * ga * x1 + ga * x1 ^2);
%         zeta_2 = (sqrt(ep) * zeta_b + ga * x2 ^ 2 * zeta_s) / ...
%             sqrt(D) / (1 + 2 * ga * x2 + ga * x2 ^2);
        zeta_1 = zeta_b + zeta_s ;
        zeta_2 = zeta_s ;
        %disp(zeta_1)
        %disp(zeta_2)
        L1 = (1 + ga * x1) / (1 + 2 * ga * x1 + ga * x1 ^2);
        L2 = (1 + ga * x2) / (1 + 2 * ga * x2 + ga * x2 ^2);
    
        vSA_1 = rsa(Tg, zeta_1, 2 * pi / omega_1, ag) * 9806;
        vSA_2 = rsa(Tg, zeta_2, 2 * pi / omega_2, ag) * 9806;
        q1 = L1 * vSA_1 / (omega_1 ^ 2);
        q2 = L2 * vSA_2 / (omega_2 ^ 2);
        xm_new = sqrt(q1 ^ 2 + q2 ^ 2);
        if abs(xm - xm_new) < 1e-1 || times > 10000
            vSA_s = rsa(Tg, zeta_s, 2 * pi / omega_s, ag) * 9806;
            vi = sqrt((x1 * q1) ^ 2 + (x2 * q2) ^2) * ks;
            disp(vi)
            vf = ms * vSA_s;
            disp(vf)
            beta = vi / vf;
            break
        end
        xm = xm_new;
        times = times + 1;
    end
end
