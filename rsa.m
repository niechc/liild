function[RSA] = rsa(T_g, zeta, T, ag)
    gamma=0.9+(0.05-zeta)/(0.3+6*zeta);
    eta1=0.02+(0.05-zeta)/(4+32*zeta);
    eta2=1+(0.05-zeta)/(0.08+1.6*zeta);
    alpha_max = 2.25;
    RSA_0 = 1.0/alpha_max;

    if T < 0.1
        RSA = RSA_0 + 10.0 * (eta2 - RSA_0) * T;
    elseif T >= 0.1 && T < T_g
        RSA = eta2;
    elseif T >= T_g && T < 5 * T_g
        RSA = eta2 * (T_g / T) ^ gamma;
    else
        RSA = eta2*0.2^gamma-eta1*(T-5*T_g);
    end

    RSA = RSA * alpha_max * ag;
end