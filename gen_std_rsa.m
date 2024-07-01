function[T, RSA] = gen_std_rsa(T_g, zeta, ag)
    T = linspace(0, 6, 500);
    RSA = zeros(size(T));
    for i = 1:length(T)
        RSA(i) = rsa(T_g, zeta, T(i), ag);
    end
end