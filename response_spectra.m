function [T, RSA] = response_spectra(ag, dt, zeta)
    N = 500;
    T = linspace(dt, 6, N);
    RSA = zeros(1,N);
    for i = 1:N
        omg = 2 * pi / T(i);
        [u,v] = SDoF(ag, zeta, omg, dt);
        a = - 2.0*zeta*omg*v - omg*omg*u;
        RSA(i) = max(abs(a));
    end
end