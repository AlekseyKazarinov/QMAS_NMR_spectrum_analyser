function [Y] = calcY(phi,mu,eta)
% calcY - расчёт члена в формуле для частоты, зависящего от eta, mu, phi
    z = cos(2*phi);
    D = 21/16 - 7/8*eta*z + 7/48*eta*eta*z*z;
    E = -9/8 + eta*eta/12 + eta*z - 7/24*eta*eta*z*z;
    F = 5/16 - 1/8*eta*z + 7/48*eta*eta*z*z;
    mu2 = mu*mu;
    Y = D*mu2*mu2 + E*mu2 + F;
end

