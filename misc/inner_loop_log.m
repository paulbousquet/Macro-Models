function [Q, SSQ, Kappa, Chi, Iota] = inner_loop_log(eta, q0, a_e, a_h, rho_e, rho_h, sigma, phi, alpha)

N = length(eta);
deta = [eta(1); diff(eta)];

% Initialize variables
Q = ones(N, 1);
SSQ = zeros(N, 1);
Kappa = zeros(N, 1);

Rho = eta * rho_e + (1 - eta) * rho_h;

kappa = 0; q_old = q0; q = q0; ssq = sigma;

for i = 1:N
    F = [kappa * (a_e - a_h) + a_h - (q - 1)/phi - q * Rho(i);
         ssq * (q - (q - q_old)/deta(i) * (alpha * kappa - eta(i))) - sigma * q;
         a_e - a_h - q * alpha * (alpha * kappa - eta(i)) / (eta(i) * (1 - eta(i))) * ssq^2];

    J = zeros(3, 3);
    J(1, :) = [-1/phi - Rho(i), a_e - a_h, 0];
    J(2, :) = [ssq * (1 - (alpha * kappa - eta(i)) / deta(i)) - sigma, ...
               -ssq * (q - q_old) / deta(i) * alpha, ...
               q - (q - q_old) / deta(i) * (alpha * kappa - eta(i))];
    J(3, :) = [-alpha * (alpha * kappa - eta(i)) / (eta(i) * (1 - eta(i))) * ssq^2, ...
               -q * alpha^2 / (eta(i) * (1 - eta(i))) * ssq^2, ...
               -2 * q * alpha * (alpha * kappa - eta(i)) / (eta(i) * (1 - eta(i))) * ssq];

    z = [q, kappa, ssq]' - J \ F;

    if z(2) >= 1
        break;
    end

    q = z(1); kappa = z(2); ssq = z(3);
    Q(i) = q; Kappa(i) = kappa; SSQ(i) = ssq;
    q_old = q;
end

n1 = i;
for i = n1:N
    q = (1 + a_e * phi) / (1 + Rho(i) * phi);
    qp = (q - q_old) / deta(i);

    Q(i) = q; Kappa(i) = 1;
    SSQ(i) = sigma / (1 - (max(alpha, eta(i)) - eta(i)) * qp / q);
    q_old = q;
end

Chi = max(alpha * Kappa, eta);
Iota = (Q - 1) / phi;

end
