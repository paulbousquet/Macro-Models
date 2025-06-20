%% Parameters and grid
a_e = 0.11; a_h = 0.03; % production rates
rho_0 = 0.04; % time preference
rho_e_d = 0.01; rho_h_d = 0.01; % death rates
rho_e = rho_0 + rho_e_d; % expert’s discount rate
rho_h = rho_0 + rho_h_d; % household’s discount rate
zeta = 0.05;
delta = 0.05; sigma = 0.1;
phi = 10; alpha = 0.5;

N = 501; % grid size
eta = linspace(0.0001, 0.999, N)';

%% Solution
% Solve for q(0)
q0 = (1 + a_h*phi)/(1 + rho_h*phi);

% Inner loop
[Q, SSQ, Kappa, Chi, Iota] = inner_loop_log(eta, q0, a_e, a_h, rho_e, rho_h, sigma, phi, alpha);

S = (Chi - eta) .* SSQ;
Sg_e = S ./ eta;
Sg_h = -S ./ (1 - eta);

VarS_e = Chi ./ eta .* SSQ;
VarS_h = (1 - Chi) ./ (1 - eta) .* SSQ;

CN_e = rho_e;
CN_h = rho_h;

MU = eta .* (1 - eta) .* ( ...
    (VarS_e - SSQ) .* (Sg_e + SSQ) ...
  - (VarS_h - SSQ) .* (Sg_h + SSQ) ...
  - (CN_e - CN_h) ...
  + (rho_h_d * zeta * (1 - eta) - rho_e_d * (1 - zeta) * eta) ./ (eta .* (1 - eta)) ...
);
