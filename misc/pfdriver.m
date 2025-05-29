%Created for Dynare 4.6.4

dynare nkpf
tic
T=500;
kstate = oo_.dr.kstate;
nstatic = M_.nstatic;
nfwrd = M_.nfwrd;
nspred = M_.nspred;
nboth = M_.nboth;
nsfwrd = M_.nsfwrd;
nd = size(kstate,1);
nz = nnz(M_.lead_lag_incidence);
dynfun = [M_.fname, '.dynamic'];
y = repmat(oo_.steady_state, M_.maximum_lag + M_.maximum_lead + 1, 1);
x = zeros(1, M_.exo_nbr + M_.exo_det_nbr);
gx = oo_.dr.gx;
[~, g1] = feval(dynfun, y, x, M_.params, oo_.steady_state, 1);
f1 = sparse(g1(:,nonzeros(M_.lead_lag_incidence(M_.maximum_endo_lag+2:end,oo_.dr.order_var))));
f0 = sparse(g1(:,nonzeros(M_.lead_lag_incidence(M_.maximum_endo_lag+1,oo_.dr.order_var))));
fudet = sparse(g1(:,nz+M_.exo_nbr+1:end));
M1 = inv(f0+[zeros(M_.endo_nbr,nstatic) f1*gx zeros(M_.endo_nbr,nsfwrd-nboth)]);
M2 = M1*f1;

C = cell(T,1);
C{1} = -M1*fudet;
for i = 2:T
    C{i} = -M2*C{i-1}(end-nsfwrd+1:end,:);
end
A = [oo_.dr.ghx zeros(4,2)];
B = oo_.dr.ghu;
Sigma = eye(3);
Sigma(1,1) = .03;
Sigma(2,2) = .0025;
Sigma(3,3) = .01; 
% Generate stochastic shocks
eps = mvnrnd(zeros(3,1), Sigma, T)';  % [m x T]
x = 1.5 * Sigma(1,1)/10000 * randn(1, T);
n = size(A,1);
m = size(B,2);
k = size(x,1);

% Preallocate
y = zeros(n, T+1);


%eps = oo_.exo_simul';
% Simulate
for t = 1:T
    deterministic = zeros(n,1);
    for i = 0:(T - t)
        deterministic = deterministic + C{i+1} * x(t + i);
    end 
    y(:,t+1) = A * y(:,t) + B * eps(:,t) + deterministic;
end
toc 
