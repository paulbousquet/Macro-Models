/*********************************************************************
  Basic 3-eq New-Keynesian model with two monetary-policy shocks
*********************************************************************/
var pi x i rn;                // output gap, inflation, output gap, rate, nat. rate
varexo  s vrn u;                      // stochastic MP shock ε_s,t
varexo_det f;                   // deterministic “target” shock f_t (fully anticipated)

parameters beta sigma kappa phi rho_s sigma_s;
beta   = 0.99;
sigma  = 1;                     // 1/σ is the IES
kappa  = 0.024;                 // Calvo(θ)=0.75 → κ ≈ (1-θ)(1-βθ)/θ
phi    = 1.5;                   // Taylor coefficient
rho_s  = 0.8;
sigma_s= 0.03;                  // stdev(ε_s)

model(linear);
  /* IS */          x = x(+1) - (1/sigma)*( i - pi(+1) - rn );
  /* NKPC */        pi = beta*pi(+1) + kappa*x+u;
  /* Taylor */      i =  .8*i(-1)+phi*pi  + f + s;
  /* Nat. rate */   rn = .9*rn(-1)+vrn;
end;
verbatim;
% Set seed for reproducibility
rng(123);

% Generate the shock sequence
f_values = 1.5 * 0.03 * randn(1, 20000);
end;
shocks;
  var s; stderr sigma_s;
  var vrn; stderr .0025;
  var u; stderr .01; 
  var f;
  periods 1:5;
  values .02;
end;

steady;
check;


verbatim;
  % assign shock values
  f_values = 0 * randn(20000,1);
  %oo_.exo_det_simul(:, find(strcmp(M_.exo_det_names, 'f'))) = f_values';
end;
/* Only compute decision rules; IRFs not needed */
stoch_simul(order=1,irf=0,periods=200);


