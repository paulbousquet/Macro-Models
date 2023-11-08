function [nx, ny,fx,fxp,fy,fyp,fypyp,fypy,fypxp,fypx,fyyp,fyy,fyxp,fyx,fxpyp,fxpy,fxpxp,fxpx,fxyp,fxy,fxxp,fxx,...
    fypypyp,fypypy,fypypxp,fypypx,fypyyp,fypyy,fypyxp,fypyx,fypxpyp,fypxpy,fypxpxp,fypxpx,fypxyp,fypxy,fypxxp,fypxx,...
    fyypyp,fyypy,fyypxp,fyypx,fyyyp,fyyy,fyyxp,fyyx,fyxpyp,fyxpy,fyxpxp,fyxpx,fyxyp,fyxy,fyxxp,fyxx,...
    fxpypyp,fxpypy,fxpypxp,fxpypx,fxpyyp,fxpyy,fxpyxp,fxpyx,fxpxpyp,fxpxpy,fxpxpxp,fxpxpx,fxpxyp,fxpxy,fxpxxp,fxpxx,...
    fxypyp,fxypy,fxypxp,fxypx,fxyyp,fxyy,fxyxp,fxyx,fxxpyp,fxxpy,fxxpxp,fxxpx,fxxyp,fxxy,fxxxp,fxxx,f] = eir_model

syms RSTAR BETTA DBAR DELTAA ALFA  PHI  RHO  STD_EPS_A ETASHOCK SIGG OMEGA PSSI

syms c cp h hp d dp k kp kfu kfup r rp a ap

%Equilibrium conditions. The symbols e1, e2, ... denote equation 1, equation2, ...

%Evolution of debt
e1 = -dp + (1+r)*d + c + kp - (1-DELTAA)*k + .5*PHI*(kp -k)^(2) - a*k^ALFA*h^(1-ALFA); 

%FOC w.r.t. h
e4 = -h^(OMEGA-1)+ (1-ALFA) * a * (k/h)^ALFA;

%FOC w.r.t. debt
e5 = -(c-h^OMEGA/OMEGA)^(-SIGG)  + BETTA * (1+rp) * (cp-hp^OMEGA/OMEGA)^(-SIGG) ;

%FOC w.r.t. capital
e6 = -(c-h^OMEGA/OMEGA)^(-SIGG)* (1+PHI*(kp-k)) + BETTA * (cp-hp^OMEGA/OMEGA)^(-SIGG) * (1-DELTAA + ALFA * ap * (kp/hp)^(ALFA-1) + PHI * (kfup-kp));

%Country premium
e7 = -rp + RSTAR + PSSI * (exp(dp-DBAR) -1);

%Evolution of TFP
e13 = -log(ap) + RHO * log(a);

%make kfu=kp
e14 = -kfu+kp;

%Create function f
f = [e1;e4;e5;e6;e7;e13;e14];

% Define the vector of controls in periods t and t+1, controlvar and controlvarp, and the vector of states in periods t and t+1, statevar and statevarp

%States to substitute from levels to logs

states_in_logs = [k a];
states_in_logsp = [kp ap];

x = [d r states_in_logs];

xp = [dp rp states_in_logsp];
 
y = [c h kfu];

yp = [cp hp kfup];

%Number of states
nx = length(x);
ny = length(y);

%Make f a function of the logarithm of the state and control vector

%variables to substitute from levels to logs
%variables_in_logs = transpose([states_in_logs, controls_in_logs, states_in_logsp, controls_in_logsp]);

%f = subs(f, variables_in_logs, exp(variables_in_logs));

% Make f a function of the logarithm of the state and control vector --  remember to change steady state and data accordingly
%f = subs(f, [x,y,xp,yp], exp([x,y,xp,yp]));

% Compute analytical derivatives of f
global approx
[fx,fxp,fy,fyp,fypyp,fypy,fypxp,fypx,fyyp,fyy,fyxp,fyx,fxpyp,fxpy,fxpxp,fxpx,fxyp,fxy,fxxp,fxx,...
    fypypyp,fypypy,fypypxp,fypypx,fypyyp,fypyy,fypyxp,fypyx,fypxpyp,fypxpy,fypxpxp,fypxpx,fypxyp,fypxy,fypxxp,fypxx,...
    fyypyp,fyypy,fyypxp,fyypx,fyyyp,fyyy,fyyxp,fyyx,fyxpyp,fyxpy,fyxpxp,fyxpx,fyxyp,fyxy,fyxxp,fyxx,...
    fxpypyp,fxpypy,fxpypxp,fxpypx,fxpyyp,fxpyy,fxpyxp,fxpyx,fxpxpyp,fxpxpy,fxpxpxp,fxpxpx,fxpxyp,fxpxy,fxpxxp,fxpxx,...
    fxypyp,fxypy,fxypxp,fxypx,fxyyp,fxyy,fxyxp,fxyx,fxxpyp,fxxpy,fxxpxp,fxxpx,fxxyp,fxxy,fxxxp,fxxx] = anal_deriv(f,x,y,xp,yp,approx);





