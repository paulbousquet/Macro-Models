clear all

global approx
global fx fxp fy fyp fypyp fypy fypxp fypx fyyp fyy fyxp fyx fxpyp fxpy fxpxp fxpx fxyp fxy fxxp fxx
global fypypyp fypypy fypypxp fypypx fypyyp fypyy fypyxp fypyx fypxpyp fypxpy fypxpxp fypxpx fypxyp fypxy fypxxp fypxx
global fyypyp fyypy fyypxp fyypx fyyyp fyyy fyyxp fyyx fyxpyp fyxpy fyxpxp fyxpx fyxyp fyxy fyxxp fyxx
global fxpypyp fxpypy fxpypxp fxpypx fxpyyp fxpyy fxpyxp fxpyx fxpxpyp fxpxpy fxpxpxp fxpxpx fxpxyp fxpxy fxpxxp fxpxx
global fxypyp fxypy fxypxp fxypx fxyyp fxyy fxyxp fxyx fxxpyp fxxpy fxxpxp fxxpx fxxyp fxxy fxxxp fxxx f
global nfx nfxp nfy nfyp nfypyp nfypy nfypxp nfypx nfyyp nfyy nfyxp nfyx nfxpyp nfxpy nfxpxp nfxpx nfxyp nfxy nfxxp nfxx
global nfypypyp nfypypy nfypypxp nfypypx nfypyyp nfypyy nfypyxp nfypyx nfypxpyp nfypxpy nfypxpxp nfypxpx nfypxyp nfypxy nfypxxp nfypxx
global nfyypyp nfyypy nfyypxp nfyypx nfyyyp nfyyy nfyyxp nfyyx nfyxpyp nfyxpy nfyxpxp nfyxpx nfyxyp nfyxy nfyxxp nfyxx
global nfxpypyp nfxpypy nfxpypxp nfxpypx nfxpyyp nfxpyy nfxpyxp nfxpyx nfxpxpyp nfxpxpy nfxpxpxp nfxpxpx nfxpxyp nfxpxy nfxpxxp nfxpxx
global nfxypyp nfxypy nfxypxp nfxypx nfxyyp nfxyy nfxyxp nfxyx nfxxpyp nfxxpy nfxxpxp nfxxpx nfxxyp nfxxy nfxxxp nfxxx

%Order of approximation desired 
approx = 3;

[nx,ny,fx,fxp,fy,fyp,fypyp,fypy,fypxp,fypx,fyyp,fyy,fyxp,fyx,fxpyp,fxpy,fxpxp,fxpx,fxyp,fxy,fxxp,fxx,...
    fypypyp,fypypy,fypypxp,fypypx,fypyyp,fypyy,fypyxp,fypyx,fypxpyp,fypxpy,fypxpxp,fypxpx,fypxyp,fypxy,fypxxp,fypxx,...
    fyypyp,fyypy,fyypxp,fyypx,fyyyp,fyyy,fyyxp,fyyx,fyxpyp,fyxpy,fyxpxp,fyxpx,fyxyp,fyxy,fyxxp,fyxx,...
    fxpypyp,fxpypy,fxpypxp,fxpypx,fxpyyp,fxpyy,fxpyxp,fxpyx,fxpxpyp,fxpxpy,fxpxpxp,fxpxpx,fxpxyp,fxpxy,fxpxxp,fxpxx,...
    fxypyp,fxypy,fxypxp,fxypx,fxyyp,fxyy,fxyxp,fxyx,fxxpyp,fxxpy,fxxpxp,fxxpx,fxxyp,fxxy,fxxxp,fxxx,f] = eir_model;

%Numerical Evaluation
%Steady State and Parameter Values
eir_ss

eta=[0; 0; 0; STD_EPS_A];
eta=reshape(eta,[4  1]);
esk=eta*0+1;

%Obtain numerical derivatives of f
num_eval

%First-order approximation
[gx,hx] = gx_hx(nfy,nfx,nfyp,nfxp);

%approx = 2;

%Second-order approximation
if approx > 1
    [gxx,hxx] = gxx_hxx(nfx,nfxp,nfy,nfyp,nfypyp,nfypy,nfypxp,nfypx,nfyyp,nfyy,nfyxp,nfyx,nfxpyp,nfxpy,nfxpxp,nfxpx,nfxyp,nfxy,nfxxp,nfxx,hx,gx);
    [gss,hss] = gss_hss(nfx,nfxp,nfy,nfyp,nfypyp,nfypy,nfypxp,nfypx,nfyyp,nfyy,nfyxp,nfyx,nfxpyp,nfxpy,nfxpxp,nfxpx,nfxyp,nfxy,nfxxp,nfxx,hx,gx,gxx,eta);
else
    gxx = zeros(ny,nx,nx);
    hxx = zeros(nx,nx,nx);
    gss = zeros(ny,1);
    hss = zeros(nx,1);
end

if approx > 2
    [gxxx,hxxx] = gxxx_hxxx(nfx,nfxp,nfy,nfyp,nfypyp,nfypy,nfypxp,nfypx,nfyyp,nfyy,nfyxp,nfyx,nfxpyp,nfxpy,nfxpxp,nfxpx,nfxyp,nfxy,nfxxp,nfxx,...
        nfypypyp,nfypypy,nfypypxp,nfypypx,nfypyyp,nfypyy,nfypyxp,nfypyx,nfypxpyp,nfypxpy,nfypxpxp,nfypxpx,nfypxyp,nfypxy,nfypxxp,nfypxx,...
        nfyypyp,nfyypy,nfyypxp,nfyypx,nfyyyp,nfyyy,nfyyxp,nfyyx,nfyxpyp,nfyxpy,nfyxpxp,nfyxpx,nfyxyp,nfyxy,nfyxxp,nfyxx,...
        nfxpypyp,nfxpypy,nfxpypxp,nfxpypx,nfxpyyp,nfxpyy,nfxpyxp,nfxpyx,nfxpxpyp,nfxpxpy,nfxpxpxp,nfxpxpx,nfxpxyp,nfxpxy,nfxpxxp,nfxpxx,...
        nfxypyp,nfxypy,nfxypxp,nfxypx,nfxyyp,nfxyy,nfxyxp,nfxyx,nfxxpyp,nfxxpy,nfxxpxp,nfxxpx,nfxxyp,nfxxy,nfxxxp,nfxxx,...
        hx,gx,hxx,gxx);
    [gxss,hxss] = gxss_hxss(nfx,nfxp,nfy,nfyp,nfypyp,nfypy,nfypxp,nfypx,nfyyp,nfyy,nfyxp,nfyx,nfxpyp,nfxpy,nfxpxp,nfxpx,nfxyp,nfxy,nfxxp,nfxx,...
        nfypypyp,nfypypy,nfypypxp,nfypypx,nfypyyp,nfypyy,nfypyxp,nfypyx,nfypxpyp,nfypxpy,nfypxpxp,nfypxpx,nfypxyp,nfypxy,nfypxxp,nfypxx,...
        nfyypyp,nfyypy,nfyypxp,nfyypx,nfyyyp,nfyyy,nfyyxp,nfyyx,nfyxpyp,nfyxpy,nfyxpxp,nfyxpx,nfyxyp,nfyxy,nfyxxp,nfyxx,...
        nfxpypyp,nfxpypy,nfxpypxp,nfxpypx,nfxpyyp,nfxpyy,nfxpyxp,nfxpyx,nfxpxpyp,nfxpxpy,nfxpxpxp,nfxpxpx,nfxpxyp,nfxpxy,nfxpxxp,nfxpxx,...
        nfxypyp,nfxypy,nfxypxp,nfxypx,nfxyyp,nfxyy,nfxyxp,nfxyx,nfxxpyp,nfxxpy,nfxxpxp,nfxxpx,nfxxyp,nfxxy,nfxxxp,nfxxx,...
        hx,gx,hxx,gxx,gxxx,gss,hss,eta);
    [gsss,hsss] = gsss_hsss(nfx,nfxp,nfy,nfyp,nfypyp,nfypy,nfypxp,nfypx,nfyyp,nfyy,nfyxp,nfyx,nfxpyp,nfxpy,nfxpxp,nfxpx,nfxyp,nfxy,nfxxp,nfxx,...
        nfypypyp,nfypypy,nfypypxp,nfypypx,nfypyyp,nfypyy,nfypyxp,nfypyx,nfypxpyp,nfypxpy,nfypxpxp,nfypxpx,nfypxyp,nfypxy,nfypxxp,nfypxx,...
        nfyypyp,nfyypy,nfyypxp,nfyypx,nfyyyp,nfyyy,nfyyxp,nfyyx,nfyxpyp,nfyxpy,nfyxpxp,nfyxpx,nfyxyp,nfyxy,nfyxxp,nfyxx,...
        nfxpypyp,nfxpypy,nfxpypxp,nfxpypx,nfxpyyp,nfxpyy,nfxpyxp,nfxpyx,nfxpxpyp,nfxpxpy,nfxpxpxp,nfxpxpx,nfxpxyp,nfxpxy,nfxpxxp,nfxpxx,...
        nfxypyp,nfxypy,nfxypxp,nfxypx,nfxyyp,nfxyy,nfxyxp,nfxyx,nfxxpyp,nfxxpy,nfxxpxp,nfxxpx,nfxxyp,nfxxy,nfxxxp,nfxxx,...
        hx,gx,hxx,gxx,gxxx,gss,hss,esk);
else
    gxxx = zeros(ny,nx,nx,nx);
    hxxx = zeros(nx,nx,nx,nx);
    gxss = zeros(ny,nx);
    hxss = zeros(nx,nx);
    gsss = zeros(ny,1);
    hsss = zeros(nx,1);
end

% simulating and plotting moments 
sim_prune 

v=STD_EPS_A^2;
io = OMEGA^(-1);
%Computing Euler equation errors
T = T-1; 
ee1 = zeros(T,1);
ee2 = zeros(T,1);
ee3 = zeros(T,1);
% 2nd order k coef
kcoef = [hxx(3,:,1);hxx(3,:,2);hxx(3,:,3);hxx(3,:,4)];
% 2nd order c coef 
ccoef = [gxx(1,:,1);gxx(1,:,2);gxx(1,:,3);gxx(1,:,4)];
% 2nd order h coef 
hcoef = [gxx(2,:,1);gxx(2,:,2);gxx(2,:,3);gxx(2,:,4)];
disp('start')
for t=1:T
    vec = X_sim(:,t)-Xbar;
    % k prime
    ktp = k+hx(3,:)*vec+0.5*(vec'*kcoef*vec+hss(3)*v);
    % cons
    ct = c+gx(1,:)*vec+0.5*(vec'*ccoef*vec+gss(1)*v);
    %hours
    ht = h+gx(2,:)*vec+0.5*(vec'*hcoef*vec+gss(2)*v);
    % tomorrows states
    % we will replace tech to do quadrature 
    sp = ktp-k;
    % calculating output to get tomorrow's debt
    yt = A*K^ALFA*ht^(1-ALFA);
    ndebt = (1+rn)*D+ct+ktp-(1-DELTA)*K+PHI/2*(ktp-K)^2-yt;
    sp3 = ndebt - d;
    rf = r+PSSI*(exp(sp3)-1);
    % setting up quadrature 
    sp4 = rf-r;
    [x, w] = GaussHermite(ord);
    w = w/sum(w);
    x = sqrt(2)*STD_EPS_A*x+RHO*X_sim(4,t);
    sp2 = exp(x)-a;
    vec = [sp3*ones(1,ord);sp4*ones(1,ord);sp3*ones(ord);sp2];
    % k double prime
    kpp = k+hx(3,:)*vec+0.5*(vec'*kcoef*vec+hss(3)*v);
    % cons prime
    ctp = c+gx(1,:)*vec+0.5*(vec'*ccoef*vec+gss(1)*v);
    % hours prime
    htp = h+gx(2,:)*vec+0.5*(vec'*hcoef*vec+gss(2)*v);
    % c MU at t
    muc = (ct-io*ht^OMEGA)^(-SIGG);
    % c MU at t+1
    mucp =(ctp-io*htp^OMEGA)^(-SIGG);
    % h MU divided by c MU (at t)
    muh = ht^(OMEGA-1);
    % debt euler
    tempe = -muc+BETTA*(1+rf)*mucp;
    ee1(i) = ee1(i)+ tempe^2/NN;
    % intra euler
    tempe = -muh+(1-ALFA)*yt/ht;
    ee2(i) = ee2(i)+ tempe^2/NN;
    % capital euler
    delk = ktp-K;
    delkp = kpp-ktp;
    tempe = -muc*(1+PHI*delk)+BETTA*mucp*(ALFA*A^RHO*(htp/ktp)^(1-ALFA)+1-DELTA + PHI*delkp);
    ee3(i) = ee3(i)+ tempe^2/NN;
end 

for i = 1:N
    for j = 1:N
        for kk = 1:N
             %states
            K = kgrid(i);
            A = agrid(j);
            D = dgrid(kk);
            st = K-k;
            st2 = A-a;
            st3 = D-d;
            rn = r+PSSI*(exp(st3)-1);
            st4 = rn -r;
            % state vec
            vec = [st3; st4;st;st2];
            % k prime 
            ktp = k+hx(3,:)*vec+0.5*(vec'*kcoef*vec+hss(3)*v);
            % cons 
            ct = c+gx(1,:)*vec+0.5*(vec'*ccoef*vec+gss(1)*v);
            %hours 
            ht = h+gx(2,:)*vec+0.5*(vec'*hcoef*vec+gss(2)*v);
            % tomorrows states
            sp = ktp-k;
            sp2 = exp(log(A)*RHO)-a;
            % calculating output to get tomorrow's debt
            yt = A*K^ALFA*ht^(1-ALFA);
            ndebt = (1+rn)*D+ct+ktp-(1-DELTA)*K+PHI/2*(ktp-K)^2-yt;
            sp3 = ndebt - d;
            rf = r+PSSI*(exp(sp3)-1);
            sp4 = rf-r;
            vec = [sp3; sp4;sp;sp2];
            % k double prime
            kpp = k+hx(3,:)*vec+0.5*(vec'*kcoef*vec+hss(3)*v);
            % cons prime 
            ctp = c+gx(1,:)*vec+0.5*(vec'*ccoef*vec+gss(1)*v);
            % hours prime 
            htp = h+gx(2,:)*vec+0.5*(vec'*hcoef*vec+gss(2)*v);
            % c MU at t 
            muc = (ct-io*ht^OMEGA)^(-SIGG);
            % c MU at t+1 
            mucp =(ctp-io*htp^OMEGA)^(-SIGG);
            % h MU divided by c MU (at t)
            muh = ht^(OMEGA-1); 
            % debt euler 
            tempe = -muc+BETTA*(1+rf)*mucp; 
            ee1(i) = ee1(i)+ tempe^2/NN; 
            % intra euler
            tempe = -muh+(1-ALFA)*yt/ht; 
            ee2(i) = ee2(i)+ tempe^2/NN; 
            % capital euler
            delk = ktp-K;
            delkp = kpp-ktp;
            tempe = -muc*(1+PHI*delk)+BETTA*mucp*(ALFA*A^RHO*(htp/ktp)^(1-ALFA)+1-DELTA + PHI*delkp);
            ee3(i) = ee3(i)+ tempe^2/NN;  
        end 
    end 
    kgrid(i)
end


figure(1)
plot(kgrid,log10(abs(ee1)),'LineWidth',2)
figure(2)
plot(kgrid,log10(abs(ee2)),'LineWidth',2)
figure(3)
plot(kgrid,log10(abs(ee3)),'LineWidth',2)
