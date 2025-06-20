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
v=STD_EPS_A^2;
sim_prune 

io = OMEGA^(-1);
%Computing Euler equation errors 
% using quarature of order "ord" to compute expectation 
ord=3;
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
%t=2;
for t=1:T
    % linearized states
    vec = X_sim(:,t)-Xbar;
    % regular state vars
    [D, rn, K, A] = deal(X_sim(1,t), X_sim(2,t), X_sim(3,t), X_sim(4,t));
    % k prime
    ktp = k+hx(3,:)*vec+0.5*(vec'*kcoef*vec+hss(3)*v);
    % cons
    ct = c+gx(1,:)*vec+0.5*(vec'*ccoef*vec+gss(1)*v);
    %hours
    ht = h+gx(2,:)*vec+0.5*(vec'*hcoef*vec+gss(2)*v);
    % tomorrows states
    % we will replace tech to do quadrature 
    sp = ktp-k;
    delk = ktp-K;
    % calculating output to get tomorrow's debt
    yt = A*K^ALFA*ht^(1-ALFA);
    ndebt = (1+rn)*D+ct+ktp-(1-DELTAA)*K+.5*PHI*delk^2-yt;
    sp3 = ndebt - d;
    rf = r+PSSI*(exp(sp3)-1);
    sp4 = rf-r;
    % setting up quadrature 
    [x, w] = GaussHermite(ord);
    w = w/sum(w);
    % technology nodes based on quadrature change of variables
    Ap = exp(sqrt(2*v)*x+RHO*log(A));
    sp2 = Ap-a;
    vec = [sp3*ones(1,ord);sp4*ones(1,ord);sp3*ones(1,ord);sp2'];
    % getting zero+first order coeficient to streamline indexing
    hxk = k+vec.'*hx(3,:).';
    gxc = c+vec.'*gx(1,:).';
    gxh = h+vec.'*gx(2,:).';
    % vector of k double prime, one for element for each tech node 
    kpp = hxk+0.5*(arrayfun(@(n) vec(:,n)' * kcoef * vec(:,n), 1:size(vec,2)).'+hss(3)*v);
    % vector of cons prime, one for element for each tech node
    ctp = gxc+0.5*(arrayfun(@(n) vec(:,n)' * ccoef * vec(:,n), 1:size(vec,2)).'+gss(1)*v);
    % vector of hours prime, one for element for each tech node
    htp = gxh+0.5*(arrayfun(@(n) vec(:,n)' * hcoef * vec(:,n), 1:size(vec,2)).'+gss(2)*v);
    % c MU at t
    muc = (ct-io*ht^OMEGA)^(-SIGG);
    % c MU at t+1
    mucp =(ctp-io*htp.^OMEGA).^(-SIGG);
    % h MU divided by c MU (at t)
    muh = ht^(OMEGA-1);
    % debt euler
    tempe = -muc+BETTA*(1+rf).*mucp;
    %raw squared error 
    %ee1(t) = sum(tempe.*w)^2
    % error normalized by magnitude (av of LHS and RHS)  
    ee1(t) = sum(tempe.*w)^2/(sum(tempe/2+muc)/ord);
    % intra euler
    tempe = -muh+(1-ALFA)*yt/ht;
    % errors here are very low, no normalization needed 
    ee2(t) = tempe^2;
    % capital euler
    delkp = kpp-ktp;
    tempe = -muc*(1+PHI*delk)+BETTA.*mucp.*(ALFA*Ap.*(htp./ktp).^(1-ALFA)+1-DELTAA + PHI*delkp);
    %raw squared error 
    %ee3(t) = sum(tempe.*w)^2;
    % error normalized by magnitude (av of LHS and RHS)  
    ee3(t) = sum(tempe.*w)^2/(sum(tempe/2+muc*(1+PHI*delk))/ord);
end 

figure(2)
histogram(ee1)
figure(3)
histogram(ee2)
figure(4)
histogram(ee3)

