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

approx = 2;

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

v=STD_EPS_A^2;
io = OMEGA^(-1);
N = 101;
NN = N^2;
%Computing Euler equation errors
kgrid = linspace(0.5*k,1.5*k,N);
agrid = linspace(0.5*a,1.5*a,N);
dgrid = linspace(0.5*d,1.5*d,N);
ee1 = zeros(N,1);
ee2 = zeros(N,1);
ee3 = zeros(N,1);
disp('start')
for i = 1:N
    for j = 1:N
        for kk = 1:N
            %states
            K = kgrid(i);
            A = agrid(j);
            D = dgrid(kk);
            st = kgrid(i)-k;
            st2 = agrid(j)-a;
            st3 = dgrid(kk)-d;
            rn = r+PSSI*(exp(st3-d));
            st4 = rn -r;
            % state vec
            vec = [st3; st4;st;st2];
            % 2nd order k coef
            kcoef = [hxx(3,:,1);hxx(3,:,2);hxx(3,:,3);hxx(3,:,4)];
            % k prime 
            ktp = k+hx(3,:)*vec+0.5*(vec'*kcoef*vec+hss(3)*v);
            % 2nd order c coef 
            ccoef = [gxx(1,:,1);gxx(1,:,2);gxx(1,:,3);gxx(1,:,4)];
            % cons 
            ct = c+gx(1,:)*vec+0.5*(vec'*ccoef*vec+gss(1)*v);
            % tomorrows states
            sp = ktp-k;
            sp2 = A*RHO-a;
            % calculating output to get tomorrow's debt 
            hcoef = [gxx(4,:,1);gxx(4,:,2);gxx(4,:,3);gxx(4,:,4)];
            ht = h+gx(4,:)*vec+0.5*(vec'*hcoef*vec+gss(4)*v);
            yt = A*K^ALFA*ht*(1-ALFA);
            ndebt = (1+rn)*D+ct+ktp-(1-DELTA)*K+PHI/2*(ktp-K)^2-yt;
            sp3 = ndebt - d;
            rf = r+PSSI*(exp(sp3-d));
            sp4 = rf-r;
            vec = [sp3; sp4;sp;sp2];
            % k double prime
            kpp = k+hx(3,:)*vec+0.5*(vec'*kcoef*vec+hss(3)*v);
            % cons prime 
            ctp = c+gx(1,:)*vec+0.5*(vec'*ccoef*vec+gss(1)*v);
            % hours prime 
            htp = h+gx(4,:)*vec+0.5*(vec'*hcoef*vec+gss(4)*v);
            % c MU at t 
            muc = (ct-io*ht*OMEGA)^(-SIGG);
            % c MU at t+1 
            mucp =(ctp-io*htp*OMEGA)^(-SIGG);
            % h MU divided by c MU (at t)
            muh = -h*(OMEGA-1); 
            % debt euler 
            tempe = -muc+BETTA*(1+rf)*mucp; 
            ee1(i) = ee1(i)+ tempe^2/NN; 
            % intra euler
            tempe = -muh+(1-ALFA)*yt/ht; 
            ee2(i) = ee2(i)+ tempe^2/NN; 
            % capital euler
            delk = ktp-K;
            delkp = kpp-ktp;
            tempe = -muc*(1+delk)+BETTA*mucp*(A*RHO*(htp/ktp)*(1-ALFA)+1-DELTA + delkp);
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
