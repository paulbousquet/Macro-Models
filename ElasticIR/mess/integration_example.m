%Gauss-Hermite quadrature (normal random variable)
ex1 = zeros(15,1); ex2 = zeros(15,1); ex3 = zeros(15,1); ex4 = zeros(15,1);
sigma = 2; mu = 1;
for n=1:15
    [x, w] = GaussHermite(n);
    w = w/sum(w);
    x = sqrt(2)*sigma*x+mu;
    ex1(n) = sum(x.*w);
    ex2(n) = sum((x-ex1(n)).*(x-ex1(n)).*w);
    ex3(n) = sum((x-ex1(n)).*(x-ex1(n)).*(x-ex1(n)).*w);
    ex4(n) = sum((x-ex1(n)).*(x-ex1(n)).*(x-ex1(n)).*(x-ex1(n)).*w);
end
figure(1)
subplot(2,2,1),plot([1:1:15],ex1,'LineWidth',2)
hold
subplot(2,2,1),plot([1:1:15],ones(15,1)*mu,'LineWidth',2)
subplot(2,2,2),plot([1:1:15],ex2,'LineWidth',2)
hold
subplot(2,2,2),plot([1:1:15],ones(15,1)*sigma^2,'LineWidth',2)
subplot(2,2,3),plot([1:1:15],ex3/sigma^3,'LineWidth',2)
hold
subplot(2,2,3),plot([1:1:15],zeros(15,1),'LineWidth',2)
subplot(2,2,4),plot([1:1:15],ex4/sigma^4-3,'LineWidth',2)
hold
subplot(2,2,4),plot([1:1:15],zeros(15,1),'LineWidth',2)

%Gauss-Hermite quadrature (log-normal random variable)
ex1 = zeros(35,1); ex2 = zeros(35,1); ex3 = zeros(35,1); ex4 = zeros(35,1);
sigma = 2; mu = 1;
for n=1:35
    [x, w] = GaussHermite(30);
    w = w/sum(w);
    x = sqrt(2)*STD_EPS_A*x+mu;
    x = exp(x);
    ex1(n) = sum(x.*w);
    ex2(n) = sum((x-ex1(n)).*(x-ex1(n)).*w);
    ex3(n) = sum((x-ex1(n)).*(x-ex1(n)).*(x-ex1(n)).*w);
    ex4(n) = sum((x-ex1(n)).*(x-ex1(n)).*(x-ex1(n)).*(x-ex1(n)).*w);
end
figure(2)
subplot(2,2,1),plot([1:1:35],ex1,'LineWidth',2)
hold
subplot(2,2,1),plot([1:1:35],ones(35,1)*exp(mu+sigma^2/2),'LineWidth',2)
subplot(2,2,2),plot([1:1:35],ex2,'LineWidth',2)
hold
subplot(2,2,2),plot([1:1:35],ones(35,1)*(exp(sigma^2)-1)*exp(2*mu+sigma^2),'LineWidth',2)
subplot(2,2,3),plot([1:1:35],ex3./ex2.^(3/2),'LineWidth',2)
hold
subplot(2,2,3),plot([1:1:35],ones(35,1)*(exp(sigma^2)+2)*sqrt(exp(sigma^2)-1),'LineWidth',2)
subplot(2,2,4),plot([1:1:35],ex4./ex2.^2-3,'LineWidth',2)
hold
subplot(2,2,4),plot([1:1:35],ones(35,1)*(exp(4*sigma^2)+2*exp(3*sigma^2)+3*exp(2*sigma^2)-6),'LineWidth',2)

%Monte Carlo
T = [100 1000 10000 100000];
ex1 = zeros(4,1); ex2 = zeros(4,1); ex3 = zeros(4,1); ex4 = zeros(4,1);
for t=1:4
    eps = sigma*randn(T(t),1)+mu;
    ex1(t) = sum(eps)/T(t);
    ex2(t) = sum((eps-ex1(t)).*(eps-ex1(t)))/T(t);
    ex3(t) = sum((eps-ex1(t)).*(eps-ex1(t)).*(eps-ex1(t)))/T(t);
    ex4(t) = sum((eps-ex1(t)).*(eps-ex1(t)).*(eps-ex1(t)).*(eps-ex1(t)))/T(t);
end
figure(3)
subplot(2,2,1),plot(log10(T),ex1,'LineWidth',2)
hold
subplot(2,2,1),plot(log10(T),ones(4,1)*mu,'LineWidth',2)
subplot(2,2,2),plot(log10(T),ex2,'LineWidth',2)
hold
subplot(2,2,2),plot(log10(T),ones(4,1)*sigma^2,'LineWidth',2)
subplot(2,2,3),plot(log10(T),ex3/sigma^3,'LineWidth',2)
hold
subplot(2,2,3),plot(log10(T),zeros(4,1),'LineWidth',2)
subplot(2,2,4),plot(log10(T),ex4/sigma^4-3,'LineWidth',2)
hold
subplot(2,2,4),plot(log10(T),zeros(4,1),'LineWidth',2)

%Monte Carlo with antithetic variates
T = [100 1000 10000 100000];
ex1 = zeros(4,1); ex2 = zeros(4,1); ex3 = zeros(4,1); ex4 = zeros(4,1);
for t=1:4
    eps1 = sigma*randn(T(t),1)+mu;
    eps2 = 2/T(t)*sum(eps1)-eps1;
    ex1(t) = sum(eps1)/T(t)/2+sum(eps2)/T(t)/2;
    ex2(t) = sum((eps1-ex1(t)).*(eps1-ex1(t)))/T(t)/2+...
        sum((eps2-ex1(t)).*(eps2-ex1(t)))/T(t)/2;
    ex3(t) = sum((eps1-ex1(t)).*(eps1-ex1(t)).*(eps1-ex1(t)))/T(t)/2+...
        sum((eps2-ex1(t)).*(eps2-ex1(t)).*(eps2-ex1(t)))/T(t)/2;
    ex4(t) = sum((eps1-ex1(t)).*(eps1-ex1(t)).*(eps1-ex1(t)).*(eps1-ex1(t)))/T(t)/2+...
        sum((eps2-ex1(t)).*(eps2-ex1(t)).*(eps2-ex1(t)).*(eps2-ex1(t)))/T(t)/2;
end
figure(4)
subplot(2,2,1),plot(log10(T),ex1,'LineWidth',2)
hold
subplot(2,2,1),plot(log10(T),ones(4,1)*mu,'LineWidth',2)
subplot(2,2,2),plot(log10(T),ex2,'LineWidth',2)
hold
subplot(2,2,2),plot(log10(T),ones(4,1)*sigma^2,'LineWidth',2)
subplot(2,2,3),plot(log10(T),ex3/sigma^3,'LineWidth',2)
hold
subplot(2,2,3),plot(log10(T),zeros(4,1),'LineWidth',2)
subplot(2,2,4),plot(log10(T),ex4/sigma^4-3,'LineWidth',2)
hold
subplot(2,2,4),plot(log10(T),zeros(4,1),'LineWidth',2)

%adaptive Romberg integration
f = @(x) x.*exp(-(x-mu).^2/2/sigma^2)/sqrt(2*pi)/sigma;
Q1 = integral(f,-Inf,Inf)
f = @(x) (x-Q1).^2.*exp(-(x-mu).^2/2/sigma^2)/sqrt(2*pi)/sigma;
Q2 = integral(f,-Inf,Inf)
f = @(x) (x-Q1).^3.*exp(-(x-mu).^2/2/sigma^2)/sqrt(2*pi)/sigma;
Q3 = integral(f,-Inf,Inf)/sigma^3
f = @(x) (x-Q1).^4.*exp(-(x-mu).^2/2/sigma^2)/sqrt(2*pi)/sigma;
Q4 = integral(f,-Inf,Inf)/sigma^4-3