% Some indices
ny      = size(gx,1);
nx      = size(hx,1);

% Allocating memory
T = 100000;
Y_sim   = zeros(ny,T);
X_sim   = zeros(nx,T+1);
xf_sim  = zeros(nx,T);
xs_sim  = zeros(nx,T);
xrd_sim = zeros(nx,T);

% Defining matrices
HHxxtil  = 1/2*reshape(hxx,nx,nx^2);
HHxxxtil = 1/6*reshape(hxxx,nx,nx^3);
GGxxtil  = 1/2*reshape(gxx,ny,nx^2);
GGxxxtil = 1/6*reshape(gxxx,ny,nx^3);

Ybar = [c;h;k];
Xbar = [d;r;k;a];

shocks = zeros(length(eta),T+1);
shocks(4,:) = randn(T+1,1);

X_sim(:,1) = Xbar;
xf_sim(:,1) = Xbar;
xs_sim(:,1) = Xbar;
xrd_sim(:,1) = Xbar;
for t=1:T
    xf = xf_sim(:,t)-Xbar;
    xs = xs_sim(:,t)-Xbar;
    xrd = xrd_sim(:,t)-Xbar;
    AA = kron(xf,xf);
    BB = kron(xf,xs);
    Y_sim(:,t) = Ybar+gx*(xf+xs+xrd)+GGxxtil*(AA+2*BB)+0.5*gss+GGxxxtil*kron(xf,AA);
    xf_sim(:,t+1) = Xbar+hx*xf+eta'*shocks(:,t+1);
    xs_sim(:,t+1) = Xbar+hx*xs+HHxxtil*AA+0.5*hss;
    xrd_sim(:,t+1) = Xbar+hx*xrd+2*HHxxtil*BB+HHxxxtil*kron(xf,AA);
    X_sim(:,t+1) = Xbar+(xf_sim(:,t+1)-Xbar)+(xs_sim(:,t+1)-Xbar)+(xrd_sim(:,t+1)-Xbar);
end


tx = ["debt", "IR", "capital", "tech"];
ty = ["cons", "hours"];
mux = round(mean(X_sim,2),2)';
muy = round(mean(Y_sim(1:2,:),2),2)';
stx = round(std(X_sim, 0, 2),2)';
sty = round(std(Y_sim(1:2,:), 0, 2),2)';
titx = arrayfun(@(y,z) sprintf('$\\mu=%.2f$, $\\sigma=%.2f$', y, z), mux, stx, 'UniformOutput', false);
tity = arrayfun(@(y,z) sprintf('$\\mu=%.2f$, $\\sigma=%.2f$', y, z), muy, sty, 'UniformOutput', false);

figure;

% superfluous code to print moments in legend 
% otherwise titles too cumbersome. thanks chatGPT

for i = 1:4
    subplot(2, 3, i); 
    histogram(X_sim(i, :)); 
    title(tx{i});
    hold on; % Keep the current plot, so we can add a dummy plot for the legend
    h = plot(NaN,NaN,'w'); % Invisible plot
    legend(h, titx{i}, 'Interpreter', 'latex');
    hold off; 
end

for j = 1:2
    subplot(2, 3, j + 4); 
    histogram(Y_sim(j, :)); 
    title(ty{j});
    hold on; % Keep the current plot
    h = plot(NaN,NaN,'w'); % Invisible plot for the legend
    legend(h, tity{j}, 'Interpreter', 'latex');
    hold off; % Release the plot hold
end

% more sims 

ne      = size(eta,2);
num_sim = size(shocks,2);

% Allocating memory
y_sim   = zeros(ny,num_sim);
xf_sim  = zeros(nx,num_sim);
xs_sim  = zeros(nx,num_sim);
xrd_sim = zeros(nx,num_sim);

% The simulation
xf    = zeros(nx,1);
xs    = zeros(nx,1);
xrd   = zeros(nx,1);


for t=1:num_sim
        % The controls
        xf_sim(:,t) = xf;
        xs_sim(:,t) = xs;
        AA          = kron(xf,xf);
        y_sim(:,t)  = gx*(xf+xs) + GGxxtil*AA + 1/2*gss;
        
        % The state variables
        xf_p  = hx*xf + eta'*shocks(:,t);
        xs_p  = hx*xs + HHxxtil*AA + 1/2*hss;
        
        % Updating xf, xs
        xf  = xf_p;
        xs  = xs_p;
end

owmega = inv(chol(cov(X_sim([1,3,4],:)'))); 

