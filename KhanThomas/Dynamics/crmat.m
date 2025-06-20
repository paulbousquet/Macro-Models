% Step 1: Run setParameters.m and save its outputs
clearvars
setParameters
vars_before = who;
save('economicParameters.mat');

% Step 2: Run setGrids.m, keeping parameters in memory
vars_before = who;    % record existing vars (from setParameters)
computeApprox
vars_after = who;
new_vars = setdiff(vars_after, vars_before);
save('approximationParameters.mat', new_vars{:});


% Step 2: Run setGrids.m, keeping parameters in memory
%vars_before = who;    % record existing vars (from setParameters)
computeGrids
vars_after = who;
new_vars = setdiff(vars_after, vars_before);
save('grids.mat', new_vars{:});

% Step 3: Run setShocks.m, again keeping prior vars
vars_before = who;
computePolynomials
vars_after = who;
new_vars = setdiff(vars_after, vars_before);
save('polynomials.mat', new_vars{:});

% Step 4: Same for steady state
