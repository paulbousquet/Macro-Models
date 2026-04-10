function [residual,vMoments,vParameters,aggregateConsumption,marginalUtility,aggregateOutput,aggregateCapital,aggregateInvestment] = ...
	computeLMCResidualPolynomials(wage,vParameters,vMoments,mGridMoments);

% Computes the residual of the labor market clearing condition, approximating the distribution
% using exponential polynomials
% 
% Inputs
%   (1) wage: candidate wage
%	(2) vParmaters: initial guess of distribution parameters
%	(3) vMoments: initial guess of distribution moments
%	(4) mGridMoments: grid of centralized moments, corresponding to vMoments
%
% Outputs
%   (1) resid: aggregate labor demand - supply
%	(2) vMoments: moments of parametric family (optional)
%	(3) vParameters: parameters of parametric family (optional)
%	(4) aggregateConsumption (optional)
%	(5) marginalUtility (optional)
%	(6) aggregateOutput (optional)
%	(7) aggregateCapital (optional)
%	(8) aggregateInvestment (optional)
% 
% Thomas Winberry, September 8th, 2016

% Declare global variables
global mStateGrid ttheta nnu ddelta mStatePoly vStatePolySquared tolerance maxIterations nState mFinePoly mFineGrid ...
	aaUpper aaLower capitalMin capitalMax nShocks aProdPrimeFinePoly nStateFine nProd nCapital bbeta ppsi2Capital dampening ...
	nProdFine nCapitalFine vShocksWeights ppsiCapital nProdQuadrature nCapitalQuadrature nStateQuadrature mQuadratureGrid ...
	mQuadraturePoly aProdPrimeQuadraturePoly nMeasure nMeasureCoefficients mProdPrimeQuadrature vQuadratureWeights nSS ssigma

% Persistent cache: return cached results if called with same wage and full output requested
persistent cachedPolyWage cachedPolyResults
if ~isempty(cachedPolyWage) && cachedPolyWage == wage && nargout > 1
	residual = cachedPolyResults.residual;
	vMoments = cachedPolyResults.vMoments;
	vParameters = cachedPolyResults.vParameters;
	aggregateConsumption = cachedPolyResults.aggregateConsumption;
	marginalUtility = cachedPolyResults.marginalUtility;
	aggregateOutput = cachedPolyResults.aggregateOutput;
	aggregateCapital = cachedPolyResults.aggregateCapital;
	aggregateInvestment = cachedPolyResults.aggregateInvestment;
	return;
end


%----------------------------------------------------------------
% Compute value function
%----------------------------------------------------------------

global vLaborDemandGrid vProfitGrid
vLaborDemandGrid = ((exp(mStateGrid(:,1)) .* (mStateGrid(:,2) .^ ttheta) * nnu) / wage) .^ (1 / (1 - nnu));
vProfitGrid = exp(mStateGrid(:,1)) .* (mStateGrid(:,2) .^ ttheta) .* (vLaborDemandGrid .^ nnu) - wage * vLaborDemandGrid;

% Initialize value function (warm-start from previous call if available)
persistent vCoefficientsWarmStart_Polynomials
if isempty(vCoefficientsWarmStart_Polynomials)
	vInitGrid = vProfitGrid + (1 - ddelta) * mStateGrid(:,2);
	vCoefficients = sum(mStatePoly' .* (ones(nState,1) * vInitGrid'),2);
	vCoefficients = vCoefficients ./ vStatePolySquared;
else
	vCoefficients = vCoefficientsWarmStart_Polynomials;
end
err = 100; iterations = 1;
t0 = tic;

% Iterate
while err > tolerance && iterations <= maxIterations

	[vCoefficientsNew,vCapitalAdjust] = updateCoefficients(vCoefficients);
	err = max(abs(vCoefficientsNew - vCoefficients));
	iterations = iterations + 1;
	vCoefficients = dampening * vCoefficients + (1 - dampening) * vCoefficientsNew;

end

% Save converged coefficients for warm-starting next call
vCoefficientsWarmStart_Polynomials = vCoefficients;


%----------------------------------------------------------------
% Compute value and policy functions over quadrature grid
%----------------------------------------------------------------

% Compute polynomial approximation of investment policy conditional on adjustment
vCapitalCoefficients = sum(mStatePoly' .* (ones(nState,1) * vCapitalAdjust'),2);
vCapitalCoefficients = vCapitalCoefficients ./ vStatePolySquared;

% Compute policy functions
[vCapitalAdjust,vCapitalConstrained,vCutoff] = computePolicies(vCoefficients,vCapitalCoefficients,wage,mQuadratureGrid,mQuadraturePoly,aProdPrimeQuadraturePoly);

% Enforce choices to be within bounds
vCapitalAdjust = min(max(capitalMin * ones(nStateQuadrature,1),vCapitalAdjust),capitalMax * ones(nStateQuadrature,1));
vCapitalConstrained = min(max(capitalMin * ones(nStateQuadrature,1),vCapitalConstrained),capitalMax * ones(nStateQuadrature,1));

% Replicate choices for each draw of future idioysncratic shock
vCapitalAdjustPrime = repmat(vCapitalAdjust', [nShocks 1]);
vCapitalConstrainedPrime = repmat(vCapitalConstrained', [nShocks 1]);
vCutoffPrime = repmat(vCutoff', [nShocks 1]);


%----------------------------------------------------------------
% Compute stationary distribution by iterating on law of motion
%----------------------------------------------------------------

% Initialize objects for the iteration
err = 100; iteration = 1; maxIterations = 200; dampening = 0;

% Iterate
while err > 1e-6 && iteration <= maxIterations

	% Compute parameters of distribution using Newton steps
	% Minimizes f(p) = w' * exp(M*p) which is strictly convex
	vParams = vParameters(2:nMeasureCoefficients+1,1);
	for iNewton = 1:20
		vExpTerm = exp(mGridMoments * vParams);
		vWeightedExp = vQuadratureWeights .* vExpTerm;
		g = mGridMoments' * vWeightedExp;                       % gradient (nMeasureCoefficients x 1)
		H = (mGridMoments .* vWeightedExp)' * mGridMoments;     % Hessian (nMeasureCoefficients x nMeasureCoefficients)
		vStep = H \ g;
		vParams = vParams - vStep;
		if norm(vStep) < 1e-10
			break;
		end
	end
	normalization = vQuadratureWeights' * exp(mGridMoments * vParams);
	vParameters = [1 / normalization; vParams];

	% Cache density evaluation (used multiple times below)
	vDensity = vParameters(1,1) .* exp(mGridMoments * vParameters(2:nMeasureCoefficients+1));

	% Compute new moments and centered moments grid
	% First moments
	vMomentsNew = zeros(2,1);
	for i = 0:1
		mIntegrand 				= (mProdPrimeQuadrature .^ (1 - i)) .* ((vCutoffPrime ./ ppsiCapital) .* (log(vCapitalAdjustPrime) .^ i) + ...
											(1 - (vCutoffPrime ./ ppsiCapital)) .* (log(vCapitalConstrainedPrime) .^ i));
		vIntegrand 				= (vShocksWeights' * mIntegrand)';
		vMomentsNew(i+1,1) 	= vQuadratureWeights' * (vIntegrand .* vDensity);
	end
	mGridMomentsNew 			= [mQuadratureGrid(:,1) - vMomentsNew(1,1) log(mQuadratureGrid(:,2)) - vMomentsNew(2,1)];

	% Higher-order moments
	for i = 2:nMeasure
		for j = 0:i
			mIntegrand 			= (mProdPrimeQuadrature .^ (i - j)) .* ((vCutoffPrime ./ ppsiCapital) .* ((log(vCapitalAdjustPrime) - vMomentsNew(2,1)) .^ j) + ...
											(1 - (vCutoffPrime ./ ppsiCapital)) .* ((log(vCapitalConstrainedPrime) - vMomentsNew(2,1)) .^ j));
			vIntegrand 			= (vShocksWeights' * mIntegrand)';
			vMomentsNew 		= [vMomentsNew; vQuadratureWeights' * (vIntegrand .* vDensity)];
			mGridMomentsNew 	= [mGridMomentsNew (((mQuadratureGrid(:,1) - vMomentsNew(1,1)) .^ (i - j)) .* ((log(mQuadratureGrid(:,2)) - ...
											vMomentsNew(2,1)) .^ j) - vQuadratureWeights' * (vIntegrand .* vDensity))];
		end
	end

	% Update iteration
	err 					= max(abs(vMomentsNew - vMoments));
	iteration 			= iteration + 1;
	vMoments 			= (1 - dampening) * vMomentsNew + dampening * vMoments;
	mGridMoments 	= (1 - dampening) * mGridMomentsNew + dampening * mGridMoments;

end


%----------------------------------------------------------------
% Compute labor market clearing residual
%----------------------------------------------------------------

% Cache density for aggregates (recompute after distribution iteration updated mGridMoments)
vDensity = vParameters(1,1) .* exp(mGridMoments * vParameters(2:nMeasureCoefficients+1));

% Compute labor demand along grid
vLaborDemand 		= ((exp(mQuadratureGrid(:,1)) .* (mQuadratureGrid(:,2) .^ ttheta) * nnu) / wage) .^ (1 / (1 - nnu)) + ...
								((vCutoff .^ 2) ./ (2 * ppsiCapital));

% Integrate
aggLaborDemand 	= vQuadratureWeights' * (vLaborDemand .* vDensity);
residual 				= aggLaborDemand - nSS;


%----------------------------------------------------------------
% Optional: compute steady state aggregates
%----------------------------------------------------------------

% Aggregate consumption (always computed for caching; cheap on 80-element vectors)
vLaborDemand 			= ((exp(mQuadratureGrid(:,1)) .* (mQuadratureGrid(:,2) .^ ttheta) * nnu) / wage) .^ (1 / (1 - nnu));
vIntegrand 				= exp(mQuadratureGrid(:,1)) .* (mQuadratureGrid(:,2) .^ ttheta) .* (vLaborDemand .^ nnu) + ...
									(vCutoff ./ ppsiCapital) .* (-(vCapitalAdjust - (1 - ddelta) * mQuadratureGrid(:,2))) + (1 - (vCutoff ./ ppsiCapital)) .* ...
									(-(vCapitalConstrained - (1 - ddelta) * mQuadratureGrid(:,2)));
aggregateConsumption = vQuadratureWeights' * (vIntegrand .* vDensity);

% Marginal utility
marginalUtility 			= aggregateConsumption ^ (-ssigma);

% Output
vIntegrand 				= exp(mQuadratureGrid(:,1)) .* (mQuadratureGrid(:,2) .^ ttheta) .* (vLaborDemand .^ nnu);
aggregateOutput 		= vQuadratureWeights' * (vIntegrand .* vDensity);

% Capital
aggregateCapital 		= vQuadratureWeights' * (mQuadratureGrid(:,2) .* vDensity);

% Investment
vIntegrand 				= (vCutoff ./ ppsiCapital) .* (vCapitalAdjust - (1 - ddelta) * mQuadratureGrid(:,2)) + (1 - (vCutoff ./ ppsiCapital)) .* ...
									(vCapitalConstrained - (1 - ddelta) * mQuadratureGrid(:,2));
aggregateInvestment 	= vQuadratureWeights' * (vIntegrand .* vDensity);

% Cache results for subsequent calls with the same wage
cachedPolyWage = wage;
cachedPolyResults = struct('residual',residual,'vMoments',vMoments,'vParameters',vParameters,...
	'aggregateConsumption',aggregateConsumption,'marginalUtility',marginalUtility,...
	'aggregateOutput',aggregateOutput,'aggregateCapital',aggregateCapital,...
	'aggregateInvestment',aggregateInvestment);
