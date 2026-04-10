% Compute polynomials for approximating various objects
%
% Thomas Winberry, February 14th, 2018

%---------------------------------------------------------------
% Polynomials for value function and capital accumulation conditional on adjusting
%		(Chebyshev polynomials evaluated over collocation nodes from computeGrids.m)
%---------------------------------------------------------------

global mStatePoly vStatePolySquared aProdPrimePoly

% Create one-dimensional polynomials
mProdPoly 		= computeChebyshev(nProd,mStateGridZeros(:,1));
mCapitalPoly 	= computeChebyshev(nCapital,mStateGridZeros(:,2));

% Create tensor product in polynomials (vectorized)
mStatePoly = repmat(mProdPoly, 1, nCapital) .* repelem(mCapitalPoly, 1, nProd);
clear mProdPoly mCapitalPoly

% Compute squared terms for interpolation formulas
mProdPoly 				= computeChebyshev(nProd,vProdZeros);
mCapitalPoly 			= computeChebyshev(nCapital,vCapitalZeros);
vProd 					= sum(mProdPoly .^ 2)';
vCapital 				= sum(mCapitalPoly .^ 2)';
vStatePolySquared 	= reshape(vProd * vCapital',nState,1);
clear mProdPoly mCapitalPoly vProd vCapital

% Compute polynomials over future productivity shocks (for computing expected value function)
aProdPrimePoly 		= reshape(computeChebyshev(nProd,reshape(mProdPrimeZeros,nShocks * nState,1)),nShocks,nState,nProd);


%---------------------------------------------------------------
% Chebyshev polynomials over fine grid
%		(useful for plotting and computing law of motion for histogram)
%---------------------------------------------------------------

global mFinePoly aProdPrimeFinePoly

% Create one-dimensional polynomials
mProdFinePoly 		= computeChebyshev(nProd,mFineGridZeros(:,1));
mCapitalFinePoly 	= computeChebyshev(nCapital,mFineGridZeros(:,2));

% Compute tensor product of polynomials (vectorized)
mFinePoly = repmat(mProdFinePoly, 1, nCapital) .* repelem(mCapitalFinePoly, 1, nProd);

% Compute polynomials over future shocks
aProdPrimeFinePoly 	= reshape(computeChebyshev(nProd,reshape(mProdPrimeFineZeros,nShocks * nStateFine,1)),nShocks,...
									nStateFine,nProd);
clear mProdFinePoly mCapitalFinePoly


%---------------------------------------------------------------
% Chebyshev polynomials over quadrature grid 
%		(useful for compting law of motion for parametric family)
%---------------------------------------------------------------

global mQuadraturePoly aProdPrimeQuadraturePoly

% Create individual polynomials
mProdQuadraturePoly 		= computeChebyshev(nProd,mQuadratureGridZeros(:,1));
mCapitalQuadraturePoly 	= computeChebyshev(nCapital,mQuadratureGridZeros(:,2));

% Compute tensor product of polynomials (vectorized)
mQuadraturePoly = repmat(mProdQuadraturePoly, 1, nCapital) .* repelem(mCapitalQuadraturePoly, 1, nProd);

% Compute polynomials over future shocks
aProdPrimeQuadraturePoly = reshape(computeChebyshev(nProd,reshape(mProdPrimeQuadratureZeros,nShocks * ...
										nStateQuadrature,1)),nShocks,nStateQuadrature,nProd);

	
%---------------------------------------------------------------
% Derivative of Chebyshev polynomials over collocation nodes
%		(used to compute first-order condition in dynamic model)	
%---------------------------------------------------------------

global mStatePoly2

% Create individual polynomials
mProdPoly 			= computeChebyshev(nProd,mStateGridZeros(:,1));
mCapitalPoly2 	= computeChebyshev2(nCapital-1,mStateGridZeros(:,2));

% Create tensor product in polynomials (vectorized)
mStatePoly2 = repmat(mProdPoly, 1, nCapital-1) .* repelem(mCapitalPoly2, 1, nProd);
clear mCapitalPoly2


%---------------------------------------------------------------
% Derivative of Chebyshev polynomials over fine grid
%		(used to plot marginal value function)	
%---------------------------------------------------------------

global mFinePoly2

% Create individual polynomials
mProdPoly 			= computeChebyshev(nProd,mFineGridZeros(:,1));
mCapitalPoly2 	= computeChebyshev2(nCapital-1,mFineGridZeros(:,2));

% Create tensor product in polynomials (vectorized)
mFinePoly2 = repmat(mProdPoly, 1, nCapital-1) .* repelem(mCapitalPoly2, 1, nProd);
clear mCapitalPoly2
