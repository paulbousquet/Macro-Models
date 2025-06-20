% Sets parameter values 
%
% Thomas Winberry,February 14th, 2018


%----------------------------------------------------------------
% Set parameters governing approximations
%----------------------------------------------------------------

global nProd nCapital nState prodMin prodMax capitalMin capitalMax nShocks nProdFine nCapitalFine nStateFine ...
	maxIterations tolerance acc dampening nMeasure nStateQuadrature nMeasureCoefficients nProdQuadrature ...
	nCapitalQuadrature kRepSS wRepSS

% Compute representative agent steady state (used in constructing the grids)
global kRepSS wRepSS
kRepSS 			= ((ttheta * (nSS ^ nnu)) / ((1 / bbeta) - (1 - ddelta))) ^ (1 / (1 - ttheta));
wRepSS 		= (kRepSS .^ ttheta) * nnu * (nSS ^ (nnu - 1));

% Order of approximation of value function
nProd 			= 3;										% order of polynomials in productivity
nCapital 		= 5;										% order of polynomials in capital
nState 			= nProd * nCapital;					% total number of coefficients

% Bounds on grid space
prodMin 		= -3 * ssigmaProd / sqrt(1 - rrhoProd ^ 2);
prodMax 		= 3 * ssigmaProd / sqrt(1 - rrhoProd ^ 2);
capitalMin 		= .1 * (exp(prodMin) ^ (1 / (1 - ttheta))) * kRepSS;
capitalMax 		= 2.5 * (exp(prodMax) ^ (1 / (1 - ttheta))) * kRepSS;

% Shocks 
nShocks 		= 3;										% order of Gauss-Hermite quadrature over idiosyncratic shocks

% Finer grid for analyzing policy functions and computing histogram
nProdFine 		= 60;
nCapitalFine 	= 40;
nStateFine 	= nProdFine * nCapitalFine;

% Iteration on value function
maxIterations	= 100;
tolerance 		= 1e-6;
acc 				= 500;									% number of iterations in "Howard improvement step"
dampening 	= 0;										% weight on old iteration in updating step

% Approximation of distribution
nMeasure 					= 2;							% order of polynomial approximating distribution
nProdQuadrature 		= 8; 							% number of quadrature points in productivity dimension
nCapitalQuadrature 	= 10;						% number of quadrature points in capital dimension
nStateQuadrature 		= nProdQuadrature * nCapitalQuadrature;
nMeasureCoefficients 	= (nMeasure * (nMeasure + 1)) / 2 + nMeasure;
