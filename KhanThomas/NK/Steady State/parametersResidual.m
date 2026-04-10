function [value vDerivs] = parametersResidual(vParameters,mGridMoments);

% Computes the objective function for the minimzation problem for approximating the distribution
% 
% Inputs
% (1) vParameters: candidate parameters
%	(2) mGridMoments: parameters of distribution
%
% Outputs
% (1) value: value of the objective function
%	(2) (optional) vDerivs: Jacobian of the objective function (to help with numerical minimizer)
% 
% Thomas Winberry, February 14th, 2018

global vQuadratureWeights nMeasureCoefficients

% Cache exp() evaluation
vExpTerm = exp(mGridMoments * vParameters);

% Value of function
value = vQuadratureWeights' * vExpTerm;

% Derivatives
if nargout > 1
	vDerivs = mGridMoments' * (vQuadratureWeights .* vExpTerm);
end