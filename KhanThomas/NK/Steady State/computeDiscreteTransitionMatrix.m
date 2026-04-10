function mTransition = computeDiscreteTransitionMatrix(vCapitalAdjust,vCapitalConstrained,vCutoff)

% Computes discrete transition matrix associated with decision rules following Young (2010)
%
% Inputs
% (1) vCapitalAdjust: investment decision, conditional on adjusting (over fine grid)
%	(2) vCapitalConstrained: investment decision, conditional on not adjusting (over fine grid)
%	(3) vCutoff: fixed cost cutoff for adjusting capital (over fine grid)
%
% Outputs
% (1) mTransition: transition matrix
%
% Thomas Winberry, Feburary 14th, 2018

% Declare global variables used in this function
global vCapitalGridFine nProdFine nCapitalFine nStateFine mProdTransition ppsiCapital


%----------------------------------------------------------------
% Compute weights for investment decision, conditional on adjusting
%----------------------------------------------------------------

% Compute nearest neighbor using discretize (no Statistics Toolbox needed)
vMidpoints = (vCapitalGridFine(1:end-1) + vCapitalGridFine(2:end)) / 2;
edges = [-Inf; vMidpoints; Inf];
vIndices 			= discretize(vCapitalAdjust, edges);
vGridIndices 		= vCapitalGridFine(vIndices);

% Find indices and gridpoints above and below
vIndicesBelow 	= vIndices; vIndicesAbove = vIndices;

vIndicesBelow(vGridIndices > vCapitalAdjust) 	= vIndicesBelow(vGridIndices > vCapitalAdjust) - 1;
vIndicesBelow(vIndicesBelow < 1) 					= 1;
vGridBelow 													= vCapitalGridFine(vIndicesBelow);

vIndicesAbove(vGridIndices <= vCapitalAdjust) 	= vIndicesAbove(vGridIndices <= vCapitalAdjust) + 1;
vIndicesAbove(vIndicesAbove > nCapitalFine)		= nCapitalFine;
vGridAbove 													= vCapitalGridFine(vIndicesAbove);

% Compute weighting matrix
vAdjustWeightBelow 											= (vGridAbove - vCapitalAdjust) ./ (vGridAbove - vGridBelow);
vAdjustWeightBelow(vCapitalAdjust < vGridBelow) 	= 1;
vAdjustWeightBelow(vCapitalAdjust > vGridAbove) 	= 0;

vAdjustWeightAbove 											= (vCapitalAdjust - vGridBelow) ./ (vGridAbove - vGridBelow);
vAdjustWeightAbove(vCapitalAdjust < vGridBelow) 	= 0;
vAdjustWeightAbove(vCapitalAdjust > vGridAbove) 	= 1;

% Rename indices to denote conditional on adjusting
vIndicesAdjustBelow = vIndicesBelow; vIndicesAdjustAbove = vIndicesAbove;
clear vIndicesBelow vIndicesAbove


%----------------------------------------------------------------
% Compute weights for investment decision, conditional on not adjusting
%----------------------------------------------------------------

% Compute nearest neighbor using discretize (no Statistics Toolbox needed)
vIndices 		= discretize(vCapitalConstrained, edges);
vGridIndices 	= vCapitalGridFine(vIndices);

% Find indices and gridpoints above and below
vIndicesBelow = vIndices; vIndicesAbove = vIndices;

vIndicesBelow(vGridIndices > vCapitalConstrained) 	= vIndicesBelow(vGridIndices > vCapitalConstrained) - 1;
vIndicesBelow(vIndicesBelow < 1) 						= 1;
vGridBelow 														= vCapitalGridFine(vIndicesBelow);

vIndicesAbove(vGridIndices <= vCapitalConstrained) 	= vIndicesAbove(vGridIndices <= vCapitalConstrained) + 1;
vIndicesAbove(vIndicesAbove > nCapitalFine) 			= nCapitalFine;
vGridAbove 															= vCapitalGridFine(vIndicesAbove);

% Compute weighting matrix
vConstrainedWeightBelow 													= (vGridAbove - vCapitalConstrained) ./ (vGridAbove - vGridBelow);
vConstrainedWeightBelow(vCapitalConstrained < vGridBelow) 	= 1;
vConstrainedWeightBelow(vCapitalConstrained > vGridAbove) 	= 0;

vConstrainedWeightAbove 													= (vCapitalConstrained - vGridBelow) ./ (vGridAbove - vGridBelow);
vConstrainedWeightAbove(vCapitalConstrained < vGridBelow) 	= 0;
vConstrainedWeightAbove(vCapitalConstrained > vGridAbove) 	= 1;

% Rename indices to denote conditional on not adjusting
vIndicesConstrainedBelow = vIndicesBelow; vIndicesConstrainedAbove = vIndicesAbove;
clear vIndicesBelow vIndicesAbove

%----------------------------------------------------------------
% Compute Transition Matrix
%----------------------------------------------------------------

nF = nStateFine;
nK = nCapitalFine;

% Build sparse nF x nK capital transition matrices (each row has at most 2 nonzeros)
vRows = (1:nF)';
mAdjustTransition_cap = sparse(vRows, vIndicesAdjustBelow, vAdjustWeightBelow, nF, nK) + ...
                         sparse(vRows, vIndicesAdjustAbove, vAdjustWeightAbove, nF, nK);
mConstrainedTransition_cap = sparse(vRows, vIndicesConstrainedBelow, vConstrainedWeightBelow, nF, nK) + ...
                              sparse(vRows, vIndicesConstrainedAbove, vConstrainedWeightAbove, nF, nK);

% Expand along productivity dimension using kron: each capital bin maps to all nProdFine states
mExpandOp = kron(speye(nK), ones(1, nProdFine));  % nK x nStateFine sparse
mAdjustTransition = mAdjustTransition_cap * mExpandOp;
mConstrainedTransition = mConstrainedTransition_cap * mExpandOp;

% Combine the adjust and constrained decision, and weight by productivity draws
mProdTransitionExpanded 	= repmat(mProdTransition,nCapitalFine);
dAdj = spdiags(vCutoff / ppsiCapital, 0, nF, nF);
dCon = spdiags(1 - vCutoff / ppsiCapital, 0, nF, nF);
mTransition = (dAdj * mAdjustTransition + dCon * mConstrainedTransition) .* ...
											mProdTransitionExpanded;
