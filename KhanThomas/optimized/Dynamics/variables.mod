// Declare variables for "dynamicModel.mod"
//
// Thomas Winberry, Feburary 15th, 2018

//----------------------------------------------------------------
// Value function coefficients
//----------------------------------------------------------------

@#for iState in 1 : nState
    var valueCoefficient_@{iState};
@#endfor


//----------------------------------------------------------------
// Adjust policy function along productivity grid
//----------------------------------------------------------------

@#for iState in 1 : nProd
    var capitalAdjust_@{iState};
@#endfor


//----------------------------------------------------------------
// Moments
//----------------------------------------------------------------

@#for iMoment in 1 : nMeasureCoefficients
    var moment_@{iMoment};
@#endfor


//----------------------------------------------------------------
// Distribution parameters
//----------------------------------------------------------------

@#for iParameter in 1 : nMeasureCoefficients
    var measureCoefficient_@{iParameter};
@#endfor


//----------------------------------------------------------------
// Prices
//----------------------------------------------------------------

var wage marginalUtility;


//----------------------------------------------------------------
// Aggregate shocks
//----------------------------------------------------------------

var aggregateTFP aggregateQ aggregateBeta;


//----------------------------------------------------------------
// Auxiliary variables we're interested in
//----------------------------------------------------------------

var aggregateConsumption aggregateHours expectedMarginalUtilityPrime realInterestRate
	logAggregateOutput logAggregateConsumption logAggregateInvestment logAggregateHours logWage
	logMarginalUtility;


//----------------------------------------------------------------
// Auxiliary variables promoted from # model-local to var
// (prevents expression inlining during Dynare symbolic differentiation)
//----------------------------------------------------------------

// capitalAdjustExpanded: expand capitalAdjust along entire grid (nState = 15)
@#for iState in 1 : nState
    var capitalAdjustExpanded_@{iState};
@#endfor

// adjustExpectedValueFunction: expected value conditional on adjusting (nState = 15)
@#for iState in 1 : nState
    var adjustExpectedValueFunction_@{iState};
@#endfor

// capitalConstrained: constrained capital choice (nState = 15)
@#for iState in 1 : nState
    var capitalConstrained_@{iState};
@#endfor

// constrainedExpectedValueFunction: expected value conditional on NOT adjusting (nState = 15)
@#for iState in 1 : nState
    var constrainedExpectedValueFunction_@{iState};
@#endfor

// cutoff: adjustment threshold (nState = 15)
@#for iState in 1 : nState
    var cutoff_@{iState};
@#endfor

// capitalAdjustCoefficient: Chebyshev interpolation coefficients (nState = 15)
@#for iState in 1 : nState
    var capitalAdjustCoefficient_@{iState};
@#endfor

// measurePDF: PDF of distribution over quadrature grid (nStateQuadrature = 80)
@#for iState in 1 : nStateQuadrature
    var measurePDF_@{iState};
@#endfor

// capitalAdjustQuadrature: capital adjust on quadrature grid (nStateQuadrature = 80)
@#for iState in 1 : nStateQuadrature
    var capitalAdjustQuadrature_@{iState};
@#endfor

// capitalConstrainedQuadrature: constrained capital on quadrature grid (nStateQuadrature = 80)
@#for iState in 1 : nStateQuadrature
    var capitalConstrainedQuadrature_@{iState};
@#endfor

// cutoffQuadrature: adjustment threshold on quadrature grid (nStateQuadrature = 80)
@#for iState in 1 : nStateQuadrature
    var cutoffQuadrature_@{iState};
@#endfor

// totalMass: total mass of distribution for normalization (1)
var totalMass;


//----------------------------------------------------------------
// Innovations to aggregate shocks
//----------------------------------------------------------------

varexo aggregateTFPShock aggregateQShock aggregateBetaShock;
