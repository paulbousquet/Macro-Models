// Specifies equations for "dynamicModel.mod"
//
// Thomas Winberry, February 15th, 2018


//----------------------------------------------------------------
// Bellman equation (#equations = nState)
//----------------------------------------------------------------

// Preliminary: expand capitalAdjust (which is defined only for productivity) along entire grid (productivity and capital)
// (promoted from # model-local to var endogenous equations)
@#for iCapital in 1 : nCapital
	@#for iProd in 1 : nProd
		@#define iState = nProd * (iCapital - 1) + iProd
		capitalAdjustExpanded_@{iState} = capitalAdjust_@{iProd};
	@#endfor
@#endfor


// Define the Bellman equation for each point in iState in the individual state space
@#for iState in 1 : nState

	// STEP 1: COMPUTE EXPECTED VALUE FUNCTION NEXT PERIOD
	// Step 1a: compute expected value function next period, conditional on adjusting
	// (promoted from # model-local to var endogenous equation)
	adjustExpectedValueFunction_@{iState} = 0
	@#for iShock in 1 : nShocks
		+ shocksWeights_@{iShock} * (
		@#for iPowerCapital in 1 : nCapital
			@#for iPowerProd in 1 : nProd
				@#define power = nProd * (iPowerCapital - 1) + iPowerProd
				+ valueCoefficient_@{power}(+1) * valueFunctionPrimePolys_@{iState}_@{iShock}_@{iPowerProd} *
					cos((@{iPowerCapital} - 1) * acos(min(max(2 * ((capitalAdjustExpanded_@{iState} - capitalMin) /
					(capitalMax - capitalMin)) - 1,-1),1)))
			@#endfor
		@#endfor
		)
	@#endfor
	;

	// Step 1b: compute constrained expected value function
	// (promoted from # model-local to var endogenous equations)
	capitalConstrained_@{iState} = min(max(capitalAdjustExpanded_@{iState},
		(1 - ddelta + aaLower) * valueGrid_@{iState}_2),(1 - ddelta + aaUpper) * valueGrid_@{iState}_2);
	constrainedExpectedValueFunction_@{iState} = 0
	@#for iShock in 1 : nShocks
		+ shocksWeights_@{iShock} * (
		@#for iPowerCapital in 1 : nCapital
			@#for iPowerProd in 1 : nProd
				@#define power = nProd * (iPowerCapital - 1) + iPowerProd
				+ valueCoefficient_@{power}(+1) * valueFunctionPrimePolys_@{iState}_@{iShock}_@{iPowerProd} *
					cos((@{iPowerCapital} - 1) * acos(min(max(2 * ((capitalConstrained_@{iState} - capitalMin) /
					(capitalMax - capitalMin)) - 1,-1),1)))
			@#endfor
		@#endfor
		)
	@#endfor
	;

	// Step 1c: compute adjustment threshold
	// (promoted from # model-local to var endogenous equation)
	cutoff_@{iState} = min(max((1 / (wage * marginalUtility)) * (-(marginalUtility * tobinQ) * (capitalAdjustExpanded_@{iState} -
        capitalConstrained_@{iState}) + bbeta * exp(aggregateBeta) * (adjustExpectedValueFunction_@{iState} - constrainedExpectedValueFunction_@{iState})),
		0),ppsiCapital);


	// STEP 2: COMPUTE ENTIRE BELLMAN EQUATION
	// Left hand side: current value function
	0
	@#for iCoefficient in 1 : nState
		+ valueCoefficient_@{iCoefficient} * valueFunctionPolys_@{iState}_@{iCoefficient}
	@#endfor
	// Right hand side: flow profits plus expected value function
	= marginalUtility * (nnu ^ (nnu / (1 - nnu)) - nnu ^ (1 / (1 - nnu))) * ((intermediateGoodsPrice * exp(aggregateTFP) *
		exp(valueGrid_@{iState}_1)) ^ (1 / (1 - nnu))) * (valueGrid_@{iState}_2 ^ (ttheta / (1 - nnu))) *
		(wage ^ (-nnu / (1 - nnu)))	+ (cutoff_@{iState} / ppsiCapital) *
        (-(marginalUtility * tobinQ) * (capitalAdjustExpanded_@{iState} - (1 - ddelta) * valueGrid_@{iState}_2) -
        (cutoff_@{iState} / 2) * wage * marginalUtility + bbeta * exp(aggregateBeta) * adjustExpectedValueFunction_@{iState}) +
        (1 - (cutoff_@{iState} / ppsiCapital)) * (-(marginalUtility * tobinQ) * (capitalConstrained_@{iState} -
        (1 - ddelta) * valueGrid_@{iState}_2) + bbeta * exp(aggregateBeta) * constrainedExpectedValueFunction_@{iState});

@#endfor


//----------------------------------------------------------------
// FOC for adjust capital decision (#equations = nProd)
//----------------------------------------------------------------

@#for iProd in 1 : nProd
	(marginalUtility * tobinQ) = bbeta * exp(aggregateBeta) * (
	@#for iShock in 1 : nShocks
		+ shocksWeights_@{iShock} * (0
		@#for iPowerCapital in 1 : nCapital - 1
			@#for iPowerProd in 1 : nProd
				@#define power = nProd * iPowerCapital + iPowerProd
					+ (2 / (capitalMax - capitalMin)) * @{iPowerCapital} * valueCoefficient_@{power}(+1) * valueFunctionPrimePolys_@{iProd}_@{iShock}_@{iPowerProd} *
						(sin((@{iPowerCapital}) * acos(min(max(2 * ((capitalAdjust_@{iProd} - capitalMin) / (capitalMax - capitalMin)) - 1, - 1), 1))) /
						sin(acos(min(max(2 * ((capitalAdjust_@{iProd} - capitalMin) / (capitalMax - capitalMin)) - 1, - 1), 1))))
			@#endfor
		@#endfor
		)
	@#endfor
	);
@#endfor


//----------------------------------------------------------------
// ASIDE: compute objects over quadrature grid for integrating distribution
//		(necessary for the remaining equilibrium conditions, so its convenient
//		 to do all at once here)
//----------------------------------------------------------------

// Compute coefficients of Chebyshev interpolation of capital accumulation decision, conditional on adjusting
// Will use this interpolation to evaluate capitalAdjust over the entire grid for integrating distribution
// (promoted from # model-local to var endogenous equations)
@#for iState in 1 : nState
	capitalAdjustCoefficient_@{iState} = 0
	@#for iCoefficient in 1 : nState
		+ capitalAdjustExpanded_@{iCoefficient} * valueFunctionPolys_@{iCoefficient}_@{iState} / valueFunctionPolySquared_@{iState}
	@#endfor
	;
@#endfor


// Other objects over the grid
// (all promoted from # model-local to var endogenous equations)
@#for iState in 1 : nStateQuadrature

	// PDF of distribution
	measurePDF_@{iState} = exp(0 + measureCoefficient_1 * (quadratureGrid_@{iState}_1 - moment_1(-1)) +
	measureCoefficient_2 * (log(quadratureGrid_@{iState}_2) - moment_2(-1))
	@#define nIndexCounter = 3
	@#for iOrder in 2 : nMeasure
		@#for iPower in 0 : iOrder
			+ measureCoefficient_@{nIndexCounter} * (((quadratureGrid_@{iState}_1 - moment_1(-1)) ^ (@{iOrder} - @{iPower})) *
			((log(quadratureGrid_@{iState}_2) - moment_2(-1)) ^ @{iPower}) - moment_@{nIndexCounter}(-1))
			@#define nIndexCounter = nIndexCounter + 1
		@#endfor
	@#endfor
	);

	// Capital accumulation decision, conditional on adjusting
	capitalAdjustQuadrature_@{iState} = min(max(0 +
	@#for iCoefficient in 1 : nState
		+ capitalAdjustCoefficient_@{iCoefficient} * quadraturePolys_@{iState}_@{iCoefficient}
	@#endfor
	,capitalMin),capitalMax);

	// Capital accumulation decision, conditional on not adjusting
	capitalConstrainedQuadrature_@{iState} = min(max(min(max(capitalAdjustQuadrature_@{iState},(1 - ddelta + aaLower) * quadratureGrid_@{iState}_2),
		(1 - ddelta + aaUpper) * quadratureGrid_@{iState}_2),capitalMin),capitalMax);

	// Fixed cost threshold for adjusting or not
	cutoffQuadrature_@{iState} = min(max((1 / (wage * marginalUtility)) * (-(marginalUtility * tobinQ) * (capitalAdjustQuadrature_@{iState} -
        capitalConstrainedQuadrature_@{iState}) + bbeta * exp(aggregateBeta) * (0
	@#for iShock in 1 : nShocks
		+ shocksWeights_@{iShock} * (
		@#for iPowerCapital in 1 : nCapital
			@#for iPowerProd in 1 : nProd
				@#define power = nProd * (iPowerCapital - 1) + iPowerProd
				+ valueCoefficient_@{power}(+1) * quadraturePrimePolys_@{iState}_@{iShock}_@{iPowerProd} * (
					cos((@{iPowerCapital} - 1) * acos(min(max(2 * ((capitalAdjustQuadrature_@{iState} - capitalMin)
					/ (capitalMax - capitalMin)) - 1, - 1), 1))) - cos((@{iPowerCapital} - 1) * acos(min(max(2 *
					((capitalConstrainedQuadrature_@{iState} - capitalMin) / (capitalMax - capitalMin)) - 1, - 1), 1))))
			@#endfor
		@#endfor
		)
	@#endfor
	)),0),ppsiCapital);

@#endfor


// Compute total mass of distribution for normalization
// (promoted from # model-local to var endogenous equation)
totalMass = 0
@#for iState in 1 : nStateQuadrature
	+ quadratureWeights_@{iState} * measurePDF_@{iState}
@#endfor
;


//----------------------------------------------------------------
// Relationship between moments of distribution and parameters
// (#equations = nMeasureCoefficients)
//----------------------------------------------------------------

// First moments (uncentered)
moment_1(-1) = (0
@#for iState in 1 : nStateQuadrature
	+ quadratureWeights_@{iState} * quadratureGrid_@{iState}_1 * measurePDF_@{iState}
@#endfor
) / totalMass;

moment_2(-1) = (0
@#for iState in 1 : nStateQuadrature
	+ quadratureWeights_@{iState} * log(quadratureGrid_@{iState}_2) * measurePDF_@{iState}
@#endfor
) / totalMass;


// Higher order moments (centered)
@#define nIndexCounter = 3
@#for iOrder in 2 : nMeasure
	@#for iPower in 0 : iOrder
		moment_@{nIndexCounter}(-1) = (0
		@#for iState in 1 : nStateQuadrature
			+ quadratureWeights_@{iState} * measurePDF_@{iState} * ((quadratureGrid_@{iState}_1 - moment_1(-1)) ^ (@{iOrder} - @{iPower})) *
			((log(quadratureGrid_@{iState}_2) - moment_2(-1)) ^ @{iPower})
		@#endfor
		) / totalMass;
		@#define nIndexCounter = nIndexCounter + 1
	@#endfor
@#endfor


//----------------------------------------------------------------
// Law of motion for distribution
// (#equations = nMeasureCoefficients)
//----------------------------------------------------------------

// First moment of productivity (uncentered)
moment_1 = (0
@#for iShock in 1 : nShocks
	+ shocksWeights_@{iShock} * (
	@#for iState in 1 : nStateQuadrature
		+ quadratureWeights_@{iState} * measurePDF_@{iState} * quadratureProdPrimeGrid_@{iState}_@{iShock}
	@#endfor
	)
@#endfor
) / totalMass;

// First moment of capital (uncentered)
moment_2 = (0
@#for iShock in 1 : nShocks
	+ shocksWeights_@{iShock} * (
	@#for iState in 1 : nStateQuadrature
		+ quadratureWeights_@{iState} * measurePDF_@{iState} * ((cutoffQuadrature_@{iState} / ppsiCapital) * log(capitalAdjustQuadrature_@{iState}) +
			(1 - (cutoffQuadrature_@{iState} / ppsiCapital)) * log(capitalConstrainedQuadrature_@{iState}))
	@#endfor
	)
@#endfor
) / totalMass;

// Higher order moments (centered)
@#define nIndexCounter = 3
@#for iOrder in 2 : nMeasure
	@#for iPower in 0 : iOrder
		moment_@{nIndexCounter} = (0
		@#for iShock in 1 : nShocks
			+ shocksWeights_@{iShock} * (
			@#for iState in 1 : nStateQuadrature
				+ quadratureWeights_@{iState} * measurePDF_@{iState} * ((quadratureProdPrimeGrid_@{iState}_@{iShock} - moment_1) ^ (@{iOrder} - @{iPower})) *
					(((cutoffQuadrature_@{iState} / ppsiCapital) * (log(capitalAdjustQuadrature_@{iState}) - moment_2) ^ @{iPower}) +
					((1 - (cutoffQuadrature_@{iState} / ppsiCapital)) * (log(capitalConstrainedQuadrature_@{iState}) - moment_2) ^ @{iPower}))
			@#endfor
			)
		@#endfor
		) / totalMass;
		@#define nIndexCounter = nIndexCounter + 1
	@#endfor
@#endfor

//----------------------------------------------------------------
// Labor Market clearing (# equations = 2)
//----------------------------------------------------------------

// Definition of aggregate hours from labor demand
aggregateHours = (0
@#for iState in 1 : nStateQuadrature
	+ quadratureWeights_@{iState} * measurePDF_@{iState} * (((intermediateGoodsPrice * nnu * exp(aggregateTFP) * exp(quadratureGrid_@{iState}_1)
		* (quadratureGrid_@{iState}_2 ^ ttheta)) / wage) ^ (1 / (1 - nnu)) + ((cutoffQuadrature_@{iState} ^ 2) / (2 * ppsiCapital)))
@#endfor
) / totalMass;

// Labor demand = labor supply
aggregateHours = ((wage * marginalUtility) / cchi) ^ (1 / pphi);


//----------------------------------------------------------------
// Aggregate output and investment from distribution (# equations = 2)
//----------------------------------------------------------------

// Aggregate output level
aggregateOutputLevel = (0
@#for iState in 1 : nStateQuadrature
    + quadratureWeights_@{iState} * measurePDF_@{iState} * (exp(aggregateTFP) * exp(quadratureGrid_@{iState}_1) * (quadratureGrid_@{iState}_2 ^ ttheta) *
        ((((intermediateGoodsPrice * nnu * exp(aggregateTFP) * exp(quadratureGrid_@{iState}_1) * (quadratureGrid_@{iState}_2 ^ ttheta)) / wage) ^
        (1 / (1 - nnu))) ^ nnu))
@#endfor
) / totalMass;

// Aggregate investment level (net)
aggregateInvestmentLevel = (0
@#for iState in 1 : nStateQuadrature
    + quadratureWeights_@{iState} * measurePDF_@{iState} * ((cutoffQuadrature_@{iState} / ppsiCapital) * capitalAdjustQuadrature_@{iState} +
        (1 - (cutoffQuadrature_@{iState} / ppsiCapital)) * capitalConstrainedQuadrature_@{iState} - (1 - ddelta) * quadratureGrid_@{iState}_2)
@#endfor
) / totalMass;


//----------------------------------------------------------------
// Capital, Investment Price, and Resource Constraint (# equations = 4)
//----------------------------------------------------------------

// Aggregate capital stock law of motion
aggregateCapitalStock = (1 - ddelta) * aggregateCapitalStock(-1) + aggregateInvestmentLevel;

// Tobin's q (capital producer FOC)
tobinQ = (aggregateInvestmentLevel / (aggregateCapitalStock(-1) * ddelta)) ^ (1 / kkappa);

// Investment resource cost
investmentResourceCost = aggregateCapitalStock(-1) * ddelta * (kkappa / (kkappa + 1) *
    (aggregateInvestmentLevel / (aggregateCapitalStock(-1) * ddelta)) ^ ((kkappa + 1) / kkappa) +
    1 / (kkappa + 1));

// Resource constraint: Y*(1-Rotemberg_cost) = C + investment_resource_cost
aggregateConsumption = aggregateOutputLevel * (1 - vvarphi / 2 * inflation ^ 2) - investmentResourceCost;

// Marginal utility = u'(aggregate consumption)
marginalUtility = aggregateConsumption ^ (-ssigma);


//----------------------------------------------------------------
// New Keynesian Equations (# equations = 4)
//----------------------------------------------------------------

// NK Phillips Curve (Rotemberg pricing)
(1 - ggamma) + ggamma * intermediateGoodsPrice
    - vvarphi * inflation * (1 + inflation)
    + bbeta * exp(aggregateBeta) * (marginalUtility(+1) / marginalUtility) * vvarphi
      * inflation(+1) * (1 + inflation(+1)) * (aggregateOutputLevel(+1) / aggregateOutputLevel) = 0;

// Taylor Rule for monetary policy
log(1 + nominalInterestRate) = rrho_r * log(1 + nominalInterestRate(-1))
    + (1 - rrho_r) * (log(1 / bbeta) + vvarphi_pi * log(1 + inflation))
    + ssigmaMonetary * monetaryPolicyShock;

// Expected inflation (for Fisher equation)
expectedInflation = inflation(+1);

// Household Euler equation (Fisher)
1 = bbeta * exp(aggregateBeta) * (1 + nominalInterestRate)
    * expectedMarginalUtilityPrime / (marginalUtility * (1 + expectedInflation));


//----------------------------------------------------------------
// Law of motion for aggregate shocks (# equations = 2)
//----------------------------------------------------------------

aggregateTFP = rrhoTFP * aggregateTFP(-1) + ssigmaTFP * aggregateTFPShock;
aggregateBeta = rrhoBeta * aggregateBeta(-1) + ssigmaBeta * aggregateBetaShock;


//----------------------------------------------------------------
// Auxiliary variables of interest (# equations = 12)
//----------------------------------------------------------------

// Expected marginal utility (used in Euler equation above)
expectedMarginalUtilityPrime = marginalUtility(+1);

// Log aggregate output (uses aggregateOutputLevel defined above)
logAggregateOutput = log(aggregateOutputLevel);

// Log aggregate investment (uses aggregateInvestmentLevel defined above)
logAggregateInvestment = log(aggregateInvestmentLevel);

// Logs of aggregate variables already defined
logAggregateHours = log(aggregateHours);
logAggregateConsumption = log(aggregateConsumption);
logWage = log(wage);
logMarginalUtility = log(marginalUtility);

// Real interest rate (from Fisher equation)
realInterestRate = 100 * ((1 + nominalInterestRate) / (1 + expectedInflation) - 1);

// NK log variables
logInflation = log(1 + inflation);
logNominalInterestRate = log(1 + nominalInterestRate);
logTobinQ = log(tobinQ);
logIntermediateGoodsPrice = log(intermediateGoodsPrice);
