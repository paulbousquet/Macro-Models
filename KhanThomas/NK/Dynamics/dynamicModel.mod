// Solve for aggregate dynamics using DYNARE
//
// Thomas Winberry, February 14th, 2018

//----------------------------------------------------------------
// Load parameters
//----------------------------------------------------------------

@#include "parameters.mod"


//----------------------------------------------------------------
// Define variables
//----------------------------------------------------------------

@#include "variables.mod"


//----------------------------------------------------------------
// Model equations
//----------------------------------------------------------------

model;

	@#include "equations.mod"

end;


//----------------------------------------------------------------
// Set options
//----------------------------------------------------------------

// Specify shock process (=1 to include shock; =0 to not include shock)
shocks;
    var aggregateTFPShock    = 1;
    var monetaryPolicyShock  = 1;
    var aggregateBetaShock   = 1;
end;


// Accept SS from the .m file (VFI approximation has small residuals by design)
options_.steadystate.nocheck = 1;

// Show residuals for diagnostic purposes (does not block execution)
steady;
resid;

// Check Blanchard-Kahn conditions
check;


// Simulate
stoch_simul(order=@{perturbation_order},irf=40,hp_filter=100,periods=1000,pruning)
    aggregateTFP aggregateBeta inflation nominalInterestRate intermediateGoodsPrice tobinQ
    logAggregateOutput logAggregateInvestment logAggregateConsumption logAggregateHours
    logWage realInterestRate moment_2 moment_4 moment_5 logMarginalUtility;

// Change directory
cd('..')
