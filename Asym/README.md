The following discussion relates to Sections 5 and 6 of ["Assessing DSGE Nonlinearities"](http://www.sciencedirect.com/science/article/pii/S0165188917301562) by Aruoba, Bocola, and Schorfheide (2017), hereafter ABS. This is a description of the changes made to the public replication files themselves (so they can run as intended), a general discussion of the model that is implied by their code, and finally an extension based on casting in Dynare that potentially uncovered some issues with how they inputted the model. 

## Code Modifications 
* There are also two versions of the original replication files floating around, though they are mostly similar. For reference, the most up to date, publicly available ones can be currently found [here](https://www.dropbox.com/s/j9cr6pn0kvmsnt0/Replication%20codes.zip?dl=0) 
* Some programs call a function model_ss that does not exist. Additionally, the "nonlinear steady state" file is ordered differently than everything else. This means the easiest thing to do is modify the linear steady state file, since it retains the same ordering. The only change is to include $\psi_p$ and $\psi_w$, which are the last two entries of the Theta vector (changed)
* To run the metropolis hastings routine, you need to comment out the line that adds the predictive check folder, which has its own "prior.m" file (changed but not an issue in some versions)
* The ARC in the original code contains an erroneous sign, but the value of the term is generally small anyway so thankfully doesn't seem to have made a huge difference (changed)
* To run the code, you are restricted to start at the linear mode, which has been precomputed for you. You can compute your own using the Dynare file `linlin.mod`. Their default method carries two important considerations
  * In the linear model, theres is no asymmetry. In the full model, there are cost scaling parameters, which interact with the asymmetry parameters. But under no asymmetry, in equillibrium this collapses to a matter of relative costs of adjusting prices and wages. In other words, only one object (relative cost) is identified. It's thereofore not advisible to try to estimate $\phi_w$ in the linear setting. (changed) 
  * For the nonlinear model, the inverse Hessian is appended with 4s on the diagonal for the asymmetry parameters and 0 on the off-diagonals. While this is an ultimately arbirtray decision, the key for the MH procedure to be valid is merely a positive-definite inverse Hessian; having the exact mode and Hessian just improves efficency in practice. But because the starting point and MH variance are selected in an ad-hoc way, this is worth noting since these estimation processes can be fragile along a myriad of lines. From what I've done so far, it seems this choice is non-trivial. (unchanged)  
* ABS makes a point that they find lower values of the asymmetry parameter than past approaches in the literature. But because of their uniform prior, they are ruling out the possibility this could occurr to begin with. While the authors may gesture to the fact that there is no clustering around the bounds and these parameter values seemed unreasonably high to begin with, there's no reason to not widen the bounds to more robustly conclude there's no meaniningful mass at more extreme values (changed)
* I noticed some odd patterns in the acceptance rates which can be partially explained by two things
  * There is no burn in period for draws (changed)
  * The seeding method is outdated and possibly was wrong to begin with. Using `i^(2)` as a seed (where `i` is the number of draws) will hit the upper bound for seed number, though perhaps this upper bound was higher in previous versions of Matlab. In any event, the syntax used is no longer reccomended and has been replaced by the `rng` function (changed)
* There is no prior file analog for the nonlinear model, which has two additional parameters. This is "accounted for" by the uniform prior --  the part of the log-liklihood related to these parameters will always be the same across draws that fall within the bounds. However, if you care about the actual values (e.g., if you want to compare to Dynare), this would need to be changed on the back end when computing liklihood in MH. (unchanged)
* There are several minor inconsistencies with respect to the listed priors in the table in the paper (unchanged)

## Model Equations

The notable feature of this model is the presence of asymmetric setting in price and wage
setting. Specifically, we assume there is an adjustment cost $\Phi_p(\pi_t)$ and
$\Phi_p(\pi_t^w)$ and it's more costly to lower prices and wages versus raise them. This
additional downward rigidity was formalized in a New Keynesian model first by Kim and
Ruge-Murcia (2009) and is based on observations first made by Keynes that have been
continually reconfirmed. ABS code this structure by treating these adjustment cost functions
(and their derivatives) as variables and then making the realizations themselves model
equations. Many of the macro fundamentals are "detrended" (normalized by technology). Log technology is a random walk with drift, meaning only a noise process ($\text{exp}(z_t$)) scaled by the dift shows up. To be consistent with the original notation, I denote this amorphous technology term by $\tilde{A}_t$.

| Description | Equation | # |
|-------------|----------|---|
| Consumption Euler Equation | $1 = \beta \left( \frac{C_{t+1}}{C_t} \right)^{- \tau} \frac{R_t}{\Pi_{t+1} \tilde{A}_{t+1}}$ | (1) |
| Real Wage ($W_t$) Inflation  | $\Delta_t^w = \frac{W_t}{W_{t-1}} \cdot \tilde{A}_t$ | (2) |
| Resource Constraint | $\frac{G_t - 1}{G_t} \cdot Y_t + C_t = Y_t (1 - \Phi_t^p) + W_t Y_t \cdot \Phi_t^w$ | (3) |
| Wage Equation, Household's problem | $\frac{\chi_h}{\lambda_w} \cdot W_t^{-1} C_t^{\tau} Y_t^{\frac{1}{\nu}} + (1 - \Phi_t^w)\left(1 - \lambda_w^{-1}\right) =$ $\Delta_t^{w_{nom}} \cdot \Phi_t^{'w} - \beta \left( \frac{C_{t+1}}{C_t} \right)^{- \tau} \Pi_{t+1}\frac{R_t}{\tilde{A}\_{t+1}} W_{t+1}^2 Y_{t+1} \cdot \Phi_{t+1}^{'w}$ | (4) |
| Price Equation, Intermediate Firms problem | $(1 - \Phi_t^p) + \beta \left( \frac{C_{t+1}}{C_t} \right)^{- \tau} \frac{\Pi_{t+1}}{ \tilde{A}\_{t+1} } \cdot \Phi_{t+1}^{'p} \cdot Y_{t+1} \tilde{A}_{t+1} = \Pi_t \cdot \Phi_t^{'p} + \frac{\mu_t}{\Lambda_t}$  | (5) |
| Hours Equation | $W_t = (1 - \Phi_t^p) - \mu_t$ | (6) |
| Adjustment Costs, Nominal Wages | $\Phi_t^w = \frac{\phi_w}{\psi_w^2} \left( e^{-\psi_w (\Delta_{t}^{w_{nom}} - \gamma \Pi^\star)} + \psi_w (\Delta_{t}^{w_{nom}} - \gamma \Pi^\star) - 1 \right)$ | (7) |
| Adjustment Costs, Prices | $\Phi_t^p = \frac{\phi_p}{\psi_p^2} \left(e^{-\psi_p (\Pi_t - \Pi^\star)} + \psi_p (\Pi_t - \Pi^\star) - 1\right)$ | (8) |
| Derivative, Adjustment Costs Nominal Wages | $\Phi_t^{'w} = \frac{\phi_p}{\psi_p} \left( 1-e^{-\psi_p (\Delta_t^{w_{nom}} - \gamma \Pi^\star)}  \right)$ | (9) |
| Derivative, Adjustment Costs to Prices | $\Phi_t^{'p} = \frac{\phi_p}{\psi_p} \left(1 - e^{-\psi_p \left(\Pi_t - \Pi^\star\right)}\right)$ | (10) |
| Taylor Rule | $R_t = \exp(r_t);\quad r_t = \rho_r r_{t-1} + (1 - \rho_r)r_t^\star + \sigma_r \varepsilon_r$ | (11) |
| TFP Growth | $\tilde{A}\_t = \exp(a_t); \quad a_t = (1 - \rho_a)  \log{\gamma} + \rho_a  a_{t-1} + \sigma_a  \varepsilon_a$ | (12) |
| Government Spending Shocks | $G_t = \exp(g_t); \quad g_t = (1 - \rho_g)  \log{g^\star} + \rho_g  g_{t-1} + \sigma_a  \varepsilon_g$ | (13) |
| Price Markup Shock | $\Lambda_t = \exp(\lambda_t); \quad\lambda_t = (1 - \rho_p)  \log{\lambda_{p_{ss}}} + \rho_p  \lambda_{t-1} + \sigma_p  \varepsilon_p$ | (14) |
| Output change | $\Delta_t^y = Y_t/Y_{t-1}$ | (15) |
| Nominal wage change | $\Delta_t^{w_{nom}} = \Delta_t^w\Pi_t$ | (16) |
| Interest rate target | $r_t^\star = \log{\left(\frac{\gamma}{\beta} \cdot \pi^\star\right)} + \psi_1 \left(\pi_t - \log{\pi^\star}\right) + \psi_2 \left(\Delta_t^y + a_t - \log{\gamma}\right)$ | (17) |

It's also important to note that their original code has everything in log levels. This has the advantage of making the exogenous processes and change variables very clean, at the expense of having to carry around a bunch of exponential terms. For visual presentation, presenting everything in pure levels is much less cumbersome, but for completeness here is how the equations appear in their code

| Description | Equation | # |
|-------------|----------|---|
| Consumption Euler Equation | 1 = $\beta \left(\frac{e^{c_{t+1}}}{e^{c_t}}\right)^{- \tau} e^{r_t-(\pi_{t+1}+a_{t+1})}$ | (1) |
| Real Wage ($w_t$) Inflation  | $e^{\Delta_t^w} = \frac{e^{w_t}}{e^{w_{t-1}}} \cdot e^{a_t}$ | (2) |
| Resource Constraint | $e^{c_t} + \frac{e^{g_t} - 1}{e^{g_t}} \cdot e^{y_t} = e^{y_t} (1 - \Phi_t^p) + e^{w_t+y_t} \cdot \Phi_t^w$ | (3) |
| Wage Equation, Household's problem | $\frac{\chi_h}{\lambda_w} \cdot e^{-w_t + \tau c_t + \frac{1}{\nu} y_t} + (1 - \Phi_t^w)\left(1 - \lambda_w^{-1}\right)=$ $e^{\Delta_t^{w_{nom}}} \cdot \Phi_t^{'w} - \beta \left(\frac{e^{c_{t+1}}}{e^{c_t}}\right)^{- \tau} e^{r_t-a_{t+1}}e^{\pi_{t+1}+2 \Delta_{t+1}^w + \Delta_{t+1}^y} \cdot \Phi_{t+1}^{'w}$ | (4)  |
| Price Equation, Intermediate Firms problem | $(1 - \Phi_t^p)+\beta \left(\frac{e^{c_{t+1}}}{e^{c_t}}\right)^{- \tau} \frac{e^{\pi_{t+1}}}{e^{a_{t+1}}}  \cdot \Phi_{t+1}^{'p} \cdot e^{\Delta_{t+1}^y+a_{t+1}} = e^{\pi_t} \cdot \Phi_t^{'p} + \frac{\mu_t}{e^{\lambda_t}} $ | (5)  |
| Hours Equation | $e^{w_t} = (1 - \Phi_t^p) - \mu_t$ | (6) |
| Adjustment Costs, Nominal Wages | $\Phi_t^w = \frac{\phi_w}{\psi_w^2} \left(e^{-\psi_w (e^{\Delta_t^{w_{nom}}} - \gamma \pi^\star)} + \psi_w (e^{\Delta_t^{w_{nom}}}  - \gamma \pi^\star) - 1\right)$ | (7) |
| Adjustment Costs, Prices | $\Phi_t^p = \frac{\phi_p}{\psi_p^2} \left(e^{-\psi_p (e^{\pi_t} - \pi^\star)} + \psi_p (e^{\pi_t} - \pi^\star) - 1\right)$ | (8) |
| Derivative, Adjustment Costs Nominal Wages | $\Phi_t^{'w} = \frac{\phi_w}{\psi_w} \left(1 - e^{-\psi_w (e^{\Delta_t^{w_{nom}}} - \gamma \pi^\star)}\right)$ | (9) |
| Derivative, Adjustment Costs to Prices | $\Phi_t^{'p} = \frac{\phi_p}{\psi_p} \left(1 - e^{-\psi_p \left(e^{\pi_t} - \pi^\star\right)}\right)$ | (10) |
| Taylor Rule | $r_t = \rho_r r_{t-1} + (1 - \rho_r)r_t^\star + \sigma_r \varepsilon_r$ | (11) |
| TFP Growth | $a_t = (1 - \rho_a)  \log{\gamma} + \rho_a  a_{t-1} + \sigma_a  \varepsilon_a$ | (12) |
| Government Spending Shocks | $g_t = (1 - \rho_g)  \log{g^\star} + \rho_g  g_{t-1} + \sigma_a  \varepsilon_g$ | (13) |
| Price Markup Shock | $\lambda_t = (1 - \rho_p)  \log{\lambda_{p_{ss}}} + \rho_p  \lambda_{t-1} + \sigma_p  \varepsilon_p$ | (14) |
| Output change | $\Delta_t^y = y_t - y_{t-1}$ | (15) |
| Nominal wage change | $\Delta_t^{w_{nom}} = \Delta_t^w + \pi_t$ | (16) |
| Interest rate target | $r_t^\star = \log{\left(\frac{\gamma}{\beta} \cdot \pi^\star\right)} + \psi_1 \left(\pi_t - \log{\pi^\star}\right) + \psi_2 \left(\Delta_t^y + a_t - \log{\gamma}\right)$ | (17) |

The only real abuse of notation across the two tables comes from the change variables, which are notated the exact same with the excpetion of price inflation. Also, their rules for capitalization in the appendix are for detrended variables, whereas I use it to signify variables in log terms. Similarly, I use $W_t$ for real wages in levels, but they use that for nominal wages. 

## Extension (Dynare) and Possible Issues

The foundation of the model solution is Matlab's symbolic toolbox (see the Elastic Interest rate folder in this repo for an extensive discussion on this solution method). There are some key differences between this structure and an analagous representation in Dynare. The codes use $y_{t-1}$ and $w_{t-1}$ as state variables (y0,w0) and code all exogenous processes using a $t+1$ timing. Dynare is sensitive to timing -- using "t-1" as the baseline for a forward looking variable is not possible at least in this application) and similarly "t+1" for a backward looking variable. In addition, we need to manually add in the noise processes for the exogenous variables. Some other considerations for trasnlating into Dynare are implicit from the discussion of the model equations above. Finally, we use a parameterization of (roughly) the posterior modes provided in the table and use the observation equations from the simulation files.  

### Bayesian Estimation 

([a previous thread about this issue](https://forum.dynare.org/t/replicatation-of-aruoba-bocola-and-schorfheide-2017/10979)) 
