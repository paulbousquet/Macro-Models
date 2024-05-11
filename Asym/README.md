The notable feature of this model is the presence of asymmetric setting in price and wage
setting. Specifically, we assume there is an adjustment cost $\Phi_p(\pi_t)$ and
$\Phi_p(\pi_t^w)$ and it's more costly to lower prices and wages versus raise them. This
additional downward rigidity was formalized in a New Keynesian model first by Kim and
Ruge-Murcia (2009) and is based on observations first made by Keynes that have been
continually reconfirmed. ABS code this structure by treating these adjustment cost functions
(and their derivatives) as variables and then making the realizations themselves model
equations. Some of their work in the appendix suggests that the equations that appear in the
code represent "detrended" variables (normalized by technology), meaning only the
fluctuation process $\log(z_t)=\rho_z \log(z_{t-1})+\varepsilon_t^z$ shows up in the equations.
What they code slightly differs from this structure (see (12) ), so to be consistent what
they coded while noting the inconsistency with their paper, I have denote this amorphous
technology term by $\tilde{A}_t$.

| Description | Equation | # |
|-------------|----------|---|
| Consumption Euler Equation | $1 = \beta \left( \frac{C_{t+1}}{C_t} \right)^{- \tau} \frac{R_t}{\Pi_{t+1} \tilde{A}_{t+1}}$ | (1) |
| Definition for Real Wages | $\Delta_t^w = \frac{W_t}{W_{t-1}} \cdot \tilde{A}_t$ | (2) |
| Resource Constraint | $\frac{G_t - 1}{G_t} \cdot Y_t + C_t = Y_t (1 - \Phi_t^p) + W_t Y_t \cdot \Phi_t^w$ | (3) |
| Wage Equation, Household's problem | $\frac{\chi_h}{\lambda_w} \cdot W_t^{-\tau} C_t^{\tau} Y_t^{\frac{1}{\nu}} + (1 - \Phi_t^w)$ $\left(1 - \frac{1}{\lambda_w}\right) = \ \Delta_t^{w_{nom}} \cdot \Phi_t^{'w} - \beta \left( \frac{C_{t+1}}{C_t} \right)^{- \tau} \frac{R_t}{\Pi_{t+1} + \tilde{A}\_{t+1}} W_{t+1}^2 Y_{t+1} \cdot \Phi_{t+1}^{'w}$ | (4) |
| Price Equation, Intermediate Firms problem | $(1 - \Phi_t^p) - \Pi_t \cdot \Phi_t^{'p} - \frac{\mu_t}{\Lambda_t} = \beta \left( \frac{C_{t+1}}{C_t} \right)^{- \tau} \frac{\Pi_{t+1}}{ \tilde{A}\_{t+1} } \cdot \Phi_{t+1}^{'p} \cdot Y_{t+1} \tilde{A}_{t+1}$ | (5) |
| Hours Equation | $W_t = (1 - \Phi_t^p) + \mu_t$ | (6) |
| Adjustment Costs, Nominal Wages | $\Phi_t^w = \frac{\phi_w}{\psi_w^2} \left( e^{-\psi_w (W_{t}^{w_{nom}} - \gamma \cdot \pi^\star)} + \psi_w (W_{t}^{w_{nom}} - \gamma \pi^\star) - 1 \right)$ | (7) |
| Adjustment Costs, Prices | $\Phi_t^p = \frac{\phi_p}{\psi_p^2} \left(e^{-\psi_p (e^{\pi_t} - \pi^\star)} + \psi_p (e^{\pi_t} - \pi^\star) - 1\right)$ | (8) |
| Derivative, Adjustment Costs Nominal Wages | $\Phi_t^p = \frac{\phi_p}{\psi_p^2} \left( e^{-\psi_p (\Pi_t - \pi^\star)} + \psi_p (\Pi_t - \pi^\star) - 1 \right)$ | (9) |
| Derivative, Adjustment Costs to Prices | $\Phi_t^{'p} = \frac{\phi_p}{\psi_p} \left(1 - e^{-\psi_p \left(\Pi_t - \pi^\star\right)}\right)$ | (10) |
| Taylor Rule | $R_t = \exp(r_t);\quad r_t = \rho_r r_{t-1} + (1 - \rho_r)r_t^\star + \sigma_r \varepsilon_r$ | (11) |
| TFP Growth | $\tilde{A}\_t = \exp(a_t); \quad a_t = (1 - \rho_a)  \log{\gamma} + \rho_a  a_{t-1} + \sigma_a  \varepsilon_a$ | (12) |
| Government Spending Shocks | $G_t = \exp(G_t); \quad g_t = (1 - \rho_g)  \log{g^\star} + \rho_g  g_{t-1} + \sigma_a  \varepsilon_g$ | (13) |
| Price Markup Shock | $\lambda_t = \exp(\lambda_t); \quad\lambda_t = (1 - \rho_p)  \log{\lambda_{p_{ss}}} + \rho_p  \lambda_{t-1} + \sigma_p  \varepsilon_p$ | (14) |
| Output change | $\Delta_t^y = Y_t/Y_{t-1}$ | (15) |
| Nominal wage change | $\Delta_t^{w_{nom}} = \Delta_t^w\Pi_t$ | (16) |
| Interest rate target | $r_t^\star = \log{\left(\frac{\gamma}{\beta} \cdot \pi^\star\right)} + \psi_1 \left(\pi_t - \log{\pi^\star}\right) + \psi_2 \left(\Delta_t^y + a_t - \log{\gamma}\right)$ | (17) |

It's also important to note that their original code has everything in log-linear terms. This has the advantage of making the exogenous processes and change variables very clean, at the expense of having to carry around a bunch of exponential terms. For visual presentation, presenting everything in pure levels is much less cumbersome, but for completeness here is how the equations appear in their code

| Description | Equation | # |
|-------------|----------|---|
| Consumption Euler Equation | 1 = $\beta \left(\frac{e^{c_{t+1}}}{e^{c_t}}\right)^{- \tau} \frac{e^{R_t}}{e^{\pi_{t+1}+a_{t+1}}}$ | (1) |
| Definition for Real Wages | $e^{\Delta_t^w} = \frac{e^{w_t}}{e^{w_{t-1}}} \cdot e^{a_t}$ | (2) |
| Resource Constraint | $e^{c_t} + \frac{e^{g_t} - 1}{e^{g_t}} \cdot e^{y_t} = e^{y_t} (1 - \Phi_t^p) + e^{w_t+y_t} \cdot \Phi_t^w$ | (3) |
| Wage Equation, Household's problem | $\frac{\chi_h}{\lambda_w} \cdot e^{-\tau w_t + \tau c_t + \frac{1}{\nu} y_t} + (1 - \Phi_t^w)$ $\left(1 - \frac{1}{\lambda_w}\right) = \ e^{\Delta_t^{w_{nom}}} \cdot \Phi_t^{'w} - \beta \left(\frac{e^{c_{t+1}}}{e^{c_t}}\right)^{- \tau} \frac{e^{R_t}}{e^{\pi_{t+1}+a_{t+1}}}e^{ 2 \Delta_{t+1}^w + \Delta_{t+1}^y} \cdot \Phi_{t+1}^{'w}$ | (4)  |
| Price Equation, Intermediate Firms problem | $(1 - \Phi_t^p) - e^{\pi_t} \cdot \Phi_t^{'p} - \frac{\mu_t}{e^{\lambda_t}} =$ | (5) |
|  | $\ \beta \left(\frac{e^{c_{t+1}}}{e^{c_t}}\right)^{- \tau} \frac{e^{\pi_{t+1}}}{e^{a_{t+1}}}  \cdot \Phi_{t+1}^{'p} \cdot e^{\Delta_{t+1}^y+a_{t+1}}$ |  |
| Hours Equation | $e^{w_t} = (1 - \Phi_t^p) + \mu_t$ | (6) |
| Adjustment Costs, Nominal Wages | $\Phi_t^w = \frac{\phi_w}{\psi_w^2} \left(e^{-\psi_w (e^{\Delta_t^{w_{nom}}} - \gamma \cdot \pi^\star)} + \psi_w (e^{\Delta_t^{w_{nom}}}  - \gamma \pi^\star) - 1\right)$ | (7) |
| Adjustment Costs, Prices | $\Phi_t^p = \frac{\phi_p}{\psi_p^2} \left(e^{-\psi_p (e^{\pi_t} - \pi^\star)} + \psi_p (e^{\pi_t} - \pi^\star) - 1\right)$ | (8) |
| Derivative, Adjustment Costs Nominal Wages | $\Phi_t^{'w} = \frac{\phi_w}{\psi_w} \left(1 - e^{-\psi_w \left(e^{\Delta_t^w} \cdot e^{\pi_t} - \gamma \cdot \pi^\star\right)}\right)$ | (9) |
| Derivative, Adjustment Costs to Prices | $\Phi_t^{'p} = \frac{\phi_p}{\psi_p} \left(1 - e^{-\psi_p \left(e^{\pi_t} - \pi^\star\right)}\right)$ | (10) |
| Taylor Rule | $R_t = \rho_r R_{t-1} + (1 - \rho_r)R_t^\star + \sigma_r \varepsilon_r$ | (11) |
| TFP Growth | $a_t = (1 - \rho_a)  \log{\gamma} + \rho_a  a_{t-1} + \sigma_a  \varepsilon_a$ | (12) |
| Government Spending Shocks | $g_t = (1 - \rho_g)  \log{g^\star} + \rho_g  g_{t-1} + \sigma_a  \varepsilon_g$ | (13) |
| Price Markup Shock | $\lambda_t = (1 - \rho_p)  \log{\lambda_{p_{ss}}} + \rho_p  \lambda_{t-1} + \sigma_p  \varepsilon_p$ | (14) |
| Output change | $\Delta_t^y = y_t - y_{t-1}$ | (15) |
| Nominal wage change | $\Delta_t^{w_{nom}} = \Delta_t^w + \pi_t$ | (16) |
| Interest rate target | $R_t^\star = \log{\left(\frac{\gamma}{\beta} \cdot \pi^\star\right)} + \psi_1 \left(\pi_t - \log{\pi^\star}\right) + \psi_2 \left(\Delta_t^y + a_t - \log{\gamma}\right)$ | (17) |

The only real abuse of notation across the two tables comes from the change variables, which are notated the exact same.
