"""
Parameters and state dynamics for the Hall of Mirrors model.

Paper: "The Natural Rate of Interest Through a Hall of Mirrors"
Authors: Rungcharoenkitkul & Winkler
References: Table 1, Equations 5.1-5.6
"""

import numpy as np

# =============================================================
# Belief hierarchy truncation
# =============================================================
N_BELIEFS = 12   # Orders of higher-order beliefs (1 through N)
N_Z = 7          # Dimension of exogenous state vector Z_t

# Derived dimensions
DIM_M = N_Z * N_BELIEFS          # 84: dimension of belief vector m_it
DIM_X = N_Z * (N_BELIEFS + 1)    # 91: dimension of state X_it = (Z_t', m_jt')'

# =============================================================
# Structural parameters (Table 1)
# =============================================================
sigma = 1.0 / 6.0     # Inverse EIS from Table 1 (testing paper's exact value)
kappa = 0.015         # Phillips curve slope
beta = 0.9941         # Discount factor
phi_pi = 1.5          # Taylor rule: inflation coefficient
phi_y = 0.125         # Taylor rule: output gap coefficient
rho_i = 0.5           # Taylor rule: interest rate smoothing

# =============================================================
# Shock autocorrelations
# =============================================================
rho_r = 1.0           # r** process (random walk)
rho_eh = 0.8          # Private signal noise autocorrelation (household)
rho_ec = 0.8          # Private signal noise autocorrelation (central bank)
rho_f = 0.8           # Public signal noise autocorrelation
rho_h = 0.8           # Demand shock autocorrelation
rho_p = 0.8           # Cost-push shock autocorrelation
rho_c = 0.5           # Monetary policy shock autocorrelation

# =============================================================
# Shock standard deviations (innovations)
# =============================================================
sigma_v = 0.05        # r** innovation SD
sigma_eps_h = 0.9     # Private signal noise innovation SD (household)
sigma_eps_c = 0.9     # Private signal noise innovation SD (central bank)
sigma_eta = 0.9       # Public signal noise innovation SD
sigma_uh = 0.17       # Demand shock innovation SD
sigma_up = 0.05       # Cost-push shock innovation SD
sigma_uc = 0.1        # Monetary policy shock innovation SD

# =============================================================
# Initial / steady-state values
# =============================================================
r_star_star_0 = 0.0235   # Initial r** (2.35% annual, quarterly model)
pi_star = 0.02            # Steady-state inflation (2%)

# =============================================================
# State vector Z_t = (r**, e_h, e_c, f, u_h, u_p, u_c)'
# Indices for accessing Z_t components
# =============================================================
IDX_RSTAR = 0     # r_t^{**}
IDX_EH = 1        # e_{ht} (household private signal noise)
IDX_EC = 2        # e_{ct} (central bank private signal noise)
IDX_F = 3         # f_t (public signal noise)
IDX_UH = 4        # u_{ht} (demand shock)
IDX_UP = 5        # u_{pt} (cost-push shock)
IDX_UC = 6        # u_{ct} (monetary policy shock)

# =============================================================
# Transition matrix R (Eq 5.6)
# Z_t = R Z_{t-1} + q_t
# =============================================================
R = np.diag([rho_r, rho_eh, rho_ec, rho_f, rho_h, rho_p, rho_c])

# =============================================================
# Innovation covariance Sigma_q (Eq 5.6)
# q_t ~ N(0, Sigma_q)
# =============================================================
Sigma_q = np.diag([
    sigma_v**2,
    sigma_eps_h**2,
    sigma_eps_c**2,
    sigma_eta**2,
    sigma_uh**2,
    sigma_up**2,
    sigma_uc**2,
])

# =============================================================
# Selection / embedding maps (Eq C.3-C.4)
# X_it = (Z_t', m_jt')' where m_jt in R^{DIM_M}
# =============================================================
def get_C_z():
    """C_z embeds Z_t into X_it: X_it = C_z Z_t + C_m m_jt."""
    C = np.zeros((DIM_X, N_Z))
    C[:N_Z, :N_Z] = np.eye(N_Z)
    return C

def get_C_m():
    """C_m embeds m_jt into X_it."""
    C = np.zeros((DIM_X, DIM_M))
    C[N_Z:, :] = np.eye(DIM_M)
    return C

def get_D_z():
    """D_z extracts Z_t from X_it."""
    D = np.zeros((N_Z, DIM_X))
    D[:N_Z, :N_Z] = np.eye(N_Z)
    return D

def get_D_m():
    """D_m extracts m_jt from X_it.

    Uses the approximation E_it^{(N+1)}[Z_t] ~ E_it^{(N)}[Z_t]
    for the truncation.
    """
    D = np.zeros((DIM_M, DIM_X))
    D[:, N_Z:] = np.eye(DIM_M)
    return D


# =============================================================
# Signal structure helpers
# =============================================================
def get_signal_h(Z):
    """Household private signal: s_ht = r** + e_ht (Eq 5.2)."""
    return Z[IDX_RSTAR] + Z[IDX_EH]

def get_signal_c(Z):
    """Central bank private signal: s_ct = r** + e_ct (Eq 5.2)."""
    return Z[IDX_RSTAR] + Z[IDX_EC]

def get_public_signal(Z):
    """Public signal: x_t = r** + f_t (Eq 5.4)."""
    return Z[IDX_RSTAR] + Z[IDX_F]
