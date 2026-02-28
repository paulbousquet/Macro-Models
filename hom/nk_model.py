"""
NK model solver for the Hall of Mirrors model.

Solves the 3-equation New Keynesian model (IS curve, Phillips curve, Taylor rule)
conditional on exogenous belief hierarchies, following Rungcharoenkitkul & Winkler.

Key equations:
  IS curve   (Eq 3.4): E_t^h[Dy_{t+1}] = (1/sigma)(i_t - E_t^h[pi_{t+1}] - r_t* - u_ht)
  Phillips   (Eq 3.5): pi_t = beta E_t^h[pi_{t+1}] + kappa y_t + u_pt
  Taylor     (Eq 3.6): i_t = rho_i i_{t-1} + (1 - rho_i)(r_hat_t* + phi_pi pi_t + phi_y y_t + u_ct)

Solution form (Eq C.1):
  (y_t, pi_t, i_t)' = Gamma i_{t-1} + Theta E_t^h[Z_t]
                       + sum_{s>=0} theta_s E_t^h[E_{t+s}^c[Z_{t+s}]]

References: Appendix C (Equations C.1, C.7, C.9)
"""

import numpy as np
from params import (
    sigma, kappa, beta, phi_pi, phi_y, rho_i,
    N_Z, IDX_RSTAR, IDX_UH, IDX_UP, IDX_UC,
    R,
)


def solve_nk_system():
    """Solve the 3-equation NK model for Gamma, Theta, and theta_s coefficients.

    Substitutes the Taylor rule into the IS curve, producing a 2x2 simultaneous
    system for (y_tilde, pi). Forward-looking household expectations (E^h[y_{t+1}],
    E^h[pi_{t+1}]) generate a Sylvester equation for Theta and a simple matrix
    recursion for the theta_s sequence.

    Returns:
        Gamma: ndarray (3,) -- response to lagged interest rate i_{t-1}
        Theta: ndarray (3, N_Z) -- response to E_t^h[Z_t]
        theta_s_func: callable(s) -> ndarray (3, N_Z) -- response to
            E_t^h[E_{t+s}^c[Z_{t+s}]] for the s-th forward CB belief term
    """
    om = 1.0 - rho_i

    # --- Gamma: response to i_{t-1} (Eq C.1) ---
    #
    # The nonlinear system arises because E^h[y_{t+1}] = Gamma_y * i_t and
    # E^h[pi_{t+1}] = Gamma_pi * i_t feed back through the IS curve:
    #   gamma_y = gamma_i * (gamma_y + gamma_pi/sigma - 1/sigma)
    #   gamma_pi = beta * gamma_pi * gamma_i + kappa * gamma_y
    #   gamma_i = rho_i + om * (phi_pi * gamma_pi + phi_y * gamma_y)
    #
    # Solve iteratively, seeded from the static (no-feedback) solution.
    A_static = np.array([
        [1.0 + om * phi_y / sigma, om * phi_pi / sigma],
        [-kappa, 1.0],
    ])
    ypi_seed = np.linalg.solve(A_static, np.array([-rho_i / sigma, 0.0]))
    gy, gpi = ypi_seed[0], ypi_seed[1]
    gi = rho_i + om * (phi_pi * gpi + phi_y * gy)

    for _ in range(500):
        gpi_new = kappa * gy / (1.0 - beta * gi)
        gi_new = rho_i + om * (phi_pi * gpi_new + phi_y * gy)
        coeff_y = 1.0 + om * phi_y / sigma - gi_new
        coeff_pi = om * phi_pi / sigma - gi_new / sigma
        gy_new = (-rho_i / sigma - coeff_pi * gpi_new) / coeff_y
        if max(abs(gy_new - gy), abs(gpi_new - gpi), abs(gi_new - gi)) < 1e-14:
            gy, gpi, gi = gy_new, gpi_new, gi_new
            break
        gy, gpi, gi = gy_new, gpi_new, gi_new

    Gamma = np.array([gy, gpi, gi])

    # --- Auxiliary quantities for forward-looking terms ---
    # g: combined coefficient from E^h[y_{t+1}] + (1/sigma)*E^h[pi_{t+1}]
    #    acting on i_t through Gamma
    g = Gamma[0] + Gamma[1] / sigma
    bgo = beta * Gamma[1] * om

    # M_left: coefficients on (theta_y, theta_pi) after absorbing the
    # Gamma*i_t feedback from forward expectations
    M_left = np.array([
        [1.0 + om * phi_y * (1.0 / sigma - g), om * phi_pi * (1.0 / sigma - g)],
        [-kappa - bgo * phi_y, 1.0 - bgo * phi_pi],
    ])

    # M_fwd: how theta_{s-1} feeds into the s-th equation through
    # E^h[y_{t+1}] and E^h[pi_{t+1}]
    M_fwd = np.array([
        [1.0, 1.0 / sigma],
        [0.0, beta],
    ])

    # --- Theta: response to E_t^h[Z_t] (Sylvester equation, Eq C.1) ---
    #
    # The equation M_left @ Theta_ypi - M_fwd @ Theta_ypi @ R = F arises because:
    # - M_left captures the simultaneous (y,pi) determination including Gamma feedback
    # - M_fwd @ Theta_ypi @ R captures E^h[Z_{t+1}] = R*E^h[Z_t] flowing through
    #   the forward-looking IS and PC expectations
    # - F contains direct forcing from household observations (r**, u_h, u_p)
    #
    # Vectorized as: (I_NZ kron M_left - R' kron M_fwd) vec(Theta_ypi) = vec(F)
    I_nz = np.eye(N_Z)
    big_A = np.kron(I_nz, M_left) - np.kron(R.T, M_fwd)

    F_theta = np.zeros((2, N_Z))
    F_theta[0, IDX_RSTAR] = 1.0 / sigma  # E^h[r**_t] enters IS
    F_theta[0, IDX_UH] = 1.0 / sigma     # u_ht enters IS
    F_theta[1, IDX_UP] = 1.0             # u_pt enters PC

    vec_Theta_ypi = np.linalg.solve(big_A, F_theta.ravel(order='F'))
    Theta_ypi = vec_Theta_ypi.reshape((2, N_Z), order='F')

    Theta_y = Theta_ypi[0, :]
    Theta_pi = Theta_ypi[1, :]
    Theta_i = om * (phi_pi * Theta_pi + phi_y * Theta_y)
    Theta = np.vstack([Theta_y, Theta_pi, Theta_i])  # (3, N_Z)

    # --- theta_s: response to E^h_t[E^c_{t+s}[Z_{t+s}]] ---
    #
    # s=0: CB beliefs (r_hat* and u_c) enter the Taylor rule directly.
    # The forcing is: IS gets -(om/sigma) on r** and u_c columns; PC gets 0.
    # No forward coupling at s=0 (no theta_{-1} exists).
    #
    # s>=1: Forward expectations E^h[y_{t+1}], E^h[pi_{t+1}] contain theta_{s-1}.
    # The recursion is: M_left @ theta_s = M_fwd @ theta_{s-1}
    # Note: NO multiplication by R -- the beliefs E^c_{t+s}[Z_{t+s}] are indexed
    # by the SAME s in both the current and forward equations.
    #
    # Recursion: theta_s_ypi = P^s @ theta0_ypi where P = M_left^{-1} @ M_fwd

    F_cb = np.zeros((2, N_Z))
    F_cb[0, IDX_RSTAR] = -om / sigma  # -(om/sigma)*E^c[r**] in IS
    F_cb[0, IDX_UC] = -om / sigma     # -(om/sigma)*u_c in IS

    theta0_ypi = np.linalg.solve(M_left, F_cb)  # (2, N_Z)

    P = np.linalg.inv(M_left) @ M_fwd  # (2, 2) propagator

    def theta_s_func(s):
        """Compute theta_s coefficient matrix (3, N_Z) for the s-th CB belief term."""
        ypi = np.linalg.matrix_power(P, s) @ theta0_ypi  # (2, N_Z)
        th_y = ypi[0, :]
        th_pi = ypi[1, :]
        th_i = om * (phi_pi * th_pi + phi_y * th_y)
        if s == 0:
            # Direct Taylor rule: i_t includes (1-rho_i)*(E_c[r**] + u_c).
            # These enter the Taylor rule directly but were not captured by
            # the (y,pi) feedback channel above.
            th_i[IDX_RSTAR] += om
            th_i[IDX_UC] += om
        return np.vstack([th_y, th_pi, th_i])  # (3, N_Z)

    # --- FI response to u_c (needed for C.21) ---
    #
    # Under full information, E^h[Z_{t+s}] = R^s Z_t, so the u_c column
    # of E^h[E^c_{t+s}[Z_{t+s}]] contributes rho_c^s * u_c_t to theta_s.
    # The geometric series sums to: (I - rho_c * P)^{-1} @ theta0_ypi[:, IDX_UC]
    # Theta[:, IDX_UC] = 0 (u_c has no direct IS/PC forcing), so FI = sum only.
    from params import rho_c
    FI_ypi_uc = np.linalg.solve(np.eye(2) - rho_c * P, theta0_ypi[:, IDX_UC])
    FI_i_uc = om * (phi_pi * FI_ypi_uc[1] + phi_y * FI_ypi_uc[0]) + om
    FI_uc = np.array([FI_ypi_uc[0], FI_ypi_uc[1], FI_i_uc])

    return Gamma, Theta, theta_s_func, FI_uc


def compute_M_h_direct(Gamma, F_m):
    """Compute M_h directly from the NK equations and belief dynamics.

    Solves the IS-PC-TR system for M_h (3, DIM_M):
      Y_t = Gamma * i_{t-1} + M_h @ m_ht

    where the household's forward expectations use:
      E_h[m_{h,t+1}] = F_m @ m_ht

    The Taylor rule uses:
      - E_h[E_c[r**]] at position N_Z + IDX_RSTAR (second-order belief)
      - E_h[u_c] at position IDX_UC (first-order belief, since CB knows u_c)

    This bypasses the theta_s power series entirely, avoiding the D_z/D_z_cb
    extraction issue that conflates belief orders.

    Args:
        Gamma: (3,) response to lagged i_{t-1} (from solve_nk_system)
        F_m: (DIM_M, DIM_M) belief transition matrix:
             F_m = Phi_h + Psi_h @ D_m_embed + Omega_h @ R @ D_z_embed

    Returns:
        M_h: (3, DIM_M) mapping from m_ht to macro outcomes (y, pi, i)
    """
    from params import DIM_M, N_Z

    om = 1.0 - rho_i
    sig = sigma
    g0, g1 = Gamma[0], Gamma[1]

    I_M = np.eye(DIM_M)

    # Unit row vectors in R^DIM_M
    def e_row(j):
        v = np.zeros((1, DIM_M))
        v[0, j] = 1.0
        return v

    # Taylor rule forcing: om * (E_h[E_c[r**]] + E_h[u_c])
    # E_h[E_c[r**]] is at position N_Z + IDX_RSTAR (second-order belief)
    # E_h[u_c] is at position IDX_UC (first-order belief, CB knows u_c directly)
    taylor_forcing = e_row(N_Z + IDX_RSTAR) + e_row(IDX_UC)

    # IS forcing: (1/sigma) * (E_h[r**] + E_h[u_h])
    is_z_forcing = e_row(IDX_RSTAR) / sig + e_row(IDX_UH) / sig

    # PC forcing: E_h[u_p]
    pc_forcing = e_row(IDX_UP)

    # Scalar: alpha = (Gamma[0] - (1 - Gamma[1])/sigma) * om
    alpha = (g0 - (1.0 - g1) / sig) * om
    bgo = beta * g1 * om

    # IS: row_y @ L_IS = row_pi @ R_IS + f_IS
    L_IS = I_M - F_m - alpha * phi_y * I_M
    R_IS = alpha * phi_pi * I_M + (1.0 / sig) * F_m

    # f_IS: Taylor rule forcing through IS + direct IS forcing
    f_IS = alpha * taylor_forcing + is_z_forcing

    # PC: row_pi @ L_PC = row_y @ R_PC + f_PC
    L_PC = I_M - beta * F_m - bgo * phi_pi * I_M
    R_PC = (kappa + bgo * phi_y) * I_M

    # f_PC: Taylor rule forcing through PC + direct PC forcing
    f_PC = bgo * taylor_forcing + pc_forcing

    # Solve: substitute PC into IS
    L_PC_inv = np.linalg.inv(L_PC)
    Xi = L_IS - R_PC @ L_PC_inv @ R_IS

    rhs_y = f_PC @ L_PC_inv @ R_IS + f_IS
    row_y = rhs_y @ np.linalg.inv(Xi)

    row_pi = (row_y @ R_PC + f_PC) @ L_PC_inv
    row_i = om * (phi_pi * row_pi + phi_y * row_y + taylor_forcing)

    return np.vstack([row_y, row_pi, row_i])  # (3, DIM_M)


def compute_M_h(Gamma, Theta, theta_s_func, A_Xh_plus_B_mh, D_z, D_z_cb=None,
                max_s=500, tol=1e-12):
    """Compute M_h matrix mapping household beliefs to macro outcomes (Eq C.7).

    M_h = Theta @ D_z + sum_{s=0}^{inf} theta_s @ D_z_cb @ (A_Xh + B_mh)^s

    The Theta term uses D_z (positions 0:N_Z) for E^h[Z_t].
    The theta_s terms use D_z_cb, a hybrid extraction map that accounts for
    which Z components the CB knows vs. estimates:
      - Components unknown to CB (e.g., r**): extract from positions N_Z:2*N_Z
        (CB's first-order Z beliefs within the HH state).
      - Components known to CB (e.g., u_c): extract from positions 0:N_Z
        (first-order Z beliefs), since E_c[u_c] = u_c identically.

    Args:
        Gamma: (3,) response to lagged i_{t-1}
        Theta: (3, N_Z) response to E_t^h[Z_t]
        theta_s_func: callable(s) -> (3, N_Z)
        A_Xh_plus_B_mh: (dim_state, dim_state) household belief transition
        D_z: (N_Z, dim_state) extraction for E^h[Z_t] (positions 0:N_Z)
        D_z_cb: (N_Z, dim_state) hybrid extraction for theta_s terms
        max_s: maximum summation terms
        tol: convergence tolerance on individual terms

    Returns:
        M_h: (3, dim_state) mapping from m_ht to macro outcomes
    """
    dim_state = A_Xh_plus_B_mh.shape[0]

    M_h = Theta @ D_z  # (3, dim_state)

    # D_z_cb extracts CB's Z beliefs (positions N_Z:2*N_Z); if not provided,
    # fall back to D_z (original paper formula, positions 0:N_Z).
    D_cb = D_z_cb if D_z_cb is not None else D_z

    power = np.eye(dim_state)
    prev_norm = np.inf
    for s in range(max_s):
        theta_s = theta_s_func(s)  # (3, N_Z)
        contribution = theta_s @ D_cb @ power  # (3, dim_state)
        cur_norm = np.max(np.abs(contribution))

        # During CK iteration, (A_Xh + B_mh) may temporarily have spectral
        # radius > 1 in the m_c block.  If the contribution norm grows for
        # two consecutive steps, the series is diverging -- truncate early
        # to keep M_h bounded.  Once the equilibrium converges the spectral
        # radius drops below 1 and the full sum converges normally.
        if s > 0 and cur_norm > prev_norm and cur_norm > 1.0:
            break

        M_h += contribution
        if cur_norm < tol:
            break
        prev_norm = cur_norm
        power = power @ A_Xh_plus_B_mh

    return M_h


def compute_expected_macro(Gamma, M_h, A_Xh_plus_B_mh, m_ht, horizon=40):
    """Compute expected future macro outcomes under household beliefs (Eq C.9).

    E_t^h(y_{t+s}, pi_{t+s}, i_{t+s})' = Gamma * E_t^h[i_{t+s-1}]
                                           + M_h @ (A_Xh + B_mh)^s @ m_ht

    The interest rate expectation evolves recursively:
      E_t^h[i_{t+s}] = Gamma[2]*E_t^h[i_{t+s-1}] + (M_h @ (A+B)^s @ m_ht)[2]

    Args:
        Gamma: (3,) response to lagged interest rate
        M_h: (3, dim_state) belief-to-macro mapping
        A_Xh_plus_B_mh: (dim_state, dim_state) household belief transition
        m_ht: (dim_state,) current household belief vector
        horizon: number of future periods

    Returns:
        y_tilde: (horizon,) expected output gap path
        pi: (horizon,) expected inflation path
        i: (horizon,) expected interest rate path
    """
    y_tilde = np.zeros(horizon)
    pi = np.zeros(horizon)
    i_path = np.zeros(horizon)

    power_m = m_ht.copy()
    E_i_prev = 0.0  # caller can add Gamma*i_{t-1} externally

    for s in range(horizon):
        macro_s = Gamma * E_i_prev + M_h @ power_m  # (3,)  Eq (C.9)
        y_tilde[s] = macro_s[0]
        pi[s] = macro_s[1]
        i_path[s] = macro_s[2]

        E_i_prev = i_path[s]
        power_m = A_Xh_plus_B_mh @ power_m

    return y_tilde, pi, i_path
