"""
Belief system and Kalman filter for the Hall of Mirrors model.

Implements the higher-order belief hierarchy, state space construction,
observation equations, and Kalman filter for the two-agent learning model.

Paper: "The Natural Rate of Interest Through a Hall of Mirrors"
Authors: Rungcharoenkitkul & Winkler
References: Section 5.2 (Eqs 5.7-5.11), Appendix C (Eqs C.2-C.26)

Dimension conventions:
  N_Z    = 7   -- exogenous state Z_t
  DIM_M  = 84  -- belief vector m_it (= N_Z * N_BELIEFS)
  DIM_X  = 91  -- full state X_it = (Z_t', m_jt')' (= N_Z * (N_BELIEFS+1))

The conjecture (5.9) uses DIM_M-dimensional objects:
  m_jt = Phi_j @ m_{jt-1} + Psi_j @ m_{it-1} + Omega_j @ Z_t
where Phi_j is (DIM_M, DIM_M), Psi_j is (DIM_M, DIM_M), Omega_j is (DIM_M, N_Z).

The Kalman filter operates in DIM_X. Selection maps C_z, C_m, D_z, D_m
convert between DIM_M and DIM_X representations. Key identity: D_m @ C_m = I_{DIM_M}.
"""

import numpy as np
import params


# ======================================================================
# Selection / embedding maps (cached at module level)
# ======================================================================
C_z = params.get_C_z()   # (DIM_X, N_Z): embeds Z_t into X_it
C_m = params.get_C_m()   # (DIM_X, DIM_M): embeds m_jt into X_it
D_z = params.get_D_z()   # (N_Z, DIM_X): extracts Z_t from X_it
D_m = params.get_D_m()   # (DIM_M, DIM_X): extracts m_jt from X_it

# m_i extraction/embedding: m_i = X̂[0:DIM_M] (first-order beliefs + higher-order)
# D_mi extracts m_i from X̂; C_mi embeds m_i into X̂ (with truncation for last block)
D_mi = np.zeros((params.DIM_M, params.DIM_X))
D_mi[:params.DIM_M, :params.DIM_M] = np.eye(params.DIM_M)

C_mi = np.zeros((params.DIM_X, params.DIM_M))
C_mi[:params.DIM_M, :params.DIM_M] = np.eye(params.DIM_M)
C_mi[params.DIM_M:, params.DIM_M - params.N_Z:] = np.eye(params.N_Z)  # truncation


# ======================================================================
# State space construction (Eq C.5)
# ======================================================================

def build_state_space(Phi_j, Psi_j, Omega_j):
    """Build state space matrices from conjectured belief law of motion.

    Given conjecture m_jt = Phi_j m_{jt-1} + Psi_j m_{it-1} + Omega_j Z_t,
    constructs the state space:
        X_it = A_Xi X_{it-1} + B_qi q_t + B_mi m_{it-1}  (Eq C.5)

    where:
        A_Xi  = (C_z + C_m Omega_j) R D_z + C_m Phi_j D_m
        B_qi  = C_z + C_m Omega_j
        B_mi  = C_m Psi_j D_m  (padded: maps DIM_X -> DIM_X)

    Args:
        Phi_j: (DIM_M, DIM_M) own-lag coefficient for agent j.
        Psi_j: (DIM_M, DIM_M) cross-lag coefficient for agent j.
        Omega_j: (DIM_M, N_Z) shock loading for agent j.

    Returns:
        A_Xi: (DIM_X, DIM_X) state transition.
        B_qi: (DIM_X, N_Z) shock loading.
        B_mi: (DIM_X, DIM_X) belief loading (padded).
    """
    R = params.R

    CzpCmOmega = C_z + C_m @ Omega_j                      # (DIM_X, N_Z)
    A_Xi = CzpCmOmega @ R @ D_z + C_m @ Phi_j @ D_m       # (DIM_X, DIM_X)
    B_qi = CzpCmOmega                                      # (DIM_X, N_Z)
    B_mi = C_m @ Psi_j @ D_mi                              # (DIM_X, DIM_X)

    return A_Xi, B_qi, B_mi


# ======================================================================
# Observation equations
# ======================================================================

def build_H_h():
    """Build household observation equation (Eq C.6).

    The household observes 5 signals:
      0) s_ht = r** + e_ht          (private signal)
      1) x_t  = r** + f_t           (public signal)
      2) u_ht                        (demand shock)
      3) u_pt                        (cost-push shock)
      4) E_t^c[r**] + u_ct          (inferred from Taylor rule)

    H_h is exogenous (does not depend on equilibrium beliefs).

    Returns:
        H_h: (5, DIM_X) observation matrix.
    """
    N_Z = params.N_Z
    N_OBS = 5

    # Z_t component: maps Z_t = (r**, e_h, e_c, f, u_h, u_p, u_c) to obs
    H_z = np.zeros((N_OBS, N_Z))
    H_z[0, params.IDX_RSTAR] = 1.0   # s_ht: r**
    H_z[0, params.IDX_EH] = 1.0      # s_ht: + e_ht
    H_z[1, params.IDX_RSTAR] = 1.0   # x_t: r**
    H_z[1, params.IDX_F] = 1.0       # x_t: + f_t
    H_z[2, params.IDX_UH] = 1.0      # u_ht
    H_z[3, params.IDX_UP] = 1.0      # u_pt
    H_z[4, params.IDX_UC] = 1.0      # u_ct (part of 5th signal)

    # m_ct (CB belief) component: only 5th signal depends on CB beliefs.
    # E_t^c[r**] is at position IDX_RSTAR in the first block of m_ct.
    H_m = np.zeros((N_OBS, params.DIM_M))
    H_m[4, params.IDX_RSTAR] = 1.0   # E_ct^(1)[r**]

    return H_z @ D_z + H_m @ D_m                           # (5, DIM_X)


def build_H_c(M_h):
    """Build central bank observation equation (Eq C.8).

    The central bank observes 5 signals:
      0) s_ct = r** + e_ct          (private signal)
      1) x_t  = r** + f_t           (public signal)
      2) u_ct                        (monetary policy shock)
      3) y_tilde_t                   (output gap, endogenous)
      4) pi_t                        (inflation, endogenous)

    The last two depend on M_h (household belief -> macro outcome mapping).

    Args:
        M_h: (3, DIM_M) matrix from Eq C.7. Rows: (y_tilde, pi, i).

    Returns:
        H_c: (5, DIM_X) observation matrix.
    """
    N_Z = params.N_Z
    N_OBS = 5

    H_z = np.zeros((N_OBS, N_Z))
    H_z[0, params.IDX_RSTAR] = 1.0   # s_ct: r**
    H_z[0, params.IDX_EC] = 1.0      # s_ct: + e_ct
    H_z[1, params.IDX_RSTAR] = 1.0   # x_t: r**
    H_z[1, params.IDX_F] = 1.0       # x_t: + f_t
    H_z[2, params.IDX_UC] = 1.0      # u_ct

    # Macro outcomes depend on HH beliefs through M_h @ D_m @ X_ct
    H_m = np.zeros((N_OBS, params.DIM_M))
    H_m[3, :] = M_h[0, :]            # y_tilde_t
    H_m[4, :] = M_h[1, :]            # pi_t

    return H_z @ D_z + H_m @ D_m                           # (5, DIM_X)


# ======================================================================
# HoM perceived observation equations (Eqs C.19-C.21)
# ======================================================================

def build_H_c_perceived_by_h():
    """CB obs matrix as perceived by household (Eq C.19).

    HoM assumption: household thinks CB ignores y_tilde and pi.
    CB only observes: s_ct, x_t, u_ct (3 exogenous signals).

    Returns:
        H: (3, DIM_X)
    """
    N_Z = params.N_Z
    H_z = np.zeros((3, N_Z))
    H_z[0, params.IDX_RSTAR] = 1.0
    H_z[0, params.IDX_EC] = 1.0
    H_z[1, params.IDX_RSTAR] = 1.0
    H_z[1, params.IDX_F] = 1.0
    H_z[2, params.IDX_UC] = 1.0
    return H_z @ D_z


def build_H_h_perceived_by_c():
    """HH obs matrix as perceived by CB (Eq C.20).

    HoM assumption: CB thinks household ignores the interest rate.
    HH only observes: s_ht, x_t, u_ht, u_pt (4 exogenous signals).

    Returns:
        H: (4, DIM_X)
    """
    N_Z = params.N_Z
    H_z = np.zeros((4, N_Z))
    H_z[0, params.IDX_RSTAR] = 1.0
    H_z[0, params.IDX_EH] = 1.0
    H_z[1, params.IDX_RSTAR] = 1.0
    H_z[1, params.IDX_F] = 1.0
    H_z[2, params.IDX_UH] = 1.0
    H_z[3, params.IDX_UP] = 1.0
    return H_z @ D_z


def build_H_c_own_hom(M_h_c, direct_uc_macro):
    """CB's own obs matrix under HoM (Eq C.21).

    The CB's perceived macro outcomes use M_h^{|c} (computed under the
    perceived HH info set without interest rate). The direct Z_t effect
    on y_tilde and pi from u_c uses the current NK-block loading on the
    household's first-order u_c belief.

    Args:
        M_h_c: (3, DIM_M) M_h under CB-perceived HH info set.
        direct_uc_macro: (2,) direct impact of a current household u_c belief
            on (y_tilde, pi) in the NK block.

    Returns:
        H: (5, DIM_X)
    """
    N_Z = params.N_Z
    N_OBS = 5

    H_z = np.zeros((N_OBS, N_Z))
    H_z[0, params.IDX_RSTAR] = 1.0
    H_z[0, params.IDX_EC] = 1.0
    H_z[1, params.IDX_RSTAR] = 1.0
    H_z[1, params.IDX_F] = 1.0
    H_z[2, params.IDX_UC] = 1.0
    H_z[3, params.IDX_UC] = direct_uc_macro[0]
    H_z[4, params.IDX_UC] = direct_uc_macro[1]

    H_m = np.zeros((N_OBS, params.DIM_M))
    H_m[3, :] = M_h_c[0, :]          # y_tilde belief-driven part
    H_m[4, :] = M_h_c[1, :]          # pi belief-driven part

    return H_z @ D_z + H_m @ D_m


# ======================================================================
# Kalman filter: iterative DARE solver (Eqs C.11-C.14)
# ======================================================================

def solve_kalman_dare(A_Xi, B_qi, H_i, Sigma_q,
                      tol=1e-10, max_iter=1000, damping=0.5):
    """Solve the discrete algebraic Riccati equation iteratively.

    Iterates Eqs C.11-C.14 to find steady-state posterior covariance P_i:
        P_i^- = A_Xi P_i A_Xi' + B_qi Sigma_q B_qi'     (C.11)
        S_i   = H_i P_i^- H_i'                           (C.12)
        G_i   = P_i^- H_i' S_i^{-1}                      (C.13)
        P_i   = P_i^- - G_i S_i G_i'                     (C.14)

    Note: The DARE uses A_Xi (not A_Xi+B_mi) because B_mi acts on the
    known own-belief, which does not contribute to estimation uncertainty.

    Args:
        A_Xi: (DIM_X, DIM_X) state transition.
        B_qi: (DIM_X, N_Z) shock loading.
        H_i: (n_obs, DIM_X) observation matrix.
        Sigma_q: (N_Z, N_Z) innovation covariance.
        tol: convergence tolerance.
        max_iter: maximum iterations.
        damping: update damping factor.

    Returns:
        G_i: (DIM_X, n_obs) Kalman gain.
        P_i: (DIM_X, DIM_X) steady-state posterior covariance.
    """
    BqSigBq = B_qi @ Sigma_q @ B_qi.T
    P_i = BqSigBq.copy()

    for _ in range(max_iter):
        # Eq C.11
        P_prior = A_Xi @ P_i @ A_Xi.T + BqSigBq
        # Eq C.12 (add small regularization for numerical stability)
        S_i = H_i @ P_prior @ H_i.T + np.eye(H_i.shape[0]) * 1e-14
        # Eq C.13 (use solve for better numerical stability)
        G_i = np.linalg.solve(S_i, H_i @ P_prior.T).T
        # Eq C.14
        P_new = P_prior - G_i @ S_i @ G_i.T
        P_new = 0.5 * (P_new + P_new.T)  # symmetrize

        diff = np.max(np.abs(P_new - P_i))
        if diff < tol:
            return G_i, P_new

        P_i = damping * P_new + (1.0 - damping) * P_i

    return G_i, P_i


# ======================================================================
# Equilibrium coefficient extraction (Eqs C.16-C.18)
# ======================================================================

def extract_coefficients(G_i, H_i, A_Xi, B_mi, Phi_j, Omega_j):
    """Extract DIM_M equilibrium coefficients from Kalman solution.

    Projects the DIM_X Kalman update onto DIM_M-dimensional conjecture
    coefficients using D_m (extraction) and C_m (embedding).

    From C.10 + C.15, the full DIM_X coefficient on m_{it-1}_full is:
        (I - G_i H_i) A_Xi + B_mi     (Eq C.16)
    Projecting to DIM_M: Phi_i = D_m @ [...] @ C_m

    The cross and shock coefficients come from the Y_it substitution:
        Psi_i = D_m @ G_i @ H_i @ C_m @ Phi_j   (Eq C.17)
        Omega_i = D_m @ G_i @ H_i @ (C_z + C_m @ Omega_j)  (Eq C.18)

    Args:
        G_i: (DIM_X, n_obs) Kalman gain.
        H_i: (n_obs, DIM_X) observation matrix.
        A_Xi: (DIM_X, DIM_X) state transition.
        B_mi: (DIM_X, DIM_X) belief loading (C_m Psi_j D_m).
        Phi_j: (DIM_M, DIM_M) other agent's own-lag.
        Omega_j: (DIM_M, N_Z) other agent's shock loading.

    Returns:
        Phi_i: (DIM_M, DIM_M) own-lag coefficient.
        Psi_i: (DIM_M, DIM_M) cross-lag coefficient.
        Omega_i: (DIM_M, N_Z) shock loading.
    """
    I_X = np.eye(params.DIM_X)
    IminGH = I_X - G_i @ H_i

    # Eq C.16: project to DIM_M using D_mi (own-agent extraction) and C_mi (own embedding)
    Phi_i = D_mi @ (IminGH @ A_Xi + B_mi) @ C_mi

    # Eq C.17: C_m stays (embeds m_j, the other agent)
    Psi_i = D_mi @ G_i @ H_i @ C_m @ Phi_j

    # Eq C.18: C_z and C_m stay (Z and m_j embedding)
    Omega_i = D_mi @ G_i @ H_i @ (C_z + C_m @ Omega_j)

    return Phi_i, Psi_i, Omega_i


# ======================================================================
# CK equilibrium solver
# ======================================================================

def solve_ck_equilibrium(M_h_func, tol=1e-8, max_iter=200, damping=0.5):
    """Solve for the Common Knowledge (CK) learning equilibrium.

    Iterates the algorithm from the end of Appendix C:
    1. Guess (Phi_i, Psi_i, Omega_i, H_i) for i=h,c
    2. Build state space from C.5
    3. Solve Kalman DARE (C.11-C.14)
    4. Extract new coefficients (C.16-C.18)
    5. Update H_c via C.8 (H_h is fixed)
    6. Iterate until convergence

    Args:
        M_h_func: Callable(Phi_h, Psi_h, Omega_h) -> M_h (3, DIM_M).
                  Computes HH belief -> macro outcome mapping directly
                  from HH Kalman filter coefficients.
        tol: convergence tolerance.
        max_iter: maximum iterations.
        damping: update damping.

    Returns:
        dict with equilibrium coefficients, Kalman gains, etc.
    """
    DIM_M = params.DIM_M
    N_Z = params.N_Z
    Sigma_q = params.Sigma_q

    # Initialize Phi with block-diagonal R to seed belief persistence.
    # With Phi=0 the equilibrium has a trivial fixed point where beliefs
    # have no memory (Phi_j=0 makes A_Xi @ C_m = 0, so extracted Phi_i=0
    # regardless of Kalman gains). Seeding with the natural persistence
    # R_block = kron(I_N, R) breaks this degeneracy.
    R_block = np.kron(np.eye(params.N_BELIEFS), params.R)  # (DIM_M, DIM_M)
    Phi_h = R_block.copy()
    Psi_h = np.zeros((DIM_M, DIM_M))
    Omega_h = 0.01 * np.tile(np.eye(N_Z), (params.N_BELIEFS, 1))
    Phi_c = R_block.copy()
    Psi_c = np.zeros((DIM_M, DIM_M))
    Omega_c = 0.01 * np.tile(np.eye(N_Z), (params.N_BELIEFS, 1))

    H_h = build_H_h()                # fixed
    M_h = np.zeros((3, DIM_M))
    H_c = build_H_c(M_h)

    diff = np.inf
    for it in range(max_iter):
        # --- Household ---
        A_Xh, B_qh, B_mh = build_state_space(Phi_c, Psi_c, Omega_c)
        G_h, P_h = solve_kalman_dare(A_Xh, B_qh, H_h, Sigma_q)
        Phi_h_new, Psi_h_new, Omega_h_new = extract_coefficients(
            G_h, H_h, A_Xh, B_mh, Phi_c, Omega_c)

        # --- Central bank ---
        A_Xc, B_qc, B_mc = build_state_space(Phi_h, Psi_h, Omega_h)
        G_c, P_c = solve_kalman_dare(A_Xc, B_qc, H_c, Sigma_q)
        Phi_c_new, Psi_c_new, Omega_c_new = extract_coefficients(
            G_c, H_c, A_Xc, B_mc, Phi_h, Omega_h)

        # --- Update H_c ---
        M_h_new = M_h_func(Phi_h_new, Psi_h_new, Omega_h_new)
        H_c_new = build_H_c(M_h_new)

        # --- Convergence ---
        diff = max(
            np.max(np.abs(Phi_h_new - Phi_h)),
            np.max(np.abs(Psi_h_new - Psi_h)),
            np.max(np.abs(Omega_h_new - Omega_h)),
            np.max(np.abs(Phi_c_new - Phi_c)),
            np.max(np.abs(Psi_c_new - Psi_c)),
            np.max(np.abs(Omega_c_new - Omega_c)),
        )

        # Damped update
        d = damping
        Phi_h = d * Phi_h_new + (1 - d) * Phi_h
        Psi_h = d * Psi_h_new + (1 - d) * Psi_h
        Omega_h = d * Omega_h_new + (1 - d) * Omega_h
        Phi_c = d * Phi_c_new + (1 - d) * Phi_c
        Psi_c = d * Psi_c_new + (1 - d) * Psi_c
        Omega_c = d * Omega_c_new + (1 - d) * Omega_c
        M_h = d * M_h_new + (1 - d) * M_h
        H_c = d * H_c_new + (1 - d) * H_c

        if diff < tol:
            break

    return {
        'Phi_h': Phi_h, 'Psi_h': Psi_h, 'Omega_h': Omega_h,
        'Phi_c': Phi_c, 'Psi_c': Psi_c, 'Omega_c': Omega_c,
        'G_h': G_h, 'G_c': G_c, 'P_h': P_h, 'P_c': P_c,
        'H_h': H_h, 'H_c': H_c,
        'A_Xh': A_Xh, 'A_Xc': A_Xc, 'B_mh': B_mh, 'B_mc': B_mc,
        'M_h': M_h,
        'converged': diff < tol,
        'iterations': it + 1,
        'final_diff': diff,
    }


# ======================================================================
# HoM equilibrium solver
# ======================================================================

def solve_hom_equilibrium(M_h_func, FI_uc, ck_solution=None,
                          tol=1e-8, max_iter=500, damping=0.3):
    """Solve for the Hall-of-Mirrors (HoM) equilibrium.

    Two-phase algorithm:
    Phase 1: Solve each agent's PERCEIVED filtering problem (with
             restricted signal sets for the other agent).
    Phase 2: Compute ACTUAL equilibrium law of motion (Eqs C.24-C.26)
             by combining perceived Kalman gains with actual signal matrices.

    Args:
        M_h_func: Callable(Phi_h, Psi_h, Omega_h) -> M_h (3, DIM_M).
        FI_uc: retained for compatibility with callers; phase 1 uses the
            current household-side u_c loading instead.
        ck_solution: Optional CK equilibrium dict to warm-start iteration.
        tol: convergence tolerance.
        max_iter: maximum iterations.
        damping: update damping.

    Returns:
        dict with perceived and actual equilibrium coefficients.
    """
    DIM_M = params.DIM_M
    DIM_X = params.DIM_X
    N_Z = params.N_Z
    Sigma_q = params.Sigma_q

    # Fixed observation matrices
    H_h_actual = build_H_h()                      # (5, DIM_X)
    H_c_perc_h = build_H_c_perceived_by_h()       # (3, DIM_X) -- Eq C.19
    H_h_perc_c = build_H_h_perceived_by_c()       # (4, DIM_X) -- Eq C.20

    # Initialize from CK solution (warm start) or block-diagonal R seed.
    z_mm = lambda: np.zeros((DIM_M, DIM_M))
    if ck_solution is not None:
        Phi_h_h = ck_solution['Phi_h'].copy()
        Psi_h_h = ck_solution['Psi_h'].copy()
        Omega_h_h = ck_solution['Omega_h'].copy()
        Phi_c_h = ck_solution['Phi_c'].copy()
        Psi_c_h = ck_solution['Psi_c'].copy()
        Omega_c_h = ck_solution['Omega_c'].copy()
        Phi_h_c = ck_solution['Phi_h'].copy()
        Psi_h_c = ck_solution['Psi_h'].copy()
        Omega_h_c = ck_solution['Omega_h'].copy()
        Phi_c_c = ck_solution['Phi_c'].copy()
        Psi_c_c = ck_solution['Psi_c'].copy()
        Omega_c_c = ck_solution['Omega_c'].copy()
        M_h_c = ck_solution['M_h'].copy()
    else:
        R_block = np.kron(np.eye(params.N_BELIEFS), params.R)
        omega_seed = 0.01 * np.tile(np.eye(N_Z), (params.N_BELIEFS, 1))
        Phi_h_h, Psi_h_h, Omega_h_h = R_block.copy(), z_mm(), omega_seed.copy()
        Phi_c_h, Psi_c_h, Omega_c_h = R_block.copy(), z_mm(), omega_seed.copy()
        Phi_h_c, Psi_h_c, Omega_h_c = R_block.copy(), z_mm(), omega_seed.copy()
        Phi_c_c, Psi_c_c, Omega_c_c = R_block.copy(), z_mm(), omega_seed.copy()
        M_h_c = np.zeros((3, DIM_M))

    diff = np.inf
    best_diff = np.inf
    best_state = None
    d = damping

    for it in range(max_iter):
        # ---- Household's perceived world ----
        # HH uses full signals; CB perceived to use restricted signals (Eq C.19)
        A_Xh_h, B_qh_h, B_mh_h = build_state_space(Phi_c_h, Psi_c_h, Omega_c_h)
        G_h_h, P_h_h = solve_kalman_dare(A_Xh_h, B_qh_h, H_h_actual, Sigma_q)
        Phi_h_h_n, Psi_h_h_n, Omega_h_h_n = extract_coefficients(
            G_h_h, H_h_actual, A_Xh_h, B_mh_h, Phi_c_h, Omega_c_h)

        A_Xc_h, B_qc_h, B_mc_h = build_state_space(Phi_h_h, Psi_h_h, Omega_h_h)
        G_c_h, _ = solve_kalman_dare(A_Xc_h, B_qc_h, H_c_perc_h, Sigma_q)
        Phi_c_h_n, Psi_c_h_n, Omega_c_h_n = extract_coefficients(
            G_c_h, H_c_perc_h, A_Xc_h, B_mc_h, Phi_h_h, Omega_h_h)

        # ---- CB's perceived world ----
        # HH perceived to use restricted signals (Eq C.20); CB uses own obs (Eq C.21)
        A_Xh_c, B_qh_c, B_mh_c = build_state_space(Phi_c_c, Psi_c_c, Omega_c_c)
        G_h_c, _ = solve_kalman_dare(A_Xh_c, B_qh_c, H_h_perc_c, Sigma_q)
        Phi_h_c_n, Psi_h_c_n, Omega_h_c_n = extract_coefficients(
            G_h_c, H_h_perc_c, A_Xh_c, B_mh_c, Phi_c_c, Omega_c_c)

        M_h_c_new = M_h_func(Phi_h_c, Psi_h_c, Omega_h_c)
        M_h_h_cur = M_h_func(Phi_h_h, Psi_h_h, Omega_h_h)
        H_c_c = build_H_c_own_hom(M_h_c_new, M_h_h_cur[:2, params.IDX_UC])

        A_Xc_c, B_qc_c, B_mc_c = build_state_space(Phi_h_c, Psi_h_c, Omega_h_c)
        G_c_c, P_c_c = solve_kalman_dare(A_Xc_c, B_qc_c, H_c_c, Sigma_q)
        Phi_c_c_n, Psi_c_c_n, Omega_c_c_n = extract_coefficients(
            G_c_c, H_c_c, A_Xc_c, B_mc_c, Phi_h_c, Omega_h_c)

        # ---- Convergence ----
        diff = max(
            np.max(np.abs(Phi_h_h_n - Phi_h_h)),
            np.max(np.abs(Phi_c_h_n - Phi_c_h)),
            np.max(np.abs(Phi_h_c_n - Phi_h_c)),
            np.max(np.abs(Phi_c_c_n - Phi_c_c)),
            np.max(np.abs(Omega_h_h_n - Omega_h_h)),
            np.max(np.abs(Omega_c_h_n - Omega_c_h)),
            np.max(np.abs(Omega_h_c_n - Omega_h_c)),
            np.max(np.abs(Omega_c_c_n - Omega_c_c)),
        )

        for old, new in [
            (Phi_h_h, Phi_h_h_n), (Psi_h_h, Psi_h_h_n), (Omega_h_h, Omega_h_h_n),
            (Phi_c_h, Phi_c_h_n), (Psi_c_h, Psi_c_h_n), (Omega_c_h, Omega_c_h_n),
            (Phi_h_c, Phi_h_c_n), (Psi_h_c, Psi_h_c_n), (Omega_h_c, Omega_h_c_n),
            (Phi_c_c, Phi_c_c_n), (Psi_c_c, Psi_c_c_n), (Omega_c_c, Omega_c_c_n),
        ]:
            old[:] = d * new + (1 - d) * old
        M_h_c = d * M_h_c_new + (1 - d) * M_h_c

        # Track best iterate and detect divergence
        if diff < best_diff:
            best_diff = diff
            best_state = {
                'h_h': (Phi_h_h.copy(), Psi_h_h.copy(), Omega_h_h.copy()),
                'c_h': (Phi_c_h.copy(), Psi_c_h.copy(), Omega_c_h.copy()),
                'h_c': (Phi_h_c.copy(), Psi_h_c.copy(), Omega_h_c.copy()),
                'c_c': (Phi_c_c.copy(), Psi_c_c.copy(), Omega_c_c.copy()),
                'M_h_c': M_h_c.copy(),
                'G_h_h': G_h_h.copy(), 'G_c_c': G_c_c.copy(),
                'P_h_h': P_h_h.copy(), 'P_c_c': P_c_c.copy(),
                'iter': it,
            }
        elif diff > best_diff * 100:
            # Diverging: restore best state and stop
            Phi_h_h, Psi_h_h, Omega_h_h = best_state['h_h']
            Phi_c_h, Psi_c_h, Omega_c_h = best_state['c_h']
            Phi_h_c, Psi_h_c, Omega_h_c = best_state['h_c']
            Phi_c_c, Psi_c_c, Omega_c_c = best_state['c_c']
            M_h_c = best_state['M_h_c']
            G_h_h = best_state['G_h_h']
            G_c_c = best_state['G_c_c']
            P_h_h = best_state['P_h_h']
            P_c_c = best_state['P_c_c']
            diff = best_diff
            # Rebuild state spaces from restored coefficients
            A_Xh_h, B_qh_h, B_mh_h = build_state_space(Phi_c_h, Psi_c_h, Omega_c_h)
            A_Xc_c, B_qc_c, B_mc_c = build_state_space(Phi_h_c, Psi_h_c, Omega_h_c)
            break

        if diff < tol:
            break

    # Recompute the phase-1 objects at the final fixed point so phase 2 uses
    # internally consistent gains and observation matrices.
    A_Xh_h, B_qh_h, B_mh_h = build_state_space(Phi_c_h, Psi_c_h, Omega_c_h)
    G_h_h, P_h_h = solve_kalman_dare(A_Xh_h, B_qh_h, H_h_actual, Sigma_q)

    A_Xc_h, B_qc_h, B_mc_h = build_state_space(Phi_h_h, Psi_h_h, Omega_h_h)
    G_c_h, _ = solve_kalman_dare(A_Xc_h, B_qc_h, H_c_perc_h, Sigma_q)

    A_Xh_c, B_qh_c, B_mh_c = build_state_space(Phi_c_c, Psi_c_c, Omega_c_c)
    G_h_c, _ = solve_kalman_dare(A_Xh_c, B_qh_c, H_h_perc_c, Sigma_q)

    M_h_actual = M_h_func(Phi_h_h, Psi_h_h, Omega_h_h)
    M_h_c = M_h_func(Phi_h_c, Psi_h_c, Omega_h_c)
    H_c_c = build_H_c_own_hom(M_h_c, M_h_actual[:2, params.IDX_UC])

    A_Xc_c, B_qc_c, B_mc_c = build_state_space(Phi_h_c, Psi_h_c, Omega_h_c)
    G_c_c, P_c_c = solve_kalman_dare(A_Xc_c, B_qc_c, H_c_c, Sigma_q)

    # ---- Phase 2: Actual equilibrium (Eqs C.24-C.26) ----
    # The paper's phase-2 recursion uses the phase-1 Kalman gains G_i^{|i},
    # own signal matrices H_i^{|i}, and cross signal matrices H_i^{|j}.
    # In this implementation:
    # - H_h^{|h} and H_h^{|c} coincide with the household's actual signal map.
    # - H_c^{|c} is the CB's own HoM observation map from Eq. C.21.
    # - H_c^{|h} uses the macro mapping implied by the household's perceived world.
    H_h_h = H_h_actual
    H_h_c = H_h_actual
    H_c_h = build_H_c(M_h_actual)

    I_X = np.eye(DIM_X)
    I_M = np.eye(DIM_M)

    own_h = Phi_h_h - D_mi @ G_h_h @ H_h_h @ C_m @ Psi_c_h
    own_c = Phi_c_c - D_mi @ G_c_c @ H_c_c @ C_m @ Psi_h_c

    cross_h = D_mi @ G_h_h @ H_h_c @ C_m
    cross_c = D_mi @ G_c_c @ H_c_h @ C_m

    inv_h = np.linalg.inv(I_M - cross_h @ cross_c)
    inv_c = np.linalg.inv(I_M - cross_c @ cross_h)

    # Eq C.24
    Phi_h_act = inv_h @ own_h
    Phi_c_act = inv_c @ own_c

    # Eq C.25
    Psi_h_act = inv_h @ cross_h @ own_c
    Psi_c_act = inv_c @ cross_c @ own_h

    # Eq C.26
    Omega_h_act = inv_h @ D_mi @ G_h_h @ H_h_c @ (
        I_X + C_m @ D_mi @ G_c_c @ H_c_h) @ C_z
    Omega_c_act = inv_c @ D_mi @ G_c_c @ H_c_h @ (
        I_X + C_m @ D_mi @ G_h_h @ H_h_c) @ C_z

    return {
        'perceived': {
            'h': {'Phi': Phi_h_h, 'Psi': Psi_h_h, 'Omega': Omega_h_h,
                   'G': G_h_h, 'P': P_h_h, 'A_X': A_Xh_h, 'B_m': B_mh_h},
            'c': {'Phi': Phi_c_c, 'Psi': Psi_c_c, 'Omega': Omega_c_c,
                   'G': G_c_c, 'P': P_c_c, 'A_X': A_Xc_c, 'B_m': B_mc_c},
            'h_sees_c': {'Phi': Phi_c_h, 'Psi': Psi_c_h, 'Omega': Omega_c_h},
            'c_sees_h': {'Phi': Phi_h_c, 'Psi': Psi_h_c, 'Omega': Omega_h_c},
        },
        'actual': {
            'Phi_h': Phi_h_act, 'Psi_h': Psi_h_act, 'Omega_h': Omega_h_act,
            'Phi_c': Phi_c_act, 'Psi_c': Psi_c_act, 'Omega_c': Omega_c_act,
        },
        'debug': {
            'H_h_actual': H_h_actual,
            'H_c_actual': H_c_h,
            'H_c_perc_h': H_c_perc_h,
            'H_h_perc_c': H_h_perc_c,
            'H_h_h': H_h_h,
            'H_h_c': H_h_c,
            'H_c_h': H_c_h,
            'H_c_c': H_c_c,
            'G_h_h': G_h_h,
            'G_c_h': G_c_h,
            'G_h_c': G_h_c,
            'G_c_c': G_c_c,
        },
        'H_h': H_h_actual, 'H_c': H_c_h, 'M_h': M_h_actual,
        'converged': diff < tol,
        'iterations': it + 1,
        'final_diff': diff,
    }
