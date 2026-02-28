"""
Driver script for the Hall of Mirrors model pipeline.

Solves the NK model, CK and HoM equilibria, computes IRFs for all shock types,
and generates plots matching Figures E.1-E.4 in Rungcharoenkitkul & Winkler.

Usage:
    cd /home/pblit/Music/Projects/cflow/test
    python3 solve_hom.py
"""

import os
import numpy as np
import matplotlib.pyplot as plt

import params
from nk_model import solve_nk_system, compute_M_h_direct
from belief_system import (
    solve_ck_equilibrium,
    solve_hom_equilibrium,
    C_z, C_m, D_z, D_m,
)

RESULTS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
PP = 1  # model variables are already in percentage points


# ======================================================================
# Main pipeline
# ======================================================================
def main():
    os.makedirs(RESULTS_DIR, exist_ok=True)

    # ==================================================================
    # Step 1: Solve NK model
    # ==================================================================
    print("Step 1: Solving NK model...")
    Gamma, Theta, theta_s_func, FI_uc = solve_nk_system()
    print(f"  Gamma = {Gamma}")
    print(f"  Theta shape = {Theta.shape}")
    print(f"  FI_uc = {FI_uc}  (should be non-zero)")

    # ==================================================================
    # Step 2: M_h callback (direct NK equation approach)
    # ==================================================================
    N_Z = params.N_Z
    N_B = params.N_BELIEFS
    DIM_M = params.DIM_M

    D_z_mat = params.get_D_z()  # (N_Z, DIM_X)
    D_m_mat = params.get_D_m()  # (DIM_M, DIM_X)

    # E_embed: (DIM_X, DIM_M) maps belief hierarchy m_ht to posterior hat{X}_{ht}
    D_z_m = np.zeros((N_Z, DIM_M))
    D_z_m[:N_Z, :N_Z] = np.eye(N_Z)

    T_shift = np.zeros((DIM_M, DIM_M))
    for b in range(N_B - 1):
        T_shift[b * N_Z:(b + 1) * N_Z, (b + 1) * N_Z:(b + 2) * N_Z] = np.eye(N_Z)
    T_shift[(N_B - 1) * N_Z:, (N_B - 1) * N_Z:] = np.eye(N_Z)

    E_embed = C_z @ D_z_m + C_m @ T_shift  # (DIM_X, DIM_M)

    # Precompute selection maps for F_m construction
    D_m_embed = D_m_mat @ E_embed  # (DIM_M, DIM_M): maps m_h to m_j via X̂_h
    D_z_embed = D_z_mat @ E_embed  # (N_Z, DIM_M): picks first-order Z beliefs

    def M_h_func(Phi_h, Psi_h, Omega_h):
        """Compute M_h (3, DIM_M) directly from the NK equations.

        Solves the IS-PC-TR system for M_h given the household's Kalman filter
        coefficients. The belief transition is:
            E_h[m_{h,t+1}] = F_m @ m_ht
        where F_m = Phi_h + Psi_h @ D_m_embed + Omega_h @ R @ D_z_embed.

        This correctly handles the mixed-order nature of the Taylor rule:
        - E_h[E_c[r**]] enters as a second-order belief (position N_Z+IDX_RSTAR)
        - E_h[u_c] enters as a first-order belief (position IDX_UC)
        """
        F_m = Phi_h + Psi_h @ D_m_embed + Omega_h @ params.R @ D_z_embed
        return compute_M_h_direct(Gamma, F_m)

    # ==================================================================
    # Step 3: Solve equilibria
    # ==================================================================
    print("\nStep 3a: Solving CK equilibrium...")
    ck = solve_ck_equilibrium(M_h_func)
    print(f"  CK converged: {ck['converged']} ({ck['iterations']} iters, "
          f"diff={ck['final_diff']:.2e})")

    print("\nStep 3b: Solving HoM equilibrium (warm-started from CK)...")
    hom = solve_hom_equilibrium(M_h_func, FI_uc, ck_solution=ck)
    print(f"  HoM converged: {hom['converged']} ({hom['iterations']} iters, "
          f"diff={hom['final_diff']:.2e})")

    # ==================================================================
    # Step 4: Compute IRFs
    # ==================================================================
    print("\nStep 4: Computing IRFs...")

    R = params.R

    # ------------------------------------------------------------------
    # IRF simulation functions
    # ------------------------------------------------------------------
    om_fi = 1.0 - params.rho_i
    sig = params.sigma

    # Block matrices for the stacked FI system (3 eqs per period)
    B_diag = np.array([
        [1.0 + om_fi * params.phi_y / sig, om_fi * params.phi_pi / sig, 0.0],
        [-params.kappa, 1.0, 0.0],
        [-om_fi * params.phi_y, -om_fi * params.phi_pi, 1.0],
    ])
    B_upper = np.array([
        [-1.0, -1.0 / sig, 0.0],
        [0.0, -params.beta, 0.0],
        [0.0, 0.0, 0.0],
    ])
    B_lower = np.array([
        [0.0, 0.0, params.rho_i / sig],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, -params.rho_i],
    ])

    def simulate_irf_fi(q0, horizon):
        """Full information IRF via stacked block-tridiagonal linear system.

        Solves the IS-PC-TR system simultaneously across all periods,
        giving the exact rational-expectations solution. Works correctly
        for both transitory (rho<1) and random-walk (rho=1) shocks.

        Returns: (rstar_true, rstar_true_copy, i_path, y_path, pi_path)
        Under FI, both agents' r* beliefs = true r**.
        """
        T = horizon
        n = 3  # (y, pi, i) per period

        # Build exogenous Z path: Z_t = R^t @ q0
        Z_path = np.zeros((T, N_Z))
        Z_path[0] = q0.copy()
        for t in range(1, T):
            Z_path[t] = R @ Z_path[t - 1]

        # Assemble 3T x 3T block-tridiagonal system
        A = np.zeros((n * T, n * T))
        b = np.zeros(n * T)

        for t in range(T):
            r = t * n
            A[r:r + n, r:r + n] = B_diag
            if t < T - 1:
                A[r:r + n, r + n:r + 2 * n] = B_upper
            if t > 0:
                A[r:r + n, r - n:r] = B_lower

            Z = Z_path[t]
            b[r] = ((params.rho_i / sig) * Z[params.IDX_RSTAR]
                    + (1.0 / sig) * Z[params.IDX_UH]
                    - (om_fi / sig) * Z[params.IDX_UC])
            b[r + 1] = Z[params.IDX_UP]
            b[r + 2] = om_fi * (Z[params.IDX_RSTAR] + Z[params.IDX_UC])

        x = np.linalg.solve(A, b)

        y_path = x[0::3]
        pi_path = x[1::3]
        i_path = x[2::3]
        rstar = Z_path[:, params.IDX_RSTAR]

        return rstar, rstar.copy(), i_path, y_path, pi_path

    def simulate_irf_filtering(regime_name, eq_data, q0, horizon):
        """Simulate IRF using the equilibrium law of motion for beliefs.

        Uses equilibrium coefficients (Phi, Psi, Omega) directly:
            m_{h,t} = Phi_h @ m_{h,t-1} + Psi_h @ m_{c,t-1} + Omega_h @ Z_t
            m_{c,t} = Phi_c @ m_{c,t-1} + Psi_c @ m_{h,t-1} + Omega_c @ Z_t
        Macro outcomes: (y, pi, i) = Gamma * i_{t-1} + M_h @ m_{h,t}

        Returns: (rstar_h, rstar_c, i_path, y_path, pi_path)
        """
        DIM_M_local = params.DIM_M

        if regime_name == "CK":
            Phi_h = eq_data['Phi_h']
            Psi_h = eq_data['Psi_h']
            Omega_h = eq_data['Omega_h']
            Phi_c = eq_data['Phi_c']
            Psi_c = eq_data['Psi_c']
            Omega_c = eq_data['Omega_c']
            M_h_m = eq_data['M_h']       # (3, DIM_M)
        else:  # HoM -- use actual equilibrium law of motion (Eqs C.24-C.26)
            act = eq_data['actual']
            Phi_h = act['Phi_h']
            Psi_h = act['Psi_h']
            Omega_h = act['Omega_h']
            Phi_c = act['Phi_c']
            Psi_c = act['Psi_c']
            Omega_c = act['Omega_c']
            M_h_m = eq_data['M_h']       # (3, DIM_M)

        rstar_h = np.zeros(horizon)
        rstar_c = np.zeros(horizon)
        i_path = np.zeros(horizon)
        y_path = np.zeros(horizon)
        pi_path = np.zeros(horizon)

        Z_t = q0.copy()
        m_h = np.zeros(DIM_M_local)
        m_c = np.zeros(DIM_M_local)
        i_prev = 0.0

        for t in range(horizon):
            # Equilibrium belief update (simultaneous)
            m_h_new = Phi_h @ m_h + Psi_h @ m_c + Omega_h @ Z_t
            m_c_new = Phi_c @ m_c + Psi_c @ m_h + Omega_c @ Z_t
            m_h, m_c = m_h_new, m_c_new

            # Macro outcomes from HH beliefs
            macro = Gamma * i_prev + M_h_m @ m_h
            y_path[t] = macro[0]
            pi_path[t] = macro[1]
            i_path[t] = macro[2]

            # Agents' first-order beliefs about r**
            rstar_h[t] = m_h[params.IDX_RSTAR]
            rstar_c[t] = m_c[params.IDX_RSTAR]

            i_prev = i_path[t]
            Z_t = R @ Z_t

        return rstar_h, rstar_c, i_path, y_path, pi_path

    def compute_all_irfs(shock_idx, shock_size, horizon):
        """Compute IRFs for all three regimes.

        Returns dict with keys 'FI', 'CK', 'HoM', each containing
        (rstar_h, rstar_c, i_path, y_path, pi_path).
        """
        q0 = np.zeros(N_Z)
        q0[shock_idx] = shock_size

        fi = simulate_irf_fi(q0, horizon)
        ck_irf = simulate_irf_filtering("CK", ck, q0, horizon)
        hom_irf = simulate_irf_filtering("HoM", hom, q0, horizon)

        return {'FI': fi, 'CK': ck_irf, 'HoM': hom_irf}

    # ------------------------------------------------------------------
    # Shock configurations matching paper's Figures E.1-E.4
    # ------------------------------------------------------------------

    # E.1: Monetary policy shock -- normalized to 25bp i_t reduction under HoM
    #   Use negative u_c (accommodative). Find scale s.t. HoM i_t[0] = -25bp.
    test = compute_all_irfs(params.IDX_UC, -1.0, 20)
    hom_i0 = test['HoM'][2][0]  # i_path[0] for unit negative u_c under HoM
    if abs(hom_i0) > 1e-15:
        # Scale so that HoM i_t[0] * PP = -0.25 pp (25bp reduction)
        mp_shock_size = 0.25 / hom_i0  # negative (accommodative)
    else:
        mp_shock_size = -params.sigma_uc  # fallback
    print(f"\n  Monetary policy: HoM i_t[0] per unit negative u_c = {hom_i0:.6f}")
    print(f"  Shock size for 25bp = {mp_shock_size:.6f}")

    mp_irfs = compute_all_irfs(params.IDX_UC, mp_shock_size, 20)

    # E.2: Demand shock -- 1 SD, negative (adverse demand)
    demand_irfs = compute_all_irfs(params.IDX_UH, -params.sigma_uh, 40)

    # E.3: Cost-push shock -- 1 SD, positive (inflationary)
    cp_irfs = compute_all_irfs(params.IDX_UP, params.sigma_up, 40)

    # E.4: r** fundamentals shock -- 1 percentage point decrease
    rstar_irfs = compute_all_irfs(params.IDX_RSTAR, -1.0, 40)

    # Verify normalization
    print(f"\n  E.1 HoM i_t[0] = {mp_irfs['HoM'][2][0] * PP:.3f} pp "
          f"(target: -0.250)")
    print(f"  E.1 FI  i_t[0] = {mp_irfs['FI'][2][0] * PP:.3f} pp")
    print(f"  E.4 FI  r** beliefs[0] = {rstar_irfs['FI'][0][0] * PP:.3f} pp "
          f"(target: -1.000)")

    all_shocks = [
        ("monetary_policy", mp_irfs, 20,
         "Figure E.1: Responses to an expansionary monetary policy shock"),
        ("demand", demand_irfs, 40,
         "Figure E.2: Responses to a demand shock"),
        ("cost_push", cp_irfs, 40,
         "Figure E.3: Responses to a cost-push shock"),
        ("rstar", rstar_irfs, 40,
         "Figure E.4: Responses to a change in r-star fundamentals"),
    ]
    print("  All IRFs computed.")

    # ==================================================================
    # Step 5: Plot IRFs (matching paper layout)
    # ==================================================================
    print("\nStep 5: Generating plots...")

    def prepend_zero(tup):
        """Prepend a zero entry to each array in a 5-tuple of IRF paths."""
        return tuple(np.concatenate([[0.0], arr]) for arr in tup)

    for shock_name, irf_data, horizon, title in all_shocks:
        fi = irf_data['FI']    # (rstar_h, rstar_c, i, y, pi)
        ck_d = irf_data['CK']
        hom_d = irf_data['HoM']

        # r** shock: paper plots pre-shock steady state at t=0
        if shock_name == "rstar":
            fi = prepend_zero(fi)
            ck_d = prepend_zero(ck_d)
            hom_d = prepend_zero(hom_d)
            quarters = np.arange(horizon + 1)
        else:
            quarters = np.arange(horizon)

        fig, axes = plt.subplots(2, 2, figsize=(12, 8))
        fig.suptitle(title, fontsize=13, fontweight='bold')

        # ---- Top-left: r* beliefs ----
        ax = axes[0, 0]
        # True r** (black dotted)
        ax.plot(quarters, fi[0] * PP, 'k:', linewidth=1.0, label=r'$r^{**}$')
        # HH beliefs r* (blue: solid=HM, dashed=CK)
        ax.plot(quarters, hom_d[0] * PP, 'b-', linewidth=1.5,
                label=r'$r^*_{HM}$')
        ax.plot(quarters, ck_d[0] * PP, 'b--', linewidth=1.5,
                label=r'$r^*_{CK}$')
        # CB beliefs r̂* (red: solid=HM, dashed=CK)
        ax.plot(quarters, hom_d[1] * PP, 'r-', linewidth=1.5,
                label=r'$\hat{r}^*_{HM}$')
        ax.plot(quarters, ck_d[1] * PP, 'r--', linewidth=1.5,
                label=r'$\hat{r}^*_{CK}$')
        ax.axhline(y=0, color='gray', linewidth=0.5)
        ax.set_title("r* beliefs")
        ax.legend(fontsize=7, loc='best')
        ax.grid(True, alpha=0.3)

        # ---- Top-right: Interest rate ----
        ax = axes[0, 1]
        ax.plot(quarters, hom_d[2] * PP, 'b-', linewidth=1.5,
                label=r'$i_{HM}$')
        ax.plot(quarters, ck_d[2] * PP, 'r--', linewidth=1.5,
                label=r'$i_{CK}$')
        ax.plot(quarters, fi[2] * PP, 'k:', linewidth=1.0,
                label=r'$i_{FI}$')
        ax.axhline(y=0, color='gray', linewidth=0.5)
        ax.set_title(r"Interest rate $i_t$")
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)

        # ---- Bottom-left: Output gap ----
        ax = axes[1, 0]
        ax.plot(quarters, hom_d[3] * PP, 'b-', linewidth=1.5,
                label=r'$\tilde{y}_{HM}$')
        ax.plot(quarters, ck_d[3] * PP, 'r--', linewidth=1.5,
                label=r'$\tilde{y}_{CK}$')
        ax.plot(quarters, fi[3] * PP, 'k:', linewidth=1.0,
                label=r'$\tilde{y}_{FI}$')
        ax.axhline(y=0, color='gray', linewidth=0.5)
        ax.set_title(r"Output gap $\tilde{y}_t$")
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)

        # ---- Bottom-right: Inflation ----
        ax = axes[1, 1]
        ax.plot(quarters, hom_d[4] * PP, 'b-', linewidth=1.5,
                label=r'$\pi_{HM}$')
        ax.plot(quarters, ck_d[4] * PP, 'r--', linewidth=1.5,
                label=r'$\pi_{CK}$')
        ax.plot(quarters, fi[4] * PP, 'k:', linewidth=1.0,
                label=r'$\pi_{FI}$')
        ax.axhline(y=0, color='gray', linewidth=0.5)
        ax.set_title(r"Inflation $\pi_t$")
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)

        for ax in axes.flat:
            ax.set_xlabel("Quarters")

        plt.tight_layout()
        fig_path = os.path.join(RESULTS_DIR, f"irf_{shock_name}.png")
        fig.savefig(fig_path, dpi=150)
        plt.close(fig)
        print(f"  Saved: {fig_path}")

    # ==================================================================
    # Step 6: Summary diagnostics
    # ==================================================================
    print("\n" + "=" * 60)
    print("Pipeline complete.")
    print(f"  CK converged: {ck['converged']} ({ck['iterations']} iters)")
    print(f"  HoM converged: {hom['converged']} ({hom['iterations']} iters)")
    print(f"  Results saved to: {RESULTS_DIR}")

    # Eigenvalue diagnostics
    ck_rho_h = np.max(np.abs(np.linalg.eigvals(ck['Phi_h'])))
    ck_rho_c = np.max(np.abs(np.linalg.eigvals(ck['Phi_c'])))
    print(f"\n  CK Phi_h spectral radius: {ck_rho_h:.4f}")
    print(f"  CK Phi_c spectral radius: {ck_rho_c:.4f}")

    hom_rho_h = np.max(np.abs(np.linalg.eigvals(hom['actual']['Phi_h'])))
    hom_rho_c = np.max(np.abs(np.linalg.eigvals(hom['actual']['Phi_c'])))
    print(f"  HoM actual Phi_h spectral radius: {hom_rho_h:.4f}")
    print(f"  HoM actual Phi_c spectral radius: {hom_rho_c:.4f}")

    # M_h diagnostics
    print(f"\n  CK M_h max |entry|: {np.max(np.abs(ck['M_h'])):.4f}")
    print(f"  HoM M_h max |entry|: {np.max(np.abs(hom['M_h'])):.4f}")

    # Quick sanity: FI long-run i for r** shock should approach -1pp (= r**)
    fi_rstar = rstar_irfs['FI']
    print(f"\n  FI i_t at t=39 for r** shock: {fi_rstar[2][-1] * PP:.3f} pp "
          f"(target: -1.0)")
    print(f"  FI y_t at t=39 for r** shock: {fi_rstar[3][-1] * PP:.4f} pp "
          f"(target: 0)")
    print(f"  FI pi_t at t=39 for r** shock: {fi_rstar[4][-1] * PP:.4f} pp "
          f"(target: 0)")
    print("=" * 60)


if __name__ == "__main__":
    main()
