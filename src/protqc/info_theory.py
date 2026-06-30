import numpy as np
import pandas as pd
from scipy.stats import entropy
from functools import wraps
from pathlib import Path
from typing import Optional, Sequence

# Global debug flag - can be set to True/False to control verbose logging
DEBUG = False

def debug_mode(func):
    """
    Decorator to control verbose print statements in functions.
    When DEBUG=True, all print statements execute normally.
    When DEBUG=False, print statements are suppressed (only returns are executed).
    
    Usage:
        @debug_mode
        def my_function():
            print("This will only show when DEBUG=True")
            return result
    """
    @wraps(func)
    def wrapper(*args, **kwargs):
        if DEBUG:
            # Debug mode ON - execute function normally with all prints
            return func(*args, **kwargs)
        else:
            # Debug mode OFF - suppress print statements
            import builtins
            original_print = builtins.print
            builtins.print = lambda *args, **kwargs: None
            try:
                result = func(*args, **kwargs)
            finally:
                builtins.print = original_print
            return result
    return wrapper


@debug_mode
def residue_information_gain(Z_samples: np.ndarray,
                             E_total: np.ndarray,
                             bins: int = 10000) -> pd.DataFrame:
    """
    Quantify per-residue information gain: how much knowing residue r's state
    changes the ensemble energy distribution.

    Args:
        Z_samples: (shots, n_res) array of ±1 or {0,1} samples
        E_total:   (shots,) array of total configuration energies
        bins:      histogram bins for energy discretisation

    Returns:
        pandas DataFrame with columns:
        ['res_index', 'p_exc', 'KL_ground', 'KL_exc', 'I_total', 'ΔE']

        `KL_ground` and `KL_exc` are normalised after computation so each
        column sums to 1 across residues.
    """
    import pandas as pd
    # Z_samples[Z_samples==-1] = 0
    # normalise to {0,1}
    if np.min(Z_samples) < 0:
        states = (Z_samples == -1).astype(int) # Indicator: Z-samples == +1
    else:
        states = Z_samples.copy().astype(int)
    shots, n_res = states.shape

    # global distribution
    counts_all, edges = np.histogram(E_total, bins=bins, density=False)
    P_all = counts_all / np.sum(counts_all) + 1e-12
    centers = 0.5 * (edges[:-1] + edges[1:])

    print(f"P_all sum is {np.sum(P_all)}, std dev is {np.std(P_all)}, \
       min is {np.min(P_all)}, P1 max is {np.max(P_all)}")
    #breakpoint()

    results = []
    for r in range(n_res):
        mask0 = (states[:, r] == 0) # Indicator:states[:,r] == 0
        mask1 = (states[:, r] == 1)
        # if np.sum(mask0) < 5 or np.sum(mask1) < 5:
        #breakpoint()

        def hist_for_mask(mask):
            c, _ = np.histogram(E_total[mask], bins=edges, density=False)
            p = c / np.sum(c) + 1e-5
            return p

        P0 = hist_for_mask(mask0)
        P1 = hist_for_mask(mask1)

        print(f"Residue {r+1}: mask0 sum is {np.sum(mask0)}, \
          mask1 sum is {np.sum(mask1)}")

        KL0 = entropy(P0, P_all, base=np.e)
        KL1 = entropy(P1, P_all, base=np.e)
        # dE = np.mean(E_total[mask1]) - np.mean(E_total[mask0])
        print(f"Residue is {r+1}: P0 is {P0}, P1 is {P1}, P_all is {P_all}, KL_ground={KL0:.3f}, KL_exc={KL1:.3f}")
        print(f"P0 sum is {np.sum(P0)}, P0 std dev is {np.std(P0)}, P0 min is {np.min(P0)}, P0 max is {np.max(P0)}")
        print(f"P1 sum is {np.sum(P1)}, P1 std dev is {np.std(P1)}, P1 min is {np.min(P1)}, P1 max is {np.max(P1)}")
        print(f"E_total sum is {np.sum(E_total)}, E_total std dev is {np.std(E_total)}, E_total min is {np.min(E_total)}, E_total max is {np.max(E_total)}")
        results.append({
            "res_index": r+1,
            "KL_ground": KL0,
            "KL_exc": KL1,
            # "ΔE": dE
        })

    df = pd.DataFrame(results)
    # Convert raw KL divergences into per-residue normalised contributions.
    df["KL_ground"] /= np.sum(df["KL_ground"])
    df["KL_exc"] /= np.sum(df["KL_exc"])
    # df.sort_values("res_index", ascending=False, inplace=True)
    return df

# def residue_information_gain(Z_samples: np.ndarray,
#                              E_total: np.ndarray,
#                              bins: int = 10000) -> pd.DataFrame:
#     """
#     Quantify per-residue information gain: how much knowing residue r's state
#     changes the ensemble energy distribution.

#     Args:
#         Z_samples: (shots, n_res) array of ±1 or {0,1} samples
#         E_total:   (shots,) array of total configuration energies
#         bins:      histogram bins for energy discretisation

#     Returns:
#         pandas DataFrame with columns:
#         ['res_index', 'p_exc', 'KL_ground', 'KL_exc', 'I_total', 'ΔE']
#     """
#     # Z_samples[Z_samples==-1] = 0
#     # normalise to {0,1}
#     if np.min(Z_samples) < 0:
#         states = (Z_samples == -1).astype(int) # Indicator: Z-samples == +1
#     else:
#         states = Z_samples.copy().astype(int)
#     shots, n_res = states.shape

#     # global distribution
#     counts_all, edges = np.histogram(E_total, bins=bins, density=False)
#     P_all = counts_all / np.sum(counts_all) + 1e-12
#     centers = 0.5 * (edges[:-1] + edges[1:])

#     results = []
#     for r in range(n_res):
#         mask0 = (states[:, r] == 0) # Indicator:states[:,r] == 0
#         mask1 = (states[:, r] == 1)
#         # if np.sum(mask0) < 5 or np.sum(mask1) < 5:
#         #     continue

#         def hist_for_mask(mask):
#             c, _ = np.histogram(E_total[mask], bins=edges, density=False)
#             ### original line has bug here: PRATIK please check below lines
#             total = np.sum(c)
#             # Handle empty histogram (avoid division by zero)
#             # if total == 0:
#             #     return np.ones_like(c) / len(c)  # Uniform distribution
#             p = c / total + 1e-5
#             return p

#         P0 = hist_for_mask(mask0)
#         P1 = hist_for_mask(mask1)

#         KL0 = entropy(P0, P_all, base=np.e)
#         KL1 = entropy(P1, P_all, base=np.e)
#         # dE = np.mean(E_total[mask1]) - np.mean(E_total[mask0])
#         print(f"Residue is {r+1}: P0 is {P0}, P1 is {P1}, P_all is {P_all}, KL_ground={KL0:.3f}, KL_exc={KL1:.3f}")
#         print(f"P0 sum is {np.sum(P0)}, P0 std dev is {np.std(P0)}, P0 min is {np.min(P0)}, P0 max is {np.max(P0)}")
#         print(f"P1 sum is {np.sum(P1)}, P1 std dev is {np.std(P1)}, P1 min is {np.min(P1)}, P1 max is {np.max(P1)}")

#         results.append({
#             "res_index": r+1,
#             "KL_ground": KL0,
#             "KL_exc": KL1,
#             # "ΔE": dE
#         })

#     df = pd.DataFrame(results)
#     breakpoint()
#     # df.sort_values("res_index", ascending=False, inplace=True)
#     return df


# @debug_mode
# def directional_excitation_information(Z_samples, threshold=0.0):
#     """
#     This function is not Mutual Information anymore - it is Directional_excitation_information
#         n_res = n_qubits or number of residues
#         n_shots = number of samples collected from circuit measurements
#         -------
#         args: Z_samples: (shots, n_res) array of ±1 or {0,1} samples
#         returns: MI: (n_res, n_res) matrix (nats)
#     """
#     n_shots, n_res = Z_samples.shape

#     eps = 1e-6
#     DEI = np.zeros((n_res, n_res))
#     ref_states = {1:-1, 0: +1} # Circuit state for |0> and |1> states

#     for i in range(n_res):
#         xi = Z_samples[:, i]
#         for j in range(n_res):
#             xj = Z_samples[:, j]
#             if i == j:
#                 continue
#             # joint probabilities for (0,0), (0,1), (1,0), (1,1)
#             joint = np.zeros((2, 2))
#             for a in [0, 1]:
#                 for b in [0, 1]:
#                     # joint[a, b] = np.mean((xi == ref_states[a]) & (xj == ref_states[b]))
#                     num_j_b = np.sum(xj == ref_states[b])
#                     x_i_given_j_b = xi[xj == ref_states[b]]
#                     x_i_a_given_j_b = np.sum(x_i_given_j_b == ref_states[a])
#                     joint[a,b] = x_i_a_given_j_b / num_j_b
#                     # print(f"i,j is {i},{j} - A,B is {a},{b}; num_j_b={num_j_b},x_i_given_j_b={x_i_given_j_b};x_i_a_given_j_b={x_i_a_given_j_b} p({i}={a}|{j}={b})={joint[a,b]}")
#                     pi = np.mean((xi == ref_states[a])) + eps
#                     pj = np.mean((xj == ref_states[b])) + eps
#                     joint[a,b] = joint[a,b]*np.log2(joint[a,b]/(pi))

#                     # print(f"i,j is {i},{j} - A,B is {a},{b}; pi={pi},pj={pj}; joint({i},{j})={joint[a,b]}")

#             # MI[i, j] = np.max(joint*np.log2(joint/(pi @ pj)))
#             DEI[i,j] = max(joint[1,1] - threshold, 0.0)
#             # MI[i,j] += 0.4
#             print(f"DEI({i},{j}) = {DEI[i,j]}")
#             # MI[j, i] = MI[i, j]

#     return DEI

@debug_mode
def directional_excitation_information(Z_samples, threshold=0.0):
    """
    Directional Excitation Information (DEI).

    This is NOT mutual information.

    For each ordered residue pair (i, j), this function computes:

        DEI[i, j] = D_KL( P(R_i | R_j = 1) || P(R_i) )

    where:
        R_k = 1 means residue k is in the excited state
        R_k = 0 means residue k is in the ground state

    Matrix convention:
        DEI[i, j] measures directional excitation information from
        source residue j to target residue i.

        In words:
            "How much does knowing that residue j is excited change the
             full two-state probability distribution of residue i?"

    Explicitly, for R_i in {0, 1}:

        DEI[i, j] =
            sum_{a in {0,1}}
                P(R_i = a | R_j = 1)
                log2[
                    P(R_i = a | R_j = 1) / P(R_i = a)
                ]

    Important:
        - This is a conditional KL divergence.
        - It is directional: DEI[i, j] != DEI[j, i] in general.
        - It is non-negative by construction.
        - It is not weighted by P(R_j = 1).
        - Units are bits because log2 is used.
        - If residue j is never observed excited, DEI[:, j] is set to 0.

    Parameters
    ----------
    Z_samples : np.ndarray
        Array of shape (shots, n_res).

        Accepted encodings:
            ±1 encoding:
                +1 = ground
                -1 = excited

            binary encoding:
                0 = ground
                1 = excited

    threshold : float, default=0.0
        Optional noise floor. Values below threshold are set to zero after
        computing the KL divergence. This does not weight the score.

    Returns
    -------
    DEI : np.ndarray
        Directional excitation information matrix of shape (n_res, n_res),
        where DEI[i, j] is the influence j -> i in bits.
    """

    Z_samples = np.asarray(Z_samples)
    if Z_samples.ndim != 2:
        raise ValueError("Z_samples must be a 2D array with shape (shots, n_res).")

    n_shots, n_res = Z_samples.shape
    if n_shots == 0:
        raise ValueError("Z_samples must contain at least one shot.")

    eps = 1e-12

    # ------------------------------------------------------------
    # Convert samples to binary excitation states:
    #     R = 1 excited
    #     R = 0 ground
    # ------------------------------------------------------------
    if np.min(Z_samples) < 0:
        # Pauli-Z convention used in the circuit:
        # +1 = ground, -1 = excited
        R = (Z_samples == -1).astype(np.uint8)
    else:
        # Binary convention:
        # 0 = ground, 1 = excited
        R = (Z_samples > 0).astype(np.uint8)

    # ------------------------------------------------------------
    # Marginal distribution P(R_i)
    # p_exc[i] = P(R_i = 1)
    # p_gnd[i] = P(R_i = 0)
    # ------------------------------------------------------------
    p_exc = R.mean(axis=0).astype(float)
    p_gnd = 1.0 - p_exc

    # Clip only for numerical safety inside the logarithm.
    p_exc_safe = np.clip(p_exc, eps, 1.0 - eps)
    p_gnd_safe = np.clip(p_gnd, eps, 1.0 - eps)

    DEI = np.zeros((n_res, n_res), dtype=float)

    # ------------------------------------------------------------
    # For each source residue j, condition on R_j = 1.
    # Then compare P(R_i | R_j = 1) against P(R_i).
    # ------------------------------------------------------------
    for j in range(n_res):
        source_excited_mask = (R[:, j] == 1)
        source_support = int(np.sum(source_excited_mask))

        # Conditional distribution is undefined if residue j is never excited.
        # We set the corresponding source column to zero.
        if source_support == 0:
            continue

        R_given_j_excited = R[source_excited_mask, :]

        # P(R_i = 1 | R_j = 1)
        p_exc_cond = R_given_j_excited.mean(axis=0).astype(float)

        # P(R_i = 0 | R_j = 1)
        p_gnd_cond = 1.0 - p_exc_cond

        # Numerical safety for log terms.
        p_exc_cond_safe = np.clip(p_exc_cond, eps, 1.0 - eps)
        p_gnd_cond_safe = np.clip(p_gnd_cond, eps, 1.0 - eps)

        # --------------------------------------------------------
        # Bernoulli KL divergence:
        #
        # D_KL( P(R_i | R_j=1) || P(R_i) )
        #
        # = q1 log2(q1 / p1) + q0 log2(q0 / p0)
        #
        # where:
        #   q1 = P(R_i=1 | R_j=1)
        #   q0 = P(R_i=0 | R_j=1)
        #   p1 = P(R_i=1)
        #   p0 = P(R_i=0)
        # --------------------------------------------------------
        kl = (
            p_exc_cond_safe * np.log2(p_exc_cond_safe / p_exc_safe)
            +
            p_gnd_cond_safe * np.log2(p_gnd_cond_safe / p_gnd_safe)
        )

        # Remove tiny numerical negatives.
        kl = np.maximum(kl, 0.0)

        # Optional noise floor.
        if threshold > 0.0:
            kl = np.maximum(kl - float(threshold), 0.0)

        DEI[:, j] = kl

    # A residue is not considered directionally informative to itself here.
    np.fill_diagonal(DEI, 0.0)

    return DEI


def mutual_information_matrix(Z_samples: np.ndarray, eps: float = 1e-12) -> np.ndarray:
    """
    Pairwise mutual information between residue excitation states.

    This is symmetric and non-directional:

        MI(i, j) = D_KL(P(R_i, R_j) || P(R_i) P(R_j))

    The KL divergence is computed with scipy.stats.entropy using base 2,
    so values are reported in bits.

    Accepted encodings match directional_excitation_information:
        +1 = ground, -1 = excited for Pauli-Z samples
        0 = ground, 1 = excited for binary samples
    """
    Z_samples = np.asarray(Z_samples)
    if Z_samples.ndim != 2:
        raise ValueError("Z_samples must be a 2D array with shape (shots, n_res).")

    n_shots, n_res = Z_samples.shape
    if n_shots == 0:
        raise ValueError("Z_samples must contain at least one shot.")

    if np.min(Z_samples) < 0:
        states = (Z_samples == -1).astype(np.uint8)
    else:
        states = (Z_samples > 0).astype(np.uint8)

    mi = np.zeros((n_res, n_res), dtype=float)

    for i in range(n_res):
        xi = states[:, i]
        px = np.array(
            [
                np.mean(xi == 0),
                np.mean(xi == 1),
            ],
            dtype=float,
        )

        for j in range(i + 1, n_res):
            yj = states[:, j]
            py = np.array(
                [
                    np.mean(yj == 0),
                    np.mean(yj == 1),
                ],
                dtype=float,
            )

            joint = np.array(
                [
                    np.mean((xi == 0) & (yj == 0)),
                    np.mean((xi == 0) & (yj == 1)),
                    np.mean((xi == 1) & (yj == 0)),
                    np.mean((xi == 1) & (yj == 1)),
                ],
                dtype=float,
            )
            independent = np.outer(px, py).reshape(-1)

            joint_safe = joint + eps
            independent_safe = independent + eps
            joint_safe = joint_safe / np.sum(joint_safe)
            independent_safe = independent_safe / np.sum(independent_safe)

            value = float(entropy(joint_safe, independent_safe, base=2))
            value = max(value, 0.0)
            mi[i, j] = value
            mi[j, i] = value

    np.fill_diagonal(mi, 0.0)
    return mi


def joint_state_matrices(Z_samples):
    """
    Count co-occurrences of (0,0), (0,1), (1,0), (1,1)
    for all residue pairs (i,j).

    Returns
    -------
    counts_00, counts_01, counts_10, counts_11 : np.ndarray
        Each of shape (n_res, n_res), containing the
        counts (or probabilities if normalised) of each joint state.
    """
    n_shots, n_res = Z_samples.shape

    counts_00 = np.zeros((n_res, n_res))
    counts_01 = np.zeros((n_res, n_res))
    counts_10 = np.zeros((n_res, n_res))
    counts_11 = np.zeros((n_res, n_res))

    for i in range(n_res):
        xi = Z_samples[:, i]
        for j in range(i+1, n_res):
            xj = Z_samples[:, j]

            # Compute Boolean masks
            mask00 = (xi == 1) & (xj == 1)
            mask01 = (xi == 1) & (xj == -1)
            mask10 = (xi == -1) & (xj == 1)
            mask11 = (xi == -1) & (xj == -1)

            # Count or normalise by n_shots
            counts_00[i, j] = counts_00[j, i] = np.sum(mask00)
            counts_01[i, j] = counts_10[j, i] = np.sum(mask01)
            counts_10[i, j] = counts_01[j, i] = np.sum(mask10)
            counts_11[i, j] = counts_11[j, i] = np.sum(mask11)
            total = counts_00[i, j] + counts_01[i, j] + counts_10[i, j] + counts_11[i, j]
            # print(f"{i},{j}; Counts -- 00:{counts_00[i,j]};10: {counts_10[i,j]};01: {counts_01[i,j]};11: {counts_11[i,j]}")
            # print(f"{i},{j}; Counts Normalised -- 00:{counts_00[i,j]/n_shots};10: {counts_10[i,j]/n_shots};01: {counts_01[i,j]/n_shots};11: {counts_11[i,j]/n_shots}")
            assert total == n_shots, f"{i},{j} -- {counts_00[i, j]},{counts_01[i, j]},{counts_10[i, j]},{counts_11[i, j]} doesn't match the n_shots"

    # Optionally normalise to probabilities:
    # counts_00 /= n_shots; counts_01 /= n_shots; counts_10 /= n_shots; counts_11 /= n_shots

    return counts_00/n_shots, counts_01/n_shots, counts_10/n_shots, counts_11/n_shots
 

def compute_static_mi_from_counts(counts_00, counts_01, counts_10, counts_11, eps=1e-12):
    """
    Static-only mutual information contribution per pair:
        I_static_ij = p00 * log(p00/(p_i0 p_j0)) + p11 * log(p11/(p_i1 p_j1))
    where all probabilities and marginals are computed from the *pair's* full
    4-state distribution. This avoids cross-pair mixing and double-counting.

    Args:
        counts_00, counts_01, counts_10, counts_11: (N,N) joint counts for states
            (0,0), (0,1), (1,0), (1,1). Raw counts are fine.

    Returns:
        I_static: (N,N) matrix (nats) with the static (00+11) contribution to MI.
    """
    N = counts_00.shape[0]
    I_static = np.zeros((N, N), dtype=float)

    # Pairwise totals and probabilities
    total = counts_00 + counts_01 + counts_10 + counts_11 + eps

    p00 = counts_00 / total
    p01 = counts_01 / total
    p10 = counts_10 / total
    p11 = counts_11 / total

    # Pairwise marginals from the full (00,01,10,11) distribution
    # p_i(0) = p00 + p01 ; p_i(1) = p10 + p11
    # p_j(0) = p00 + p10 ; p_j(1) = p01 + p11
    p_i0 = p00 + p01
    p_i1 = p10 + p11
    p_j0 = p00 + p10
    p_j1 = p01 + p11

    # Static-only MI contribution (sum over 00 and 11 terms)
    with np.errstate(divide='ignore', invalid='ignore'):
        term00 = np.where(p00 > 0, p00 * (np.log(p00 + eps) - np.log((p_i0 * p_j0) + eps)), 0.0)
        term11 = np.where(p11 > 0, p11 * (np.log(p11 + eps) - np.log((p_i1 * p_j1) + eps)), 0.0)
        I_static = term00 + term11

    # Symmetrise (numerically) and clear diagonal
    I_static = 0.5 * (I_static + I_static.T)
    np.fill_diagonal(I_static, 0.0)
    return I_static


def _to_excited_binary(Z_samples: np.ndarray) -> np.ndarray:
    """Convert samples to {0,1}, where 1 means excited."""
    samples = np.asarray(Z_samples)
    if samples.ndim != 2:
        raise ValueError("Z_samples must be a 2D array with shape (shots, n_res).")

    if np.min(samples) < 0:
        # In this codebase, -1 corresponds to the excited state.
        return (samples == -1).astype(np.uint8)
    return (samples > 0).astype(np.uint8)


def _residue_coordinates_from_pdb(
    pdb_path: str,
    chain: Optional[str] = "A",
    protein_only: bool = True,
    anchor_atom: str = "CA",
) -> tuple[np.ndarray, list[str]]:
    """
    Extract per-residue coordinates from a PDB file in MDAnalysis residue order.
    """
    def _residue_chain_label(residue) -> str:
        chain_label = str(getattr(residue, "chainID", "") or "").strip()
        if chain_label:
            return chain_label
        try:
            chain_ids = getattr(residue.atoms, "chainIDs", None)
            if chain_ids is not None and len(chain_ids) > 0:
                chain_label = str(chain_ids[0] or "").strip()
                if chain_label:
                    return chain_label
        except Exception:
            pass
        return str(getattr(residue, "segid", "") or "").strip()

    try:
        import MDAnalysis as mda
    except Exception as exc:
        raise ImportError(
            "MDAnalysis is required for PDB-based voxelisation. "
            "Install it or pass residue_coords directly."
        ) from exc

    try:
        from protqc.helper import AA3
    except Exception:
        AA3 = None

    universe = mda.Universe(str(pdb_path))
    coords: list[np.ndarray] = []
    labels: list[str] = []

    for residue in universe.residues:
        chain_label = _residue_chain_label(residue)
        if chain is not None and chain_label != chain:
            continue

        resname = residue.resname.strip()
        if protein_only and AA3 is not None and resname not in AA3:
            continue

        atoms = residue.atoms
        if atoms.n_atoms == 0:
            continue

        anchor = atoms.select_atoms(f"name {anchor_atom}") if anchor_atom else atoms
        if anchor.n_atoms > 0:
            position = anchor.positions[0]
        else:
            position = atoms.positions.mean(axis=0)

        coords.append(np.asarray(position, dtype=float))
        labels.append(f"{int(residue.resid)}-{resname}")

    if not coords:
        raise ValueError(
            "No residues found for voxelisation. Check pdb_path, chain, and protein_only options."
        )

    return np.asarray(coords, dtype=float), labels


def _kmeans_lloyd(X: np.ndarray, n_clusters: int, random_state: int = 0, max_iter: int = 100) -> np.ndarray:
    """Minimal k-means implementation used by spectral clustering."""
    X = np.asarray(X, dtype=float)
    n_samples = X.shape[0]
    if n_samples == 0:
        return np.array([], dtype=int)
    if n_samples == 1 or n_clusters <= 1:
        return np.zeros(n_samples, dtype=int)

    n_clusters = int(max(1, min(n_clusters, n_samples)))
    rng = np.random.default_rng(random_state)

    unique_rows = np.unique(X, axis=0)
    if unique_rows.shape[0] < n_clusters:
        n_clusters = unique_rows.shape[0]
    centroids = unique_rows[rng.choice(unique_rows.shape[0], size=n_clusters, replace=False)]

    labels = np.zeros(n_samples, dtype=int)
    for _ in range(max_iter):
        sq_dist = np.sum((X[:, None, :] - centroids[None, :, :]) ** 2, axis=2)
        new_labels = np.argmin(sq_dist, axis=1)
        if np.array_equal(new_labels, labels):
            break
        labels = new_labels

        for k in range(n_clusters):
            mask = labels == k
            if np.any(mask):
                centroids[k] = X[mask].mean(axis=0)
            else:
                centroids[k] = X[rng.integers(0, n_samples)]

    return labels


def _spectral_cluster_from_affinity(
    affinity: np.ndarray,
    n_clusters: int,
    random_state: int = 0,
) -> np.ndarray:
    """Spectral clustering from a symmetric non-negative affinity matrix."""
    A = np.asarray(affinity, dtype=float)
    n = A.shape[0]
    if A.ndim != 2 or A.shape[0] != A.shape[1]:
        raise ValueError("affinity must be a square matrix.")
    if n == 0:
        return np.array([], dtype=int)
    if n == 1 or n_clusters <= 1:
        return np.zeros(n, dtype=int)

    n_clusters = int(max(1, min(n_clusters, n)))
    A = 0.5 * (A + A.T)
    A = np.clip(A, 0.0, None)
    np.fill_diagonal(A, 0.0)

    degree = A.sum(axis=1)
    if np.all(degree <= 1e-12):
        return np.zeros(n, dtype=int)

    inv_sqrt_degree = 1.0 / np.sqrt(np.maximum(degree, 1e-12))
    laplacian = np.eye(n) - (inv_sqrt_degree[:, None] * A) * inv_sqrt_degree[None, :]
    _, eigvecs = np.linalg.eigh(laplacian)

    # embedding = eigvecs[:, :n_clusters]
    # Skip the trivial first eigenvector
    embedding = eigvecs[:, 1:n_clusters+1]
    norms = np.linalg.norm(embedding, axis=1, keepdims=True)
    norms[norms == 0.0] = 1.0
    embedding /= norms

    return _kmeans_lloyd(embedding, n_clusters=n_clusters, random_state=random_state)


def _grid_voxelise(
    coords: np.ndarray,
    voxel_size: float,
    eps: float = 1e-12,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, list[np.ndarray]]:
    """Hard voxelization on a Cartesian grid."""
    X = np.asarray(coords, dtype=float)
    n_res = X.shape[0]
    if n_res == 0:
        return (
            np.array([], dtype=int),
            np.zeros((0, 3), dtype=int),
            np.zeros((0, 3), dtype=float),
            np.array([], dtype=int),
            [],
        )
    if voxel_size <= 0.0:
        raise ValueError("voxel_size must be > 0.")

    origin = np.min(X, axis=0)
    grid_idx = np.floor((X - origin) / max(float(voxel_size), eps)).astype(int)

    key_to_id: dict[tuple[int, int, int], int] = {}
    voxel_ids = np.empty(n_res, dtype=int)
    for i, key_arr in enumerate(grid_idx):
        key = (int(key_arr[0]), int(key_arr[1]), int(key_arr[2]))
        if key not in key_to_id:
            key_to_id[key] = len(key_to_id)
        voxel_ids[i] = key_to_id[key]

    n_voxels = len(key_to_id)
    voxel_grid = np.zeros((n_voxels, 3), dtype=int)
    voxel_members: list[np.ndarray] = []
    voxel_centroids = np.zeros((n_voxels, 3), dtype=float)
    voxel_sizes = np.zeros(n_voxels, dtype=int)

    for key, vid in key_to_id.items():
        voxel_grid[vid] = np.asarray(key, dtype=int)
    for vid in range(n_voxels):
        members = np.flatnonzero(voxel_ids == vid)
        voxel_members.append(members)
        voxel_sizes[vid] = int(members.size)
        voxel_centroids[vid] = X[members].mean(axis=0)

    return voxel_ids, voxel_grid, voxel_centroids, voxel_sizes, voxel_members


def _overlapping_grid_voxelise(
    coords: np.ndarray,
    voxel_size: float,
    overlap_radius: Optional[float] = None,
    kernel_sigma: Optional[float] = None,
    membership_threshold: float = 0.05,
    max_memberships: Optional[int] = 8,
    eps: float = 1e-12,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, list[np.ndarray], np.ndarray]:
    """Soft Gaussian voxelization on an overlapping Cartesian grid."""
    X = np.asarray(coords, dtype=float)
    n_res = X.shape[0]
    if n_res == 0:
        return (
            np.array([], dtype=int),
            np.zeros((0, 3), dtype=int),
            np.zeros((0, 3), dtype=float),
            np.array([], dtype=float),
            np.array([], dtype=int),
            [],
            np.zeros((0, 0), dtype=float),
        )
    if voxel_size <= 0.0:
        raise ValueError("voxel_size must be > 0.")

    radius = float(voxel_size if overlap_radius is None else overlap_radius)
    sigma = float(0.75 * voxel_size if kernel_sigma is None else kernel_sigma)
    if radius <= 0.0:
        raise ValueError("voxel_overlap_radius must be > 0.")
    if sigma <= 0.0:
        raise ValueError("voxel_kernel_sigma must be > 0.")
    if not (0.0 <= float(membership_threshold) < 1.0):
        raise ValueError("voxel_membership_threshold must be in [0, 1).")
    if max_memberships is not None and int(max_memberships) <= 0:
        raise ValueError("max_voxel_memberships must be positive or None.")

    size = float(voxel_size)
    origin = np.min(X, axis=0) - 0.5 * size
    radius_cells = int(max(0, np.ceil(radius / size)))

    residue_candidates: list[list[tuple[tuple[int, int, int], float]]] = []
    key_to_id: dict[tuple[int, int, int], int] = {}

    for point in X:
        base = np.floor((point - origin) / max(size, eps)).astype(int)
        candidates: list[tuple[tuple[int, int, int], float, float]] = []
        for dx in range(-radius_cells, radius_cells + 1):
            for dy in range(-radius_cells, radius_cells + 1):
                for dz in range(-radius_cells, radius_cells + 1):
                    key_arr = base + np.array([dx, dy, dz], dtype=int)
                    key = (int(key_arr[0]), int(key_arr[1]), int(key_arr[2]))
                    center = origin + (np.asarray(key, dtype=float) + 0.5) * size
                    distance = float(np.linalg.norm(point - center))
                    if distance <= radius + eps:
                        weight = float(np.exp(-0.5 * (distance / sigma) ** 2))
                        candidates.append((key, weight, distance))

        if not candidates:
            key_arr = base.astype(int)
            key = (int(key_arr[0]), int(key_arr[1]), int(key_arr[2]))
            center = origin + (np.asarray(key, dtype=float) + 0.5) * size
            distance = float(np.linalg.norm(point - center))
            candidates.append((key, float(np.exp(-0.5 * (distance / sigma) ** 2)), distance))

        max_weight = max(weight for _, weight, _ in candidates)
        threshold = float(membership_threshold) * max_weight
        candidates = [item for item in candidates if item[1] >= threshold]
        candidates.sort(key=lambda item: (item[1], -item[2]), reverse=True)
        if max_memberships is not None:
            candidates = candidates[: int(max_memberships)]

        total = sum(weight for _, weight, _ in candidates)
        if total <= eps:
            total = eps
        normalized = [(key, float(weight / total)) for key, weight, _ in candidates]
        residue_candidates.append(normalized)

        for key, _ in normalized:
            if key not in key_to_id:
                key_to_id[key] = len(key_to_id)

    n_voxels = len(key_to_id)
    voxel_grid = np.zeros((n_voxels, 3), dtype=int)
    for key, vid in key_to_id.items():
        voxel_grid[vid] = np.asarray(key, dtype=int)

    weights = np.zeros((n_res, n_voxels), dtype=float)
    for res_idx, candidates in enumerate(residue_candidates):
        for key, weight in candidates:
            weights[res_idx, key_to_id[key]] = weight

    row_sums = weights.sum(axis=1, keepdims=True)
    weights = np.divide(weights, np.maximum(row_sums, eps))

    dominant_voxel_ids = np.argmax(weights, axis=1).astype(int)
    effective_sizes = weights.sum(axis=0)
    voxel_members = [np.flatnonzero(weights[:, vid] > 0.0) for vid in range(n_voxels)]
    voxel_member_counts = np.asarray([len(members) for members in voxel_members], dtype=int)

    voxel_centroids = np.zeros((n_voxels, 3), dtype=float)
    grid_centers = origin + (voxel_grid.astype(float) + 0.5) * size
    for vid in range(n_voxels):
        if effective_sizes[vid] > eps:
            voxel_centroids[vid] = (weights[:, vid][:, None] * X).sum(axis=0) / effective_sizes[vid]
        else:
            voxel_centroids[vid] = grid_centers[vid]

    return (
        dominant_voxel_ids,
        voxel_grid,
        voxel_centroids,
        effective_sizes,
        voxel_member_counts,
        voxel_members,
        weights,
    )


def _residue_directional_excitation_matrix(
    Z_samples: np.ndarray,
    min_support_res: int,
    eps: float = 1e-12,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Residue-level directional excitation matrix:
        D_res[i,j] = p(i excited | j excited) * log((p(i excited | j excited)+eps)/(p(i excited)+eps))
    where excited means Z == -1.
    """
    X = np.asarray(Z_samples)
    if X.ndim != 2:
        raise ValueError("Z_samples must be a 2D array with shape (shots, n_res).")

    if np.min(X) < 0:
        excited = (X == -1).astype(np.uint8)
    else:
        excited = (X > 0).astype(np.uint8)

    n_shots, n_res = excited.shape
    if n_shots <= 0:
        raise ValueError("Z_samples must contain at least one shot.")
    if min_support_res <= 0:
        raise ValueError("min_support_res must be > 0.")

    p = excited.mean(axis=0)
    D_res = np.zeros((n_res, n_res), dtype=float)
    source_support = np.zeros(n_res, dtype=int)

    for j in range(n_res):
        mask = excited[:, j] == 1
        source_support[j] = int(np.sum(mask))
        if source_support[j] < min_support_res:
            continue
        p_given = excited[mask].mean(axis=0)
        score = p_given * np.log((p_given + eps) / (p + eps))
        score[~np.isfinite(score)] = 0.0
        D_res[:, j] = score

    np.fill_diagonal(D_res, 0.0)
    return D_res, p, source_support


def _aggregate_block(
    block: np.ndarray,
    mode: str,
    topk: int,
) -> float:
    """Aggregate a DEI block for voxel-level directional coupling."""
    B = np.asarray(block, dtype=float)
    if B.size == 0:
        return 0.0

    mode_norm = str(mode).strip().lower()
    if mode_norm == "mean":
        return float(np.mean(B))
    if mode_norm == "sum_abs":
        return float(np.sum(np.abs(B)))
    if mode_norm == "topk_mean":
        flat = B.ravel()
        k = int(max(1, min(topk, flat.size)))
        idx = np.argpartition(np.abs(flat), -k)[-k:]
        return float(np.mean(flat[idx]))

    raise ValueError("voxel_agg_mode must be one of {'mean', 'sum_abs', 'topk_mean'}.")


def _voxel_coupling_from_residue_dei(
    D_res: np.ndarray,
    voxel_members: Sequence[np.ndarray],
    mode: str = "mean",
    topk: int = 5,
) -> np.ndarray:
    """
    Voxel directional coupling by DEI block aggregation:
        C_vox[V,W] = aggregate( D_res[np.ix_(RV, RW)] )
    """
    D = np.asarray(D_res, dtype=float)
    n_vox = len(voxel_members)
    C_vox = np.zeros((n_vox, n_vox), dtype=float)

    for v in range(n_vox):
        rv = np.asarray(voxel_members[v], dtype=int)
        if rv.size == 0:
            continue
        for w in range(n_vox):
            rw = np.asarray(voxel_members[w], dtype=int)
            if rw.size == 0:
                continue
            block = D[np.ix_(rv, rw)]
            C_vox[v, w] = _aggregate_block(block=block, mode=mode, topk=topk)

    np.fill_diagonal(C_vox, 0.0)
    C_vox[~np.isfinite(C_vox)] = 0.0
    return C_vox


def _weighted_aggregate_block(
    block: np.ndarray,
    weights: np.ndarray,
    mode: str,
    topk: int,
    eps: float = 1e-12,
) -> float:
    """Aggregate a DEI block using target/source voxel membership weights."""
    B = np.asarray(block, dtype=float)
    W = np.asarray(weights, dtype=float)
    valid = W > 0.0
    if B.size == 0 or not np.any(valid):
        return 0.0

    mode_norm = str(mode).strip().lower()
    if mode_norm == "mean":
        return float(np.sum(B * W) / (np.sum(W) + eps))
    if mode_norm == "sum_abs":
        return float(np.sum(np.abs(B) * W))
    if mode_norm == "topk_mean":
        flat_values = B.ravel()
        flat_weights = W.ravel()
        valid_idx = np.flatnonzero(flat_weights > 0.0)
        if valid_idx.size == 0:
            return 0.0
        k = int(max(1, min(topk, valid_idx.size)))
        ranked = valid_idx[np.argsort(np.abs(flat_values[valid_idx]))[::-1][:k]]
        return float(
            np.sum(flat_values[ranked] * flat_weights[ranked])
            / (np.sum(flat_weights[ranked]) + eps)
        )

    raise ValueError("voxel_agg_mode must be one of {'mean', 'sum_abs', 'topk_mean'}.")


def _weighted_voxel_coupling_from_residue_dei(
    D_res: np.ndarray,
    voxel_weights: np.ndarray,
    mode: str = "mean",
    topk: int = 5,
    eps: float = 1e-12,
) -> np.ndarray:
    """
    Voxel directional coupling from soft memberships:
        C_vox[V,W] = aggregate(D_res * membership[:,V] * membership[:,W])
    """
    D = np.asarray(D_res, dtype=float)
    W = np.asarray(voxel_weights, dtype=float)
    n_vox = W.shape[1]
    C_vox = np.zeros((n_vox, n_vox), dtype=float)

    for v in range(n_vox):
        target_weights = W[:, v]
        if not np.any(target_weights > 0.0):
            continue
        for w in range(n_vox):
            source_weights = W[:, w]
            if not np.any(source_weights > 0.0):
                continue
            block_weights = np.outer(target_weights, source_weights)
            C_vox[v, w] = _weighted_aggregate_block(
                block=D,
                weights=block_weights,
                mode=mode,
                topk=topk,
                eps=eps,
            )

    np.fill_diagonal(C_vox, 0.0)
    C_vox[~np.isfinite(C_vox)] = 0.0
    return C_vox


def _voxel_to_residue_from_dei(
    D_res: np.ndarray,
    voxel_members: Sequence[np.ndarray],
) -> np.ndarray:
    """Aggregate source-residue columns in D_res to voxel->residue directional effects."""
    D = np.asarray(D_res, dtype=float)
    n_vox = len(voxel_members)
    n_res = D.shape[0]
    out = np.zeros((n_vox, n_res), dtype=float)

    for v in range(n_vox):
        rv = np.asarray(voxel_members[v], dtype=int)
        if rv.size == 0:
            continue
        out[v] = np.mean(D[:, rv], axis=1)

    out[~np.isfinite(out)] = 0.0
    return out


def _weighted_voxel_to_residue_from_dei(
    D_res: np.ndarray,
    voxel_weights: np.ndarray,
    eps: float = 1e-12,
) -> np.ndarray:
    """Aggregate source-residue columns in D_res using soft voxel memberships."""
    D = np.asarray(D_res, dtype=float)
    W = np.asarray(voxel_weights, dtype=float)
    n_vox = W.shape[1]
    n_res = D.shape[0]
    out = np.zeros((n_vox, n_res), dtype=float)

    for v in range(n_vox):
        weights = W[:, v]
        total = float(np.sum(weights))
        if total <= eps:
            continue
        out[v] = D @ weights / total

    out[~np.isfinite(out)] = 0.0
    return out


def _symmetric_affinity_from_directional_coupling(
    coupling: np.ndarray,
    eps: float = 1e-12,
    scale_percentile: float = 95.0,
) -> tuple[np.ndarray, float]:
    """
    Symmetric affinity from directional coupling:
        Affinity[V,W] = |C[V,W]| + |C[W,V]|
    with robust normalization by percentile scale.
    """
    C = np.asarray(coupling, dtype=float)
    affinity = np.abs(C) + np.abs(C.T)
    np.fill_diagonal(affinity, 0.0)

    positive = affinity[affinity > 0.0]
    if positive.size > 0:
        scale = float(np.percentile(positive, float(scale_percentile)))
        if (not np.isfinite(scale)) or scale <= eps:
            scale = float(np.max(positive))
    else:
        scale = float(np.max(affinity))

    affinity = np.clip(affinity / (scale + eps), 0.0, 1.0)
    np.fill_diagonal(affinity, 0.0)
    affinity[~np.isfinite(affinity)] = 0.0
    return affinity, scale


def _top_block_contributors(
    D_res: np.ndarray,
    rv: np.ndarray,
    rw: np.ndarray,
    residue_labels: Sequence[str],
    topk: int,
) -> list[tuple[float, str, str]]:
    """Top residue pairs inside a voxel block by directional score D_res[i,j]."""
    block = np.asarray(D_res[np.ix_(rv, rw)], dtype=float)
    if block.size == 0:
        return []

    flat = block.ravel()
    k = int(max(1, min(topk, flat.size)))
    idx_sorted = np.argsort(flat)[::-1][:k]
    n_cols = block.shape[1]

    # Tuple is (score, target_residue_label_i, source_residue_label_j).
    pairs: list[tuple[float, str, str]] = []
    for idx in idx_sorted:
        i_loc = int(idx // n_cols)
        j_loc = int(idx % n_cols)
        i = int(rv[i_loc])
        j = int(rw[j_loc])
        pairs.append((float(block[i_loc, j_loc]), str(residue_labels[i]), str(residue_labels[j])))
    return pairs


def _build_allostery_summary_text(
    D_res: np.ndarray,
    residue_labels: Sequence[str],
    voxel_coupling: np.ndarray,
    voxel_affinity: np.ndarray,
    voxel_members: Sequence[np.ndarray],
    residue_assignments: pd.DataFrame,
    voxel_summary: pd.DataFrame,
    min_support_res: int,
    voxel_agg_mode: str,
    voxel_agg_topk: int,
    affinity_scale: float,
    top_edges: int = 20,
    top_residue_pairs: int = 5,
) -> str:
    """Build a human-readable summary with voxel definitions and DEI-derived couplings."""
    voxel_defs = {}
    max_voxel = 0
    for row in voxel_summary.itertuples(index=False):
        row_data = row._asdict()
        voxel_label = int(row.voxel_id) + 1
        max_voxel = max(max_voxel, voxel_label)
        if {"voxel_ix", "voxel_iy", "voxel_iz"}.issubset(row_data):
            location = (
                f"grid=({int(row_data['voxel_ix'])},"
                f"{int(row_data['voxel_iy'])},"
                f"{int(row_data['voxel_iz'])})"
            )
        else:
            location = (
                f"centroid=({float(row_data.get('centroid_x', np.nan)):.3f},"
                f"{float(row_data.get('centroid_y', np.nan)):.3f},"
                f"{float(row_data.get('centroid_z', np.nan)):.3f})"
            )
        voxel_defs[voxel_label] = {
            "location": location,
            "n_residues": int(row.n_residues),
            "supported_sources": int(row.supported_source_residues),
            "members": [],
        }

    for row in residue_assignments.itertuples(index=False):
        voxel_label = int(row.voxel_id) + 1
        voxel_defs.setdefault(
            voxel_label,
            {"location": "location=(?,?,?)", "n_residues": 0, "supported_sources": 0, "members": []},
        )["members"].append(str(row.res_label))
        max_voxel = max(max_voxel, voxel_label)

    lines = []
    lines.append("Allostery Summary")
    lines.append("================")
    lines.append("")
    lines.append("Excitation convention:")
    lines.append("Z[t,i] = -1 (excited), Z[t,i] = +1 (not excited)")
    lines.append("")
    lines.append("Definitions:")
    lines.append(f"Residue source support cutoff: #{{t : Z[t,j] == -1}} >= {int(min_support_res)}")
    if "n_voxel_memberships" in residue_assignments.columns:
        mean_memberships = float(np.mean(residue_assignments["n_voxel_memberships"]))
        lines.append(
            "Voxel membership mode: overlapping soft voxels "
            f"(mean memberships per residue={mean_memberships:.3f})"
        )
    else:
        lines.append("Voxel membership mode: hard one-residue-to-one-voxel assignment")
    lines.append("Residue DEI:")
    lines.append("D_res[i,j] = p(i excited | j excited) * log((p(i excited|j excited)+eps)/(p(i excited)+eps))")
    lines.append(f"Voxel coupling aggregation mode: {voxel_agg_mode} (topk={int(voxel_agg_topk)})")
    lines.append("Directional voxel coupling: C_vox[V,W] = aggregate(D_res block for V x W)")
    lines.append("Symmetric affinity:")
    lines.append("Affinity[V,W] = |C[V,W]| + |C[W,V]|, then robustly scaled and clipped to [0,1]")
    lines.append(f"Affinity scale value used: {affinity_scale:.12f}")
    lines.append("")
    lines.append("Voxel Definitions")
    lines.append("-----------------")

    for voxel_label in range(1, max_voxel + 1):
        meta = voxel_defs.get(
            voxel_label,
            {"location": "location=(?,?,?)", "n_residues": 0, "supported_sources": 0, "members": []},
        )
        members = ", ".join(meta["members"]) if meta["members"] else "(none)"
        lines.append(
            f"Voxel {voxel_label}: {meta['location']}, n_residues={meta['n_residues']}, "
            f"supported_sources={meta['supported_sources']}, dominant_residues=[{members}]"
        )

    lines.append("")
    lines.append("Directional Voxel Coupling")
    lines.append("--------------------------")
    lines.append("Format: Voxel i -> Voxel j => C_vox(i,j)")

    C = np.asarray(voxel_coupling, dtype=float)
    ranked_pairs: list[tuple[float, float, int, int]] = []
    for v in range(C.shape[0]):
        for w in range(C.shape[1]):
            if v == w:
                continue
            score = float(C[v, w])
            if abs(score) <= 0.0 or not np.isfinite(score):
                continue
            ranked_pairs.append((abs(score), score, v, w))

    ranked_pairs.sort(key=lambda x: x[0], reverse=True)
    if not ranked_pairs:
        lines.append("No non-zero voxel-to-voxel directional couplings were found.")
    else:
        for _, score, v, w in ranked_pairs[:top_edges]:
            src_voxel = int(v) + 1
            tgt_voxel = int(w) + 1
            lines.append(
                f"Voxel {src_voxel} -> Voxel {tgt_voxel} => {score:.12f}"
            )
            rv = np.asarray(voxel_members[v], dtype=int)
            rw = np.asarray(voxel_members[w], dtype=int)
            top_pairs = _top_block_contributors(
                D_res=D_res,
                rv=rv,
                rw=rw,
                residue_labels=residue_labels,
                topk=top_residue_pairs,
            )
            if not top_pairs:
                lines.append("  contributors: none")
                continue
            contributors = "; ".join(
                [f"({tgt_res}<-{src_res})={pair_score:.6f}" for pair_score, tgt_res, src_res in top_pairs]
            )
            lines.append(f"  contributors: {contributors}")

    lines.append("")
    lines.append("Strongest Symmetric Affinities")
    lines.append("------------------------------")
    lines.append("Format: (Voxel i, Voxel j) => Affinity(i,j)")

    A = np.asarray(voxel_affinity, dtype=float)
    affinity_pairs: list[tuple[float, int, int]] = []
    for v in range(A.shape[0]):
        for w in range(v + 1, A.shape[1]):
            score = float(A[v, w])
            if score <= 0.0 or not np.isfinite(score):
                continue
            affinity_pairs.append((score, v + 1, w + 1))

    affinity_pairs.sort(key=lambda x: x[0], reverse=True)
    if not affinity_pairs:
        lines.append("No positive symmetric affinities were found.")
    else:
        for score, v, w in affinity_pairs[:top_edges]:
            lines.append(f"(Voxel {v}, Voxel {w}) => {score:.12f}")

    return "\n".join(lines) + "\n"


def _save_voxel_affinity_heatmap(
    affinity_df: pd.DataFrame,
    output_png: Path,
    output_pdf: Optional[Path] = None,
) -> tuple[Path, Optional[Path]]:
    """Save a publication-style voxel-vs-voxel affinity heatmap."""
    import matplotlib.pyplot as plt
    import seaborn as sns

    n_vox = int(affinity_df.shape[0])
    side = min(max(0.32 * max(n_vox, 1) + 3.2, 6.0), 14.0)
    fig, ax = plt.subplots(figsize=(side, side))
    cmap = sns.color_palette("light:darkblue", as_cmap=True)
    sns.heatmap(
        affinity_df,
        ax=ax,
        cmap=cmap,
        vmin=0.0,
        vmax=1.0,
        square=True,
        linewidths=0.35,
        linecolor="#f0f0f0",
        cbar_kws={"label": "Voxel affinity", "shrink": 0.78, "pad": 0.02},
    )
    ax.set_title("Voxel-vs-voxel affinity matrix", fontsize=13, fontweight="bold", pad=14)
    ax.set_xlabel("Voxel")
    ax.set_ylabel("Voxel")
    ax.tick_params(axis="x", rotation=90, labelsize=8)
    ax.tick_params(axis="y", rotation=0, labelsize=8)
    fig.tight_layout()

    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_png, dpi=600, bbox_inches="tight")
    saved_pdf = None
    if output_pdf is not None:
        fig.savefig(output_pdf, format="pdf", bbox_inches="tight")
        saved_pdf = output_pdf
    plt.close(fig)
    return output_png, saved_pdf


def _save_region_affinity_heatmap(
    affinity_df: pd.DataFrame,
    output_png: Path,
    output_pdf: Optional[Path] = None,
    title: str = "Region-vs-region affinity matrix",
    colorbar_label: str = "Region affinity",
) -> tuple[Path, Optional[Path]]:
    """Save a publication-style region-vs-region affinity heatmap."""
    import matplotlib.pyplot as plt
    import seaborn as sns

    n_regions = int(affinity_df.shape[0])
    side = min(max(0.32 * max(n_regions, 1) + 3.2, 6.0), 14.0)
    fig, ax = plt.subplots(figsize=(side, side))
    cmap = sns.color_palette("light:darkblue", as_cmap=True)
    sns.heatmap(
        affinity_df,
        ax=ax,
        cmap=cmap,
        vmin=0.0,
        vmax=1.0,
        square=True,
        linewidths=0.35,
        linecolor="#f0f0f0",
        cbar_kws={"label": colorbar_label, "shrink": 0.78, "pad": 0.02},
    )
    ax.set_title(title, fontsize=13, fontweight="bold", pad=14)
    ax.set_xlabel("Motif region")
    ax.set_ylabel("Motif region")
    ax.tick_params(axis="x", rotation=90, labelsize=8)
    ax.tick_params(axis="y", rotation=0, labelsize=8)
    fig.tight_layout()

    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_png, dpi=600, bbox_inches="tight")
    saved_pdf = None
    if output_pdf is not None:
        fig.savefig(output_pdf, format="pdf", bbox_inches="tight")
        saved_pdf = output_pdf
    plt.close(fig)
    return output_png, saved_pdf


def _edge_key(i: int, j: int) -> tuple[int, int]:
    return (int(i), int(j)) if int(i) < int(j) else (int(j), int(i))


def _edge_term_value(
    terms: dict[tuple[int, int], dict[str, float]],
    i: int,
    j: int,
    name: str,
) -> float:
    return float(terms.get(_edge_key(i, j), {}).get(name, 0.0))


def _noncovalent_strength(
    terms: dict[tuple[int, int], dict[str, float]],
    i: int,
    j: int,
) -> float:
    pair_terms = terms.get(_edge_key(i, j), {})
    return float(
        max(
            pair_terms.get("kappa_HB", 0.0),
            pair_terms.get("kappa_salt", 0.0),
            pair_terms.get("kappa_vdW", 0.0),
        )
    )


def _build_motif_regions_from_structure(
    pdb_path: str,
    residue_labels: Sequence[str],
    chain: Optional[str] = "A",
    protein_only: bool = True,
    peptide_threshold: float = 0.5,
    noncovalent_threshold: float = 0.0,
    include_dimers: bool = True,
    include_peptide_trimers: bool = True,
    include_triangulated_triads: bool = True,
    eps: float = 1e-12,
) -> tuple[pd.DataFrame, pd.DataFrame, np.ndarray, list[np.ndarray]]:
    """
    Build motif regions from peptide dimers, peptide trimers, and
    triangulated peptide/non-covalent triads.
    """
    import itertools
    import MDAnalysis as mda
    from protqc.coupling_factor import build_Kij_pairs

    labels = [str(label) for label in residue_labels]
    n_res = len(labels)
    universe = mda.Universe(str(pdb_path))
    _, kijs, terms = build_Kij_pairs(
        universe,
        protein_only=protein_only,
        chain=chain,
        return_terms=True,
    )

    peptide_edges: set[tuple[int, int]] = set()
    noncovalent_edges: set[tuple[int, int]] = set()
    for i in range(n_res):
        for j in range(i + 1, n_res):
            if _edge_term_value(terms, i, j, "kappa_BB") >= peptide_threshold:
                peptide_edges.add((i, j))
            if _noncovalent_strength(terms, i, j) > noncovalent_threshold:
                noncovalent_edges.add((i, j))

    regions: list[dict] = []
    seen: set[tuple[str, tuple[int, ...]]] = set()

    def add_region(region_type: str, members: Sequence[int]) -> None:
        members_tuple = tuple(sorted(int(x) for x in members))
        key = (region_type, members_tuple)
        if key in seen:
            return
        seen.add(key)
        region_edges = [_edge_key(i, j) for i, j in itertools.combinations(members_tuple, 2)]
        peptide = [edge for edge in region_edges if edge in peptide_edges]
        noncov = [edge for edge in region_edges if edge in noncovalent_edges]
        kij_values = [float(kijs.get(edge, 0.0)) for edge in region_edges]
        noncov_values = [_noncovalent_strength(terms, edge[0], edge[1]) for edge in region_edges]
        regions.append(
            {
                "region_type": region_type,
                "members": members_tuple,
                "motif_size": len(members_tuple),
                "peptide_edges": peptide,
                "noncovalent_edges": noncov,
                "mean_Kij": float(np.mean(kij_values)) if kij_values else 0.0,
                "max_noncovalent_strength": float(np.max(noncov_values)) if noncov_values else 0.0,
            }
        )

    if include_dimers:
        for edge in sorted(peptide_edges):
            add_region("peptide_dimer", edge)

    if include_peptide_trimers:
        for i in range(max(0, n_res - 2)):
            if _edge_key(i, i + 1) in peptide_edges and _edge_key(i + 1, i + 2) in peptide_edges:
                add_region("peptide_trimer", (i, i + 1, i + 2))

    if include_triangulated_triads:
        combined_edges = peptide_edges | noncovalent_edges
        for triad in itertools.combinations(range(n_res), 3):
            triad_edges = [_edge_key(i, j) for i, j in itertools.combinations(triad, 2)]
            has_peptide = any(edge in peptide_edges for edge in triad_edges)
            has_noncovalent = any(edge in noncovalent_edges for edge in triad_edges)
            if not (has_peptide and has_noncovalent):
                continue
            # Require connected support across peptide/non-covalent edges.
            connected_nodes: set[int] = set()
            for edge in triad_edges:
                if edge in combined_edges:
                    connected_nodes.update(edge)
            if len(connected_nodes) == 3:
                add_region("peptide_noncovalent_triad", triad)

    if not regions:
        for edge in sorted(peptide_edges):
            add_region("peptide_dimer", edge)
    if not regions:
        raise ValueError("No peptide dimer or triad motif regions were found.")

    weights = np.zeros((n_res, len(regions)), dtype=float)
    summary_rows = []
    membership_rows = []
    region_members: list[np.ndarray] = []

    for region_id, region in enumerate(regions):
        members = tuple(region["members"])
        member_weight = 1.0 / max(len(members), 1)
        weights[list(members), region_id] = member_weight
        region_members.append(np.asarray(members, dtype=int))
        region_label = f"R{region_id + 1}"
        member_labels = [labels[i] for i in members]

        for i in members:
            membership_rows.append(
                {
                    "region_id": region_id,
                    "region_label": region_label,
                    "region_type": region["region_type"],
                    "res_index": i + 1,
                    "res_label": labels[i],
                    "membership_weight": member_weight,
                }
            )

        summary_rows.append(
            {
                "region_id": region_id,
                "region_label": region_label,
                "region_type": region["region_type"],
                "motif_size": region["motif_size"],
                "residue_indices": ";".join(str(i + 1) for i in members),
                "residue_labels": ";".join(member_labels),
                "peptide_edges": ";".join(f"{i + 1}-{j + 1}" for i, j in region["peptide_edges"]),
                "noncovalent_edges": ";".join(f"{i + 1}-{j + 1}" for i, j in region["noncovalent_edges"]),
                "mean_Kij": region["mean_Kij"],
                "max_noncovalent_strength": region["max_noncovalent_strength"],
            }
        )

    row_sums = weights.sum(axis=1, keepdims=True)
    # A residue can be part of multiple motif regions. Normalize each residue's
    # memberships so region-level aggregations remain comparable.
    weights = weights / np.maximum(row_sums, eps)
    for row in membership_rows:
        res_idx = int(row["res_index"]) - 1
        region_id = int(row["region_id"])
        row["membership_weight"] = float(weights[res_idx, region_id])

    return (
        pd.DataFrame(summary_rows),
        pd.DataFrame(membership_rows),
        weights,
        region_members,
    )


def _top_region_block_contributors(
    D_res: np.ndarray,
    target_members: np.ndarray,
    source_members: np.ndarray,
    residue_labels: Sequence[str],
    topk: int,
) -> list[tuple[float, str, str]]:
    return _top_block_contributors(
        D_res=D_res,
        rv=target_members,
        rw=source_members,
        residue_labels=residue_labels,
        topk=topk,
    )


def _build_motif_allostery_summary_text(
    region_summary: pd.DataFrame,
    region_affinity: np.ndarray,
    region_coupling: np.ndarray,
    region_members: Sequence[np.ndarray],
    residue_labels: Sequence[str],
    D_res: np.ndarray,
    min_support_res: int,
    agg_mode: str,
    agg_topk: int,
    affinity_scale: float,
    top_edges: int = 20,
    top_residue_pairs: int = 5,
) -> str:
    lines = []
    lines.append("Motif Allostery Summary")
    lines.append("=======================")
    lines.append("")
    lines.append("Excitation convention:")
    lines.append("Z[t,i] = -1 (excited), Z[t,i] = +1 (not excited)")
    lines.append("")
    lines.append("Definitions:")
    lines.append(f"Residue source support cutoff: #{{t : Z[t,j] == -1}} >= {int(min_support_res)}")
    lines.append("Motif regions are structural dimer/trimer regions.")
    lines.append("Dimers are peptide-bonded residue pairs.")
    lines.append("Peptide trimers are consecutive peptide-bonded triples.")
    lines.append("Triangulated triads contain at least one peptide edge and at least one non-covalent edge.")
    lines.append(f"Region coupling aggregation mode: {agg_mode} (topk={int(agg_topk)})")
    lines.append("Directional region coupling: C_region[A,B] = aggregate(D_res weighted by motif memberships)")
    lines.append("Symmetric affinity:")
    lines.append("Affinity[A,B] = |C[A,B]| + |C[B,A]|, then robustly scaled and clipped to [0,1]")
    lines.append(f"Affinity scale value used: {affinity_scale:.12f}")
    lines.append("")
    lines.append("Motif Region Definitions")
    lines.append("------------------------")

    for row in region_summary.itertuples(index=False):
        lines.append(
            f"{row.region_label}: type={row.region_type}, size={int(row.motif_size)}, "
            f"residues=[{row.residue_labels}], peptide_edges=[{row.peptide_edges}], "
            f"noncovalent_edges=[{row.noncovalent_edges}], supported_sources={int(row.supported_source_residues)}"
        )

    lines.append("")
    lines.append("Directional Motif Coupling")
    lines.append("--------------------------")
    lines.append("Format: Region A -> Region B => C_region(A,B)")

    C = np.asarray(region_coupling, dtype=float)
    ranked_pairs: list[tuple[float, float, int, int]] = []
    for a in range(C.shape[0]):
        for b in range(C.shape[1]):
            if a == b:
                continue
            score = float(C[a, b])
            if abs(score) <= 0.0 or not np.isfinite(score):
                continue
            ranked_pairs.append((abs(score), score, a, b))

    ranked_pairs.sort(key=lambda item: item[0], reverse=True)
    if not ranked_pairs:
        lines.append("No non-zero motif-to-motif directional couplings were found.")
    else:
        for _, score, a, b in ranked_pairs[:top_edges]:
            src_label = str(region_summary.loc[a, "region_label"])
            tgt_label = str(region_summary.loc[b, "region_label"])
            lines.append(f"{src_label} -> {tgt_label} => {score:.12f}")
            top_pairs = _top_region_block_contributors(
                D_res=D_res,
                target_members=np.asarray(region_members[a], dtype=int),
                source_members=np.asarray(region_members[b], dtype=int),
                residue_labels=residue_labels,
                topk=top_residue_pairs,
            )
            if not top_pairs:
                lines.append("  contributors: none")
                continue
            contributors = "; ".join(
                [f"({target}<-{source})={pair_score:.6f}" for pair_score, target, source in top_pairs]
            )
            lines.append(f"  contributors: {contributors}")

    lines.append("")
    lines.append("Strongest Symmetric Motif Affinities")
    lines.append("------------------------------------")
    lines.append("Format: (Region A, Region B) => Affinity(A,B)")

    A = np.asarray(region_affinity, dtype=float)
    affinity_pairs: list[tuple[float, int, int]] = []
    for a in range(A.shape[0]):
        for b in range(a + 1, A.shape[1]):
            score = float(A[a, b])
            if score <= 0.0 or not np.isfinite(score):
                continue
            affinity_pairs.append((score, a, b))

    affinity_pairs.sort(key=lambda item: item[0], reverse=True)
    if not affinity_pairs:
        lines.append("No positive symmetric motif affinities were found.")
    else:
        for score, a, b in affinity_pairs[:top_edges]:
            label_a = str(region_summary.loc[a, "region_label"])
            label_b = str(region_summary.loc[b, "region_label"])
            lines.append(f"({label_a}, {label_b}) => {score:.12f}")

    return "\n".join(lines) + "\n"


def analyse_motif_allostery(
    Z_samples: np.ndarray,
    residue_labels: Optional[Sequence[str]] = None,
    save_dir: Optional[str] = None,
    save_prefix: str = "motif_allostery",
    random_state: int = 0,
    pdb_path: Optional[str] = None,
    chain: Optional[str] = "A",
    protein_only: bool = True,
    min_support: int = 20,
    min_support_res: Optional[int] = None,
    motif_agg_mode: str = "mean",
    motif_agg_topk: int = 5,
    affinity_scale_percentile: float = 95.0,
    n_clusters: Optional[int] = None,
    peptide_threshold: float = 0.5,
    noncovalent_threshold: float = 0.0,
    include_dimers: bool = True,
    include_peptide_trimers: bool = True,
    include_triangulated_triads: bool = True,
    save_affinity_plot: bool = True,
    eps: float = 1e-12,
) -> dict:
    """
    Motif-based directional allostery using peptide dimers and triads.

    Motif regions replace spatial voxels. Region membership is based on
    peptide connectivity and non-covalent structural interactions from the
    Kij decomposition.
    """
    if pdb_path is None:
        raise ValueError("pdb_path is required for motif allostery.")

    excited = _to_excited_binary(Z_samples)
    _, n_res = excited.shape

    if residue_labels is None:
        residue_labels = [f"R{i + 1}" for i in range(n_res)]
    else:
        if len(residue_labels) != n_res:
            raise ValueError("residue_labels length must match number of residues.")
        residue_labels = [str(label) for label in residue_labels]

    min_support_res_eff = int(min_support if min_support_res is None else min_support_res)
    if min_support_res_eff <= 0:
        raise ValueError("min_support_res must be > 0.")
    if motif_agg_topk <= 0:
        raise ValueError("motif_agg_topk must be > 0.")
    if not (0.0 < float(affinity_scale_percentile) <= 100.0):
        raise ValueError("affinity_scale_percentile must be in (0, 100].")
    if eps <= 0.0:
        raise ValueError("eps must be > 0.")

    region_summary, residue_memberships, region_weights, region_members = _build_motif_regions_from_structure(
        pdb_path=str(pdb_path),
        residue_labels=residue_labels,
        chain=chain,
        protein_only=protein_only,
        peptide_threshold=peptide_threshold,
        noncovalent_threshold=noncovalent_threshold,
        include_dimers=include_dimers,
        include_peptide_trimers=include_peptide_trimers,
        include_triangulated_triads=include_triangulated_triads,
        eps=eps,
    )

    residue_allostery, p_exc, residue_source_support = _residue_directional_excitation_matrix(
        Z_samples=Z_samples,
        min_support_res=min_support_res_eff,
        eps=eps,
    )
    region_to_residue = _weighted_voxel_to_residue_from_dei(
        D_res=residue_allostery,
        voxel_weights=region_weights,
        eps=eps,
    )
    directional_coupling = _weighted_voxel_coupling_from_residue_dei(
        D_res=residue_allostery,
        voxel_weights=region_weights,
        mode=motif_agg_mode,
        topk=motif_agg_topk,
        eps=eps,
    )
    affinity, affinity_scale = _symmetric_affinity_from_directional_coupling(
        coupling=directional_coupling,
        eps=eps,
        scale_percentile=affinity_scale_percentile,
    )

    n_regions = int(region_weights.shape[1])
    if n_regions <= 1:
        n_clusters_eff = 1
    elif n_clusters is None:
        n_clusters_eff = max(2, min(n_regions, 8))
    else:
        n_clusters_eff = int(max(1, min(int(n_clusters), n_regions)))

    region_cluster_labels = _spectral_cluster_from_affinity(
        affinity=affinity,
        n_clusters=n_clusters_eff,
        random_state=random_state,
    )

    source_supported_mask = residue_source_support >= min_support_res_eff
    supported_source_counts = np.asarray(
        [int(np.sum(source_supported_mask[members])) for members in region_members],
        dtype=int,
    )
    supported_source_weights = region_weights[source_supported_mask, :].sum(axis=0)
    active_region_mask = supported_source_weights > 0.0

    region_summary = region_summary.copy()
    region_summary["supported_source_residues"] = supported_source_counts
    region_summary["supported_source_weight"] = supported_source_weights
    region_summary["active"] = active_region_mask
    region_summary["cluster_id"] = region_cluster_labels

    residue_memberships = residue_memberships.copy()
    residue_memberships["cluster_id"] = residue_memberships["region_id"].map(
        {i: int(label) for i, label in enumerate(region_cluster_labels)}
    )

    region_labels = region_summary["region_label"].astype(str).tolist()
    residue_matrix_df = pd.DataFrame(
        residue_allostery,
        index=residue_labels,
        columns=residue_labels,
    )
    region_matrix_df = pd.DataFrame(
        region_to_residue,
        index=region_labels,
        columns=residue_labels,
    )
    membership_df = pd.DataFrame(
        region_weights,
        index=residue_labels,
        columns=region_labels,
    )
    affinity_df = pd.DataFrame(
        affinity,
        index=region_labels,
        columns=region_labels,
    )
    directional_coupling_df = pd.DataFrame(
        directional_coupling,
        index=region_labels,
        columns=region_labels,
    )
    cluster_summary = (
        region_summary.groupby("cluster_id", as_index=False)
        .agg(
            n_regions=("region_id", "count"),
            n_active_regions=("active", "sum"),
            mean_supported_source_weight=("supported_source_weight", "mean"),
        )
        .sort_values("cluster_id")
        .reset_index(drop=True)
    )

    saved_paths = {}
    if save_dir:
        out_dir = Path(save_dir)
        out_dir.mkdir(parents=True, exist_ok=True)

        residue_matrix_path = out_dir / f"{save_prefix}_residue_allostery_matrix.csv"
        region_matrix_path = out_dir / f"{save_prefix}_region_to_residue_matrix.csv"
        region_summary_path = out_dir / f"{save_prefix}_region_summary.csv"
        memberships_path = out_dir / f"{save_prefix}_residue_memberships.csv"
        membership_matrix_path = out_dir / f"{save_prefix}_membership_matrix.csv"
        cluster_summary_path = out_dir / f"{save_prefix}_cluster_summary.csv"
        affinity_path = out_dir / f"{save_prefix}_affinity_matrix.csv"
        affinity_heatmap_png_path = out_dir / f"{save_prefix}_affinity_matrix_heatmap.png"
        affinity_heatmap_pdf_path = out_dir / f"{save_prefix}_affinity_matrix_heatmap.pdf"
        directional_coupling_path = out_dir / f"{save_prefix}_directional_region_coupling.csv"
        summary_path = out_dir / "allostery_summary.txt"

        residue_matrix_df.to_csv(residue_matrix_path)
        region_matrix_df.to_csv(region_matrix_path)
        region_summary.to_csv(region_summary_path, index=False)
        residue_memberships.to_csv(memberships_path, index=False)
        membership_df.to_csv(membership_matrix_path)
        cluster_summary.to_csv(cluster_summary_path, index=False)
        affinity_df.to_csv(affinity_path)
        directional_coupling_df.to_csv(directional_coupling_path)

        saved_affinity_heatmap_png = None
        saved_affinity_heatmap_pdf = None
        if save_affinity_plot:
            saved_affinity_heatmap_png, saved_affinity_heatmap_pdf = _save_region_affinity_heatmap(
                affinity_df=affinity_df,
                output_png=affinity_heatmap_png_path,
                output_pdf=affinity_heatmap_pdf_path,
                title="Motif-region affinity matrix",
                colorbar_label="Motif affinity",
            )

        split_saved_paths = {}
        for group_name, group_mask, group_title in [
            (
                "dimer",
                region_summary["region_type"].astype(str).str.contains("dimer", case=False, regex=False),
                "Dimer-region affinity matrix",
            ),
            (
                "trimer",
                region_summary["region_type"].astype(str).str.contains("trimer", case=False, regex=False)
                | region_summary["region_type"].astype(str).str.contains("triad", case=False, regex=False),
                "Trimer/triad-region affinity matrix",
            ),
        ]:
            group_labels = region_summary.loc[group_mask, "region_label"].astype(str).tolist()
            if not group_labels:
                continue

            group_affinity_df = affinity_df.loc[group_labels, group_labels]
            group_coupling_df = directional_coupling_df.loc[group_labels, group_labels]
            group_region_summary = region_summary.loc[group_mask].copy()

            group_affinity_path = out_dir / f"{save_prefix}_{group_name}_affinity_matrix.csv"
            group_coupling_path = out_dir / f"{save_prefix}_{group_name}_directional_region_coupling.csv"
            group_summary_path = out_dir / f"{save_prefix}_{group_name}_region_summary.csv"
            group_heatmap_png_path = out_dir / f"{save_prefix}_{group_name}_affinity_matrix_heatmap.png"
            group_heatmap_pdf_path = out_dir / f"{save_prefix}_{group_name}_affinity_matrix_heatmap.pdf"

            group_affinity_df.to_csv(group_affinity_path)
            group_coupling_df.to_csv(group_coupling_path)
            group_region_summary.to_csv(group_summary_path, index=False)

            split_saved_paths[f"{group_name}_affinity_matrix"] = str(group_affinity_path)
            split_saved_paths[f"{group_name}_directional_region_coupling"] = str(group_coupling_path)
            split_saved_paths[f"{group_name}_region_summary"] = str(group_summary_path)

            if save_affinity_plot:
                group_png, group_pdf = _save_region_affinity_heatmap(
                    affinity_df=group_affinity_df,
                    output_png=group_heatmap_png_path,
                    output_pdf=group_heatmap_pdf_path,
                    title=group_title,
                    colorbar_label=f"{group_name.capitalize()} affinity",
                )
                split_saved_paths[f"{group_name}_affinity_matrix_heatmap_png"] = str(group_png)
                if group_pdf is not None:
                    split_saved_paths[f"{group_name}_affinity_matrix_heatmap_pdf"] = str(group_pdf)

        summary_text = _build_motif_allostery_summary_text(
            region_summary=region_summary,
            region_affinity=affinity,
            region_coupling=directional_coupling,
            region_members=region_members,
            residue_labels=residue_labels,
            D_res=residue_allostery,
            min_support_res=min_support_res_eff,
            agg_mode=motif_agg_mode,
            agg_topk=motif_agg_topk,
            affinity_scale=affinity_scale,
        )
        summary_path.write_text(summary_text, encoding="utf-8")

        saved_paths = {
            "residue_allostery_matrix": str(residue_matrix_path),
            "region_to_residue_matrix": str(region_matrix_path),
            "region_summary": str(region_summary_path),
            "residue_memberships": str(memberships_path),
            "membership_matrix": str(membership_matrix_path),
            "cluster_summary": str(cluster_summary_path),
            "affinity_matrix": str(affinity_path),
            "directional_region_coupling": str(directional_coupling_path),
            "allostery_summary": str(summary_path),
        }
        if saved_affinity_heatmap_png is not None:
            saved_paths["affinity_matrix_heatmap_png"] = str(saved_affinity_heatmap_png)
        if saved_affinity_heatmap_pdf is not None:
            saved_paths["affinity_matrix_heatmap_pdf"] = str(saved_affinity_heatmap_pdf)
        saved_paths.update(split_saved_paths)

    return {
        "residue_allostery_matrix": residue_allostery,
        "region_to_residue_matrix": region_to_residue,
        "membership_matrix": region_weights,
        "affinity_matrix": affinity,
        "cluster_labels": region_cluster_labels,
        "directional_region_coupling_matrix": directional_coupling,
        "residue_source_support": residue_source_support,
        "active_region_mask": active_region_mask,
        "activation_support": supported_source_weights,
        "affinity_scale": affinity_scale,
        "residue_matrix_df": residue_matrix_df,
        "region_matrix_df": region_matrix_df,
        "membership_df": membership_df,
        "affinity_df": affinity_df,
        "directional_coupling_df": directional_coupling_df,
        "region_summary_df": region_summary,
        "residue_memberships_df": residue_memberships,
        "cluster_summary_df": cluster_summary,
        "saved_paths": saved_paths,
    }


def analyse_allostery(
    Z_samples: np.ndarray,
    residue_coords: Optional[np.ndarray] = None,
    voxel_size: float = 6.0,
    min_local_excited: int = 2,
    distal_distance_cutoff: float = 10.0,
    n_clusters: Optional[int] = None,
    spatial_weight: float = 0.6,
    residue_labels: Optional[Sequence[str]] = None,
    save_dir: Optional[str] = None,
    save_prefix: str = "allostery",
    random_state: int = 0,
    pdb_path: Optional[str] = None,
    chain: Optional[str] = "A",
    protein_only: bool = True,
    anchor_atom: str = "CA",
    min_support: int = 20,
    min_support_res: Optional[int] = None,
    voxel_agg_mode: str = "mean",
    voxel_agg_topk: int = 5,
    affinity_scale_percentile: float = 95.0,
    overlapping_voxels: bool = True,
    voxel_overlap_radius: Optional[float] = None,
    voxel_kernel_sigma: Optional[float] = None,
    voxel_membership_threshold: float = 0.05,
    max_voxel_memberships: Optional[int] = 8,
    save_affinity_plot: bool = True,
    eps: float = 1e-12,
) -> dict:
    """
    Voxel-based directional allostery where voxel->voxel influence is aggregated
    from residue-level directional excitation information (DEI) blocks.

    Parameters
    ----------
    Z_samples : np.ndarray
        Shape (shots, n_res), with values in {-1,+1} or {0,1}.
        Excitation convention is preserved: -1 means excited.
    residue_coords : Optional[np.ndarray]
        Shape (n_res, 3) Cartesian coordinates per residue.
        If None, coordinates are extracted from `pdb_path`.
    voxel_size : float
        Cartesian voxel edge length (angstrom).
    min_local_excited : int
        Legacy parameter retained for compatibility; not used in this method.
    distal_distance_cutoff : float
        Legacy parameter retained for compatibility; not used in this method.
    n_clusters : Optional[int]
        Number of voxel clusters for spectral clustering.
    spatial_weight : float
        Legacy parameter retained for compatibility; not used in this method.
    residue_labels : Optional[Sequence[str]]
        Optional labels (length n_res). Defaults to R1, R2, ...
    save_dir : Optional[str]
        If provided, save CSV outputs to this directory.
    save_prefix : str
        Prefix used for saved CSV file names.
    random_state : int
        Random seed used for spectral clustering k-means step.
    pdb_path : Optional[str]
        PDB path used to extract residue coordinates when `residue_coords` is None.
    chain : Optional[str]
        Chain ID filter (uses `chainID`, with fallback to atom `chainIDs` and `segid`).
    protein_only : bool
        If True, keep only standard amino-acid residues.
    anchor_atom : str
        Atom name used as residue coordinate anchor (default: CA).
    min_support : int
        Legacy alias for `min_support_res` when `min_support_res` is None.
    min_support_res : Optional[int]
        Minimum source-residue excitation support for DEI column updates.
    voxel_agg_mode : str
        DEI block aggregation mode for C_vox: {"mean", "sum_abs", "topk_mean"}.
    voxel_agg_topk : int
        Top-k parameter for "topk_mean" aggregation.
    affinity_scale_percentile : float
        Percentile used for robust affinity scaling.
    overlapping_voxels : bool
        If True, residues receive Gaussian memberships in multiple nearby
        voxels instead of a single hard voxel assignment.
    voxel_overlap_radius : Optional[float]
        Maximum distance from a residue coordinate to a voxel center for
        non-zero membership. Defaults to voxel_size.
    voxel_kernel_sigma : Optional[float]
        Gaussian length scale for soft voxel weights. Defaults to
        0.75 * voxel_size.
    voxel_membership_threshold : float
        Relative threshold for keeping non-dominant memberships. For example,
        0.05 keeps memberships at least 5% of the strongest membership for
        that residue.
    max_voxel_memberships : Optional[int]
        Maximum number of voxel memberships retained per residue. Use None
        to keep all memberships above threshold.
    save_affinity_plot : bool
        If True and save_dir is provided, save a voxel-vs-voxel affinity
        heatmap as PNG and PDF.
    eps : float
        Numerical stability guard added inside log ratios and normalization.

    Returns
    -------
    dict
        Dictionary with matrices, assignment tables, and optional saved paths.
    """
    excited = _to_excited_binary(Z_samples)
    if residue_coords is None:
        if pdb_path is None:
            raise ValueError("Provide either residue_coords or pdb_path.")
        coords, pdb_labels = _residue_coordinates_from_pdb(
            pdb_path=pdb_path,
            chain=chain,
            protein_only=protein_only,
            anchor_atom=anchor_atom,
        )
        if residue_labels is None:
            residue_labels = pdb_labels
    else:
        coords = np.asarray(residue_coords, dtype=float)

    if coords.ndim != 2 or coords.shape[1] != 3:
        raise ValueError("residue_coords must have shape (n_res, 3).")

    _, n_res = excited.shape
    if coords.shape[0] != n_res:
        raise ValueError("residue_coords first dimension must match Z_samples residues.")
    if voxel_size <= 0.0:
        raise ValueError("voxel_size must be > 0.")
    min_support_res_eff = int(min_support if min_support_res is None else min_support_res)
    if min_support_res_eff <= 0:
        raise ValueError("min_support_res must be > 0.")
    if voxel_agg_topk <= 0:
        raise ValueError("voxel_agg_topk must be > 0.")
    if not (0.0 < float(affinity_scale_percentile) <= 100.0):
        raise ValueError("affinity_scale_percentile must be in (0, 100].")
    if eps <= 0.0:
        raise ValueError("eps must be > 0.")

    if residue_labels is None:
        residue_labels = [f"R{i + 1}" for i in range(n_res)]
    else:
        if len(residue_labels) != n_res:
            raise ValueError("residue_labels length must match number of residues.")
        residue_labels = [str(x) for x in residue_labels]

    if overlapping_voxels:
        (
            dominant_voxel_ids,
            voxel_grid,
            voxel_centroids,
            voxel_effective_sizes,
            voxel_member_counts,
            voxel_members,
            voxel_prob_weights,
        ) = _overlapping_grid_voxelise(
            coords=coords,
            voxel_size=voxel_size,
            overlap_radius=voxel_overlap_radius,
            kernel_sigma=voxel_kernel_sigma,
            membership_threshold=voxel_membership_threshold,
            max_memberships=max_voxel_memberships,
            eps=eps,
        )
    else:
        dominant_voxel_ids, voxel_grid, voxel_centroids, voxel_sizes, voxel_members = _grid_voxelise(
            coords=coords,
            voxel_size=voxel_size,
            eps=eps,
        )
        voxel_effective_sizes = voxel_sizes.astype(float)
        voxel_member_counts = voxel_sizes.astype(int)
        voxel_prob_weights = np.zeros((n_res, len(voxel_members)), dtype=float)
        voxel_prob_weights[np.arange(n_res), dominant_voxel_ids] = 1.0
    n_voxels = int(voxel_centroids.shape[0])

    residue_allostery, p_exc, residue_source_support = _residue_directional_excitation_matrix(
        Z_samples=Z_samples,
        min_support_res=min_support_res_eff,
        eps=eps,
    )
    voxel_to_residue = _weighted_voxel_to_residue_from_dei(
        D_res=residue_allostery,
        voxel_weights=voxel_prob_weights,
        eps=eps,
    )
    directional_coupling = _weighted_voxel_coupling_from_residue_dei(
        D_res=residue_allostery,
        voxel_weights=voxel_prob_weights,
        mode=voxel_agg_mode,
        topk=voxel_agg_topk,
        eps=eps,
    )
    affinity, affinity_scale = _symmetric_affinity_from_directional_coupling(
        coupling=directional_coupling,
        eps=eps,
        scale_percentile=affinity_scale_percentile,
    )

    if n_voxels <= 1:
        n_clusters_eff = 1
    elif n_clusters is None:
        n_clusters_eff = max(2, min(n_voxels, 8))
    else:
        n_clusters_eff = int(max(1, min(int(n_clusters), n_voxels)))

    voxel_cluster_labels = _spectral_cluster_from_affinity(
        affinity=affinity,
        n_clusters=n_clusters_eff,
        random_state=random_state,
    )

    cluster_labels = voxel_cluster_labels[dominant_voxel_ids]
    source_supported_mask = residue_source_support >= min_support_res_eff
    voxel_supported_source_counts = np.array(
        [int(np.sum(source_supported_mask[members])) for members in voxel_members],
        dtype=int,
    )
    voxel_supported_source_weight = voxel_prob_weights[source_supported_mask, :].sum(axis=0)
    active_voxel_mask = voxel_supported_source_weight > 0.0

    residue_matrix_df = pd.DataFrame(
        residue_allostery,
        index=residue_labels,
        columns=residue_labels,
    )
    voxel_labels = [f"V{i + 1}" for i in range(n_voxels)]
    voxel_matrix_df = pd.DataFrame(
        voxel_to_residue,
        index=voxel_labels,
        columns=residue_labels,
    )
    soft_membership_df = pd.DataFrame(
        voxel_prob_weights,
        index=residue_labels,
        columns=voxel_labels,
    )
    affinity_df = pd.DataFrame(
        affinity,
        index=voxel_labels,
        columns=voxel_labels,
    )
    directional_coupling_df = pd.DataFrame(
        directional_coupling,
        index=voxel_labels,
        columns=voxel_labels,
    )

    residue_assignments = pd.DataFrame(
        {
            "res_index": np.arange(1, n_res + 1),
            "res_label": residue_labels,
            "x": coords[:, 0],
            "y": coords[:, 1],
            "z": coords[:, 2],
            "voxel_id": dominant_voxel_ids,
            "voxel_membership": voxel_prob_weights[np.arange(n_res), dominant_voxel_ids],
            "n_voxel_memberships": np.sum(voxel_prob_weights > 0.0, axis=1),
            "region_centroid_x": voxel_centroids[dominant_voxel_ids, 0],
            "region_centroid_y": voxel_centroids[dominant_voxel_ids, 1],
            "region_centroid_z": voxel_centroids[dominant_voxel_ids, 2],
            "region_size": voxel_member_counts[dominant_voxel_ids],
            "region_effective_size": voxel_effective_sizes[dominant_voxel_ids],
            "source_support": residue_source_support,
            "source_supported": residue_source_support >= min_support_res_eff,
            "activation_support": voxel_supported_source_counts[dominant_voxel_ids],
            "activation_support_weight": voxel_supported_source_weight[dominant_voxel_ids],
            "voxel_active": active_voxel_mask[dominant_voxel_ids],
            "cluster_id": cluster_labels,
            "p_excited": p_exc,
        }
    )

    membership_res_idx, membership_vox_idx = np.nonzero(voxel_prob_weights > 0.0)
    residue_voxel_memberships = pd.DataFrame(
        {
            "res_index": membership_res_idx + 1,
            "res_label": [residue_labels[int(i)] for i in membership_res_idx],
            "voxel_id": membership_vox_idx,
            "voxel_label": [f"V{int(v) + 1}" for v in membership_vox_idx],
            "membership_weight": voxel_prob_weights[membership_res_idx, membership_vox_idx],
            "dominant_voxel": dominant_voxel_ids[membership_res_idx] == membership_vox_idx,
            "voxel_centroid_x": voxel_centroids[membership_vox_idx, 0],
            "voxel_centroid_y": voxel_centroids[membership_vox_idx, 1],
            "voxel_centroid_z": voxel_centroids[membership_vox_idx, 2],
            "cluster_id": voxel_cluster_labels[membership_vox_idx],
        }
    )

    voxel_summary = pd.DataFrame(
        {
            "voxel_id": np.arange(n_voxels, dtype=int),
            "centroid_x": voxel_centroids[:, 0],
            "centroid_y": voxel_centroids[:, 1],
            "centroid_z": voxel_centroids[:, 2],
            "voxel_ix": voxel_grid[:, 0],
            "voxel_iy": voxel_grid[:, 1],
            "voxel_iz": voxel_grid[:, 2],
            "n_residues": voxel_member_counts,
            "effective_size": voxel_effective_sizes.astype(float),
            "supported_source_residues": voxel_supported_source_counts,
            "supported_source_weight": voxel_supported_source_weight,
            "active": active_voxel_mask,
            "cluster_id": voxel_cluster_labels,
        }
    )

    cluster_summary = (
        residue_assignments.groupby("cluster_id", as_index=False)
        .agg(
            n_residues=("res_index", "count"),
            n_voxels=("voxel_id", "nunique"),
            mean_p_excited=("p_excited", "mean"),
        )
        .sort_values("cluster_id")
        .reset_index(drop=True)
    )

    saved_paths = {}
    if save_dir:
        out_dir = Path(save_dir)
        out_dir.mkdir(parents=True, exist_ok=True)

        residue_matrix_path = out_dir / f"{save_prefix}_residue_allostery_matrix.csv"
        voxel_matrix_path = out_dir / f"{save_prefix}_voxel_to_residue_matrix.csv"
        assignments_path = out_dir / f"{save_prefix}_residue_assignments.csv"
        memberships_path = out_dir / f"{save_prefix}_residue_voxel_memberships.csv"
        voxel_summary_path = out_dir / f"{save_prefix}_voxel_summary.csv"
        soft_membership_path = out_dir / f"{save_prefix}_soft_membership.csv"
        cluster_summary_path = out_dir / f"{save_prefix}_cluster_summary.csv"
        affinity_path = out_dir / f"{save_prefix}_affinity_matrix.csv"
        affinity_heatmap_png_path = out_dir / f"{save_prefix}_affinity_matrix_heatmap.png"
        affinity_heatmap_pdf_path = out_dir / f"{save_prefix}_affinity_matrix_heatmap.pdf"
        directional_coupling_path = out_dir / f"{save_prefix}_directional_voxel_coupling.csv"
        summary_path = out_dir / "allostery_summary.txt"

        residue_matrix_df.to_csv(residue_matrix_path)
        voxel_matrix_df.to_csv(voxel_matrix_path)
        residue_assignments.to_csv(assignments_path, index=False)
        residue_voxel_memberships.to_csv(memberships_path, index=False)
        voxel_summary.to_csv(voxel_summary_path, index=False)
        soft_membership_df.to_csv(soft_membership_path)
        cluster_summary.to_csv(cluster_summary_path, index=False)
        affinity_df.to_csv(affinity_path)
        saved_affinity_heatmap_png = None
        saved_affinity_heatmap_pdf = None
        if save_affinity_plot:
            saved_affinity_heatmap_png, saved_affinity_heatmap_pdf = _save_voxel_affinity_heatmap(
                affinity_df=affinity_df,
                output_png=affinity_heatmap_png_path,
                output_pdf=affinity_heatmap_pdf_path,
            )
        directional_coupling_df.to_csv(directional_coupling_path)
        summary_text = _build_allostery_summary_text(
            D_res=residue_allostery,
            residue_labels=residue_labels,
            voxel_coupling=directional_coupling,
            voxel_affinity=affinity,
            voxel_members=voxel_members,
            residue_assignments=residue_assignments,
            voxel_summary=voxel_summary,
            min_support_res=min_support_res_eff,
            voxel_agg_mode=voxel_agg_mode,
            voxel_agg_topk=voxel_agg_topk,
            affinity_scale=affinity_scale,
        )
        summary_path.write_text(summary_text, encoding="utf-8")

        saved_paths = {
            "residue_allostery_matrix": str(residue_matrix_path),
            "voxel_to_residue_matrix": str(voxel_matrix_path),
            "residue_assignments": str(assignments_path),
            "residue_voxel_memberships": str(memberships_path),
            "voxel_summary": str(voxel_summary_path),
            "soft_membership": str(soft_membership_path),
            "cluster_summary": str(cluster_summary_path),
            "affinity_matrix": str(affinity_path),
            "directional_voxel_coupling": str(directional_coupling_path),
            "allostery_summary": str(summary_path),
        }
        if saved_affinity_heatmap_png is not None:
            saved_paths["affinity_matrix_heatmap_png"] = str(saved_affinity_heatmap_png)
        if saved_affinity_heatmap_pdf is not None:
            saved_paths["affinity_matrix_heatmap_pdf"] = str(saved_affinity_heatmap_pdf)

    return {
        "residue_allostery_matrix": residue_allostery,
        "voxel_to_residue_matrix": voxel_to_residue,
        "soft_membership_matrix": voxel_prob_weights,
        "affinity_matrix": affinity,
        "cluster_labels": cluster_labels,
        "voxel_cluster_labels": voxel_cluster_labels,
        "voxel_ids": dominant_voxel_ids,
        "directional_voxel_coupling_matrix": directional_coupling,
        "residue_source_support": residue_source_support,
        "active_voxel_mask": active_voxel_mask,
        "activation_support": voxel_supported_source_counts,
        "affinity_scale": affinity_scale,
        "residue_matrix_df": residue_matrix_df,
        "voxel_matrix_df": voxel_matrix_df,
        "soft_membership_df": soft_membership_df,
        "affinity_df": affinity_df,
        "directional_coupling_df": directional_coupling_df,
        "residue_assignments_df": residue_assignments,
        "residue_voxel_memberships_df": residue_voxel_memberships,
        "voxel_summary_df": voxel_summary,
        "cluster_summary_df": cluster_summary,
        "saved_paths": saved_paths,
    }


def analyze_allostery(*args, **kwargs):
    """US spelling alias for analyse_allostery."""
    return analyse_allostery(*args, **kwargs)


def run_basin(args):
    import json
    import matplotlib.pyplot as plt
    from protqc.plots import _save_figure_with_pdf, set_academic_style

    z_df = pd.read_csv(args.z_samples)
    if z_df.shape[1] < 2:
        raise ValueError("z_samples CSV must contain 'shot' plus at least one residue column.")

    first_col = str(z_df.columns[0]).strip().lower()
    if first_col != "shot":
        raise ValueError("First column in z_samples CSV must be 'shot'.")

    residue_labels = list(z_df.columns[1:])
    Z = z_df.iloc[:, 1:].to_numpy(dtype=float)

    valid_mask = np.isin(Z, (-1.0, 1.0))
    if not np.all(valid_mask):
        bad_values = np.unique(Z[~valid_mask])
        raise ValueError(f"z_samples contains values outside ±1 basis: {bad_values.tolist()}")

    if (args.E0 is None) ^ (args.E1 is None):
        raise ValueError("Provide both --E0 and --E1 together, or omit both to use energies_ref defaults.")

    if args.E0 is not None and args.E1 is not None:
        E0 = np.asarray(np.load(args.E0, allow_pickle=False), dtype=float).reshape(-1)
        E1 = np.asarray(np.load(args.E1, allow_pickle=False), dtype=float).reshape(-1)
    else:
        from protqc.circuit import energies_ref as circuit_energies_ref

        residue_names = []
        for label in residue_labels:
            token = str(label).strip()
            if "-" in token:
                token = token.split("-", 1)[1].strip()
            resname = token.upper()
            if not resname:
                raise ValueError(f"Unable to parse residue name from label: {label!r}")
            residue_names.append(resname)

        ref_energy = circuit_energies_ref()
        missing = sorted({res for res in residue_names if res not in ref_energy})
        if missing:
            raise ValueError(
                f"Residues not found in energies_ref table: {missing}. "
                "Provide --E0 and --E1 explicitly for custom residue mappings."
            )

        E0 = np.asarray([ref_energy[res]["e0"] for res in residue_names], dtype=float)
        E1 = np.asarray([ref_energy[res]["e1"] for res in residue_names], dtype=float)

    n_residues = Z.shape[1]
    if E0.shape[0] != n_residues or E1.shape[0] != n_residues:
        raise ValueError(
            f"Length mismatch: residues={n_residues}, len(E0)={E0.shape[0]}, len(E1)={E1.shape[0]}"
        )

    per_residue_energy = np.where(Z == 1.0, E0[np.newaxis, :], E1[np.newaxis, :])
    E = per_residue_energy.sum(axis=1)

    if len(args.basins) % 2 != 0:
        raise ValueError("--basins must contain an even number of floats.")
    basin_pairs = [
        (float(args.basins[i]), float(args.basins[i + 1]))
        for i in range(0, len(args.basins), 2)
    ]

    out_path = Path(args.out)
    if out_path.suffix.lower() == ".json":
        summary_path = out_path
        output_dir = summary_path.parent
        output_dir.mkdir(parents=True, exist_ok=True)
    else:
        output_dir = out_path
        output_dir.mkdir(parents=True, exist_ok=True)
        summary_path = output_dir / "summary.json"

    residue_indices = []
    for idx, label in enumerate(residue_labels, start=1):
        token = str(label).split("-", 1)[0].strip()
        residue_indices.append(int(token) if token.isdigit() else idx)
    residue_indices = np.asarray(residue_indices, dtype=int)

    def _mutual_information_matrix_pm1(z_block: np.ndarray, epsilon: float = 1e-12) -> np.ndarray:
        n_samples, n_res = z_block.shape
        if n_samples == 0:
            return np.zeros((n_res, n_res), dtype=float)

        state_indicators = np.stack((z_block == -1.0, z_block == 1.0), axis=2).astype(float)
        p_i = state_indicators.mean(axis=0)  # (n_res, 2)
        p_ij = (
            np.einsum("sik,sjl->ijkl", state_indicators, state_indicators, optimize=True)
            / float(n_samples)
        )  # (n_res, n_res, 2, 2)

        denom = p_i[:, None, :, None] * p_i[None, :, None, :]
        mi = np.sum(
            p_ij * np.log((p_ij + epsilon) / (denom + epsilon)),
            axis=(2, 3),
        )
        mi = 0.5 * (mi + mi.T)
        return mi

    def _cross_excitation_probability_matrix_pm1(z_block: np.ndarray) -> np.ndarray:
        n_samples, n_res = z_block.shape
        if n_samples == 0:
            return np.full((n_res, n_res), np.nan, dtype=float)
        excited = (z_block == -1.0).astype(float)
        p_exc_cross = (excited.T @ excited) / float(n_samples)
        p_exc_cross = 0.5 * (p_exc_cross + p_exc_cross.T)
        return p_exc_cross

    per_basin_statistics = []
    basin_stabilities = []
    basin_sample_index_arrays = []
    basin_mi_matrices = []
    basin_cross_excitation_matrices = []
    for basin_index, (low, high) in enumerate(basin_pairs, start=1):
        mask = (E >= low) & (E <= high)
        sample_indices = np.where(mask)[0].astype(int).tolist()
        basin_sample_indices = np.asarray(sample_indices, dtype=int)
        count = len(sample_indices)

        if count > 0:
            mean_energy = float(np.mean(E[mask]))
            std_energy = float(np.std(E[mask]))
            stability = np.mean(Z[mask] == 1.0, axis=0)
        else:
            mean_energy = None
            std_energy = None
            stability = np.full(n_residues, np.nan, dtype=float)

        stability_csv_path = output_dir / f"basin_{basin_index:02d}_stability.csv"
        stability_plot_path = output_dir / f"basin_{basin_index:02d}_stability.png"

        stability_df = pd.DataFrame(
            {
                "residue": residue_labels,
                "ground_state_probability": stability.astype(float),
            }
        )
        stability_df.to_csv(stability_csv_path, index=False)

        mi_matrix = _mutual_information_matrix_pm1(Z[mask], epsilon=1e-12)
        mi_npy_path = output_dir / f"basin_{basin_index:02d}_MI.npy"
        mi_plot_path = output_dir / f"basin_{basin_index:02d}_MI_heatmap.png"
        np.save(mi_npy_path, mi_matrix)

        with plt.rc_context(rc=dict(plt.rcParamsDefault)):
            fig_mi, ax_mi = plt.subplots(figsize=(6, 5))
            im_mi = ax_mi.imshow(
                mi_matrix,
                interpolation="nearest",
                aspect="auto",
                cmap="magma",
            )
            cbar_mi = fig_mi.colorbar(im_mi, ax=ax_mi)
            cbar_mi.set_label("Mutual Information")
            ax_mi.set_xlabel("Residue")
            ax_mi.set_ylabel("Residue")
            ax_mi.set_xticks(np.arange(n_residues))
            ax_mi.set_yticks(np.arange(n_residues))
            ax_mi.set_xticklabels(residue_labels, rotation=90, fontsize=7)
            ax_mi.set_yticklabels(residue_labels, fontsize=7)
            ax_mi.set_title(f"Basin {basin_index} MI [{low}, {high}]")
            fig_mi.tight_layout()
            fig_mi.savefig(mi_plot_path, dpi=600, bbox_inches="tight")
            plt.close(fig_mi)

        set_academic_style()
        fig, ax = plt.subplots(figsize=(10, 4))
        ax.plot(
            residue_indices,
            stability,
            marker="o",
            linewidth=1.4,
            color="#E74C3C",
            markerfacecolor="#E74C3C",
            markeredgecolor="black",
            markersize=4.0,
        )
        ax.set_xlabel("Residue")
        ax.set_xticks(residue_indices)
        ax.set_xticklabels(residue_labels, rotation=90, fontsize=7)
        ax.set_ylabel("Ground-state probability")
        ax.set_ylim(0.0, 1.05)
        ax.set_title(f"Basin {basin_index} [{low}, {high}]")
        fig.tight_layout()
        _save_figure_with_pdf(fig, str(stability_plot_path), dpi=600)
        plt.close(fig)

        per_basin_statistics.append(
            {
                "basin_index": basin_index,
                "low": low,
                "high": high,
                "sample_indices": sample_indices,
                "count": count,
                "mean_energy": mean_energy,
                "std_energy": std_energy,
                "stability_csv": str(stability_csv_path),
                "stability_plot": str(stability_plot_path),
                "stability_plot_pdf": str(stability_plot_path.with_suffix(".pdf")),
                "mi_npy": str(mi_npy_path),
                "mi_heatmap": str(mi_plot_path),
            }
        )
        basin_stabilities.append(stability.astype(float))
        basin_sample_index_arrays.append(basin_sample_indices)
        basin_mi_matrices.append(mi_matrix)
        basin_cross_excitation_matrices.append(_cross_excitation_probability_matrix_pm1(Z[mask]))

    # Combined basin stability plot for direct comparison on a single axis.
    basin_compare_plot_path = output_dir / "basin_stability_comparison.png"
    set_academic_style()
    fig_cmp, ax_cmp = plt.subplots(figsize=(12, 5))
    for basin_index, ((low, high), stability) in enumerate(
        zip(basin_pairs, basin_stabilities), start=1
    ):
        ax_cmp.plot(
            residue_indices,
            stability,
            marker="o",
            linewidth=1.4,
            markersize=3.5,
            markeredgecolor="black",
            label=f"Basin {basin_index} [{low}, {high}]",
        )
    ax_cmp.set_xlabel("Residue")
    ax_cmp.set_xticks(residue_indices)
    ax_cmp.set_xticklabels(residue_labels, rotation=90, fontsize=7)
    ax_cmp.set_ylabel("Ground-state probability")
    ax_cmp.set_ylim(0.0, 1.05)
    ax_cmp.set_title("Basin Ground-State Probability Comparison")
    ax_cmp.legend(loc="upper left", fontsize=8, frameon=True, ncol=1)
    fig_cmp.tight_layout()
    _save_figure_with_pdf(fig_cmp, str(basin_compare_plot_path), dpi=600)
    plt.close(fig_cmp)

    between_basin_delta_stability = []
    between_basin_flip_frequency = []
    between_basin_flip_frequency_cross_residue = []
    between_basin_delta_mi = []
    between_basin_delta_p_excitation = []
    delta_stability_series = []
    flip_frequency_series = []
    n_basins = len(basin_pairs)
    max_flip_samples = int(args.flip_m)
    if max_flip_samples < 0:
        raise ValueError("--flip_m must be non-negative.")
    rng = np.random.default_rng()

    for basin_a_idx in range(n_basins):
        for basin_b_idx in range(n_basins):
            if basin_a_idx == basin_b_idx:
                continue

            a_low, a_high = basin_pairs[basin_a_idx]
            b_low, b_high = basin_pairs[basin_b_idx]
            delta_stability = basin_stabilities[basin_b_idx] - basin_stabilities[basin_a_idx]

            delta_csv_path = output_dir / (
                f"basin_{basin_a_idx + 1:02d}_to_{basin_b_idx + 1:02d}_delta_stability.csv"
            )
            delta_plot_path = output_dir / (
                f"basin_{basin_a_idx + 1:02d}_to_{basin_b_idx + 1:02d}_delta_stability.png"
            )

            delta_df = pd.DataFrame(
                {"residue": residue_labels, "delta_probability": delta_stability.astype(float)}
            )
            delta_df.to_csv(delta_csv_path, index=False)

            set_academic_style()
            fig, ax = plt.subplots(figsize=(10, 4))
            ax.plot(
                residue_indices,
                delta_stability,
                marker="o",
                linewidth=1.4,
                color="#E74C3C",
                markerfacecolor="#E74C3C",
                markeredgecolor="black",
                markersize=4.0,
            )
            ax.axhline(0.0, color="black", linewidth=0.8, linestyle="--")
            ax.set_xlabel("Residue")
            ax.set_xticks(residue_indices)
            ax.set_xticklabels(residue_labels, rotation=90, fontsize=7)
            ax.set_ylabel("Delta probability")
            ax.set_title(
                f"Delta Stability Basin {basin_a_idx + 1} -> {basin_b_idx + 1} "
                f"([{a_low}, {a_high}] -> [{b_low}, {b_high}])"
            )
            fig.tight_layout()
            _save_figure_with_pdf(fig, str(delta_plot_path), dpi=600)
            plt.close(fig)

            between_basin_delta_stability.append(
                {
                    "from_basin_index": basin_a_idx + 1,
                    "to_basin_index": basin_b_idx + 1,
                    "from_low": a_low,
                    "from_high": a_high,
                    "to_low": b_low,
                    "to_high": b_high,
                    "delta_stability_csv": str(delta_csv_path),
                    "delta_stability_plot": str(delta_plot_path),
                    "delta_stability_plot_pdf": str(delta_plot_path.with_suffix(".pdf")),
                }
            )
            delta_stability_series.append(
                (
                    f"Basin {basin_a_idx + 1} -> {basin_b_idx + 1}",
                    delta_stability.astype(float),
                )
            )

            delta_mi = basin_mi_matrices[basin_b_idx] - basin_mi_matrices[basin_a_idx]
            delta_mi_npy_path = output_dir / (
                f"basin_{basin_a_idx + 1:02d}_to_{basin_b_idx + 1:02d}_delta_MI.npy"
            )
            delta_mi_plot_path = output_dir / (
                f"basin_{basin_a_idx + 1:02d}_to_{basin_b_idx + 1:02d}_delta_MI_heatmap.png"
            )
            np.save(delta_mi_npy_path, delta_mi)

            with plt.rc_context(rc=dict(plt.rcParamsDefault)):
                fig_dmi, ax_dmi = plt.subplots(figsize=(6, 5))
                im_dmi = ax_dmi.imshow(
                    delta_mi,
                    interpolation="nearest",
                    aspect="auto",
                    cmap="magma",
                )
                cbar_dmi = fig_dmi.colorbar(im_dmi, ax=ax_dmi)
                cbar_dmi.set_label("Delta MI")
                ax_dmi.set_xlabel("Residue")
                ax_dmi.set_ylabel("Residue")
                ax_dmi.set_xticks(np.arange(n_residues))
                ax_dmi.set_yticks(np.arange(n_residues))
                ax_dmi.set_xticklabels(residue_labels, rotation=90, fontsize=7)
                ax_dmi.set_yticklabels(residue_labels, fontsize=7)
                ax_dmi.set_title(
                    f"Delta MI Basin {basin_a_idx + 1} -> {basin_b_idx + 1}"
                )
                fig_dmi.tight_layout()
                fig_dmi.savefig(delta_mi_plot_path, dpi=600, bbox_inches="tight")
                plt.close(fig_dmi)

            between_basin_delta_mi.append(
                {
                    "from_basin_index": basin_a_idx + 1,
                    "to_basin_index": basin_b_idx + 1,
                    "from_low": a_low,
                    "from_high": a_high,
                    "to_low": b_low,
                    "to_high": b_high,
                    "delta_mi_npy": str(delta_mi_npy_path),
                    "delta_mi_heatmap": str(delta_mi_plot_path),
                }
            )

            delta_p_excitation = (
                basin_cross_excitation_matrices[basin_b_idx]
                - basin_cross_excitation_matrices[basin_a_idx]
            )
            delta_p_excitation_npy_path = output_dir / (
                f"basin_{basin_a_idx + 1:02d}_to_{basin_b_idx + 1:02d}_delta_P_excitation.npy"
            )
            delta_p_excitation_plot_path = output_dir / (
                f"basin_{basin_a_idx + 1:02d}_to_{basin_b_idx + 1:02d}_delta_P_excitation_heatmap.png"
            )
            np.save(delta_p_excitation_npy_path, delta_p_excitation)

            with plt.rc_context(rc=dict(plt.rcParamsDefault)):
                fig_dpe, ax_dpe = plt.subplots(figsize=(6, 5))
                im_dpe = ax_dpe.imshow(
                    delta_p_excitation,
                    interpolation="nearest",
                    aspect="auto",
                    cmap="magma",
                )
                cbar_dpe = fig_dpe.colorbar(im_dpe, ax=ax_dpe)
                cbar_dpe.set_label("Delta P(excitation)")
                ax_dpe.set_xlabel("Residue")
                ax_dpe.set_ylabel("Residue")
                ax_dpe.set_xticks(np.arange(n_residues))
                ax_dpe.set_yticks(np.arange(n_residues))
                ax_dpe.set_xticklabels(residue_labels, rotation=90, fontsize=7)
                ax_dpe.set_yticklabels(residue_labels, fontsize=7)
                ax_dpe.set_title(
                    f"Delta P(excitation) Basin {basin_a_idx + 1} -> {basin_b_idx + 1}"
                )
                fig_dpe.tight_layout()
                fig_dpe.savefig(delta_p_excitation_plot_path, dpi=600, bbox_inches="tight")
                plt.close(fig_dpe)

            between_basin_delta_p_excitation.append(
                {
                    "from_basin_index": basin_a_idx + 1,
                    "to_basin_index": basin_b_idx + 1,
                    "from_low": a_low,
                    "from_high": a_high,
                    "to_low": b_low,
                    "to_high": b_high,
                    "delta_p_excitation_npy": str(delta_p_excitation_npy_path),
                    "delta_p_excitation_heatmap": str(delta_p_excitation_plot_path),
                }
            )

            indices_a = basin_sample_index_arrays[basin_a_idx]
            indices_b = basin_sample_index_arrays[basin_b_idx]
            sample_size = min(max_flip_samples, indices_a.size, indices_b.size)

            if sample_size > 0:
                sampled_a = rng.choice(indices_a, size=sample_size, replace=False)
                sampled_b = rng.choice(indices_b, size=sample_size, replace=False)
                flip_events = (Z[sampled_a] != Z[sampled_b]).astype(float)
                flip_frequency = np.mean(flip_events, axis=0).astype(float)
                flip_frequency_cross_residue = (flip_events.T @ flip_events) / float(sample_size)
                flip_frequency_cross_residue = 0.5 * (
                    flip_frequency_cross_residue + flip_frequency_cross_residue.T
                )
            else:
                flip_frequency = np.full(n_residues, np.nan, dtype=float)
                flip_frequency_cross_residue = np.full(
                    (n_residues, n_residues), np.nan, dtype=float
                )

            flip_csv_path = output_dir / (
                f"basin_{basin_a_idx + 1:02d}_to_{basin_b_idx + 1:02d}_flip_frequency.csv"
            )
            flip_plot_path = output_dir / (
                f"basin_{basin_a_idx + 1:02d}_to_{basin_b_idx + 1:02d}_flip_frequency.png"
            )
            flip_cross_npy_path = output_dir / (
                f"basin_{basin_a_idx + 1:02d}_to_{basin_b_idx + 1:02d}_flip_frequency_cross_residue.npy"
            )
            flip_cross_plot_path = output_dir / (
                f"basin_{basin_a_idx + 1:02d}_to_{basin_b_idx + 1:02d}_flip_frequency_cross_residue_heatmap.png"
            )
            np.save(flip_cross_npy_path, flip_frequency_cross_residue)

            flip_df = pd.DataFrame(
                {
                    "residue": residue_labels,
                    "flip_frequency": flip_frequency.astype(float),
                }
            )
            flip_df.to_csv(flip_csv_path, index=False)

            set_academic_style()
            fig, ax = plt.subplots(figsize=(10, 4))
            ax.plot(
                residue_indices,
                flip_frequency,
                marker="o",
                linewidth=1.4,
                color="#16A085",
                markerfacecolor="#16A085",
                markeredgecolor="black",
                markersize=4.0,
            )
            ax.set_xlabel("Residue")
            ax.set_xticks(residue_indices)
            ax.set_xticklabels(residue_labels, rotation=90, fontsize=7)
            ax.set_ylabel("Flip frequency")
            ax.set_ylim(0.0, 1.0)
            ax.set_title(
                f"Flip Frequency Basin {basin_a_idx + 1} -> {basin_b_idx + 1} "
                f"(m={sample_size})"
            )
            fig.tight_layout()
            _save_figure_with_pdf(fig, str(flip_plot_path), dpi=600)
            plt.close(fig)

            with plt.rc_context(rc=dict(plt.rcParamsDefault)):
                fig_fcr, ax_fcr = plt.subplots(figsize=(6, 5))
                im_fcr = ax_fcr.imshow(
                    flip_frequency_cross_residue,
                    interpolation="nearest",
                    aspect="auto",
                    cmap="magma",
                )
                cbar_fcr = fig_fcr.colorbar(im_fcr, ax=ax_fcr)
                cbar_fcr.set_label("Flip frequency probability")
                ax_fcr.set_xlabel("Residue")
                ax_fcr.set_ylabel("Residue")
                ax_fcr.set_xticks(np.arange(n_residues))
                ax_fcr.set_yticks(np.arange(n_residues))
                ax_fcr.set_xticklabels(residue_labels, rotation=90, fontsize=7)
                ax_fcr.set_yticklabels(residue_labels, fontsize=7)
                ax_fcr.set_title(
                    f"Flip Frequency (ri,rj) Basin {basin_a_idx + 1} -> {basin_b_idx + 1}"
                )
                fig_fcr.tight_layout()
                fig_fcr.savefig(flip_cross_plot_path, dpi=600, bbox_inches="tight")
                plt.close(fig_fcr)

            between_basin_flip_frequency.append(
                {
                    "from_basin_index": basin_a_idx + 1,
                    "to_basin_index": basin_b_idx + 1,
                    "from_low": a_low,
                    "from_high": a_high,
                    "to_low": b_low,
                    "to_high": b_high,
                    "sample_size": int(sample_size),
                    "flip_frequency_csv": str(flip_csv_path),
                    "flip_frequency_plot": str(flip_plot_path),
                    "flip_frequency_plot_pdf": str(flip_plot_path.with_suffix(".pdf")),
                }
            )
            between_basin_flip_frequency_cross_residue.append(
                {
                    "from_basin_index": basin_a_idx + 1,
                    "to_basin_index": basin_b_idx + 1,
                    "from_low": a_low,
                    "from_high": a_high,
                    "to_low": b_low,
                    "to_high": b_high,
                    "sample_size": int(sample_size),
                    "flip_frequency_cross_residue_npy": str(flip_cross_npy_path),
                    "flip_frequency_cross_residue_heatmap": str(flip_cross_plot_path),
                }
            )
            flip_frequency_series.append(
                (
                    f"Basin {basin_a_idx + 1} -> {basin_b_idx + 1}",
                    flip_frequency.astype(float),
                )
            )

    delta_stability_compare_plot_path = output_dir / "delta_stability_comparison.png"
    set_academic_style()
    fig_delta_cmp, ax_delta_cmp = plt.subplots(figsize=(12, 5))
    for series_label, series_values in delta_stability_series:
        ax_delta_cmp.plot(
            residue_indices,
            series_values,
            marker="o",
            linewidth=1.3,
            markersize=3.2,
            markeredgecolor="black",
            label=series_label,
        )
    ax_delta_cmp.axhline(0.0, color="black", linewidth=0.8, linestyle="--")
    ax_delta_cmp.set_xlabel("Residue")
    ax_delta_cmp.set_xticks(residue_indices)
    ax_delta_cmp.set_xticklabels(residue_labels, rotation=90, fontsize=7)
    ax_delta_cmp.set_ylabel("Delta probability")
    ax_delta_cmp.set_title("Delta Stability Comparison (All Basin Pairs)")
    ax_delta_cmp.legend(loc="upper left", fontsize=7, frameon=True, ncol=2)
    fig_delta_cmp.tight_layout()
    _save_figure_with_pdf(fig_delta_cmp, str(delta_stability_compare_plot_path), dpi=600)
    plt.close(fig_delta_cmp)

    flip_frequency_compare_plot_path = output_dir / "flip_frequency_comparison.png"
    set_academic_style()
    fig_flip_cmp, ax_flip_cmp = plt.subplots(figsize=(12, 5))
    for series_label, series_values in flip_frequency_series:
        ax_flip_cmp.plot(
            residue_indices,
            series_values,
            marker="o",
            linewidth=1.3,
            markersize=3.2,
            markeredgecolor="black",
            label=series_label,
        )
    ax_flip_cmp.set_xlabel("Residue")
    ax_flip_cmp.set_xticks(residue_indices)
    ax_flip_cmp.set_xticklabels(residue_labels, rotation=90, fontsize=7)
    ax_flip_cmp.set_ylabel("Flip frequency")
    ax_flip_cmp.set_ylim(0.0, 1.0)
    ax_flip_cmp.set_title("Flip Frequency Comparison (All Basin Pairs)")
    ax_flip_cmp.legend(loc="upper left", fontsize=7, frameon=True, ncol=2)
    fig_flip_cmp.tight_layout()
    _save_figure_with_pdf(fig_flip_cmp, str(flip_frequency_compare_plot_path), dpi=600)
    plt.close(fig_flip_cmp)

    summary = {
        "total_shots": int(Z.shape[0]),
        "number_of_residues": int(n_residues),
        "per_basin_statistics": per_basin_statistics,
        "basin_stability_comparison_plot": str(basin_compare_plot_path),
        "basin_stability_comparison_plot_pdf": str(basin_compare_plot_path.with_suffix(".pdf")),
        "delta_stability_comparison_plot": str(delta_stability_compare_plot_path),
        "delta_stability_comparison_plot_pdf": str(delta_stability_compare_plot_path.with_suffix(".pdf")),
        "flip_frequency_comparison_plot": str(flip_frequency_compare_plot_path),
        "flip_frequency_comparison_plot_pdf": str(flip_frequency_compare_plot_path.with_suffix(".pdf")),
        "between_basin_delta_stability": between_basin_delta_stability,
        "between_basin_delta_mi": between_basin_delta_mi,
        "between_basin_delta_p_excitation": between_basin_delta_p_excitation,
        "between_basin_flip_frequency": between_basin_flip_frequency,
        "between_basin_flip_frequency_cross_residue": between_basin_flip_frequency_cross_residue,
    }

    with summary_path.open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)

    print(f"[saved] summary -> {summary_path}")
    return summary
