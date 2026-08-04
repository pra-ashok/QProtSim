from __future__ import annotations

from typing import Dict, List, Tuple, Sequence, Mapping, Any
import numpy as np
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array


DEFAULT_KAPPA_WEIGHTS: Dict[str, float] = {
    "BB": 1.00,
    "disulfide": 1.00,
    "HB": 0.50,
    "salt": 0.40,
    "vdW": 0.05,
    # "wat": 0.40,
}


AA_MASS_DA: Dict[str, float] = {
    "GLY": 75.07,
    "ALA": 89.09,
    "SER": 105.09,
    "PRO": 115.13,
    "VAL": 117.15,
    "THR": 119.12,
    "CYS": 121.16,
    "CYX": 121.16,
    "LEU": 131.17,
    "ILE": 131.17,
    "ASN": 132.12,
    "ASP": 133.10,
    "GLN": 146.15,
    "LYS": 146.19,
    "GLU": 147.13,
    "MET": 149.21,
    "HIS": 155.16,
    "HID": 155.16,
    "HIE": 155.16,
    "HSD": 155.16,
    "HSE": 155.16,
    "HIP": 156.17,
    "HSP": 156.17,
    "PHE": 165.19,
    "ARG": 174.20,
    "TYR": 181.19,
    "TRP": 204.23,
}

WATER = 18.02
# PEPTIDE_MASS = sum(AA_MASS_DA[res] for res in sequence) - (len(sequence) - 1) * WATER

POSITIVE_SIDECHAIN_ATOMS: Dict[str, set[str]] = {
    "LYS": {"NZ"},
    "ARG": {"NE", "NH1", "NH2"},
    "HIP": {"ND1", "NE2"},
    "HSP": {"ND1", "NE2"},
}

NEGATIVE_SIDECHAIN_ATOMS: Dict[str, set[str]] = {
    "ASP": {"OD1", "OD2"},
    "GLU": {"OE1", "OE2"},
}


def _residue_chain_label(residue: Any) -> str:
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


def _selected_residues(
    u: mda.Universe,
    protein_only: bool = True,
    chain: str | None = "A",
) -> list[Any]:
    residues = []
    protein_resnames = set(AA_MASS_DA)

    for res in u.residues:
        chain_label = _residue_chain_label(res)
        if chain is not None and chain_label != chain:
            continue
        if protein_only and str(res.resname).strip().upper() not in protein_resnames:
            continue
        residues.append(res)

    return residues


def _safe_element(atom: Any) -> str:
    """
    Robustly infer atom element.

    MDAnalysis topologies sometimes have missing atom.element fields,
    so we fall back to the first alphabetical character in atom.name.
    """
    element = getattr(atom, "element", "")

    if element is not None:
        element = str(element).strip().upper()
        if element:
            return element

    name = str(getattr(atom, "name", "")).strip().upper()

    for ch in name:
        if ch.isalpha():
            return ch

    return ""


def _heavy_atom_positions(res: Any) -> np.ndarray:
    """
    Return positions of non-hydrogen atoms in a residue.
    """
    atoms = [
        atom for atom in res.atoms
        if _safe_element(atom) != "H"
    ]

    if len(atoms) == 0:
        return np.empty((0, 3), dtype=np.float64)

    return np.asarray(
        [atom.position for atom in atoms],
        dtype=np.float64,
    )


def _polar_atom_positions(res: Any) -> np.ndarray:
    """
    Return positions of polar heavy atoms used for approximate H-bond detection.

    This is a heavy-atom approximation:
        N/O/S ... N/O/S distance <= H-bond cutoff.

    It does not require explicit hydrogens.
    """
    atoms = [
        atom for atom in res.atoms
        if _safe_element(atom) in {"N", "O", "S"}
    ]

    if len(atoms) == 0:
        return np.empty((0, 3), dtype=np.float64)

    return np.asarray(
        [atom.position for atom in atoms],
        dtype=np.float64,
    )


def _named_atom_positions(res: Any, atom_names: set[str]) -> np.ndarray:
    """
    Return positions of atoms whose names are in atom_names.
    """
    atoms = [
        atom for atom in res.atoms
        if str(atom.name).strip().upper() in atom_names
    ]

    if len(atoms) == 0:
        return np.empty((0, 3), dtype=np.float64)

    return np.asarray(
        [atom.position for atom in atoms],
        dtype=np.float64,
    )


def _residue_mass(res: Any) -> float:
    """
    Return residue mass.

    Uses topology masses if available. If unavailable or zero,
    falls back to standard residue molecular masses. If residue
    name is unknown, falls back to number of atoms.
    Adjusted for peptide (- WATER)
    """
    try:
        mass = float(np.sum(res.atoms.masses))
        if np.isfinite(mass) and mass > 0.0:
            return mass
    except Exception:
        pass

    resname = str(res.resname).strip().upper()

    if resname in AA_MASS_DA:
        return float(AA_MASS_DA[resname]) - WATER

    return float(max(len(res.atoms), 1))


def _mass_mobility_factor(m_i: float, m_j: float) -> float:
    """
    Bounded mass compatibility factor:

        Lambda_ij = 2 sqrt(M_i * M_j) / (M_i + M_j)

    This is 1 for identical masses and decreases when masses differ.
    """
    if m_i <= 0.0 or m_j <= 0.0:
        return 1.0

    return float(
        (2.0 * np.sqrt(m_i * m_j)) / (m_i + m_j)
    )


def _kappa_backbone(
    res_i: Any,
    res_j: Any,
    i: int,
    j: int,
    box: np.ndarray | None = None,
    peptide_cutoff: float = 1.8,
) -> float:
    """
    Backbone peptide-bond connectivity term.

    Returns 1 if residues are adjacent and appear peptide-bonded,
    otherwise 0.

    This is a discrete structural term:
        kappa_BB = 1 for peptide neighbours
        kappa_BB = 0 otherwise
    """
    if abs(j - i) != 1:
        return 0.0

    # Determine the earlier and later residue in sequence index.
    res_left, res_right = (res_i, res_j) if i < j else (res_j, res_i)

    c_pos = _named_atom_positions(res_left, {"C"})
    n_pos = _named_atom_positions(res_right, {"N"})

    if c_pos.size > 0 and n_pos.size > 0:
        d = distance_array(c_pos, n_pos, box=box)

        if float(np.min(d)) <= peptide_cutoff:
            return 1.0

        return 0.0

    # Fallback for coarse-grained or incomplete structures.
    # If C/N atoms are absent, sequence adjacency is used.
    return 1.0


def _kappa_disulfide(
    res_i: Any,
    res_j: Any,
    box: np.ndarray | None = None,
    disulfide_cutoff: float = 2.3,
) -> float:
    """
    Disulfide covalent crosslink term.

    Returns 1 when both residues are cysteine-like residues and their SG atoms
    are within disulfide_cutoff, otherwise 0. The default cutoff covers the
    typical protein disulfide S-S bond length (~2.03 A) with modest coordinate
    tolerance.
    """
    resname_i = str(res_i.resname).strip().upper()
    resname_j = str(res_j.resname).strip().upper()

    if resname_i not in {"CYS", "CYX"} or resname_j not in {"CYS", "CYX"}:
        return 0.0

    sg_i = _named_atom_positions(res_i, {"SG"})
    sg_j = _named_atom_positions(res_j, {"SG"})

    if sg_i.size == 0 or sg_j.size == 0:
        return 0.0

    d = distance_array(sg_i, sg_j, box=box)

    if float(np.min(d)) <= disulfide_cutoff:
        return 1.0

    return 0.0


def _kappa_hbond(
    res_i: Any,
    res_j: Any,
    box: np.ndarray | None = None,
    hb_cutoff: float = 3.5,
    hb_power: float = 2.0,
    hb_norm: float = 1.0,
) -> float:
    """
    Approximate hydrogen-bond coupling term.

    Heavy-atom approximation:

        polar atom pair = N/O/S ... N/O/S
        contact if r_ab <= hb_cutoff

    Raw score:

        raw_HB = sum_{polar contacts} (hb_cutoff / r_ab)^hb_power / hb_norm

    Bounded score:

        kappa_HB = 1 - exp(-raw_HB)

    This keeps kappa_HB in [0, 1].
    """
    xyz_i = _polar_atom_positions(res_i)
    xyz_j = _polar_atom_positions(res_j)

    if xyz_i.size == 0 or xyz_j.size == 0:
        return 0.0

    d = distance_array(xyz_i, xyz_j, box=box)
    contacts = d[d <= hb_cutoff]

    if contacts.size == 0:
        return 0.0

    contacts = np.maximum(contacts, 1.0e-6)

    raw = np.sum((hb_cutoff / contacts) ** hb_power) / hb_norm

    return float(1.0 - np.exp(-raw))


def _charged_atom_positions(res: Any) -> tuple[np.ndarray, str | None]:
    """
    Return charged side-chain atom positions and charge type.

    Returns
    -------
    xyz : np.ndarray
        Charged atom positions.

    charge_type : {"positive", "negative", None}
        Charge category.
    """
    resname = str(res.resname).strip().upper()

    if resname in POSITIVE_SIDECHAIN_ATOMS:
        xyz = _named_atom_positions(
            res,
            POSITIVE_SIDECHAIN_ATOMS[resname],
        )
        return xyz, "positive"

    if resname in NEGATIVE_SIDECHAIN_ATOMS:
        xyz = _named_atom_positions(
            res,
            NEGATIVE_SIDECHAIN_ATOMS[resname],
        )
        return xyz, "negative"

    return np.empty((0, 3), dtype=np.float64), None


def _kappa_salt_bridge(
    res_i: Any,
    res_j: Any,
    box: np.ndarray | None = None,
    salt_cutoff: float = 4.0,
    salt_power: float = 2.0,
    salt_norm: float = 1.0,
) -> float:
    """
    Salt-bridge ionic charged-contact coupling term.

    This term captures the electrostatic side-chain contact component of a
    salt bridge. The hydrogen-bond component is represented separately by
    kappa_HB in S_ij, so a geometrically compatible salt bridge can contribute
    through both w_salt * kappa_salt and w_HB * kappa_HB.

    A charged contact is counted between positive and negative side-chain
    atoms if r_ab <= salt_cutoff. Positive groups are Lys NZ, Arg guanidinium
    nitrogens, and explicitly protonated histidine names HIP/HSP. Negative
    groups are Asp/Glu carboxylate oxygens.

    Raw score:

        raw_salt = sum_{charged contacts} (salt_cutoff / r_ab)^salt_power / salt_norm

    Bounded score:

        kappa_salt = 1 - exp(-raw_salt)
    """
    xyz_i, charge_i = _charged_atom_positions(res_i)
    xyz_j, charge_j = _charged_atom_positions(res_j)

    if charge_i is None or charge_j is None:
        return 0.0

    if charge_i == charge_j:
        return 0.0

    if xyz_i.size == 0 or xyz_j.size == 0:
        return 0.0

    d = distance_array(xyz_i, xyz_j, box=box)
    contacts = d[d <= salt_cutoff]

    if contacts.size == 0:
        return 0.0

    contacts = np.maximum(contacts, 1.0e-6)

    raw = np.sum((salt_cutoff / contacts) ** salt_power) / salt_norm

    return float(1.0 - np.exp(-raw))


def _kappa_vdw(
    res_i: Any,
    res_j: Any,
    box: np.ndarray | None = None,
    vdw_cutoff: float = 4.0,
    vdw_r0: float = 3.5,
    vdw_power: float = 2.0,
    vdw_norm: float = 8.0,
) -> float:
    """
    van der Waals / packing coupling term.

    Uses heavy-atom contacts.

    Raw score:

        raw_vdW =
            sum_{heavy atom contacts}
            (vdw_r0 / r_ab)^vdw_power / vdw_norm

    Bounded score:

        kappa_vdW = 1 - exp(-raw_vdW)

    The vdw_norm parameter controls how quickly the score saturates.
    Larger vdw_norm means more atom contacts are needed to approach 1.
    """
    xyz_i = _heavy_atom_positions(res_i)
    xyz_j = _heavy_atom_positions(res_j)

    if xyz_i.size == 0 or xyz_j.size == 0:
        return 0.0

    d = distance_array(xyz_i, xyz_j, box=box)
    contacts = d[d <= vdw_cutoff]

    if contacts.size == 0:
        return 0.0

    contacts = np.maximum(contacts, 1.0e-6)

    raw = np.sum((vdw_r0 / contacts) ** vdw_power) / vdw_norm

    return float(1.0 - np.exp(-raw))


def _precompute_residue_water_sets(
    residues: Sequence[Any],
    water_xyz: np.ndarray,
    box: np.ndarray | None = None,
    water_cutoff: float = 3.5,
) -> List[set[int]]:
    """
    For each residue, find nearby water oxygen indices.

    A water oxygen is considered associated with a residue if it lies
    within water_cutoff of any residue atom.
    """
    residue_water_sets: List[set[int]] = []

    if water_xyz.size == 0:
        return [set() for _ in residues]

    for res in residues:
        xyz = np.asarray(res.atoms.positions, dtype=np.float64)

        if xyz.size == 0:
            residue_water_sets.append(set())
            continue

        d = distance_array(xyz, water_xyz, box=box)
        close_water_indices = np.where(np.any(d <= water_cutoff, axis=0))[0]

        residue_water_sets.append(set(close_water_indices.tolist()))

    return residue_water_sets


def _kappa_water_shared(
    water_set_i: set[int],
    water_set_j: set[int],
    wat_norm: float = 2.0,
) -> float:
    """
    Shared hydration shell / solvent channel coupling.

    Let n_shared be the number of water oxygens close to both residues.

        kappa_wat = 1 - exp(-n_shared / wat_norm)

    This is 0 if the residues do not share hydration waters.
    """
    if not water_set_i or not water_set_j:
        return 0.0

    n_shared = len(water_set_i.intersection(water_set_j))

    if n_shared == 0:
        return 0.0

    return float(1.0 - np.exp(-float(n_shared) / wat_norm))


def _frequency_matching_factor(
    omega_i: float,
    omega_j: float,
    sigma_omega: float,
) -> float:
    """
    Spectral overlap / frequency matching term:

        Phi_ij = exp[-(omega_i - omega_j)^2 / (2 sigma_omega^2)]
    """
    if sigma_omega <= 0.0:
        return 1.0

    delta = omega_i - omega_j

    return float(
        np.exp(
            -((delta * delta) / (2.0 * sigma_omega * sigma_omega))
        )
    )


def build_Kij_pairs(
    u: mda.Universe,
    protein_only: bool = True,
    chain: str | None = "A",

    # Important:
    # If you want kappa_BB to contribute, keep exclude_bonded=False.
    exclude_bonded: bool = False,
    min_seq_separation: int = 1,
    use_pbc: bool = False,

    weights: Mapping[str, float] | None = None,

    omega: Sequence[float] | Mapping[int, float] | None = None,
    sigma_omega: float | None = None,
    estimate_omega_from_local_stiffness: bool = True,

    # water_selection: str = (
    #     "(resname HOH WAT SOL TIP3 TIP3P TIP4P SPC OPC) "
    #     "and (name O OH2 OW)"
    # ),

    normalize: bool = True,
    return_terms: bool = True,
) -> Tuple[List[Tuple[int, int]], Dict[Tuple[int, int], float]] | Tuple[
    List[Tuple[int, int]],
    Dict[Tuple[int, int], float],
    Dict[Tuple[int, int], Dict[str, float]],
]:
    """
    Build residue-residue vibrational coupling factors:

        K_ij^vib = S_ij * Phi_ij * Lambda_ij^mass

    with

        S_ij =
            w_BB    * kappa_BB
          + w_disulfide * kappa_disulfide
          + w_HB    * kappa_HB
          + w_salt  * kappa_salt
          + w_vdW   * kappa_vdW
        #   + w_wat   * kappa_wat

        Phi_ij =
            exp[-(omega_i - omega_j)^2 / (2 sigma_omega^2)]

        Lambda_ij^mass =
            2 sqrt(M_i M_j) / (M_i + M_j)

    Notes
    -----
    The old geometric contact term C_ij is folded into kappa_vdW.
    This avoids double-counting residue contacts.

    If omega is not provided and estimate_omega_from_local_stiffness=True,
    the function estimates relative residue frequencies using:

        omega_i = sqrt(k_local_i / M_i)

    where

        k_local_i = sum_j S_ij

    This gives relative, not absolute, vibrational frequencies.

    Parameters
    ----------
    u : MDAnalysis.Universe
        MDAnalysis universe.

    protein_only : bool
        Whether to select protein residues only.

    chain : str | None
        Chain identifier. If None, all selected residues are used.

    exclude_bonded : bool
        If True, sequence-neighbour residue pairs are skipped.
        Keep this False if you want the backbone term kappa_BB.

    min_seq_separation : int
        If exclude_bonded=True, skip pairs with |i-j| <= min_seq_separation.

    use_pbc : bool
        Whether to use periodic boundary conditions.

    weights : mapping or None
        Weights for the coupling components.
        Expected keys are:
            "BB", "disulfide", "HB", "salt", "vdW", "wat"

    omega : sequence, mapping, or None
        Optional residue frequencies.
        If sequence, omega[i] is used.
        If mapping, omega[i] is used.
        If None, frequencies are estimated from local stiffness if
        estimate_omega_from_local_stiffness=True.

    sigma_omega : float or None
        Width of the frequency-matching Gaussian.
        If None, an automatic value is estimated from the omega distribution.

    estimate_omega_from_local_stiffness : bool
        Whether to estimate omega_i = sqrt(k_local_i / M_i) if omega is None.

    # water_selection : str
    #     MDAnalysis selection string for water oxygen atoms.

    normalize : bool
        If True, normalize final K_ij values by max(K_ij).

    return_terms : bool
        If True, also return the individual kappa terms and factors.

    Returns
    -------
    pairs : list[tuple[int, int]]
        Residue-pair indices.

    Kijs : dict[tuple[int, int], float]
        Dictionary mapping (i, j) to K_ij^vib.

    terms : dict, optional
        Returned only when return_terms=True.
        Contains the decomposed terms for each pair.
    """

    weights_used = dict(DEFAULT_KAPPA_WEIGHTS)

    if weights is not None:
        for key, value in weights.items():
            key_norm = str(key)
            if key_norm == "SS":
                key_norm = "disulfide"
            weights_used[key_norm] = float(value)

    residues = _selected_residues(
        u,
        protein_only=protein_only,
        chain=chain,
    )

    n_res = len(residues)

    box = u.dimensions if use_pbc else None

    masses = np.asarray(
        [_residue_mass(res) for res in residues],
        dtype=np.float64,
    )

    # Precompute water oxygen positions.
    # try:
    #     water_atoms = u.select_atoms(water_selection)
    #     water_xyz = np.asarray(water_atoms.positions, dtype=np.float64)
    # except Exception:
    #     water_xyz = np.empty((0, 3), dtype=np.float64)

    # residue_water_sets = _precompute_residue_water_sets(
    #     residues=residues,
    #     water_xyz=water_xyz,
    #     box=box,
    # )

    # ------------------------------------------------------------
    # First pass:
    # compute kappa terms and S_ij for each residue pair.
    # ------------------------------------------------------------
    pairs: List[Tuple[int, int]] = []
    Sijs: Dict[Tuple[int, int], float] = {}
    terms: Dict[Tuple[int, int], Dict[str, float]] = {}

    for i in range(n_res):
        for j in range(i + 1, n_res):

            if exclude_bonded and (j - i) <= min_seq_separation:
                continue

            kappa_BB = _kappa_backbone(
                residues[i],
                residues[j],
                i,
                j,
                box=box,
            )

            kappa_disulfide = _kappa_disulfide(
                residues[i],
                residues[j],
                box=box,
            )

            kappa_HB = _kappa_hbond(
                residues[i],
                residues[j],
                box=box,
            )

            kappa_salt = _kappa_salt_bridge(
                residues[i],
                residues[j],
                box=box,
            )

            kappa_vdW = _kappa_vdw(
                residues[i],
                residues[j],
                box=box,
            )

            # kappa_wat = _kappa_water_shared(
            #     residue_water_sets[i],
            #     residue_water_sets[j],
            # )

            Sij = (
                weights_used.get("BB", 0.0) * kappa_BB
                + weights_used.get("disulfide", 0.0) * kappa_disulfide
                + weights_used.get("HB", 0.0) * kappa_HB
                + weights_used.get("salt", 0.0) * kappa_salt
                + weights_used.get("vdW", 0.0) * kappa_vdW
                # + weights_used.get("wat", 0.0) * kappa_wat
            )

            Sij = float(max(0.0, Sij))

            if Sij <= 0.0:
                continue

            pair = (i, j)

            pairs.append(pair)
            Sijs[pair] = Sij

            terms[pair] = {
                "kappa_BB": float(kappa_BB),
                "kappa_disulfide": float(kappa_disulfide),
                "kappa_HB": float(kappa_HB),
                "kappa_salt": float(kappa_salt),
                "kappa_vdW": float(kappa_vdW),
                # "kappa_wat": float(kappa_wat),
                "Sij": float(Sij),
            }

    # ------------------------------------------------------------
    # Frequencies omega_i.
    # ------------------------------------------------------------
    omega_values: np.ndarray | None = None

    if omega is not None:
        omega_values = np.zeros(n_res, dtype=np.float64)

        for i in range(n_res):
            if isinstance(omega, Mapping):
                omega_values[i] = float(omega[i])
            else:
                omega_values[i] = float(omega[i])

    elif estimate_omega_from_local_stiffness:
        local_stiffness = np.zeros(n_res, dtype=np.float64)

        for (i, j), Sij in Sijs.items():
            local_stiffness[i] += Sij
            local_stiffness[j] += Sij

        omega_values = np.sqrt(
            np.maximum(local_stiffness, 1.0e-12)
            / np.maximum(masses, 1.0e-12)
        )

    # If we have omega values but no sigma_omega, estimate sigma automatically.
    if omega_values is not None and sigma_omega is None:
        finite_omega = omega_values[np.isfinite(omega_values)]

        if finite_omega.size > 1 and float(np.std(finite_omega)) > 0.0:
            sigma_omega_used = float(np.std(finite_omega))
        else:
            sigma_omega_used = 1.0
    elif sigma_omega is not None:
        sigma_omega_used = float(sigma_omega)
    else:
        sigma_omega_used = 1.0

    # ------------------------------------------------------------
    # Second pass:
    # K_ij = S_ij * Phi_ij * Lambda_mass
    # ------------------------------------------------------------
    Kijs: Dict[Tuple[int, int], float] = {}

    for pair in pairs:
        i, j = pair

        Sij = Sijs[pair]

        if omega_values is not None:
            Phiij = _frequency_matching_factor(
                omega_i=float(omega_values[i]),
                omega_j=float(omega_values[j]),
                sigma_omega=sigma_omega_used,
            )
        else:
            Phiij = 1.0

        Lambda_mass = _mass_mobility_factor(
            float(masses[i]),
            float(masses[j]),
        )

        Kij = float(Sij * Phiij * Lambda_mass)

        if Kij <= 0.0:
            continue

        Kijs[pair] = Kij

        terms[pair]["Phiij"] = float(Phiij)
        terms[pair]["Lambda_mass"] = float(Lambda_mass)
        terms[pair]["Kij_raw"] = float(Kij)

    # ------------------------------------------------------------
    # Optional normalization.
    # ------------------------------------------------------------
    if normalize and len(Kijs) > 0:
        Kmax = max(Kijs.values())

        if Kmax > 0.0:
            for pair in Kijs:
                Kijs[pair] = float(Kijs[pair] / Kmax)
                terms[pair]["Kij"] = float(Kijs[pair])
        else:
            for pair in Kijs:
                terms[pair]["Kij"] = float(Kijs[pair])
    else:
        for pair in Kijs:
            terms[pair]["Kij"] = float(Kijs[pair])

    if return_terms:
        return pairs, Kijs, terms

    return pairs, Kijs
