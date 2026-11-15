from __future__ import annotations

import hashlib
import os
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import warnings

import numpy as np
import MDAnalysis as mda
import freesasa

try:
    import requests
    HAVE_REQUESTS = True
except ImportError:
    HAVE_REQUESTS = False
    warnings.warn("requests module not found. PDB download and checksum validation will be disabled.")

try:
    import cupy as cp
    HAVE_CUPY = True
except Exception:
    HAVE_CUPY = False


def compute_file_checksum(file_path: str, algorithm: str = "sha256") -> str:
    """
    Compute checksum of a file.
    
    Parameters
    ----------
    file_path : str
        Path to the file
    algorithm : str
        Hash algorithm to use (md5, sha1, sha256)
        
    Returns
    -------
    str
        Hexadecimal checksum string
    """
    hash_obj = hashlib.new(algorithm)
    with open(file_path, 'rb') as f:
        for chunk in iter(lambda: f.read(8192), b""):
            hash_obj.update(chunk)
    return hash_obj.hexdigest()


def get_pdb_checksum(pdb_code: str) -> Optional[str]:
    """
    Retrieve ETag (checksum) from RCSB PDB for a given PDB code.
    
    The ETag header serves as a checksum/version identifier for the file.
    
    Parameters
    ----------
    pdb_code : str
        4-character PDB code
        
    Returns
    -------
    Optional[str]
        ETag string if available, None otherwise
    """
    if not HAVE_REQUESTS:
        return None
        
    pdb_code = pdb_code.upper()
    
    # Get ETag from RCSB file download service
    try:
        url = f"https://files.rcsb.org/download/{pdb_code}.pdb"
        response = requests.head(url, timeout=10)
        
        if response.status_code == 200:
            # RCSB provides ETag which can be used for validation
            etag = response.headers.get('ETag', '').strip('"')
            if etag:
                return etag
        elif response.status_code == 404:
            # 404 is expected for non-existent PDB codes - don't warn, just return None
            return None
        else:
            warnings.warn(f"Unexpected status code {response.status_code} for {pdb_code}")
    except requests.RequestException as e:
        warnings.warn(f"Could not retrieve checksum for {pdb_code}: {e}")
    except Exception as e:
        warnings.warn(f"Unexpected error retrieving checksum for {pdb_code}: {e}")
    
    return None


def verify_pdb_file(file_path: str, pdb_code: str, expected_checksum: Optional[str] = None) -> bool:
    """
    Verify integrity of downloaded PDB file.
    
    Parameters
    ----------
    file_path : str
        Path to the PDB file
    pdb_code : str
        PDB code
    expected_checksum : Optional[str]
        Expected checksum (if None, will skip checksum verification)
        
    Returns
    -------
    bool
        True if file is valid, False otherwise
    """
    if not Path(file_path).exists():
        return False
        
    # Basic file size check (PDB files should be at least a few KB)
    file_size = Path(file_path).stat().st_size
    if file_size < 1000:  # Less than 1KB is likely corrupt
        warnings.warn(f"PDB file {file_path} is suspiciously small ({file_size} bytes)")
        return False
    
    # If checksum verification is requested
    if expected_checksum:
        try:
            # For ETag comparison, we can use MD5
            file_checksum = compute_file_checksum(file_path, "md5")
            if file_checksum.lower() != expected_checksum.lower():
                warnings.warn(f"Checksum mismatch for {pdb_code}")
                return False
        except Exception as e:
            warnings.warn(f"Checksum verification failed: {e}")
            return False
    
    return True


def fetch_pdb(
    pdb_code: str,
    output_dir: Optional[str] = None,
    force_download: bool = False,
    verify_checksum: bool = True
) -> Path:
    """
    Fetch PDB file from RCSB PDB database using the official file download API.
    
    This function:
    1. Checks if file already exists locally (unless force_download=True)
    2. Downloads the file from https://files.rcsb.org/download/{pdb_code}.pdb
    3. Verifies checksum (ETag) if available (when verify_checksum=True)
    
    Parameters
    ----------
    pdb_code : str
        4-character PDB code (e.g., '1L2Y')
    output_dir : Optional[str]
        Directory to save PDB files. If None, uses current working directory.
    force_download : bool
        If True, download even if file exists locally
    verify_checksum : bool
        If True, verify file integrity using ETag checksum
        
    Returns
    -------
    Path
        Path to the PDB file
        
    Raises
    ------
    RuntimeError
        If requests module is not available
    FileNotFoundError
        If PDB code not found (404) or download fails
    ValueError
        If checksum verification fails
    """
    if not HAVE_REQUESTS:
        raise RuntimeError(
            "requests module is required for PDB downloads. "
            "Install with: pip install requests"
        )
    
    pdb_code = pdb_code.upper()
    
    # Setup output directory
    if output_dir is None:
        output_dir = os.getcwd()
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    final_path = output_path / f"{pdb_code}.pdb"
    
    # Check if file already exists
    if final_path.exists() and not force_download:
        print(f"[INFO] PDB file {pdb_code} already exists at {final_path}")
        
        # Verify existing file
        if verify_checksum:
            expected_checksum = get_pdb_checksum(pdb_code)
            if expected_checksum and not verify_pdb_file(str(final_path), pdb_code, expected_checksum):
                print(f"[WARNING] Existing file failed verification. Re-downloading...")
                force_download = True
            else:
                return final_path
        else:
            return final_path
    
    # Download the file
    print(f"[INFO] Downloading PDB {pdb_code} from RCSB...")
    
    # Get checksum before download if verification is enabled
    expected_checksum = None
    if verify_checksum:
        expected_checksum = get_pdb_checksum(pdb_code)
    
    try:
        # Use RCSB file download service API
        url = f"https://files.rcsb.org/download/{pdb_code}.pdb"
        response = requests.get(url, timeout=30)
        
        if response.status_code == 404:
            raise FileNotFoundError(f"PDB code {pdb_code} not found on RCSB server")
        elif response.status_code != 200:
            raise RuntimeError(f"Download failed with status code {response.status_code}")
        
        # Write file
        final_path.write_bytes(response.content)
        
        print(f"[SUCCESS] Downloaded {pdb_code} to {final_path}")
        
        # Verify the downloaded file
        if verify_checksum and expected_checksum:
            if not verify_pdb_file(str(final_path), pdb_code, expected_checksum):
                raise ValueError(f"Downloaded file failed integrity check for {pdb_code}")
            print(f"[SUCCESS] File integrity verified for {pdb_code}")
        
        return final_path
        
    except requests.RequestException as e:
        print(f"[ERROR] Network error while fetching PDB {pdb_code}: {e}")
        raise
    except Exception as e:
        print(f"[ERROR] Failed to fetch PDB {pdb_code}: {e}")
        raise


def batch_fetch_pdb(
    pdb_codes: List[str],
    output_dir: Optional[str] = None,
    force_download: bool = False,
    verify_checksum: bool = True,
    delay: float = 0.5
) -> Dict[str, Path]:
    """
    Fetch multiple PDB files synchronously with status tracking.
    
    This function downloads PDB files one at a time and waits for each
    download to complete before proceeding to the next. This ensures
    all files are available before returning.
    
    Parameters
    ----------
    pdb_codes : List[str]
        List of PDB codes to download
    output_dir : Optional[str]
        Directory to save PDB files
    force_download : bool
        If True, re-download even if files exist
    verify_checksum : bool
        If True, verify file integrity
    delay : float
        Delay in seconds between downloads (to be polite to RCSB servers)
        
    Returns
    -------
    Dict[str, Path]
        Dictionary mapping PDB codes to their file paths
        
    Raises
    ------
    RuntimeError
        If any download fails or if pdb_codes list is empty
    """
    if not pdb_codes:
        raise RuntimeError("PDB codes list cannot be empty")
    
    results = {}
    failed = []
    
    print(f"\n[INFO] Starting batch download of {len(pdb_codes)} PDB files...")
    print("=" * 60)
    
    for i, pdb_code in enumerate(pdb_codes, 1):
        print(f"\n[{i}/{len(pdb_codes)}] Processing {pdb_code}...")
        
        try:
            path = fetch_pdb(
                pdb_code=pdb_code,
                output_dir=output_dir,
                force_download=force_download,
                verify_checksum=verify_checksum
            )
            results[pdb_code.upper()] = path
            
            # Be polite to RCSB servers - add delay between downloads
            if i < len(pdb_codes) and delay > 0:
                time.sleep(delay)
                
        except Exception as e:
            print(f"[ERROR] Failed to fetch {pdb_code}: {e}")
            failed.append(pdb_code)
    
    print("\n" + "=" * 60)
    print(f"[SUMMARY] Successfully downloaded: {len(results)}/{len(pdb_codes)}")
    
    if failed:
        print(f"[WARNING] Failed downloads: {', '.join(failed)}")
        raise RuntimeError(f"Failed to download {len(failed)} PDB files: {failed}")
    
    print("[SUCCESS] All PDB files downloaded and verified!")
    return results

pdb_list = ["1L2Y", "1AKG", "1B03", "1B45", "1A13", "1C98", "1D9L", "1DPK", "1DU1", "1DVW", "1E0Q", "1EGT", "1EQX", "1LE1", "1A7F"]

# A minimal set of standard amino acids (3-letter codes)
AA3 = {
    "ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY",
    "HIS","ILE","LEU","LYS","MET","PHE","PRO","SER",
    "THR","TRP","TYR","VAL"
}

# Tien et al. 2013 theoretical MaxASA (A^2) for Gly–X–Gly tripeptides
# Source: Tien MZ et al., PLoS ONE 2013; also summarized on Wikipedia.
MAXASA_TIEN13_THEOR = {
    "ALA": 129.0, "ARG": 274.0, "ASN": 195.0, "ASP": 193.0, "CYS": 167.0,
    "GLN": 225.0, "GLU": 223.0, "GLY": 104.0, "HIS": 224.0, "ILE": 197.0,
    "LEU": 201.0, "LYS": 236.0, "MET": 224.0, "PHE": 240.0, "PRO": 159.0,
    "SER": 155.0, "THR": 172.0, "TRP": 285.0, "TYR": 263.0, "VAL": 174.0
}

def per_residue_sasa(
    pdb_path: str,
    chain: Optional[str] = "A",
    protein_only: bool = True
) -> Dict[Tuple[str, int, str], float]:
    """
    Compute normalized per-residue solvent exposure using FreeSASA.

    Parameters
    ----------
    pdb_path : str
        Path to the PDB file.
    chain : str or None
        Restrict results to this chain ID (e.g. 'A'). If None, include all chains.
    protein_only : bool
        If True, only include standard protein residues.

    Returns
    -------
    res_sasa : dict
        Mapping (chain_label, resid, resname) to dimensionless SASA divided
        by the residue's Tien et al. theoretical maximum, capped at 0.999.
    """
    structure = freesasa.Structure(pdb_path)
    result = freesasa.calc(structure)
    if chain is not None and not chain.isalnum():
        chain = ''
    res_sasa: Dict[Tuple[str, int, str], float] = {}
    residue_areas = result.residueAreas()
    for chain_label, residues in residue_areas.items():
        if chain is not None and chain_label.strip() != chain:
            continue
        for resid, rarea in residues.items():
            resid_int = int(resid) if str(resid).isdigit() else resid
            resname = rarea.residueType  # ResidueArea has residueType (3-letter name)
            if protein_only and resname not in AA3:
                continue
            # / MAXASA_TIEN13_THEOR[resname]
            res_sasa[(chain_label, resid_int, resname)] = min(0.999, float(rarea.total) / MAXASA_TIEN13_THEOR[resname])
    print("res SASA",res_sasa)
    return res_sasa

def map_sasa_to_mda_resorder(
    u: mda.Universe,
    sasa_dict: Dict[Tuple[str, int, str], float],
    protein_only: bool = True,
    chain: Optional[str] = "A",
) -> Tuple[np.ndarray, List[str]]:
    """
    Align SASA dict to MDA residue order. If FreeSASA lacks names, we match on (chain,resid,*).
    """
    if chain is not None and not chain.isalnum():
        chain = 'SYSTEM'

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

    vals, names = [], []
    # build chain-agnostic fallback (resid+name)
    by_cr = {(c, r): area for (c, r, _n), area in sasa_dict.items()}
    for res in u.residues:
        chain_label = _residue_chain_label(res)
        if chain is not None and chain_label.strip() != chain:
            continue
        if chain == "SYSTEM":
            chain_label = ' '
        rname = res.resname.strip()
        if protein_only and rname not in AA3:
            continue
        key_full = (chain_label, int(res.resid), rname)
        key_cr   = (chain_label, int(res.resid))
        val = sasa_dict.get(key_full, by_cr.get(key_cr, 0.0))
        vals.append(val)
        names.append(rname)
    return np.asarray(vals, float), names


# Useless
def per_residue_ke(
    u: mda.Universe,
    T: float = 300.0,
    protein_only: bool = True,
    seed: int | None = None,
) -> np.ndarray:
    """
    Per-residue kinetic energy (KE).

    If velocities present in trajectory:
        KE = 0.5 * sum_i m_i ||v_i||^2
    Else:
        Generate random velocities from Maxwell–Boltzmann distribution
        at temperature T to get realistic KE distribution.

    Parameters
    ----------
    u : MDAnalysis.Universe
        Universe with residues
    T : float
        Temperature in Kelvin (default: 300)
    seed : int or None
        Random seed for reproducibility
    """
    rng = np.random.default_rng(seed)
    kB = 1.380649e-23  # J/K
    amu_to_kg = 1.66053906660e-27  # kg

    has_vel = hasattr(u.trajectory.ts, "velocities") and (u.trajectory.ts.velocities is not None)

    KE = []
    for res in u.residues:
        if protein_only and res.resname.strip() not in AA3:
            continue
        masses = res.atoms.masses  # amu
        if has_vel and (res.atoms.velocities is not None):
            v = res.atoms.velocities  # (n_atoms, 3), assumed in Å/ps
            v_m_s = v * 100.0  # convert Å/ps → m/s
            v2 = np.sum(v_m_s * v_m_s, axis=1)
            m_kg = masses * amu_to_kg
            ke_r = 0.5 * float(np.sum(m_kg * v2))
        else:
            # Maxwell–Boltzmann velocities per atom
            m_kg = masses * amu_to_kg
            sigma_v = np.sqrt(kB * T / m_kg)  # std of velocity distribution (m/s)
            v = rng.normal(0.0, sigma_v[:, None], size=(len(masses), 3))
            v2 = np.sum(v * v, axis=1)
            ke_r = 0.5 * float(np.sum(m_kg * v2))
        KE.append(ke_r)

    return np.asarray(KE, float)
