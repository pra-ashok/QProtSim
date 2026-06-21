"""
PDB Structure Visualization Utilities

This module provides functions to generate or download structure images for reports:
1. Download official images from RCSB PDB
2. Render custom 3D structures using PyMOL
3. Generate interactive HTML with Mol* (Molstar) viewer
4. Create simple structure representations with matplotlib + MDAnalysis
"""

from pathlib import Path
from typing import Optional
import requests
import warnings

try:
    import pymol
    from pymol import cmd
    HAVE_PYMOL = True
except ImportError:
    HAVE_PYMOL = False

try:
    from PIL import Image as PILImage
    HAVE_PIL = True
except ImportError:
    HAVE_PIL = False


def _export_pdf_from_png(png_path: Path) -> Optional[Path]:
    if not HAVE_PIL:
        warnings.warn("PIL/Pillow not installed. Cannot export PDF from PNG.")
        return None
    try:
        img = PILImage.open(png_path)
        if img.mode == "RGBA":
            bg = PILImage.new("RGB", img.size, (255, 255, 255))
            bg.paste(img, mask=img.split()[3])
            img = bg
        pdf_path = png_path.with_suffix(".pdf")
        img.save(pdf_path, "PDF")
        print(f"[INFO] Exported PDF figure: {pdf_path}")
        return pdf_path
    except Exception as e:
        warnings.warn(f"PDF export failed for {png_path}: {e}")
        return None


def download_rcsb_image(
    pdb_code: str,
    output_path: Path,
    image_type: str = "cartoon",
    size: str = "800",
) -> Optional[Path]:
    """
    Download structure image from RCSB PDB website.
    
    RCSB provides pre-rendered images via their image API:
    - https://cdn.rcsb.org/images/structures/{pdb_code}_assembly-1.jpeg
    - https://cdn.rcsb.org/images/structures/{pdb_code}_model-1.jpeg
    
    Parameters
    ----------
    pdb_code : str
        PDB code (e.g., "1L2Y")
    output_path : Path
        Where to save the downloaded image
    image_type : str
        Type of image: "cartoon", "assembly", "model"
    size : str
        Image size (not directly controllable via API)
        
    Returns
    -------
    Path or None
        Path to downloaded image, or None if download failed
    """
    pdb_code = pdb_code.upper()
    
    # Try different RCSB image URLs
    urls = [
        f"https://cdn.rcsb.org/images/structures/{pdb_code.lower()}_assembly-1.jpeg",
        f"https://cdn.rcsb.org/images/structures/{pdb_code.lower()}_model-1.jpeg",
        f"https://cdn.rcsb.org/images/structures/1{pdb_code.lower()}_models.jpeg",
        f"https://files.rcsb.org/view/{pdb_code}.png",
    ]
    
    for url in urls:
        try:
            response = requests.get(url, timeout=10)
            if response.status_code == 200:
                output_path.parent.mkdir(parents=True, exist_ok=True)
                
                # Save with appropriate extension
                if url.endswith('.png'):
                    save_path = output_path.with_suffix('.png')
                else:
                    save_path = output_path.with_suffix('.jpeg')
                
                with save_path.open('wb') as f:
                    f.write(response.content)
                
                print(f"[INFO] Downloaded structure image from RCSB: {save_path}")
                if save_path.suffix.lower() == ".png":
                    _export_pdf_from_png(save_path)
                return save_path
                
        except Exception as e:
            warnings.warn(f"Failed to download from {url}: {e}")
            continue
    
    print(f"[WARNING] Could not download structure image from RCSB for {pdb_code}")
    return None


def render_structure_pymol(
    pdb_path: Path,
    output_path: Path,
    width: int = 800,
    height: int = 600,
    style: str = "cartoon",
    color_scheme: str = "spectrum",
    ray_trace: bool = True,
) -> Optional[Path]:
    """
    Render PDB structure using PyMOL for high-quality static images.
    
    PyMOL is a professional molecular visualization tool that can run
    in headless mode without GUI, perfect for generating publication-quality
    images programmatically.
    
    Parameters
    ----------
    pdb_path : Path
        Path to PDB file
    output_path : Path
        Where to save rendered image (PNG format)
    width : int
        Image width in pixels
    height : int
        Image height in pixels
    style : str
        Representation style: "cartoon", "sticks", "spheres", "ribbon", "surface"
    color_scheme : str
        Color scheme: "spectrum", "chain", "secondary", "element"
    ray_trace : bool
        Use ray tracing for high-quality rendering (slower but better)
        
    Returns
    -------
    Path or None
        Path to rendered image, or None if rendering failed
    """
    if not HAVE_PYMOL:
        warnings.warn("PyMOL not installed. Cannot render structure. Install with: pip install pymol-open-source")
        return None
    
    try:
        # Initialize PyMOL in headless mode (no GUI)
        pymol.finish_launching(['pymol', '-c'])
        
        # Load structure
        cmd.load(str(pdb_path), "protein")
        
        # Set background
        cmd.bg_color("white")
        
        # Apply representation style
        cmd.hide("everything", "protein")
        if style == "cartoon":
            cmd.show("cartoon", "protein")
        elif style == "sticks":
            cmd.show("sticks", "protein")
        elif style == "spheres":
            cmd.show("spheres", "protein")
        elif style == "ribbon":
            cmd.show("ribbon", "protein")
        elif style == "surface":
            cmd.show("surface", "protein")
        else:
            cmd.show("cartoon", "protein")  # default
        
        # Apply color scheme
        if color_scheme == "spectrum":
            cmd.spectrum("count", "rainbow", "protein")
        elif color_scheme == "chain":
            cmd.util.cbc("protein")  # color by chain
        elif color_scheme == "secondary":
            cmd.util.cbss("protein")  # color by secondary structure
        elif color_scheme == "element":
            cmd.util.cbag("protein")  # color by element
        else:
            cmd.spectrum("count", "rainbow", "protein")
        
        # Center and zoom
        cmd.center("protein")
        cmd.zoom("protein", buffer=2)
        
        # Set image size
        cmd.viewport(width, height)
        
        # Render and save
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_png = output_path.with_suffix('.png')
        
        if ray_trace:
            cmd.ray(width, height)
        
        cmd.png(str(output_png), width, height, dpi=150)
        
        # Clean up
        cmd.delete("protein")
        cmd.reinitialize()
        
        print(f"[INFO] Rendered structure with PyMOL: {output_png}")
        _export_pdf_from_png(output_png)
        return output_png
        
    except Exception as e:
        warnings.warn(f"PyMOL rendering failed: {e}")
        try:
            cmd.reinitialize()
        except:
            pass
        return None


def generate_molstar_html(
    pdb_path: Path,
    output_path: Path,
    title: str = "Protein Structure",
) -> Optional[Path]:
    """
    Generate standalone HTML file with embedded Mol* (Molstar) viewer.
    
    Molstar is the successor to LiteMol and PDB Component Library viewer,
    developed by PDBe and RCSB. It's a modern, WebGL-based molecular viewer
    that works in any browser without plugins.
    
    Parameters
    ----------
    pdb_path : Path
        Path to PDB file
    output_path : Path
        Where to save HTML file
    title : str
        Title for the HTML page
        
    Returns
    -------
    Path or None
        Path to generated HTML file, or None if generation failed
    """
    try:
        # Read PDB content
        with open(pdb_path, 'r') as f:
            pdb_content = f.read()
        
        # Escape for JavaScript
        pdb_content_escaped = pdb_content.replace('\\', '\\\\').replace('`', '\\`').replace('$', '\\$')
        
        # Generate HTML with embedded Molstar viewer
        html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{title}</title>
    <style>
        body {{
            margin: 0;
            padding: 0;
            width: 100%;
            height: 100vh;
            overflow: hidden;
            font-family: Arial, sans-serif;
        }}
        #viewer {{
            position: absolute;
            top: 0;
            left: 0;
            right: 0;
            bottom: 0;
        }}
        .title {{
            position: absolute;
            top: 10px;
            left: 10px;
            z-index: 1000;
            background: rgba(255, 255, 255, 0.9);
            padding: 10px;
            border-radius: 4px;
            font-size: 18px;
            font-weight: bold;
        }}
    </style>
    <link rel="stylesheet" type="text/css" href="https://www.ebi.ac.uk/pdbe/pdb-component-library/css/pdbe-molstar-3.1.0.css">
</head>
<body>
    <div class="title">{title}</div>
    <div id="viewer"></div>
    
    <script src="https://www.ebi.ac.uk/pdbe/pdb-component-library/js/pdbe-molstar-plugin-3.1.0.js"></script>
    <script>
        const pdbData = `{pdb_content_escaped}`;
        
        const viewerInstance = new PDBeMolstarPlugin();
        const options = {{
            customData: {{
                url: '',
                format: 'pdb',
                binary: false
            }},
            alphafoldView: false,
            bgColor: {{ r: 255, g: 255, b: 255 }},
            hideControls: false,
            selectInteraction: true,
            lighting: 'matte',
            validationAnnotation: false,
            domainAnnotation: false,
            lowPrecisionCoords: true,
            expanded: false,
            hideCanvasControls: ['selection', 'animation', 'controlToggle']
        }};
        
        const viewerContainer = document.getElementById('viewer');
        viewerInstance.render(viewerContainer, options);
        
        // Load PDB data
        viewerInstance.events.loadComplete.subscribe(() => {{
            console.log('Structure loaded successfully');
        }});
        
        // Load the PDB data directly
        setTimeout(() => {{
            viewerInstance.loadStructureFromData(pdbData, 'pdb');
        }}, 100);
    </script>
</body>
</html>"""
        
        # Save HTML file
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_html = output_path.with_suffix('.html')
        
        with open(output_html, 'w') as f:
            f.write(html_content)
        
        print(f"[INFO] Generated Molstar HTML viewer: {output_html}")
        print(f"[INFO] Open {output_html} in a web browser to view the structure")
        return output_html
        
    except Exception as e:
        warnings.warn(f"Molstar HTML generation failed: {e}")
        return None


def render_structure_py3dmol(
    pdb_path: Path,
    output_path: Path,
    width: int = 800,
    height: int = 600,
    style: str = "cartoon",
    color_scheme: str = "spectrum",
) -> Optional[Path]:
    """
    Render PDB structure using py3Dmol (DEPRECATED - use PyMOL or Molstar instead).
    
    Note: py3Dmol is designed for Jupyter notebooks. For static image
    generation, use render_structure_pymol() instead.
    For interactive viewing, use generate_molstar_html().
    
    Parameters
    ----------
    pdb_path : Path
        Path to PDB file
    output_path : Path
        Where to save rendered image
    width : int
        Image width in pixels
    height : int
        Image height in pixels
    style : str
        Representation style: "cartoon", "stick", "sphere", "line"
    color_scheme : str
        Color scheme: "spectrum", "chainHetatm", "ssPyMOL"
        
    Returns
    -------
    Path or None
        Always returns None - this function is deprecated
    """
    warnings.warn(
        "py3Dmol requires Jupyter environment and is not suitable for static image export. "
        "Use render_structure_pymol() for static images or generate_molstar_html() for interactive viewing.",
        DeprecationWarning
    )
    return None


def render_structure_matplotlib(
    pdb_path: Path,
    output_path: Path,
    chain: Optional[str] = "A",
    width: int = 800,
    height: int = 600,
    dpi: int = 100,
) -> Optional[Path]:
    """
    Create a simple 2D projection of the protein structure using matplotlib.
    
    This creates a basic CA trace visualization projected onto 2D.
    Not as visually appealing as 3D renders, but requires no external tools.
    
    Parameters
    ----------
    pdb_path : Path
        Path to PDB file
    output_path : Path
        Where to save rendered image
    chain : Optional[str]
        Protein chain to visualize. If None, include all chains.
    width : int
        Image width in pixels
    height : int
        Image height in pixels
    dpi : int
        Resolution (dots per inch)
        
    Returns
    -------
    Path or None
        Path to rendered image, or None if rendering failed
    """
    try:
        import MDAnalysis as mda
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D
        import numpy as np
        
        # Load structure
        u = mda.Universe(str(pdb_path))
        
        if chain is None:
            ca_atoms = u.select_atoms("protein and name CA")
        else:
            # Prefer chainID; fall back for topologies where chain IDs are absent.
            ca_atoms = u.select_atoms(f"protein and name CA and chainID {chain}")
            if len(ca_atoms) == 0:
                ca_atoms = u.select_atoms(f"protein and name CA and chainid {chain}")
            if len(ca_atoms) == 0:
                ca_atoms = u.select_atoms(f"protein and name CA and segid {chain}")
        
        if len(ca_atoms) == 0:
            warnings.warn(f"No CA atoms found for chain {chain}")
            return None
        
        positions = ca_atoms.positions
        
        # Create 3D plot
        fig = plt.figure(figsize=(width/dpi, height/dpi), dpi=dpi)
        ax = fig.add_subplot(111, projection='3d')
        
        # Plot CA trace
        ax.plot(positions[:, 0], positions[:, 1], positions[:, 2], 
                'o-', linewidth=2, markersize=6, 
                color='#2C3E50', alpha=0.8)
        
        # Color by residue index (gradient)
        n_res = len(ca_atoms)
        colors = plt.cm.viridis(np.linspace(0, 1, n_res))
        ax.scatter(positions[:, 0], positions[:, 1], positions[:, 2],
                   c=colors, s=50, alpha=0.9, edgecolors='black', linewidth=0.5)
        
        # Labels
        ax.set_xlabel('X (Å)', fontsize=10, fontweight='bold')
        ax.set_ylabel('Y (Å)', fontsize=10, fontweight='bold')
        ax.set_zlabel('Z (Å)', fontsize=10, fontweight='bold')
        chain_title = chain if chain is not None else "All"
        ax.set_title(f'Protein Structure (Chain {chain_title})',
                     fontsize=12, fontweight='bold', pad=20)
        
        # Set viewing angle
        ax.view_init(elev=20, azim=45)
        
        # Make axes equal
        max_range = np.ptp(positions, axis=0).max() / 2.0
        mid_x = (positions[:, 0].max() + positions[:, 0].min()) * 0.5
        mid_y = (positions[:, 1].max() + positions[:, 1].min()) * 0.5
        mid_z = (positions[:, 2].max() + positions[:, 2].min()) * 0.5
        ax.set_xlim(mid_x - max_range, mid_x + max_range)
        ax.set_ylim(mid_y - max_range, mid_y + max_range)
        ax.set_zlim(mid_z - max_range, mid_z + max_range)
        
        plt.tight_layout()
        
        output_path.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
        if output_path.suffix.lower() == ".png":
            pdf_path = output_path.with_suffix(".pdf")
            fig.savefig(pdf_path, format="pdf", bbox_inches="tight")
            print(f"[INFO] Exported PDF figure: {pdf_path}")
        plt.close(fig)
        
        print(f"[INFO] Rendered structure using matplotlib: {output_path}")
        return output_path
        
    except Exception as e:
        warnings.warn(f"matplotlib rendering failed: {e}")
        return None


def get_structure_image(
    pdb_path: Path,
    output_dir: Path,
    pdb_code: Optional[str] = None,
    chain: Optional[str] = "A",
    prefer_download: bool = True,
    output_path: Optional[Path] = None,
) -> Optional[Path]:
    """
    Get a structure image by downloading from RCSB or rendering locally.
    
    Strategy:
    1. If pdb_code is provided and prefer_download=True, try RCSB download
    2. Try rendering with PyMOL (best quality)
    3. Fallback to matplotlib rendering
    4. Return None if all methods fail
    
    Parameters
    ----------
    pdb_path : Path
        Path to local PDB file
    output_dir : Path
        Directory to save the image
    pdb_code : str, optional
        PDB code for downloading from RCSB
    chain : Optional[str]
        Protein chain identifier. If None, include all chains.
    prefer_download : bool
        If True, try downloading from RCSB first
        
    Returns
    -------
    Path or None
        Path to structure image, or None if generation failed
    """
    if output_path is None:
        output_path = output_dir / "structure_visualization.png"
    
    # Try RCSB download first if requested
    if pdb_code and prefer_download:
        print("[INFO] Attempting to download structure image from RCSB...")
        downloaded = download_rcsb_image(pdb_code, output_path)
        if downloaded:
            return downloaded
    
    # Try PyMOL rendering (best quality for static images)
    if HAVE_PYMOL:
        print("[INFO] Rendering structure with PyMOL...")
        rendered = render_structure_pymol(
            pdb_path=pdb_path,
            output_path=output_path,
            width=800,
            height=600,
            style="cartoon",
            color_scheme="spectrum",
            ray_trace=True,
        )
        if rendered:
            return rendered
    else:
        print("[INFO] PyMOL not available (install with: pip install pymol-open-source)")
    
    # Fallback to matplotlib rendering
    print("[INFO] Rendering structure with matplotlib (basic rendering)...")
    rendered = render_structure_matplotlib(
        pdb_path=pdb_path,
        output_path=output_path,
        chain=chain,
        width=800,
        height=600,
        dpi=100,
    )
    
    if rendered:
        return rendered
    
    # All methods failed
    warnings.warn("Could not generate structure image (tried RCSB download, PyMOL, and matplotlib rendering)")
    return None


def optimize_image_for_report(
    image_path: Path,
    max_width: int = 800,
    max_height: int = 600,
    quality: int = 85,
) -> Optional[Path]:
    """
    Optimize image size and quality for PDF reports.
    
    Parameters
    ----------
    image_path : Path
        Path to image file
    max_width : int
        Maximum width in pixels
    max_height : int
        Maximum height in pixels
    quality : int
        JPEG quality (1-100)
        
    Returns
    -------
    Path or None
        Path to optimized image
    """
    if not HAVE_PIL:
        warnings.warn("PIL/Pillow not installed. Cannot optimize image.")
        return image_path
    
    try:
        img = PILImage.open(image_path)
        
        # Resize if needed
        if img.width > max_width or img.height > max_height:
            img.thumbnail((max_width, max_height), PILImage.Resampling.LANCZOS)
        
        # Convert RGBA to RGB if needed
        if img.mode == 'RGBA':
            bg = PILImage.new('RGB', img.size, (255, 255, 255))
            bg.paste(img, mask=img.split()[3])
            img = bg
        
        # Save optimized version
        optimized_path = image_path.with_stem(image_path.stem + "_optimized")
        img.save(optimized_path, 'PNG', optimize=True, quality=quality)
        
        print(f"[INFO] Optimized image: {optimized_path}")
        if optimized_path.suffix.lower() == ".png":
            _export_pdf_from_png(optimized_path)
        return optimized_path
        
    except Exception as e:
        warnings.warn(f"Image optimization failed: {e}")
        return image_path
