#!/usr/bin/env python3
"""
CIF Clustering Analysis Tool

This script performs feature extraction and hierarchical clustering on CIF files.
It calculates d/d_min (normalized interatomic distances) for each atomic site,
converts them to histograms, and clusters CIFs based on their structural similarity.

Usage:
    python cif_clustering.py --root /path/to/data --cif-dir CIFs-cleaned-AB2

Example CSV format (if using --csv):
    ,cif,formula,stype,As,Bs
    0,528761.cif,GdPt2,"MgCu2,cF24,227",Gd,Pt
    1,528748.cif,Pt3U,"Mg3Cd,hP8,194",U,Pt
    ...
"""

import argparse
import os
import sys
import tempfile
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy.cluster import hierarchy
from scipy.cluster.hierarchy import dendrogram
from sklearn.cluster import AgglomerativeClustering
from cifkit import Cif
import multiprocessing as mp
from datetime import datetime


# Embedded paper.mplstyle
PAPER_MPLSTYLE = """
text.usetex: False
# text.latex.preamble: \\usepackage{amsmath}\\usepackage{amssymb}
font.family: serif
# font.serif: Computer Modern
savefig.bbox: tight
savefig.format: pdf

lines.linewidth: .5

axes.linewidth:     0.5
axes.labelsize:     10
axes.labelpad:      3.0

grid.linewidth: 0.2

xtick.top:           True
xtick.major.width:   .3
xtick.labelsize:     8
xtick.direction:     in
ytick.right:         True
ytick.labelsize:     8
ytick.direction:     in

legend.loc:           upper right
legend.frameon:       False
legend.fontsize:      6

figure.dpi:         300
"""


def get_first_n_neighbors(points_wd, n=100):
    """
    Calculate normalized distances to the first n neighbors.
    
    Parameters:
        points_wd: Connection data from cifkit
        n: Number of neighbors to include
    
    Returns:
        Normalized distance array (d / d_min)
    """
    point = np.array(points_wd[0][2])
    sites = [np.array(p[-1]) for p in points_wd]
    points_wd = np.vstack(sites)
    d = np.linalg.norm(points_wd - point, axis=1)
    ind = np.argsort(d)[:n]
    d = d[ind]
    d /= d.min()
    return d


def get_features(cif_path, sites=None, n_neighbors=100):
    """
    Extract d/d_min features from a CIF file.
    
    Parameters:
        cif_path: Path to CIF file
        sites: Optional list of element symbols to include (e.g., ["A", "B"])
        n_neighbors: Number of neighbors for ddmin calculation
    
    Returns:
        Dictionary mapping site names to ddmin arrays
        Returns empty dict on error
    """
    try:
        cif = Cif(cif_path)
        cif.compute_connections()
        site_symbol_map = {site: vals['element'] for site, vals in cif.atom_site_info.items()}
        
        desc = {}
        for site, points_wd in cif.connections.items():
            if sites is not None and site_symbol_map[site] not in sites:
                continue
            ddmin_hist = get_first_n_neighbors(points_wd=points_wd, n=n_neighbors)
            desc[site] = ddmin_hist
        return desc
    except Exception as e:
        print(f"Error processing {cif_path}: {e}")
        return {}


def write_features(args):
    """
    Extract features and save to .npz file.
    
    Parameters:
        args: Tuple of (cif_path, sites, n_neighbors, npdescs_dir, 
              plot_individual, bin_min, bin_max, n_bins, plots_dir, mpl_style_path, root)
    
    Returns:
        Tuple of (cif_id, success, error_message)
    """
    (cif_path, sites, n_neighbors, npdescs_dir, plot_individual, 
     bin_min, bin_max, n_bins, plots_dir, mpl_style_path, root) = args
    
    try:
        # Extract features
        desc = get_features(cif_path, sites=sites, n_neighbors=n_neighbors)
        
        if not desc:
            return (os.path.basename(cif_path), False, "No features extracted")
        
        # Get CIF ID
        cif_name = os.path.basename(cif_path)
        if cif_name.lower().endswith('.cif'):
            cif_id = cif_name[:-4]
        else:
            cif_id = cif_name
        
        # Save features
        npz_path = os.path.join(npdescs_dir, f"{cif_id}.npz")
        np.savez(npz_path, **desc)
        
        # Optional: create individual plot
        if plot_individual:
            plot_individual_cif(desc, cif_id, cif_path, 
                              bin_min, bin_max, n_bins, 
                              plots_dir, mpl_style_path, root)
        
        return (cif_id, True, None)
    
    except Exception as e:
        return (os.path.basename(cif_path), False, str(e))


def plot_individual_cif(desc, cif_id, cif_path, bin_min, bin_max, n_bins, 
                       plots_dir, mpl_style_path, root):
    """
    Create histogram plot for individual CIF.
    
    Parameters:
        desc: Feature dictionary from get_features
        cif_id: CIF identifier
        cif_path: Path to original CIF file
        bin_min, bin_max, n_bins: Histogram parameters
        plots_dir: Output directory for plots
        mpl_style_path: Path to MPL style file
        root: Root directory
    """
    try:
        # Load structure type
        cif = Cif(cif_path)
        structure_type = cif.structure
        
        # Create bins
        bins = np.linspace(bin_min, bin_max, n_bins)
        
        # Flatten all values
        _desc = []
        for k, v in desc.items():
            _desc.extend([float(_v) for _v in v])
        
        # Create histogram
        counts, bin_edges = np.histogram(_desc, bins=bins)
        
        # Plot
        plt.close('all')
        if mpl_style_path and os.path.exists(mpl_style_path):
            plt.style.use(mpl_style_path)
        
        plt.bar(
            bins[:-1],
            counts,
            width=np.diff(bins),
            align='edge',
            edgecolor='black'
        )
        plt.title(structure_type)
        plt.xlim(bin_min, bin_max)
        plt.xlabel(r"d/d\textsubscript{min} value (a. u.)")
        plt.ylabel('Count')
        
        # Save plot
        individual_dir = os.path.join(plots_dir, "individual")
        os.makedirs(individual_dir, exist_ok=True)
        
        safe_title = structure_type.replace(',', '-').replace('/', '_')
        plot_path = os.path.join(individual_dir, f"{safe_title}-{cif_id}.png")
        plt.savefig(plot_path, dpi=300)
        plt.close('all')
        
    except Exception as e:
        print(f"Error plotting {cif_id}: {e}")


def preprocess_desc(npz_path, max_sites=None, bin_min=1.0, bin_max=2.5, n_bins=100):
    """
    Convert .npz features to histogram.
    
    Parameters:
        npz_path: Path to .npz file
        max_sites: Maximum number of sites (None for no limit)
        bin_min, bin_max, n_bins: Histogram parameters
    
    Returns:
        Histogram counts array or None if site count exceeds max_sites
    """
    npdesc = np.load(npz_path)
    
    if max_sites is not None and len(list(npdesc.keys())) > max_sites:
        return None
    
    # Flatten all values
    _desc = []
    for k, v in npdesc.items():
        _desc.extend([float(_val) for _val in v])
    
    # Create histogram
    bins = np.linspace(bin_min, bin_max, n_bins)
    counts, bin_edges = np.histogram(_desc, bins=bins)
    
    return counts


def get_df_data(cif_folder):
    """
    Scan CIF folder and return DataFrame with filenames.
    
    Parameters:
        cif_folder: Path to folder containing CIF files
    
    Returns:
        DataFrame with 'cif' column containing filenames (.cif extension)
    """
    cif_files = []
    for f in os.listdir(cif_folder):
        if f.lower().endswith('.cif'):
            cif_files.append(f)
    
    df = pd.DataFrame({'cif': sorted(cif_files)})
    return df


def validate_cifs(df, cif_folder):
    """
    Validate that all CIFs in DataFrame exist in folder.
    
    Parameters:
        df: DataFrame with 'cif' column
        cif_folder: Path to CIF folder
    
    Returns:
        List of missing CIF identifiers
    """
    folder_files = set(f for f in os.listdir(cif_folder) if f.lower().endswith('.cif'))
    missing = []
    
    for cif in df['cif'].tolist():
        # Normalize: strip .cif/.CIF suffix
        if cif.lower().endswith('.cif'):
            normalized = cif
        else:
            normalized = cif + '.cif'
        
        if normalized not in folder_files:
            missing.append(cif)
    
    return missing


def load_data(cif_ids, npdescs_dir, root, cif_dir, max_sites=None,
              bin_min=1.0, bin_max=2.5, n_bins=100):
    """
    Load features and structure types for selected CIFs.
    
    Parameters:
        cif_ids: List of CIF identifiers (without .cif extension)
        npdescs_dir: Directory with .npz files
        root: Root directory
        cif_dir: Subdirectory with CIF files
        max_sites: Maximum sites per CIF
        bin_min, bin_max, n_bins: Histogram parameters
    
    Returns:
        Tuple of (data_array, structure_types_array)
    """
    data = []
    names = []
    
    for cif_id in cif_ids:
        npz_path = os.path.join(npdescs_dir, f"{cif_id}.npz")
        
        if not os.path.exists(npz_path):
            print(f"Warning: {npz_path} not found, skipping")
            continue
        
        desc = preprocess_desc(npz_path, max_sites=max_sites,
                              bin_min=bin_min, bin_max=bin_max, n_bins=n_bins)
        
        if desc is None:
            continue
        
        # Load structure type from CIF
        cif_path = os.path.join(root, cif_dir, f"{cif_id}.cif")
        try:
            cif = Cif(cif_path)
            structure_type = cif.structure
        except Exception as e:
            print(f"Warning: Could not load structure type for {cif_id}: {e}")
            structure_type = f"unknown-{cif_id}"
        
        data.append(desc)
        names.append(structure_type)
    
    return np.vstack(data), np.array(names)


def perform_clustering(data, linkage="ward"):
    """
    Perform hierarchical clustering.
    
    Parameters:
        data: Feature matrix (n_samples, n_features)
        linkage: Linkage method (ward, complete, average)
    
    Returns:
        Fitted AgglomerativeClustering model
    """
    model = AgglomerativeClustering(
        distance_threshold=1e-5,
        n_clusters=None,
        linkage=linkage,
        compute_full_tree=True
    )
    model.fit(data)
    return model


def plot_dendrogram(model, names, highlight=None, truncate_mode="level", p=30, **kwargs):
    """
    Plot dendrogram for hierarchical clustering.
    
    Parameters:
        model: Fitted AgglomerativeClustering model
        names: Array of labels for each sample
        highlight: List of labels to highlight in red
        truncate_mode, p: Dendrogram truncation parameters
        **kwargs: Additional arguments for dendrogram
    """
    if highlight is None:
        highlight = []
    
    # Build linkage matrix
    counts = np.zeros(model.children_.shape[0])
    n_samples = len(model.labels_)
    for i, merge in enumerate(model.children_):
        current_count = 0
        for child_idx in merge:
            if child_idx < n_samples:
                current_count += 1
            else:
                current_count += counts[child_idx - n_samples]
        counts[i] = current_count
    
    linkage_matrix = np.column_stack(
        [model.children_, model.distances_, counts]
    ).astype(float)
    
    # Plot dendrogram
    dendrogram(linkage_matrix, labels=names, leaf_rotation=90,
              truncate_mode=truncate_mode, p=p, **kwargs)
    
    # Highlight specific labels
    ax = plt.gca()
    xticklabels = ax.get_xticklabels()
    for lbl in xticklabels:
        text = lbl.get_text()
        if text in highlight:
            lbl.set_color("red")


def save_results(cif_ids, structure_types, labels, output_path):
    """
    Save clustering results to CSV.
    
    Parameters:
        cif_ids: List of CIF identifiers
        structure_types: List of structure types
        labels: Cluster labels from model
        output_path: Output file path
    """
    df = pd.DataFrame({
        'cif': cif_ids,
        'structure_type': structure_types,
        'cluster': labels
    })
    df.to_csv(output_path, index=False)


def log_error(error_log_path, cif_id, error_msg):
    """
    Log error to file.
    
    Parameters:
        error_log_path: Path to error log file
        cif_id: CIF identifier
        error_msg: Error message
    """
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    with open(error_log_path, 'a') as f:
        f.write(f"{timestamp} | {cif_id} | {error_msg}\n")


def check_missing_npz(cif_ids, npdescs_dir):
    """
    Check which CIFs are missing .npz files.
    
    Parameters:
        cif_ids: List of CIF identifiers
        npdescs_dir: Directory with .npz files
    
    Returns:
        List of CIF IDs without .npz files
    """
    return [cid for cid in cif_ids 
            if not os.path.exists(os.path.join(npdescs_dir, f"{cid}.npz"))]


def setup_mpl_style(user_style_path):
    """
    Setup MPL style, using embedded if user path not provided.
    
    Parameters:
        user_style_path: User-provided style path or None
    
    Returns:
        Path to style file to use
    """
    if user_style_path and os.path.exists(user_style_path):
        return user_style_path
    
    # Save embedded style to temp file
    tmp = tempfile.NamedTemporaryFile(mode='w', suffix='.mplstyle', delete=False)
    tmp.write(PAPER_MPLSTYLE)
    tmp.close()
    return tmp.name


def create_output_dirs(root, npdescs_dir, plots_dir):
    """
    Create output directories if they don't exist.
    
    Parameters:
        root: Root directory
        npdescs_dir: Features directory
        plots_dir: Plots directory
    """
    os.makedirs(os.path.join(root, npdescs_dir), exist_ok=True)
    os.makedirs(os.path.join(root, plots_dir), exist_ok=True)


def main():
    parser = argparse.ArgumentParser(
        description="CIF Clustering Analysis Tool",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
    
    # Extract features and cluster
    python cif_clustering.py --root /path/to/data --cif-dir CIFs-cleaned-AB2 
    
    # Process specific CIFs from CSV
    python cif_clustering.py --root /path/to/data --cif-dir CIFs-cleaned-AB2
    
    # Filter by sites (e.g., only A and B sites)
    python cif_clustering.py --root /path/to/data --cif-dir CIFs-cleaned-AB2 --sites "A B"

CSV format (if using --csv):
    ,cif,formula,stype,As,Bs
    0,528761.cif,GdPt2,"MgCu2,cF24,227",Gd,Pt
    1,528748.cif,Pt3U,"Mg3Cd,hP8,194",U,Pt

The 'cif' column must contain CIF filenames (with or without .cif extension).
        """
    )
    
    # Required arguments
    parser.add_argument('--cif-dir', required=True, type=str,
                       help='CIF files subdirectory (e.g., CIFs-cleaned-AB2)')
    
    # CSV input
    parser.add_argument('--csv', type=str, default=None,
                       help='Optional CSV file with \'cif\' column. If not provided, all CIFs in folder are processed.')
    
    # Feature extraction options
    parser.add_argument('--sites', type=str, default=None,
                       help='Site filter: space-separated element symbols (e.g., "A B"). If None, all sites are included.')
    parser.add_argument('--n-neighbors', type=int, default=100,
                       help='Number of neighbors for ddmin calculation (default: 100)')
    parser.add_argument('--force-recalc', action='store_true',
                       help='Recalculate features even if .npz files exist')
    parser.add_argument('--plot-individual', action='store_true',
                       help='Generate individual CIF histogram plots during feature extraction')
    parser.add_argument('--workers', type=int, default=16,
                       help='Number of multiprocessing workers (default: 16)')
    
    # Preprocessing options
    parser.add_argument('--max-sites', type=int, default=30,
                       help='Maximum sites per CIF (default: 30, use 0 for no limit)')
    parser.add_argument('--bin-min', type=float, default=1.0,
                       help='Histogram minimum value (default: 1.0)')
    parser.add_argument('--bin-max', type=float, default=2.5,
                       help='Histogram maximum value (default: 2.5)')
    parser.add_argument('--n-bins', type=int, default=100,
                       help='Number of histogram bins (default: 100)')
    
    # Clustering options
    parser.add_argument('--linkage', type=str, default='ward',
                       choices=['ward', 'complete', 'average'],
                       help='Clustering linkage method (default: ward)')
    parser.add_argument('--truncate-level', type=int, default=30,
                       help='Dendrogram truncation level (default: 30)')
    parser.add_argument('--highlight', type=str, default=None,
                       help='Comma-separated list of structure types to highlight in dendrogram')
    
    # Plotting
    parser.add_argument('--mpl-style', type=str, default=None,
                       help='Path to MPL style file')
    
    # Output
    parser.add_argument('--output', type=str, default='cluster_results.csv',
                       help='Output file for clustering results (default: cluster_results.csv)')
    
    args = parser.parse_args()
    
    # Parse sites filter
    sites = args.sites.split() if args.sites else None
    
    # Parse highlight list
    highlight = [h.strip() for h in args.highlight.split(',')] if args.highlight else []
    
    # Setup paths
    cif_folder = args.cif_dir
    npdescs_dir = 'outputs/descriptors'
    plots_dir = 'outputs/plots'
    error_log_path = 'cifs_encountered_errors.txt'
    
    # Max sites: 0 means no limit
    max_sites = None if args.max_sites == 0 else args.max_sites
    
    # Validate cif-dir exists
    if not os.path.exists(cif_folder):
        print(f"Error: CIF directory does not exist: {cif_folder}")
        sys.exit(1)
    
    # Determine CIF list
    if args.csv:
        if not os.path.exists(args.csv):
            print(f"Error: CSV file does not exist: {args.csv}")
            sys.exit(1)
        
        df = pd.read_csv(args.csv)
        if 'cif' not in df.columns:
            print(f"Error: CSV must have a 'cif' column")
            sys.exit(1)
        
        # Validate CIFs exist in folder
        missing = validate_cifs(df, cif_folder)
        if missing:
            print(f"Error: {len(missing)} CIF(s) not found in folder:")
            for m in missing[:10]:
                print(f"  - {m}")
            if len(missing) > 10:
                print(f"  ... and {len(missing) - 10} more")
            sys.exit(1)
        
        # Get CIF IDs (normalize to have .cif extension for loading)
        cif_files = df['cif'].tolist()
        cif_ids = [os.path.basename(f).replace('.cif', '').replace('.CIF', '') for f in cif_files]
    
    else:
        # Scan folder
        df = get_df_data(cif_folder)
        print(f"Found {len(df)} CIF files in {cif_folder}")
        cif_ids = [os.path.basename(f).replace('.cif', '').replace('.CIF', '') for f in df['cif']]
    
    print(f"Processing {len(cif_ids)} CIF(s)")
    
    # Setup MPL style
    mpl_style_path = setup_mpl_style(args.mpl_style)
    
    # Create output directories
    create_output_dirs('.', npdescs_dir, plots_dir)
    
    # Feature Extraction Phase
    print("\n=== Feature Extraction ===")
    
    # Determine which CIFs need feature extraction
    if args.force_recalc:
        to_process = cif_ids.copy()
        print("Force recalculation enabled, processing all CIFs")
    else:
        to_process = check_missing_npz(cif_ids, npdescs_dir)
        print(f"Missing .npz files: {len(to_process)}")
    
    if to_process:
        # Prepare tasks
        tasks = []
        for cif_id in to_process:
            cif_path = os.path.join(cif_folder, f"{cif_id}.cif")
            tasks.append((
                cif_path,
                sites,
                args.n_neighbors,
                npdescs_dir,
                args.plot_individual,
                args.bin_min,
                args.bin_max,
                args.n_bins,
                plots_dir,
                mpl_style_path,
                '.'
            ))
        
        # Process in parallel
        print(f"Processing {len(tasks)} CIF(s) with {args.workers} workers...")
        with mp.Pool(args.workers) as pool:
            results = pool.map(write_features, tasks)
        
        # Log errors
        success_count = 0
        error_count = 0
        for cif_id, success, err_msg in results:
            if success:
                success_count += 1
            else:
                error_count += 1
                log_error(error_log_path, cif_id, err_msg)
        
        print(f"Feature extraction complete: {success_count} succeeded, {error_count} failed")
        if error_count > 0:
            print(f"Errors logged to: {error_log_path}")
    
    print("Feature extraction phase complete\n")
    
    # Clustering Phase
    print("=== Clustering ===")
    
    # Check for npdescs-dir
    if not os.path.exists(npdescs_dir):
        print(f"Error: Features directory does not exist: {npdescs_dir}")
        print("No features were extracted. Check error log for details.")
        sys.exit(1)
    
    # Load data
    data, names = load_data(
        cif_ids,
        npdescs_dir,
        '.',
        args.cif_dir,
        max_sites=max_sites,
        bin_min=args.bin_min,
        bin_max=args.bin_max,
        n_bins=args.n_bins
    )
    
    if len(data) == 0:
        print("Error: No data loaded for clustering")
        sys.exit(1)
    
    print(f"Loaded {len(data)} CIF(s) with {data.shape[1]} features each")
    
    # Filter cif_ids to match loaded data (already filtered in load_data)
    cif_ids_filtered = cif_ids[:len(data)]
    
    # Perform clustering
    model = perform_clustering(data, linkage=args.linkage)
    labels = model.labels_
    n_clusters = len(np.unique(labels))
    
    print(f"Clustering complete: {n_clusters} cluster(s) identified")
    print(f"Labels: {labels}")
    print(f"Structure types: {names.tolist()}")
    
    # Plot dendrogram
    plt.figure(figsize=(20, 10))
    plt.title("Hierarchical Clustering Dendrogram")
    plot_dendrogram(
        model, names,
        highlight=highlight,
        truncate_mode="level",
        p=args.truncate_level
    )
    plt.xlabel("Structure type")
    plt.tight_layout()
    
    dendro_path = os.path.join(plots_dir, "cluster_dendrogram.png")
    os.makedirs(plots_dir, exist_ok=True)
    plt.savefig(dendro_path, dpi=300)
    plt.close('all')
    
    print(f"Dendrogram saved to: {dendro_path}")
    
    # Save results
    output_path = f"outputs/{args.output}"
    save_results(cif_ids_filtered, names, labels, output_path)
    print(f"Cluster assignments saved to: {output_path}")
    
    print("\nClustering complete")


if __name__ == "__main__":
    main()
