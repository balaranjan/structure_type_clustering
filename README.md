# CIF Clustering Analysis Tool

A unified Python script for performing structural analysis and hierarchical clustering on CIF (Crystallographic Information File) files. This tool extracts d/d_min (normalized interatomic distances) features from atomic structures and clusters similar structures together.

## Features

- **Feature Extraction**: Calculates d/d_min values for atomic site neighborhoods
- **Histogram Generation**: Converts raw features into comparable histograms
- **Hierarchical Clustering**: Groups similar structures using agglomerative clustering
- **Visualization**: Generates dendrograms and individual structure plots
- **Parallel Processing**: Multi-process feature extraction for large datasets
- **Site Filtering**: Optionally analyze only specific element types
- **Incremental Updates**: Skip already-processed files for efficiency

## Installation

### Dependencies

Install the required Python packages:

```bash
pip install numpy pandas matplotlib scipy scikit-learn cifkit
```

### Requirements

- Python 3.7+
- cifkit library for CIF file parsing
- Matplotlib (LaTeX support optional - embedded style works without LaTeX)

## Usage

### Basic Usage

```bash
# Process all CIFs in a folder (extract features + cluster)
python cif_clustering.py --cif-dir test_data

# Process CIFs listed in a CSV file
python cif_clustering.py --cif-dir test_data --csv selected.csv

# Generate individual PDF plots
python cif_clustering.py --cif-dir test_data --csv selected.csv --plot-individual
```

### Processing Selected CIFs

```bash
# Process CIFs listed in a CSV file
python cif_clustering.py --cif-dir test_data --csv selected.csv
```

### Site-Specific Analysis

```bash
# Analyze only A and B element sites
python cif_clustering.py --cif-dir test_data --sites "A B"

# Analyze all sites (default)
python cif_clustering.py --cif-dir test_data
```

## CLI Arguments

### Required Arguments

| Argument | Description |
|-|----|
| `--cif-dir` | CIF files subdirectory (e.g., `test_data`) |

### Input Arguments

| Argument | Default | Description |
|-|-|-----------|
| `--csv` | `None` | Optional CSV file with `'cif'` column. If not provided, all CIFs in folder are processed |


### Feature Extraction Options

| Argument | Default | Description |
|----------|---------|-------------|
| `--sites` | `None` | Space-separated element symbols (e.g., `"A B"`). If `None`, all sites are included |
| `--n-neighbors` | `100` | Number of neighbors for d/d_min calculation |
| `--force-recalc` | `False` | Recalculate features even if `.npz` files exist |
| `--plot-individual` | `False` | Generate individual CIF histogram plots during feature extraction |
| `--workers` | `16` | Number of multiprocessing workers |

### Preprocessing Options

| Argument | Default | Description |
|----------|---------|-------------|
| `--max-sites` | `30` | Maximum sites per CIF (use `0` for no limit) |
| `--bin-min` | `1.0` | Histogram minimum value |
| `--bin-max` | `2.5` | Histogram maximum value |
| `--n-bins` | `100` | Number of histogram bins |

### Clustering Options

| Argument | Default | Description |
|----------|---------|-------------|
| `--linkage` | `ward` | Clustering linkage method (`ward`, `complete`, `average`) |
| `--truncate-level` | `30` | Dendrogram truncation level |
| `--highlight` | `None` | Comma-separated list of structure types to highlight in red in dendrogram |

### Plotting Options

| Argument | Default | Description |
|----------|---------|-------------|
| `--mpl-style` | `None` | Path to custom MPL style file (uses embedded paper.mplstyle if not provided) |

### Output Options

| Argument | Default | Description |
|----------|---------|-------------|
| `--output` | `cluster_results.csv` | Output file for clustering results |

## CSV Format

When using `--csv`, the file must have a `'cif'` column containing CIF filenames:

```csv
,cif,formula,stype,As,Bs
0,528761.cif,GdPt2,"MgCu2,cF24,227",Gd,Pt
1,528748.cif,Pt3U,"Mg3Cd,hP8,194",U,Pt
2,250754.cif,Au3Er,"Cu3Ti,oP8,59",Er,Au
```

**Important Notes:**
- The `cif` column values can include or exclude the `.cif` extension
- All CIFs listed in the CSV must exist in the `--cif-dir` folder
- The CSV is validated before processing begins

## Output Structure

After running the script, the output will be organized as follows:

```
./
├── outputs/
│   ├── descriptors/            # Feature .npz files
│   │   ├── 528761.npz
│   │   ├── 528748.npz
│   │   └── ...
│   └── plots/                  # Plots (created only when needed)
│       ├── cluster_dendrogram.png  # Clustering dendrogram (PNG, 300 DPI)
│       └── individual/         # Per-CIF histogram plots (only with --plot-individual)
│           ├── MgCu2-528761.png
│           └── ...
├── cifs_encountered_errors.txt  # Error log for failed CIFs
└── cluster_results.csv        # Cluster assignments
```

**Note:** The `plots/individual/` folder is only created when using the `--plot-individual` flag.

## Cluster Results Format

The output CSV (`cluster_results.csv`) contains:

```csv
cif,structure_type,cluster
528761,"MgCu2,cF24,227",0
528748,"Mg3Cd,hP8,194",1
...
```

## Error Handling

### Error Log

Failed CIFs are logged to `cifs_encountered_errors.txt`:

```
2024-03-19 10:30:45 | 528761.cif | ValueError: Invalid CIF format
2024-03-19 10:30:46 | 528748.cif | No features extracted
```

### Common Issues

1. **Missing CIFs**: If using `--csv`, all listed CIFs must exist in the folder
2. **No .npz files**: Run `--mode extract_features` before clustering
3. **Empty cluster**: Check if `--max-sites` is too restrictive
4. **LaTeX rendering errors**: Set `text.usetex: False` in MPL style if LaTeX is unavailable

## Advanced Usage

### Custom Histogram Range

```bash
# Change d/d_min range to 1.0-3.0
python cif_clustering.py --root /path --cif-dir CIFs --bin-min 1.0 --bin-max 3.0 --n-bins 200 --mode both
```

### Highlight Specific Structures in Dendrogram

```bash
# Highlight specific structure types in red
python cif_clustering.py --root /path --cif-dir CIFs --highlight "PdSn2,oS24,41;CoGe2,oS24,64" --mode cluster
```

### Use Custom MPL Style

```bash
# Use custom style file
python cif_clustering.py --root /path --cif-dir CIFs --mpl-style /path/to/custom.mplstyle --mode both
```

### Full Pipeline with Custom Settings

```bash
python cif_clustering.py \
    --cif-dir CIFs-cleaned-AB2 \
    --csv selected.csv \
    --sites "A B" \
    --n-neighbors 150 \
    --max-sites 50 \
    --bin-min 1.0 \
    --bin-max 2.5 \
    --n-bins 100 \
    --linkage ward \
    --workers 32 \
    --plot-individual \
    --output results.csv
```

## Technical Details

### d/d_min Calculation

For each atomic site in a CIF:
1. Extract all atomic positions
2. Compute distances from a reference point to all other sites
3. Sort distances and select the first `n` neighbors
4. Normalize all distances by dividing by the minimum distance: `d /= d_min`
5. Result: normalized distance array where all values ≥ 1.0

### Feature Aggregation

Individual site features are aggregated into a single histogram per CIF:
- All site d/d_min values are combined
- Histogram with configurable bin count is created
- Histogram counts are used as the feature vector for clustering

### Hierarchical Clustering

Uses scikit-learn's `AgglomerativeClustering`:
- **Linkage methods**: ward (minimizes variance), complete (max distance), average (mean distance)
- **Distance threshold**: Ensures full tree computation
- **Dendrogram**: Visualizes the clustering hierarchy with structure type labels

## Workflow Diagram

```
┌─────────────────────────────────────────────────────────────┐
│                      Input                                  │
│  --cif-dir (folder)  OR  --csv (selected files)            │
└───────────────────────┬─────────────────────────────────────┘
                        │
                        ▼
            ┌───────────────────────┐
            │  Validate CIF Files   │
            └───────────┬───────────┘
                        │
        ┌───────────────┴───────────────┐
        │                               │
        ▼                               ▼
┌───────────────┐            ┌──────────────────┐
│ extract_features│            │   cluster      │
│   (or both)    │            │  (or both)      │
└───────┬───────┘            └────────┬─────────┘
        │                             │
        ▼                             ▼
┌───────────────┐            ┌──────────────────┐
│ .npz files   │            │ dendrogram      │
│ error log    │            │ results CSV     │
└───────────────┘            └──────────────────┘
```

## Troubleshooting

### Issue: "CIF directory does not exist"

**Solution**: Verify the `--cif-dir` path is correct and accessible.

### Issue: "CSV must have a 'cif' column"

**Solution**: Ensure your CSV has a column named exactly `cif` (case-sensitive).

### Issue: "0 CIF(s) passed validation"

**Solution**: Check that CIFs in the CSV actually exist in the `--cif-dir` folder.

### Issue: MemoryError during feature extraction

**Solutions**:
- Reduce `--workers` value
- Process CIFs in batches using separate CSV files
- Increase system memory

### Issue: Empty cluster results

**Solutions**:
- Check `--max-sites` isn't filtering out all CIFs
- Verify `.npz` files were created successfully
- Check error log for processing failures

## Examples

### Example 1: Quick Start

```bash
python cif_clustering.py --cif-dir test_data
```

### Example 2: Process Selected Structures with Plots

```bash
python cif_clustering.py --cif-dir test_data --csv selected.csv --plot-individual
```

### Example 3: Process All CIFs in Folder

```bash
python cif_clustering.py --cif-dir CIFs-cleaned-AB2 --workers 32
```

### Example 4: Force Recalculation

```bash
# Forces recalculation even for existing .npz files
python cif_clustering.py --cif-dir test_data --force-recalc
```

## License

This tool is provided as-is for research purposes.

## Contributing

Feel free to submit issues and enhancement requests.

## References

- **Agglomerative Clustering**: [scikit-learn documentation](https://scikit-learn.org/stable/modules/clustering.html#hierarchical-clustering)
- **CIF Format**: [International Union of Crystallography](https://www.iucr.org/resources/cif)
