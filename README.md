# NetworkOverlap-Task-RS-fMRI
# Spatial Overlap Analysis: Resting-State Connectivity and Task-Evoked Activation

This repository contains the custom Python scripts for the spatial overlap analysis presented in the manuscript:
**"[Your Manuscript Title]"** (Preprint/Submitted to [Journal Name]).

This code quantifies the spatial correspondence between intrinsic functional connectivity (FC) maps derived from resting-state fMRI and task-evoked activation maps. The pipeline is applied to two distinct functional networks: a visuo-spatial memory network (rOTC seed & OLM task) and a language learning network (lvIFG seed & APPL task).


# Spatial Overlap Analysis for tDCS-fMRI Studies

A Python tool for quantifying spatial overlap between resting-state functional connectivity networks and task-evoked activation patterns. Supports target selection for focal tDCS studies.

## Quick Start

### 1. Install Dependencies
```bash
pip install nibabel nilearn scipy numpy matplotlib pyyaml
2. Create Configuration File
Create config/config.yaml with your paths:

yaml
paths:
  task_map: "/path/to/your/task_statistical_map.nii.gz"
  rest_map: "/path/to/your/resting_state_connectivity_map.nii.gz"
  mask: "/path/to/your/brain_mask.nii.gz"  # Optional
  output_dir: "/path/to/your/results/"

analysis:
  task_df: 17
  rest_df: 36
  task_fwe_threshold: 3.6459
  rest_fwe_threshold: 3.3328
  binary_threshold_alpha: 0.001
3. Run Analysis
bash
python scripts/spatial_overlap_analysis.py
Repository Structure
text
tDCS-fMRI-SpatialOverlap/
├── scripts/
│   └── spatial_overlap_analysis.py     # Main analysis script
├── config/
│   └── config.yaml                     # Your configuration
└── README.md
Configuration Parameters
task_map: Task statistical map (t-map)

rest_map: Resting-state connectivity map

mask: Brain mask (optional)

output_dir: Results directory

task_df/rest_df: Degrees of freedom for t-tests

task_fwe_threshold/rest_fwe_threshold: FWE-corrected t-values

binary_threshold_alpha: Alpha for binary overlap (default: 0.001)

Output Files
spatial_overlap_metrics.txt: Dice coefficient, correlation, voxel counts

task_specific.nii.gz: Regions only significant in task

rest_specific.nii.gz: Regions only in resting-state network

overlap.nii.gz: Overlapping regions

combined_map.nii.gz: Combined visualization (1=task, 2=rest, 3=overlap)

threshold_sensitivity.csv: Dice values across thresholds

threshold_sensitivity.png: Sensitivity plot

Command Line Options
bash
# Use default config
python scripts/spatial_overlap_analysis.py

# Custom config file
python scripts/spatial_overlap_analysis.py --config config/my_config.yaml

# Custom output directory
python scripts/spatial_overlap_analysis.py --output-dir /path/to/results

# Skip plot generation
python scripts/spatial_overlap_analysis.py --no-plot
Method
Calculates spatial overlap using:

Dice Similarity Coefficient: Overlap of thresholded binary maps

Pearson Correlation: Voxel-wise spatial correlation

Threshold Sensitivity: Dice values across statistical thresholds

Citation
If using this tool, please cite:

bibtex
@software{spatialoverlap2024,
  title = {Spatial Overlap Analysis for tDCS-fMRI Studies},
  author = {Your Name},
  year = {2024},
  url = {https://github.com/yourusername/tDCS-fMRI-SpatialOverlap}
}
License
MIT License

text
