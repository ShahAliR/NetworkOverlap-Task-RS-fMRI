# NetworkOverlap-Task-RS-fMRI
# Spatial Overlap Analysis: Resting-State Connectivity and Task-Evoked Activation

This repository contains the custom Python scripts for the spatial overlap analysis presented in the manuscript:
**"[Cross-modal imaging-based targeting approach for network level brain stimulation]"** (Preprint/Submitted to [Journal Name]).

This code quantifies the spatial correspondence between intrinsic functional connectivity (FC) maps derived from resting-state fMRI and task-evoked activation maps. The pipeline is applied to two distinct functional networks: a visuo-spatial memory network (rOTC seed & OLM task) and a language learning network (lvIFG seed & APPL task).



## Quick Start

### 1. Install Dependencies

pip install nibabel nilearn scipy numpy matplotlib pyyaml
### 2. Create Configuration File
update config/config.yaml with your paths:

config.yaml:

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
### 3. Run Analysis

python scripts/spatial_overlap_analysis.py

# Repository Structure

"tDCS-fMRI-SpatialOverlap/
├── scripts/
│   └── spatial_overlap_analysis.py     # Main analysis script
├── config/
│   └── config.yaml                     # Your configuration
└── README.md
"

 

# Citation
If using this tool, please cite:

# License
MIT License

 
