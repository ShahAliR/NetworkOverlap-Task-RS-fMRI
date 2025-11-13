#!/usr/bin/env python3
"""
Spatial Overlap Analysis for fMRI Data

This script calculates spatial overlap between resting-state functional connectivity
maps and task-evoked activation maps using Dice coefficient and Pearson correlation.

Author: Dr. Alireza Shahbabaie
Organization: Greifswald Medical University
License: MIT
Date: 2025
"""

import numpy as np
import nibabel as nib
from nilearn import image
from scipy.stats import pearsonr, t
import os
import matplotlib.pyplot as plt
import argparse
import yaml
import sys
from pathlib import Path

# Suppress warnings
import warnings
warnings.filterwarnings('ignore', category=RuntimeWarning)


def load_config(config_path):
    """Load configuration from YAML file"""
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    return config


def setup_argparse():
    """Set up command line argument parsing"""
    parser = argparse.ArgumentParser(
        description='Spatial overlap analysis between task and resting-state fMRI maps'
    )
    parser.add_argument('--config', type=str, required=True,
                       help='Path to configuration YAML file')
    parser.add_argument('--output-dir', type=str,
                       help='Override output directory from config')
    parser.add_argument('--task-thresh', type=float,
                       help='Override task FWE threshold')
    parser.add_argument('--rest-thresh', type=float,
                       help='Override rest FWE threshold')
    parser.add_argument('--no-plot', action='store_true',
                       help='Skip generating plots')
    return parser


def load_and_align_maps(task_path, rest_path, mask_path=None):
    """Load and align maps to task space with optional masking"""
    # Load images
    task_img = nib.load(task_path)
    rest_img = nib.load(rest_path)
    
    # Resample resting-state to task space
    rest_img = image.resample_to_img(rest_img, task_img, interpolation='continuous')
    
    # Load or create mask
    if mask_path and os.path.exists(mask_path):
        mask = nib.load(mask_path).get_fdata() > 0
        mask_img = nib.Nifti1Image(
            mask.astype(np.int16), 
            task_img.affine, 
            header=task_img.header
        )
        mask = image.resample_to_img(
            mask_img,
            task_img,
            interpolation='nearest'
        ).get_fdata() > 0
    else:
        mask = (task_img.get_fdata() != 0) | (rest_img.get_fdata() != 0)
        if mask_path:
            print(f"Warning: Mask path {mask_path} not found. Using auto-generated mask.")
    
    # Validate mask
    mask_size = np.sum(mask)
    print(f"Mask contains {mask_size} voxels")
    if mask_size < 10:
        raise ValueError("Mask is too small - check mask alignment")
    
    return task_img, rest_img, mask, task_img.affine


def compute_spatial_metrics(task_data, rest_data, mask, task_thresh, rest_thresh):
    """Calculate spatial similarity metrics with proper t-thresholds"""
    # Clean data: replace NaNs and Infs
    task_data = np.nan_to_num(task_data, nan=0.0, posinf=0.0, neginf=0.0)
    rest_data = np.nan_to_num(rest_data, nan=0.0, posinf=0.0, neginf=0.0)
    
    # Extract masked data
    task_masked = task_data[mask]
    rest_masked = rest_data[mask]
    
    # Initialize metrics with fallback values
    dice, corr_coef, p_value = 0.0, np.nan, np.nan
    task_bin = np.zeros_like(task_data, dtype=bool)
    rest_bin = np.zeros_like(rest_data, dtype=bool)
    
    # Check variance before computing correlation
    if len(task_masked) > 10:  # Need sufficient voxels
        task_std = np.std(task_masked)
        rest_std = np.std(rest_masked)
        
        if task_std > 1e-6 and rest_std > 1e-6:
            corr_coef, p_value = pearsonr(task_masked, rest_masked)
        else:
            print(f"Low variance detected: task_std={task_std:.4f}, rest_std={rest_std:.4f}")
    
    # Binary masks at specified t-thresholds
    task_bin = (task_data >= task_thresh) & mask
    rest_bin = (rest_data >= rest_thresh) & mask
    
    # Overlap calculations
    intersection = np.sum(task_bin & rest_bin)
    
    # Calculate Dice only if there are voxels
    total = np.sum(task_bin) + np.sum(rest_bin)
    if total > 0:
        dice = (2 * intersection) / total
    else:
        print("No significant voxels for Dice calculation")
    
    return dice, corr_coef, p_value, task_bin, rest_bin


def threshold_sensitivity_analysis(task_data, rest_data, mask, min_t=2.0, max_t=8.0, step=0.25):
    """Compute Dice coefficients across threshold range"""
    thresholds = np.arange(min_t, max_t + step, step)
    dice_results = []
    
    # Clean data
    task_data = np.nan_to_num(task_data, nan=0.0, posinf=0.0, neginf=0.0)
    rest_data = np.nan_to_num(rest_data, nan=0.0, posinf=0.0, neginf=0.0)
    
    for t_val in thresholds:
        task_bin = (task_data >= t_val) & mask
        rest_bin = (rest_data >= t_val) & mask
        
        intersection = np.sum(task_bin & rest_bin)
        total = np.sum(task_bin) + np.sum(rest_bin)
        
        dice = (2 * intersection) / total if total > 0 else 0
        dice_results.append(dice)
    
    return thresholds, dice_results


def save_results(results, config, generate_plot=True):
    """Save all results to output directory"""
    output_dir = Path(config['paths']['output_dir'])
    output_dir.mkdir(parents=True, exist_ok=True)
    
    output_config = config['output']
    
    # Save statistical maps
    for name, data, affine in results['maps']:
        output_filename = output_config['map_names'].get(name, name)
        img = nib.Nifti1Image(data, affine)
        nib.save(img, output_dir / f"{output_filename}.nii.gz")
    
    # Save metrics report
    report_path = output_dir / output_config['report_name']
    with open(report_path, "w") as f:
        f.write("SPATIAL OVERLAP ANALYSIS REPORT\n")
        f.write("=" * 50 + "\n")
        f.write(f"Analysis: {config['analysis'].get('description', 'Task vs Resting-State')}\n")
        f.write(f"Task DF: {config['analysis']['task_df']}, Rest DF: {config['analysis']['rest_df']}\n")
        f.write(f"FWE-corrected thresholds: task t > {config['analysis']['task_fwe_threshold']:.3f}, "
                f"rest t > {config['analysis']['rest_fwe_threshold']:.3f}\n")
        f.write(f"Binary overlap thresholds: task t > {results['task_binary_thresh']:.3f}, "
                f"rest t > {results['rest_binary_thresh']:.3f}\n")
        f.write("=" * 50 + "\n")
        f.write(f"Dice Coefficient: {results['dice']:.4f}\n")
        
        # Handle NaN in correlation
        if np.isnan(results['correlation']):
            f.write("Spatial Correlation: Not computable (constant input)\n")
        else:
            f.write(f"Spatial Correlation: {results['correlation']:.4f} (p={results['corr_p']:.6f})\n")
        
        f.write(f"Significant Voxels (task): {results['task_voxels']}\n")
        f.write(f"Significant Voxels (rest): {results['rest_voxels']}\n")
        f.write(f"Overlapping Voxels: {results['overlap_voxels']}\n")
        f.write(f"Analysis completed: {results['timestamp']}\n")
    
    # Save threshold sensitivity data
    if 'thresholds' in results and 'dice_values' in results:
        sensitivity_path = output_dir / output_config['sensitivity_csv']
        np.savetxt(
            sensitivity_path,
            np.column_stack([results['thresholds'], results['dice_values']]),
            delimiter=",",
            header="Threshold,Dice"
        )
    
    # Plot threshold sensitivity
    if generate_plot and 'thresholds' in results and 'dice_values' in results:
        plt.figure(figsize=(10, 6))
        plt.plot(results['thresholds'], results['dice_values'], 'b-o', linewidth=2, markersize=4)
        plt.axvline(x=config['analysis']['task_fwe_threshold'], color='r', linestyle='--', 
                   label=f'FWE Task Threshold ({config["analysis"]["task_fwe_threshold"]:.3f})')
        plt.axvline(x=config['analysis']['rest_fwe_threshold'], color='g', linestyle='--', 
                   label=f'FWE Rest Threshold ({config["analysis"]["rest_fwe_threshold"]:.3f})')
        plt.xlabel("t-threshold")
        plt.ylabel("Dice Coefficient")
        plt.title("Threshold Sensitivity Analysis")
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        
        plot_path = output_dir / output_config['sensitivity_plot']
        plt.savefig(plot_path, dpi=300, bbox_inches='tight')
        plt.close()


def run_analysis(config):
    """Run the complete spatial overlap analysis"""
    paths = config['paths']
    analysis_config = config['analysis']
    sensitivity_config = config.get('sensitivity_analysis', {})
    
    # Calculate binary thresholds
    task_binary_thresh = t.ppf(1 - analysis_config['binary_threshold_alpha']/2, analysis_config['task_df'])
    rest_binary_thresh = t.ppf(1 - analysis_config['binary_threshold_alpha']/2, analysis_config['rest_df'])
    
    # Initialize results container
    results = {
        'dice': 0.0, 
        'correlation': np.nan, 
        'corr_p': np.nan,
        'thresholds': None, 
        'dice_values': None,
        'task_voxels': 0, 
        'rest_voxels': 0, 
        'overlap_voxels': 0,
        'task_binary_thresh': task_binary_thresh,
        'rest_binary_thresh': rest_binary_thresh,
        'maps': [],
        'timestamp': np.datetime64('now')
    }
    
    print("=" * 60)
    print("SPATIAL OVERLAP ANALYSIS")
    print("=" * 60)
    
    # 1. Load and preprocess data
    print("Loading and aligning data...")
    task_img, rest_img, mask, affine = load_and_align_maps(
        paths['task_map'], paths['rest_map'], paths.get('mask')
    )
    task_data = task_img.get_fdata()
    rest_data = rest_img.get_fdata()
    
    print(f"Data shape: {task_data.shape}")
    print(f"Task FWE threshold: t > {analysis_config['task_fwe_threshold']:.3f}")
    print(f"Rest FWE threshold: t > {analysis_config['rest_fwe_threshold']:.3f}")
    print(f"Binary thresholds: task t > {task_binary_thresh:.3f}, rest t > {rest_binary_thresh:.3f}")
    
    # Clean data immediately
    task_data = np.nan_to_num(task_data, nan=0.0, posinf=0.0, neginf=0.0)
    rest_data = np.nan_to_num(rest_data, nan=0.0, posinf=0.0, neginf=0.0)
    
    # 2. Compute spatial metrics
    print("Calculating spatial metrics...")
    (results['dice'], 
     results['correlation'], results['corr_p'],
     task_bin, rest_bin) = compute_spatial_metrics(
         task_data, rest_data, mask, task_binary_thresh, rest_binary_thresh
     )
    
    # Count significant voxels
    results['task_voxels'] = np.sum(task_bin)
    results['rest_voxels'] = np.sum(rest_bin)
    results['overlap_voxels'] = np.sum(task_bin & rest_bin)
    
    # 3. Threshold sensitivity analysis
    print("Running threshold sensitivity analysis...")
    results['thresholds'], results['dice_values'] = threshold_sensitivity_analysis(
        task_data, rest_data, mask, 
        min_t=sensitivity_config.get('min_t', 2.0), 
        max_t=sensitivity_config.get('max_t', 8.0), 
        step=sensitivity_config.get('step', 0.25)
    )
    
    # 4. Create output maps
    print("Generating output maps...")
    task_only = (task_bin == 1) & (rest_bin == 0)
    rest_only = (rest_bin == 1) & (task_bin == 0)
    overlap = (task_bin == 1) & (rest_bin == 1)
    
    # Combined visualization map
    combined_map = np.zeros_like(task_bin, dtype=np.int16)
    combined_map[task_only] = 1    # Task only
    combined_map[rest_only] = 2    # Rest only
    combined_map[overlap] = 3      # Both
    
    # Save all components
    results['maps'] = [
        ('task_specific', task_only.astype(np.int16), affine),
        ('rest_specific', rest_only.astype(np.int16), affine),
        ('overlap', overlap.astype(np.int16), affine),
        ('combined_map', combined_map, affine),
        ('continuous_overlap', (task_data * rest_data).astype(np.float32), affine)
    ]
    
    return results


def main():
    """Main analysis pipeline"""
    try:
        # Parse arguments
        args = setup_argparse().parse_args()
        config = load_config(args.config)
        
        # Override config with command line arguments if provided
        if args.output_dir:
            config['paths']['output_dir'] = args.output_dir
        if args.task_thresh:
            config['analysis']['task_fwe_threshold'] = args.task_thresh
        if args.rest_thresh:
            config['analysis']['rest_fwe_threshold'] = args.rest_thresh
            
        # Run analysis
        results = run_analysis(config)
        
        # Save results
        save_results(results, config, generate_plot=not args.no_plot)
        
        # Print summary
        print("\n" + "=" * 50)
        print("ANALYSIS SUMMARY")
        print("=" * 50)
        print(f"Dice Coefficient: {results['dice']:.4f}")
        if not np.isnan(results['correlation']):
            print(f"Spatial Correlation: {results['correlation']:.4f} (p={results['corr_p']:.6f})")
        print(f"Overlapping Voxels: {results['overlap_voxels']}")
        print(f"Results saved to: {config['paths']['output_dir']}")
        
    except FileNotFoundError as e:
        print(f"Error: File not found - {e}")
        sys.exit(1)
    except Exception as e:
        print(f"Error: {str(e)}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
