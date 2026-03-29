# Image Deblurring Benchmark with Ishikawa-Type Method

This repository contains MATLAB code for image deblurring experiments comparing a proposed Ishikawa-type restoration method with several classical and optimization-based approaches.

## Methods

The following methods are implemented and compared:

* Proposed Ishikawa-type restoration
* Wiener deconvolution
* Lucy–Richardson deconvolution
* TV-ADMM deblurring
* FISTA-based deblurring

## Features

* Motion blur + Gaussian noise simulation
* PSNR and SSIM evaluation
* Visual comparison figures
* Excel result tables
* Sensitivity analysis

## Requirements

* MATLAB
* Image Processing Toolbox

## How to Run

1. Place all `.m` files in the same folder.
2. Set input/output paths in:
   ishikawa_main_selected_clean.m
3. Run:
   ishikawa_main_selected_clean

## Input Images

The script expects the following images:

* 04.png (Starfish)
* 06.png (Plane)
* 09.png (Woman)
* 10.png (Boats)
* 11.png (Pirate)
* 12.png (Couple)

## Output

* Comparison figures (300 DPI)
* PSNR & SSIM tables (Excel)
* Bar graphs
* Sensitivity analysis results

## Reproducibility

* Random seed is fixed (rng(0))
* Images are normalized to [0,1]
* All parameters are defined in the main script
  
### Signal Processing Extension

The proposed Ishikawa-type iterative framework is also implemented for one-dimensional signal enhancement and restoration. This extension demonstrates the flexibility and general applicability of the method beyond image deblurring tasks.

## Notes

AI-assisted tools (ChatGPT) were used to support debugging and code refinement. The author takes full responsibility for the correctness and integrity of the implementation.

## Citation

If you use this code, please cite:

Esra Yolaçan, "A two-step iterative framework for signal and image deblurring using G-I-Nonexpansive Mappings", PLOS ONE (under review).

## License

This project is licensed under the MIT License. See the LICENSE file for details.
