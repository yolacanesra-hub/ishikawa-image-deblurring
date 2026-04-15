# Ishikawa-Type Deblurring Framework for Image and Signal Processing

This repository provides a complete and reproducible MATLAB implementation of an Ishikawa-type iterative framework for image deblurring and signal enhancement.

---

## 🚀 Methods Implemented

The following restoration methods are implemented and compared:

* Proposed Ishikawa-type iterative method
* Wiener deconvolution
* Lucy–Richardson deconvolution
* TV-ADMM deblurring
* FISTA-based deblurring

---

## 📂 Project Structure

* `ishikawa_main_selected_clean.m` → main script
* Helper `.m` files → algorithm components (13+ files)
* `signal_processing/` → 1D signal experiments
* `data/` → **numerical results used in the manuscript (IMPORTANT)**

---

## 📊 Minimal Dataset (PLOS Requirement)

This repository includes the **minimal dataset required to reproduce all results reported in the manuscript**, including:

* PSNR values used in Table 2
* SSIM values used in Table 3
* Iteration data used in convergence analysis (Fig. 4)
* Numerical values used to generate all graphs and figures

All these files are located in:

```
/data
```

These data fully support the reproducibility of the results presented in the article.

---

## 📥 Input Data

### Image Dataset

The benchmark dataset (Set12) used in this study is publicly available at:

https://www.kaggle.com/datasets/leweihua/set12-231008

Expected image files:

* `04.png` (Starfish)
* `06.png` (Plane)
* `09.png` (Woman)
* `10.png` (Boats)
* `11.png` (Pirate)
* `12.png` (Couple)

---

## ⚙️ Features

* Motion blur + Gaussian noise simulation
* Image and signal restoration
* PSNR & SSIM evaluation
* High-resolution (300 DPI) figure generation
* Excel-based result reporting
* Sensitivity analysis

---

## ▶️ How to Run

1. Clone the repository
2. Open MATLAB
3. Set input/output paths in:

```
ishikawa_main_selected_clean.m
```

4. Run:

```matlab
ishikawa_main_selected_clean
```

---

## 🔁 Reproducibility

* Fixed random seed: `rng(0)`
* Input normalization: `[0,1]`
* All parameters explicitly defined in the main script
* All numerical outputs stored and shared in `/data`

---

## 🔬 Signal Processing Extension

The framework is extended to 1D signal restoration, demonstrating robustness across different data modalities.

---

## 📌 Data Availability Statement (PLOS)

All data and code underlying the findings of this study are fully available without restriction.

* Source code and numerical results:
  https://github.com/yolacanesra-hub/ishikawa-deblurring-framework

* The repository contains all numerical data required to reproduce the results (PSNR, SSIM, iteration data).

* The image dataset (Set12) is publicly available from Kaggle:
  https://www.kaggle.com/datasets/leweihua/set12-231008

---

## 📝 Notes

Parts of the development process were supported by AI-assisted tools. The author is responsible for the final implementation and results.

---

## 📌 Citation

Esra Yolaçan
*A two-step iterative framework for signal and image deblurring using G-I-Nonexpansive Mappings*
PLOS ONE (under review)

---

## 📄 License

MIT License
