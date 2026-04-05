# Model-free Feature Screening and False Discovery Control for High-dimensional Quantile Regressions
This repository provides the code and data required to reproduce the simulation studies and real-data analysis presented in the paper *Model-free Feature Screening and False Discovery Control for High-dimensional Quantile Regressions*.


## Package Environments
The code in this repository was developed and tested under **R 4.5.2** with the following package versions:

1. ```splines``` 
2. ```glmnet 4.1-10```
3. ```KRLS 1.0-0```
4. ```mosum 1.2.7```
5. ```MASS 7.3-65```
6. ```cluster 2.1.8.2```
7. ```fossil 0.4.0```
8. ```arrow 23.0.1.2```
9. ```dplyr 1.2.0```
10. ```ggplot2 4.0.2```
11. ```doParallel 1.0.17```
12. ```foreach 1.5.2```
13. ```Rfast 2.1.5.2```
14. ```dcov 0.1.1```
15. ```quantreg 6.1```
16. ```QCSIS 0.1```
17. ```transport 0.15.4```
18. ```fbi 0.7.0``` 
19. ```quantdr 1.3.2```

## Simulation-code
This folder contains the R scripts used to illustrate and reproduce the simulation results reported in **Sections 5.1–5.3** and **Appendices B.1, B.2, and B.4** of the paper. Generated outputs (`.csv`, `.rds`, figures, etc.) are stored under **`Results/Simulation-code`**, following the same `Part 1` / `Part 2` / `Part 3` / `Part 4` structure. When rerunning the scripts, please configure the working directory and output paths accordingly.

- **Part 1** reproduces the simulation results in **Section 5.1** (Table 1-3) and the additional simulation results in **Appendix B.1** (Table B.1-B.3).  
- **Part 2** reproduces the simulation results in **Section 5.2** (Table 4).  
- **Part 3** reproduces the simulation results in **Section 5.3** (Tables 5–6) and the additional simulation results in **Appendix B.2** (Table B.4).  
- **Part 4** reproduces the additional simulation results in **Appendix B.4** (Figure B.1).

1. Under `Part 1`, there are three subfolders corresponding to different models:
   - `Model 1a` contains scripts such as `model_1a_n_400_p_1000.R` and `model_1a_n_400_p_5000.R` for Model 1(a) with \(n = 400\) and \(p = 1000, 5000\). Each script runs multiple signal strengths \(\alpha = 0.5, 0.75\).
   - `Model 1b` contains scripts such as `model_1b_n_400_p_1000.R`, `model_1b_n_400_p_5000_alpha_05.R`, and `model_1b_n_400_p_5000_alpha_075.R` for Model 1(b) under different dimensionalities and signal strengths \(\alpha = 0.5, 0.75\).
   - `Model 1c` contains scripts such as `model_1c_n_400_p_1000_alpha_05.R`, `model_1c_n_400_p_1000_alpha_075.R`, `model_1c_n_400_p_5000_alpha_05.R`, and `model_1c_n_400_p_5000_alpha_075.R` for Model 1(c) with \(\alpha = 0.5, 0.75\).
   - In the Part 1 scripts and saved results, the key **`fs`** denotes the proposed method (e.g., `qc_screen_with_min_model_size` in the code); other keys correspond to benchmark methods.
   - For each script, the corresponding output files (e.g., `model1_time_n400_p1000 2.csv`, `model1_results_n400_p1000 2.csv`, `model2_time_n400_p5000_alpha0_75.csv`, `model1c_results_n400_p5000_alpha0_5.csv`) record the runtime and performance metrics used in the tables and figures for **Section 5.1** and **Appendix B.1**. In this repository, these files are stored under **`Results/Simulation-code/Part 1/...`**, using the same `Model 1a` / `Model 1b` / `Model 1c` subfolder structure.

2. Under `Part 2`, the file `Table4.R` implements the simulation setup for **Section 5.2**.
   - Running `Table4.R` generates the summary `.csv` files `Table4_p1000.csv` and `Table4_p2000.csv`, which contain the results reported in **Table 4**. These files are stored under **`Results/Simulation-code/Part 2/`** when using the repository’s output layout.

3. Under `Part 3`, the files `model_3a_alpha_05.R`, `model_3a_alpha_075.R`, `model_3b_alpha_05.R`, `model_3b_alpha_075.R`, `model_3c_alpha_05.R`, and `model_3c_alpha_075.R` implement the simulation setups for the three models in **Section 5.3** under two signal strengths, \(\alpha = 0.5, 0.75\).
   - In these scripts, objects and summaries with the **`_fs`** suffix (e.g., `median_size_fs`, `Pj_fs`, `emp_fdr_fs`) correspond to the proposed method; analogous labels for benchmark methods use their own suffixes.
   - Running these scripts generates the summary `.rds` files used to produce **Tables 5–6** and the corresponding entries in **Table B.4** of **Appendix B.2**. Save these files under **`Results/Simulation-code/Part 3/`** to match the repository layout.

4. Under `Part 4`, the file `FigureB1.R` implements the simulation setup for the additional experiment reported in **Appendix B.4**.
   - Running `FigureB1.R` generates the summary `.csv` file `mean_res1_results.csv`, which is used to produce **Figure B.1**. It is stored under **`Results/Simulation-code/Part 4/`** in this repository.

## Simulation-demo
The folder **`Simualtion-demo`** (repository root) mirrors **`Simulation-code`** only for **Part 1–Part 3**—there is **no Part 4** demo. The layout is more compact than **`Simulation-code`**: Part 1 uses one script per model family instead of many \((n,p,\alpha)\) files; Part 2 and Part 3 use a single driver each.

The folder **`Simulation-demo`** (repository root) provides a simplified demonstration version of **`Simulation-code`** for **Parts 1–3**. Compared with **`Simulation-code`**, the layout is more compact: **Part 1** uses one script per model family rather than multiple scripts for different \((n, p, \alpha)\) combinations, while **Parts 2** and **3** each use a single driver script.

By default, each demo uses **2 replications** for illustration:
- `rep <- 2` in the Part 1 scripts,
- `B <- 2` in `Table4.R`,
- `nrep_total <- 2` in `model_3.R`.
Users may edit these values to change the number of replications.

1. **`Part 1/`** — contains three subfolders, **`Model 1a`**, **`Model 1b`**, and **`Model 1c`**, with scripts **`model_1a.R`**, **`model_1b.R`**, and **`model_1c.R`**, respectively. These correspond to **`Simulation-code/Part 1`**.
2. **`Part 2/`** — contains **`Table4.R`**, corresponding to **`Simulation-code/Part 2`**.
3. **`Part 3/`** — contains **`model_3.R`**, corresponding to **`Simulation-code/Part 3`**, and serving as a unified replacement for the separate **`model_3a_*`**, **`model_3b_*`**, and **`model_3c_*`** scripts.


## Real-data-code
The folder `Real-data-code` contains the data file and R scripts used to reproduce the **real-data analysis** in **Section 6** and **Appendix C** of the paper. Figures generated by `real_data.R` should be written to **`Results/Real-data-code`** rather than saved alongside the scripts, for consistency with the rest of the repository.

1. `2025-08-MD.csv` is the monthly FRED-MD macroeconomic dataset. If needed, an updated release can be obtained from the [FRED-MD database](https://www.stlouisfed.org/research/economists/mccracken/fred-databases); note that the some file names or vintages may change due to database updates.

2. `real_data.R` reads the CSV file via `fbi::fredmd` (with optional transformations and outlier handling), constructs the design matrix and response, and runs the proposed screening procedure together with benchmark methods. The script corresponds to the empirical results and comparisons reported in **Section 6** and **Appendix C**.

3. **Generated figures (PNG).** Running the plotting sections of `real_data.R` (e.g., via `png()` or `ggsave()`) produces the following files under **`Results/Real-data-code`**:
   - `Figure1.png` — figures in **Section 6**;
   - `FigureC1.png` to `FigureC5.png` — figures in **Appendix C**.

Run `real_data.R` from the `Real-data-code` directory (or set paths accordingly) so that `2025-08-MD.csv` can be located correctly.


## Results
The `Results` folder stores **generated outputs only**. Its structure mirrors the code layout for easier navigation.

1. **`Results/Simulation-code`** — outputs generated by `Simulation-code`, with the same `Part 1` / `Part 2` / `Part 3` / `Part 4` structure:
   - **Part 1:** runtime and performance `.csv` files under `Results/Simulation-code/Part 1/Model 1a`, `Model 1b`, and `Model 1c` (see **Section 5.1** and **Appendix B.1**).
   - **Part 2:** `Table4_p1000.csv` and `Table4_p2000.csv` under `Results/Simulation-code/Part 2/` (**Table 4**, Section 5.2).
   - **Part 3:** summary `.rds` files under `Results/Simulation-code/Part 3/` (**Tables 5–6**, **Table B.4**).
   - **Part 4:** `FigureB1.png` under `Results/Simulation-code/Part 4/` (**Figure B.1**, Appendix B.4).

2. **`Results/Real-data-code`** — PNG figures (`Figure1.png`, `FigureC*.png`), exported tables, and any saved intermediate objects from the real-data analysis (**Section 6**, **Appendix C**). The source data and `real_data.R` remain in `Real-data-code` at the repository root.

## wrapper-code
This folder contains `rep_main.R` and follows the same `Part 1`–`Part 4` structure as `Simulation-code` (see above for the tables and figures associated with each part).

Run `rep_main.R` with the R working directory set to `wrapper-code` (for example, `setwd(".../wrapper-code")` and then `source("rep_main.R")`). At the top of the script, Boolean flags allow each block to be turned on or off so that slower jobs can be skipped. Each script is then `source()`-d from its corresponding subfolder, with outputs written next to the code.
