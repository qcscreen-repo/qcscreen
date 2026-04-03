# Model-free Feature Screening and False Discovery Control for High-dimensional Quantile Regressions
This GitHub repository contains code and data to reproduce the simulation results and the real data analysis reported in the paper *Model-free Feature Screening and False Discovery Control for High-dimensional Quantile Regressions*.

## Package Environments
The codes in this repository requires **R 4.5.2** and the following R packages (versions as used in this rebuttal):
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

## simulation-code
This folder contains the R scripts used in the rebuttal to illustrate and partially reproduce the simulation studies reported in **Section 5.1–5.3** and **Appendices B.1, B.2, and B.4** of the paper. Generated outputs (`.csv`, `.rds`, figures, etc.) are kept under **`Results/Simulation-code`**, using the same `Part 1` / `Part 2` / `Part 3` / `Part 4` layout as here—configure the scripts’ working directory and save paths accordingly when you rerun them.

- **Part 1** reproduces the simulation results in **Section 5.1** (Table 1-3) and contributes to **Appendix B.1** (Table B.1-B.3).  
- **Part 2** reproduces the simulation results in **Section 5.2** (Table 4).  
- **Part 3** reproduces the simulation results in **Section 5.3** (Tables 5–6) and contributes to the remaining results in **Appendix B.2** (Table B.4).  
- **Part 4** reproduces the additional simulation study reported in **Appendix B.4** (Figure B.1).

1. Under the ```Part 1``` folder, there are three subfolders corresponding to different models:
   - ```Model 1a``` contains scripts such as ```model_1a_n_400_p_1000.R``` and ```model_1a_n_400_p_5000.R``` for Model 1(a) with \(n = 400\) and \(p = 1000, 5000\). Each script runs multiple signal strengths \(\alpha = 0.5, 0.75\).
   - ```Model 1b``` contains scripts such as ```model_1b_n_400_p_1000.R``` and ```model_1b_n_400_p_5000_alpha_05.R```, ```model_1b_n_400_p_5000_alpha_075.R``` for Model 1(b) with different dimensionalities and signal strengths \(\alpha = 0.5, 0.75\).
   - ```Model 1c``` contains scripts such as ```model_1c_n_400_p_1000_alpha_05.R```, ```model_1c_n_400_p_1000_alpha_075.R```, ```model_1c_n_400_p_5000_alpha_05.R``` and ```model_1c_n_400_p_5000_alpha_075.R``` for Model 1(c) with \(\alpha = 0.5, 0.75\).
   - In the Part 1 scripts and their saved results, the key **`fs`** labels **the proposed method** (e.g. ```qc_screen_with_min_model_size``` in the code); other keys refer to benchmark methods.
   - For each script, the corresponding output files (e.g. ```model1_time_n400_p1000 2.csv```, ```model1_results_n400_p1000 2.csv```, ```model2_time_n400_p5000_alpha0_75.csv```, ```model1c_results_n400_p5000_alpha0_5.csv```, etc.) record the runtime and performance metrics used in the tables and figures associated with **Section 5.1** and **Appendix B.1**; in this repository they live under **`Results/Simulation-code/Part 1/...`** (same `Model 1a` / `Model 1b` / `Model 1c` subfolders as the scripts).

2. Under the ```Part 2``` folder, the file ```Table4.R``` implements the simulation setup for Section 5.2.  
   - Running ```Table4.R``` generates the summary `.csv` files ```Table4_p1000.csv``` and ```Table4_p2000.csv```, which contain the numbers reported in **Table 4** of the paper; these are stored under **`Results/Simulation-code/Part 2/`** when using the repo’s output layout.

3. Under the ```Part 3``` folder, the files ```model_3a_alpha_05.R```, ```model_3a_alpha_075.R```, ```model_3b_alpha_05.R```, ```model_3b_alpha_075.R```, ```model_3c_alpha_05.R``` and ```model_3c_alpha_075.R``` implement the simulation setups for the three models in Section 5.3 under two signal strengths \(\alpha = 0.5, 0.75\).  
   - In these scripts, objects and summaries with the **`_fs`** suffix (e.g. ```median_size_fs```, ```Pj_fs```, ```emp_fdr_fs```) correspond to **the proposed method**; analogous labels for other methods use their own suffixes.
   - Running these scripts generates the summary `.rds` files used to produce **Tables 5–6**, and the corresponding entries **Tables B.4** in **Appendix B.2**; save them under **`Results/Simulation-code/Part 3/`** to match this repository.

4. Under the ```Part 4``` folder, the file ```FigureB1.R``` implements the simulation setup for the additional experiment in Appendix B.1.  
   - Running ```FigureB1.R``` generates the summary `.csv` file ```mean_res1_results.csv```, which is used to produce **Figure B.1**; it is kept under **`Results/Simulation-code/Part 4/`** in this layout.

## simulation-demo
The folder **`Simualtion-demo`** (repository root) mirrors **`Simulation-code`** only for **Part 1–Part 3**—there is **no Part 4** demo. The layout is more compact than **`Simulation-code`**: Part 1 uses one script per model family instead of many \((n,p,\alpha)\) files; Part 2 and Part 3 use a single driver each.

By default **each part uses 2 replications** for illustration: **`rep <- 2`** in the Part 1 scripts, **`B <- 2`** in **`Table4.R`**, and **`nrep_total <- 2`** in **`model_3.R`**. Edit those constants to run more replications.

1. **`Part 1/`** — subfolders **`Model 1a`**, **`Model 1b`**, **`Model 1c`** with **`model_1a.R`**, **`model_1b.R`**, **`model_1c.R`** respectively (corresponding to **`Simulation-code/Part 1`**).
2. **`Part 2/`** — **`Table4.R`** (corresponding to **`Simulation-code/Part 2`**).
3. **`Part 3/`** — **`model_3.R`** (corresponding to **`Simulation-code/Part 3`**, replacing separate **`model_3a_*`**, **`model_3b_*`**, **`model_3c_*`** scripts).

## Real-data-code
The folder ```Real-data-code``` contains the data file and R scripts used to reproduce the **real data analysis** in **Section 6** and **Appendix C** of the paper. Figures produced by ```real_data.R``` should be written under **`Results/Real-data-code`** (not alongside the scripts), consistent with the rest of this repository.

1. ```2025-08-MD.csv``` is the FRED-MD macroeconomic dataset (monthly). Obtain an up-to-date release from the [FRED-MD database](https://www.stlouisfed.org/research/economists/mccracken/fred-databases) if needed; the file name or vintage may change with updates.

2. ```real_data.R``` reads the CSV via ```fbi::fredmd``` (with optional transformations and outlier handling), constructs the design matrix and response, and runs the proposed screening procedure together with benchmark methods. The script corresponds to the empirical results and comparisons reported in **Section 6** and **Appendix C**.

3. **Generated figures (PNG).** When you run the plotting sections of the real-data code (e.g. ```png()``` or ```ggsave()```), the following PNG files are produced under **`Results/Real-data-code`** (set paths in ```real_data.R``` to match):
   - ```Figure1.png``` — figures in **Section 6**;
   - ```FigureC1.png``` to ```FigureC5.png``` — figures in **Appendix C**.

Run ```real_data.R``` from the ```Real-data-code``` directory (or set paths accordingly) so that ```2025-08-MD.csv``` is found.

## Results
The ```Results``` folder holds **outputs only**; it mirrors the code layout so paths are easy to find.

1. **`Results/Simulation-code`** — outputs from ```Simulation-code```, with the same ```Part 1``` / ```Part 2``` / ```Part 3``` / ```Part 4``` structure:
   - **Part 1:** runtime and performance `.csv` files under ```Results/Simulation-code/Part 1/Model 1a```, ```Model 1b```, and ```Model 1c``` (see **Section 5.1** and **Appendix B.1**).
   - **Part 2:** ```Table4_p1000.csv``` and ```Table4_p2000.csv``` under ```Results/Simulation-code/Part 2/``` (**Table 4**, Section 5.2).
   - **Part 3:** summary ```.rds``` files under ```Results/Simulation-code/Part 3/``` when generated (**Tables 5–6**, **Table B.4**).
   - **Part 4:** ```FigureB1.png``` under ```Results/Simulation-code/Part 4/``` (**Figure B.1**, Appendix B.4).

2. **`Results/Real-data-code`** — PNG figures (```Figure1.png```, ```FigureC_*.png```), exported tables, and any saved intermediate objects from the real-data analysis (**Section 6**, **Appendix C**). Source data and ```real_data.R``` remain in ```Real-data-code``` at the repo root.

## wrapper-code
This folder holds ```rep_main.R``` and the same ```Part 1```–```Part 4``` layout as ```Simulation-code``` (see above for which tables and figures each part corresponds to).

Run ```rep_main.R``` with the R working directory set to ```wrapper-code``` (e.g. ```setwd(".../wrapper-code")``` then ```source("rep_main.R")```). At the top of the script, Boolean flags turn each block on or off so you can skip slow jobs. Each script is ```source```-d from its own subfolder so outputs are written next to the code.
