# MLEVar_OSLR
Accounting for variability in reference curves from MLE estimates in one-sample log-rank tests

This repository contains the following R scripts:
- `t1e_sim.R`: Simulation study to analyse type I error rates of corrected and uncorrected one-sample log-rank tests for which the reference curve is determined by MLE estimates. One single scenario with fixed shape parameter of a Weibull distribution and fixed sample sizes in both groups is considered. The scale parameter of the Weibull distribution is chosen in such a way that the 1-year survival probability is given by 50%. Raw study results are produced and saved in the folder `results/single_scenarios`.
- `analyse_t1e.R`: Aggregation of raw simulation results from `t1e_sim.R` and production of plots for type I error rates. Plots are saved in the folder `results/plots`.
- `case_study.R`: Case study applying corrected and uncorrected one-sample log-rank tests to a real dataset. Underlying data is saved in the folder `data`. Due to copyright reasons, we cannot provide our original reconstructed data. Hence, we provide simulated data that is inspired by the recosntructed data. The resulting plot is saved in the folder `results/case_study`.

This repository contains the following data:
- In the folder `data`: simulated data inspired by reconstructed dataset is provided for the case study in `case_study.R`.

This repository contains the following results (all saved in the folder `results`):
- In the folder `single_scenarios`: Raw results from the simulation study in `t1e_sim.R`. These results will be overwritten if the simulation study is re-run.
- In the folder `plots`: Plots for type I error rates produced in `analyse_t1e.R`.
- In the folder `case_study`: Plot produced in `case_study.R`.


