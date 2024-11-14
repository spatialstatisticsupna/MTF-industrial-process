# MTF-industrial-process
Data and R code to reproduce the results described in the paper "An Adaptive Learning Approach to Multivariate Time Forecasting in Industrial Processes" (Miguelez et al., 2024).  

# Data
This folder contains the dataset `oee_data_4weeks.Rsd`, a subset of the complete dataset comprising data from week 21 to week 24 of the year 2020, with the following variables:
 - **lsp**: observation period id
 - **day**: starting date
 - **hour**: starting hour
 - **sh.id**: work shift id
 - **new.sh**: 1=first observation of the work shift, 0=otherwise
 - **wday**: weekday
 - **tday**: shift (**M**orning, **A**fternoon, **N**ight)
 - **week.id**: week id
 - **inistop**: 1=the machine is inactive at the beginning of the period, 0=the machine is active
 - **of.id**: production order id
 - **new.of**: 1=first observation of the production order, 0=otherwise
 - **ref**: reference
 - **ics**: ideal cycle speed, units per minute
 - **TU**: total units
 - **DU**: defective units
 - **TgU**: target units, `OpT`$\times$`ics`
 - **OT**: period length, minutes
 - **SBT**: stand-by time, minutes
 - **LT**: loading time = `OT`-`SBT`
 - **rcs**: real cycle speed = `TU`/`LT`
 - **lo**: loading rate = `LT`/`OT`
 - **DT**: downtime, minutes
 - **OpT**: operating time = `LT`-`DT`
 - **av**: availability rate = `OpT`/`LT`
 - **PLT**: performance losses time, minutes
 - **NOpT**: net operating time = `OpT`-`PLT`
 - **pf**: performance rate = `NOpT`/`OpT`
 - **QLT**: quality losses time, minutes
 - **VT**: valuable time = `NOpT`-`QLT`
 - **qu**: quality rate = `VT`/`NOpT`
 - **oee**: `av` $\times$ `pf` $\times$ `qu`
 - **nstops**: number of stops
 - **hum**: % of humidity
 - **temp**: temperature in ºC
 - **wPT**: work shift time
 - **av.level**: availability level ([80-100\%]: very good, [60-80\%): acceptable, [40-60\%): improvable, [0-40\%): very poor)
 - **pf.level**: performance level
 - **qu.level**: quality level
 - **oee.level**: oee level

# MVTF
This folder contains all the code required to run the multivariate and univariate models and reproduce figures similar to Figures 5.1 and 5.2 of the paper.  
 - [functions](mvtf/functions.R): `R` script with some auxiliary functions that will be used to run the model, including functions for the clustering step (section 4.3 in the paper).
 - [update_model](mvtf/update_model.R): `R` code for parameter estimation and response prediction, including Algorithm 1 and Algorithm 2 of the paper.
 - [theme_mtf](mvtf/theme_mtf.R): customized theme for figures.
 - [run_mv_model](mvtf/run_mv_model.R): `R` code to run the multivariate version of the model using a subset of the whole dataset comprising 4 consecutive weeks of data. Responses $\mathbf y_n$, covariates $\mathbf x_n$ and classification variables $\mathbf t_n$ are chosen as stated in Section 5. Different model configurations can be tested by the user changing either of them. The performance of the prediction method is measured using a 4-fold cross-validation technique alternatively using one week as the test set, whereas the remaining three weeks are used to train the model.  
 - [run_uv_model](mvtf/run_uv_model.R): `R` code to run the univariate version of the multivariate code above.
 - [mvtf_figures](mvtf/mvtf_figures.R): `R` code to reproduce figures similar to Figures 5.1 and 5.2 in the paper.

To ensure the proper working of the code please run the scripts in the following order: [run_mv_model](mvtf/run_mv_model.R) - [run_uv_model](mvtf/run_uv_model.R) - [mvtf_figures](mvtf/mvtf_figures.R).

Computations were run using R-4.2.1.

 
# References

Miguelez, F., Doncel, J. and Ugarte, M.D. (2024). An Adaptive Learning Approach to Multivariate Time Forecasting in Industrial Processes. _Submitted_. (ArXiv: https://arxiv.org/abs/2403.07554)
 





