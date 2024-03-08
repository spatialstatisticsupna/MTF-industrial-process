# MTF-industrial-process
Data and R code to reproduce the results described in the paper "An Adaptive Learning Approach to Multivariate Time Forecasting in Industrial Processes" (Miguelez et al., 2024).

# Data
This folder contains 3 datasets:
 1) `oee_data.Rds`: complete dataset from 2019-12-16 to 2021-08-25, including holidays, weekends and weeks with missing data
 2) `oee_data_complete_weeks.Rds`: subset of `oee_data.Rds` containing complete weeks without holidays and weekends
 3) `oee_data_4weeks.Rsd`: subset of `oee_data_complete_weeks.Rds` containing data from week 21 to week 24 of the year 2020

All of the datasets contain the following variables:
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
 - **DT**: down time, minutes
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
This folder contains all the code required to run the multivariate and univariate models and reproduce Figures 5.1 and 5.2 of the paper.
 - [functions](mvtf/functions.R): `R` script with some auxiliary functions that will be used to run the model
 - [update_model](mvtf/update_model.R): `R` code for parameter estimation and response prediction
 - [theme_mtf](mvtf/theme_mtf.R): customized theme for figures
 - [run_mv_model](mvtf/run_mv_model.R): `R` code to run the multivariate version of the model using a subset of the whole dataset comprising 4 consecutive weeks of data. Responses $\mathbf y_n$  
 - [run_uv_model](mvtf/run_uv_model.R)
 - [mvtf_figures](mvtf/mvtf_figures.R)

Computations were run using R-4.2.1.
