# 📡 GNSS Time-Series Analysis

This repository contains Python scripts to process GNSS time-series data using Least Squares estimation.

## 📑 Description

Based on an exercise for the Space Geodesy and InSAR course (Prof. Alessandra Borghi, ICTP, Trieste).
The script loads GNSS station data files (tenv) downloaded from https://geodesy.unr.edu/NGLStationPages/GlobalStationList. 
The functions applies a linear trend with periodic (annual & semi-annual) signals,
handles known discontinuities, and plots residuals.

## 📚 Mathematical Model

This project models GNSS time-series data using a deterministic function that accounts for:
- A linear trend (tectonic velocity)
- Annual and semi-annual periodic signals (seasonal effects)
- Known discontinuities (e.g. receiver or antenna changes)

The general model for each displacement component (*East*, *North*, *Up*) is:


Y(ti) = A + B * ti
+ C * sin(ω * ti) + D * cos(ω * ti)
+ E * sin(2ω * ti) + F * cos(2ω * ti)
+ Σ Gk * H(ti - Tk)

**Where:**
- *Y(ti)*: displacement at time *ti*
- *A*: offset
- *B*: linear velocity (mm/year)
- *C*, *D*: annual amplitude coefficients
- *E*, *F*: semi-annual amplitude coefficients
- *Gk*: coefficient for known discontinuities at time *Tk*
- *H()*: Heaviside step function
- *ω = 2π* (annual frequency in decimal years)

🧮 Least Squares Estimation

The model parameters (*A, B, C, D, E, F, Gk*) are estimated using the **Least Squares Method (LSQ)**:

X = (Aᵀ Q⁻¹ A)⁻¹ Aᵀ Q⁻¹ Yo


**Where:**
- *A*: design matrix
- *Q*: cofactor matrix (diagonal, assuming white noise)
- *Yo*: observed displacements

---

## 🔎 Goal

By removing the estimated linear and periodic trends, residuals reveal signals not explained by the model (e.g., local effects, remaining noise, or unmodelled seasonal variations).

---

📖 *Based on the exercise for the Space Geodesy and InSAR course (Prof. Alessandra Borghi).*  



## 📂 Structure

- `scripts/` : Python scripts for LSQ fitting.
- `data/` : GNSS time-series data files.
- `figs/` : Generated plots.

## 🚀 How to Run

1. Place your `.tenv` files in `data/`.
2. Edit `input_file_name` in the script 'GNSS_analysis_multiple_station.py' or 'GNSS_analysis_single_station.py'.
3. Activate a deticated conda enviroment 
4. Run the script:
    ```bash
    python scripts/gnss_timeseries_lsq.py
    ```
5. Check the output in the directory 'figs'.

## 📜 License

MIT License (add `LICENSE` file if you want).

---

Created by Dr. Luigi Sante Zampa.