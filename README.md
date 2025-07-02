# 📡 GNSS Time-Series Analysis

This repository contains Python scripts to process GNSS time-series data using Least Squares estimation.

## 📑 Description

Based on an exercise for the Space Geodesy and InSAR course (Prof. Alessandra Borghi).
The script loads GNSS station data, applies a linear trend with periodic (annual & semi-annual) signals,
handles known discontinuities, and plots residuals.

## 📂 Structure

- `scripts/` : Python scripts for LSQ fitting.
- `data/` : GNSS time-series data files.
- `figs/` : Generated plots.

## 🚀 How to Run

1. Place your `.tenv` files in `data/`.
2. Edit `input_file_name` in the script.
3. Run the script:
    ```bash
    python scripts/gnss_timeseries_lsq.py
    ```

## 📜 License

MIT License (add `LICENSE` file if you want).

---

Created by Luigi Sante Zampa.