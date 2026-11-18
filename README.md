# That's Odd! - Spot the outliers

**Hollandse Luchten Hackathon 2026 · Challenge 1 "Spot the outliers"**
Team: Piet Hoogkamer,Aliaksandr Kuznetsov, Shirin Janghorban
Submission: `thats_odd.R`

---

## Introduction

Low-cost PM2.5 sensors produce a lot of values that look strange. Some are sensor errors (frozen readings, impossible peaks), but some are real pollution, such as fireworks, construction or a nearby fire. Our script tells these apart in two stages. First it finds weird sensors over a 3-month window: sensors with too little data, frozen sensors, sensors that don't follow their neighbours, sensors that read structurally too high or too low, and sensors that are outliers in a PCA. Then it finds weird hours in the remaining sensors. For each 2-week sliding window, three independent detectors vote, and each detector judges its residual against *local* limits. Apart from a few hard rules (frozen, out of range, impossible value), a value is only called an error when the comparison with nearby sensors supports it. Every weird sensor gets a one-line reason in `step4_weird_sensors.csv`. Every hour gets a cause label, a score, a traffic-light status and, when flagged, a confidence level in `step7_sensor_hours.csv.gz`. We also test the method itself: we plant known faults in sensors that passed step 4 and measure what share of them we find.

---

## Quick start

```r
install.packages(c("data.table", "cluster", "arrow"))   # MASS ships with R
source("thats_odd.R")
main()
```

Or from the command line:

```bash
Rscript thats_odd.R
```

Results are written to `output/`. The R and package versions used are saved to `output/session_info.txt` on every run.

`arrow` (or `nanoparquet` as an alternative) is only needed for `.parquet` input. It is called through `requireNamespace()`, so it does not show as attached in RStudio's Packages tab.

### Input files

Put these in a folder `data/` next to the script. Each file may be `.parquet`, `.csv.gz` or `.csv`; the script uses whichever it finds.

| File | Columns used |
|---|---|
| `data_sensors` | `kit_id, timestamp_from, pm25_kal, pm25_raw, rh, temp, latitude, longitude` |
| `metadata_sensors` | `kit_id, pm25_sensor` (sensor type) |
| `reference_measurements` | `meetlocatie_id, component, begindatumtijd, waarde` (PM2.5 rows) |
| `reference_station_metadata.csv` | `meetlocatie_id, breedtegraad, lengtegraad` |
| `knmi_weather_data` | `time, wsi, lat, lon, T, U, DD, FH` |


---
The data may be obtained from this Hollandse Luchten repository: **[Hollandse Luchten Hackathon 2026](https://github.com/provincieNH/HackathonHollandseLuchten/tree/main) · Challenge 1 "Spot the outliers"**
## The workflow


```
Step 1  Load & clean ──► Step 2  Pick area ──► Step 3  3-month window
                                                        │
              ┌─────────────────────────────────────────┘
              ▼
Step 4  WEIRD SENSORS (whole window)
        a. > 25 % missing          → "too little data"
        b. PCA per sensor type     → far from its group
        c. frozen / low correlation / structurally too high or low
        d. map   e. table with reasons
              │  kept sensors only
              ▼
Step 6  Slide a 2-week window, 1 week at a time
              ▼
Step 7  WEIRD HOURS (inside each window)
        a. STL of the sensor's own series      ─┐
        b. local PCA with nearby sensors        ├─ each judged with local limits → vote
        c. difference with median of neighbours ─┘
           + hard rules (frozen, out of range, impossible, drift)
           + context (area-wide event, moisture)
        d. time-series plots with a symbol per cause   e. table with reasons
              ▼
Step 5  Detail plots for one selected sensor        Extra: validation with planted faults
Step 8  Change the settings and run again
```

### Step 1 — Load and clean
- **One clock (UTC).** Sensor timestamps are UTC. Official stations are reported in UTC+1, so 1 hour is subtracted. KNMI timestamps mark the *end* of the hour, so 1 hour is subtracted there too.
- GPS fixes near (0, 0) are set to missing. Each sensor gets one location: the median of its fixes. The share of hours in which a sensor was more than 200 m from that spot is computed as `moved_share`. It is not yet used in any rule and is not written to an output file.
- Impossible humidity (≤ 0 or > 100 %) becomes `NA`. Duplicate sensor-hours are removed and counted.
- Data are put in hour × sensor matrices. **Missing stays `NA`, never 0**, and negative values are kept.

### Step 2 — Areas
The script may use an alternative way of clustering the sensors clusters for more freeom. Simply set the `CFG$cluster_file` file to your way of clustering. Otherwise it runs k-medoids (PAM, k = 8) on sensor location in km. Each cluster is named after the town where *most* of its sensors are, not the town nearest its centre.

### Step 3 — Larger window
`2025-10-01` to `2026-01-01` (3 months, UTC). The sensor archive ends 2025-12-31 23:00.

### Step 4 — Weird sensors
| Rule | Limit | Why |
|---|---|---|
| 4a Too little data | coverage < 75 % | Reported separately: these sensors aren't odd, they just have no data |
| 4c Frozen | same raw value ≥ 6 h in a row, in ≥ 10 % of hours (~9 days) | A stuck sensor is malfunctioning |
| 4c Doesn't follow the others | correlation with median of the *other* sensors < 0.60 (log scale) | It is not measuring the same air |
| 4c Structurally off | median ratio to the others ≥ 1.75× or ≤ 1/1.75× | Calibration or placement problem |
| 4b PCA outlier | distance in the PC1/PC2 plane to the median of the *other* sensors > 3× the typical distance | Separate PCA per sensor type (≥ 5 sensors), because the two types form two groups |

Each weird sensor gets one reason: the first rule that fires, in the order coverage → frozen → correlation → ratio → PCA. The correlation, ratio and PCA checks need at least 4 sensors left after 4a and the frozen check.

Only kept sensors go on to step 7, and only kept sensors act as neighbours, so one broken sensor cannot make its neighbours look wrong.

### Steps 6 & 7 — Weird hours
For every kept sensor, in every 2-week window:

1. **Expected value.** This is the median of the nearby sensors at the same hour. "Nearby" means within 1 km, widened in 1 km steps up to 3 km until at least 3 neighbours are found, and the whole area if there still aren't enough. When fewer than 3 neighbours have data in a given hour, the nearest official station within 25 km is used instead.
2. **Three detectors**, each producing a residual:
   - **7a STL.** The sensor's own series is split into trend + 24-hour pattern + remainder. Gaps are filled by interpolation for the fit, but filled hours are never judged. STL is skipped if more than 25 % of the window is missing.
   - **7b Local PCA.** The sensor and its neighbours are standardised (median/MAD), PC1 is removed (what they all do together), and what remains for the sensor is its own deviation.
   - **7c Neighbour difference.** The sensor minus the expected value.
3. **Local limits.** Each residual is compared with its own moving median ± 4 robust SD (1.4826 × MAD) over a 13-hour window. The limits widen when the data are restless and narrow when they are calm. A flag also requires a difference of at least 5 µg/m³, so tiny differences at very low concentrations are ignored.
4. **Log scale.** Residuals are computed on `log(1 + x)` by default, so "2× too high" counts the same at 5 and at 50 µg/m³.
5. **Hard rules:** frozen (≥ 6 h identical raw value), out of range (< 0 or > 200 µg/m³), impossible (> 100 µg/m³ *and* > 3× expected), drift (48-hour rolling median of the ratio to the neighbours ≥ 1.75× and ≥ 4 µg/m³).
6. **Context:**
   - **Area-wide event.** If the area median is ≥ 2× its window median and ≥ 20 µg/m³, the pollution is real (fireworks, smog). During such an event only frozen and out-of-range values still count as errors.
   - **Moisture.** A sensor reading too high while sensor RH ≥ 90 % or KNMI RH ≥ 95 %. The sensor's own RH is ignored when its temperature is more than 10 °C off KNMI, because then its RH is wrong too.
7. **Sliding windows.** Each hour falls in two windows (up to three near the end of the period, where an extra window is added so the last one ends on the last hour). The final verdict comes from the window in which the hour is most central, which reduces edge effects. `windows_flagged / windows` shows how stable the verdict is.

### Scoring and colours

| Situation | Score |
|---|---|
| Supported by neighbours (7b or 7c) | 30 + 10·(z − 4) + 10·(votes − 1), max 85 |
| Only unusual for itself (STL only) | 20, information rather than an error |
| Drift | 45, or 60 if ≥ 2.5× |
| Frozen / out of range / impossible | 95 |
| **No hard evidence** | capped at 65, so it can never turn red |

"Hard evidence" means a hard rule fired, or the value is ≥ 20 µg/m³ **and** ≥ 2.5× off its neighbours, or there is drift of ≥ 2.5×.

| Colour | Meaning |
|---|---|
| 🔴 red (score ≥ 70) | certainly wrong |
| 🟠 yellow (score ≥ 30) | doubtful |
| 🟢 green | ok |
| 🔵 blue | area-wide event: real pollution, not an error |
| ⚪ grey | missing |

**Confidence:** High when all 3 detectors agree or a hard rule fired, Medium for 2, Low for 1.

**Causes**, in order of priority: missing, frozen, out of range, impossible value, moisture, drift, sudden drop, spike, unusual for this sensor only, area-wide event, ok.

---

## Validation (planted faults)

This lives in `plant_faults()` and `validate()` and is called from `main()`. It runs only for the areas set in `CFG$train_area` (Haarlem) and `CFG$test_area` (Zaandam). "Train" and "test" are labels in the output: nothing is fitted or tuned automatically on the train area. In each of these areas the script plants known faults in at most 1 in 4 kept sensors (max 8). Keeping it to a quarter means the neighbours stay honest. Per sensor it plants:

- 4 spikes (+25 to +60 µg/m³)
- 2 impossible values (300–900 µg/m³)
- 12 frozen hours
- 4 days of ×2.2 drift

It then reruns steps 6–7 for every combination of **limit 3 or 4 SD** × **log or original scale**. For each detector, for the rules alone and for the final ensemble, it reports the detection rate per fault type and the false-alarm rate. Hours within 24 h of a planted fault are not counted as false alarms. The false-alarm rate is measured on the clean hours of the **planted sensors only**, not on every sensor in the area. It also reports the **false RED** rate and compares against a **naive mean ± 3 SD** baseline. The results are in `output/validation.csv`.

---

## Outputs (`output/`)

| File | What it is |
|---|---|
| `step4_weird_sensors.csv` | Every sensor: group (kept / weird / too little data), reason, coverage, frozen share, correlation, slope, ratio, PCA coordinates and distance |
| `step4b_pca_<area>.png` | PCA plot per sensor type; triangles are weird sensors, the cross is the group median |
| `step4c_cross_<area>.png` + `_correlation.csv`, `_slope.csv` | Sensor-by-sensor correlation and slope matrices |
| `step4_scatter_<id>.png` | Each weird sensor against the average of the kept sensors, with a robust (`MASS::rlm`) fit line |
| `step4d_map.png` | Map that displays points on a longitude and latitude plane: kept (green), weird (red triangle), too little data (grey) |
| `step6_all_windows_<area>.csv.gz` | Verdict for every hour in every window |
| `step7_sensor_hours.csv.gz` | Final verdict per sensor-hour: value, expected, z-scores, votes, cause, score, confidence, status, weather |
| `step7_episodes.csv` | Consecutive flagged hours grouped into episodes: start, end, cause, peak, expected, agreement |
| `step7_overview_<area>.png` | Whole period for the most-flagged sensor in each area, with a symbol per cause |
| `step5_selected_<id>_overview.png` | Same view for the selected sensor |
| `step7_residuals_<id>_window<n>.png` | The three detector residuals with their local limits |
| `step7a_stl_decomposition_<id>_window<n>.png` | STL breakdown: measured, trend, daily pattern, remainder |
| `step7a_log_vs_original_<id>_window<n>.png` | Effect of the log transform on what gets flagged |
| `report_card.csv` | One row per sensor: share green / yellow / red, unstable hours, main issue, mounting note (temperature often > 10 °C off KNMI) |
| `validation.csv` | Planted-fault detection and false-alarm rates |
| `session_info.txt` | R and package versions |

---

## Settings (Step 8: change and run again)

Every number in the analysis lives in the `CFG` list at the top of the script, each with a comment saying why it was chosen. The settings you are most likely to change:

| Setting | Default | |
|---|---|---|
| `areas` | `c("Haarlem", "Zaandam")` | `NULL` analyses all clusters |
| `cluster_file` | `NA` | Path to alternative clustering file; `NA` uses k-medoids |
| `period_start` / `period_end` | `2025-10-01` / `2026-01-01` | The larger window (end date is exclusive) |
| `selected_sensor` / `selected_window` | `NA` | `NA` picks the most-flagged sensor and window |
| `window_days` / `step_days` | 14 / 7 | Set `step_days = 1` to slide day by day |
| `limit_sd` | 4 | Local limit in robust SD |
| `use_log` | `TRUE` | Residuals on the log or the original scale |
| `nb_radius_km` / `nb_radius_max` / `nb_min` | 1 / 3 / 3 | Neighbour search |

---


## Usage of AI
The team mostly brainstormed the methodology together. The majority of the code used in here is made using AI, but sketched using code instructions. The code has been checked for validity by the members.
