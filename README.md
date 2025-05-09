# PolyPulse_NCOMMS -- Signal Processing Pipeline
This repository contains example signal processing outputs from the **Polypulse** project.

## Data Description

The saved data falls into two broad categories:

### 1. **Wearable Sensor Ground Truth**
- **Signals**: Heart-IMU, Head-IMU, wrist-PPG and neck-IMU vibrations.
- **Files**:
  - `GT_example.mat`: Contains peak annotations (e.g., `scg_peaks_gt`, `ppg_peaks_gt`, etc.) extracted after calibration and alignment using `wearable_sensor_sync.m` illustrate in the plot `plot_GT_example.png`.
  - These peaks serve as reference points (in red dashed lines) for evaluating radar measurements.

### 2. **Radar Bin-Selected Signals**
- **Procedure**: After syncronized the radar with wearable sensors by `radar_sensor_sync.m`, `radar_find_bin.m` identifies bins (range-angle pairs) that best match cardiac or vascular pulsations at different body locations (e.g., heart, wrist, head, neck).
- **Saved files**:
  - `heart_example.mat`, `wrist_example.mat`, `head_example.mat`, `neck_example.mat`: Contain radar identified cardiac-related bins.
  - Corresponding figures: `plot_Heart_example.png`, etc., which visualize the radar signal, its derivative, and synchronized wearable signals, with red dotted lines marking cardiac arrival times.


