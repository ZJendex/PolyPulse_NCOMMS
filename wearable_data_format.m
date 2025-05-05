%% ground truth
cut_index(1,1) = cali_time*fs;
cut_index(1,2) = end_cutting_time*fs;
i = 1;

cut_time_axis_scg = time_axis_scg_s(time_axis_scg_s > cut_index(i,1)/fs);
cut_time_axis_scg = cut_time_axis_scg(cut_time_axis_scg < cut_index(i,2)/fs);
cut_time_axis_scg = cut_time_axis_scg - cut_time_axis_scg(1);
cut_scg_data = scg_data(time_axis_scg_s > cut_index(i,1)/fs);
cut_scg_data = cut_scg_data(time_axis_scg_s(time_axis_scg_s > cut_index(i,1)/fs) < cut_index(i,2)/fs);
cut_peak_time_scg = peak_time_scg(peak_time_scg > cut_index(i, 1)/fs);
cut_peak_time_scg = cut_peak_time_scg(cut_peak_time_scg < cut_index(i,2)/fs);

cut_time_axis_ppg = time_axis_ppg_s(time_axis_ppg_s > cut_index(i,1)/fs);
cut_time_axis_ppg = cut_time_axis_ppg(cut_time_axis_ppg < cut_index(i,2)/fs);
cut_ppg_data = ppg_data(time_axis_ppg_s > cut_index(i,1)/fs);
cut_ppg_data = cut_ppg_data(time_axis_ppg_s(time_axis_ppg_s > cut_index(i,1)/fs) < cut_index(i,2)/fs);
cut_time_axis_ppg = cut_time_axis_ppg - cut_time_axis_ppg(1);
cut_peak_time_ppg = peak_time_ppg(peak_time_ppg > cut_index(i, 1)/fs);
cut_peak_time_ppg = cut_peak_time_ppg(cut_peak_time_ppg < cut_index(i,2)/fs);

cut_time_axis_bcg = time_axis_bcg_s(time_axis_bcg_s > cut_index(i,1)/fs);
cut_time_axis_bcg = cut_time_axis_bcg(cut_time_axis_bcg < cut_index(i,2)/fs);
cut_bcg_data = bcg_data(time_axis_bcg_s > cut_index(i,1)/fs);
cut_bcg_data = cut_bcg_data(time_axis_bcg_s(time_axis_bcg_s > cut_index(i,1)/fs) < cut_index(i,2)/fs);
cut_time_axis_bcg = cut_time_axis_bcg - cut_time_axis_bcg(1);
cut_peak_time_bcg = peak_time_bcg(peak_time_bcg > cut_index(i, 1)/fs);
cut_peak_time_bcg = cut_peak_time_bcg(cut_peak_time_bcg < cut_index(i,2)/fs);

cut_time_axis_neck = time_axis_neck_s(time_axis_neck_s > cut_index(i,1)/fs);
cut_time_axis_neck = cut_time_axis_neck(cut_time_axis_neck < cut_index(i,2)/fs);
cut_time_axis_neck = cut_time_axis_neck - cut_time_axis_neck(1);
cut_neck_data = neck_data(time_axis_neck_s > cut_index(i,1)/fs);
cut_neck_data = cut_neck_data(time_axis_neck_s(time_axis_neck_s > cut_index(i,1)/fs) < cut_index(i,2)/fs);
cut_peak_time_neck = peak_time_neck(peak_time_neck > cut_index(i, 1)/fs);
cut_peak_time_neck = cut_peak_time_neck(cut_peak_time_neck < cut_index(i,2)/fs);

% save the data
data_name = sprintf("GT_%s", "example");
full_file_path = fullfile("output", data_name);
scg_peaks_gt = round((cut_peak_time_scg-cali_time) * 500);
ppg_peaks_gt = round((cut_peak_time_ppg-cali_time) * 500);
bcg_peaks_gt = round((cut_peak_time_bcg-cali_time) * 500);
neck_peaks_gt = round((cut_peak_time_neck-cali_time) * 500);
% Save all the variables in one .mat file
save(full_file_path, "scg_peaks_gt", "ppg_peaks_gt", "bcg_peaks_gt", "neck_peaks_gt");