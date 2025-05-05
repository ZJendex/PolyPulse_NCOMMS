function [t_i,imu,outlier_idx] = outlinerRemoval4DataPairs(t_i,imu)
num_theshold = 100;
% --- Calculate Differences Between Neighbors ---
% For t_i
diff_before_t = [0, abs(diff(t_i))]; % Difference with the previous point
diff_after_t = [abs(diff(t_i)), 0];  % Difference with the next point

% For imu
diff_before_imu = [0, abs(diff(imu))]; % Difference with the previous point
diff_after_imu = [abs(diff(imu)), 0];  % Difference with the next point

% --- Calculate Standard Deviations ---
std_diff_t = std(abs(diff(t_i)));
std_diff_imu = std(abs(diff(imu)));

% --- Define Outlier Thresholds ---
threshold_t = num_theshold * std_diff_t; % Threshold for t_i
threshold_imu = num_theshold * std_diff_imu; % Threshold for imu

% --- Identify Outliers ---
outlier_idx_t = find((diff_before_t > threshold_t) | (diff_after_t > threshold_t));
outlier_idx_imu = find((diff_before_imu > threshold_imu) | (diff_after_imu > threshold_imu));

% Combine outlier indices from both datasets
outlier_idx = unique([outlier_idx_t, outlier_idx_imu]);

% --- Remove Outliers ---
t_i(outlier_idx) = [];
imu(outlier_idx) = [];

end

