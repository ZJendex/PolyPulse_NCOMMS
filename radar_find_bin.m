%% Heart -- find the breath
fr = [0.5, 10];
heart_distance = 0.45; % meter
target_angle = -40;
target_bpm = data_BPM;
thresholdOnBpm = 10;
angleThresholdDegree = 20;
secondDerivativeFlag = 1;
[max_bin_var_idx_history, max_bin_var_idx, fs] = find_bin(imgs, heart_distance, cut_index(1, 1), cut_index(1, 2), ...
    'thresholdOnBpm', thresholdOnBpm, 'targetBpm', target_bpm, 'secondDerivative', secondDerivativeFlag, "filter_range", fr, "targetAngleDegree", target_angle, "angleThresholdDegree", angleThresholdDegree);

max_bin_var_idx_history_c = max_bin_var_idx_history(max_bin_var_idx_history(:, 5) ~= 0, :);
history_list = max_bin_var_idx_history_c;
N = 40;
if N > size(history_list, 1)
    history_list = sortrows(history_list, 2, 'descend');
else
    history_list = sortrows(history_list, 2, 'descend');
    history_list = history_list(1:N, :);
end
max_bin_var_idx = [history_list(1,3),history_list(1,4)];

found_bin = [history_list(1, 3), history_list(1, 4)]


% save heart data & log & find bin middle plot
data_name = sprintf("heart_%s", "example");
heart_bin_threshold = 10;
closest_index = 1;
% if find bin find at bottom, get all 21 bottom bins
if min(128,history_list(closest_index, 4)+heart_bin_threshold) == 128
    data = imgs(cut_index(1, 1):cut_index(1, 2),history_list(closest_index, 3), ...
    end-heart_bin_threshold*2:end);
else
    data = imgs(cut_index(1, 1):cut_index(1, 2),history_list(closest_index, 3), ...
    history_list(closest_index, 4)-heart_bin_threshold:min(128,history_list(closest_index, 4)+heart_bin_threshold));
end
save(fullfile("output", data_name), 'data');

xl = [40,60];
fig = figure(355);clf;
fig.WindowState = 'maximized';
max_bin_var_idx = found_bin;
subplot(3,1,1)
[dist, t] = plot_waveformOnBin(max_bin_var_idx, imgs, cut_index(1, 1), cut_index(1, 2), ...
    "flipped", 0, "phase", 1, "filter", 0, "filter_range", fr);
hold on
indices = cut_peak_time_scg - cut_index(1, 1)/fs;
for k = 1:length(indices)
    line([indices(k) indices(k)], ylim, 'Color', 'red', 'LineStyle', '--', 'LineWidth', 0.5);
end
hold off;
xlim(xl) 

subplot(3,1,2)
[dist, t] = plot_waveformOnBin(max_bin_var_idx, imgs, cut_index(1, 1), cut_index(1, 2), ...
    "flipped", 0, "phase", 1, "filter", 1, "filter_range", fr);
d_dist = computeSecondDerivative(dist, 1/500);
plot(t, d_dist)
title("Second Derivative of Heart Phase Waveform")
% show the gt lines
hold on
indices = cut_peak_time_scg - cut_index(1, 1)/fs;
for k = 1:length(indices)
    line([indices(k) indices(k)], ylim, 'Color', 'red', 'LineStyle', '--', 'LineWidth', 0.5);
end
hold off;
xlim(xl)

subplot(3,1,3)
plot(cut_time_axis_scg, cut_scg_data)
title("SCG Groundtruth")
hold on
indices = cut_peak_time_scg - cut_index(1, 1)/fs;
for k = 1:length(indices)
    line([indices(k) indices(k)], ylim, 'Color', 'red', 'LineStyle', '--', 'LineWidth', 0.5);
end
hold off;
xlim(xl)

%% Wrist -- find the heartbeat
target_bpm = data_BPM;
wrist_distance = 0.10; % meter
angleThresholdDegree = 20;
[max_bin_var_idx_history, max_bin_var_idx, fs] = find_bin(imgs, wrist_distance, cut_index(i, 1), cut_index(i, 2), ...
                                                        'angleThresholdDegree', 20, 'thresholdOnBpm', thresholdOnBpm, 'targetBpm', target_bpm);
max_bin_var_idx_history_c = max_bin_var_idx_history(max_bin_var_idx_history(:, 5) ~= 0, :);
most_frequent_lags = mode(max_bin_var_idx_history_c(:, 5));
threshold = 4; % the most likely range of CDF of bins candidates
sorted_max_bin_var_idx_history = sort(max_bin_var_idx_history_c(:, 5));
[most_frequent_lags, max_count] = most_appear_number(sorted_max_bin_var_idx_history, threshold);
n = size(max_bin_var_idx_history_c, 1);
% Retrieved the bins with close lag numbers
most_frequent_lags_bin = max_bin_var_idx_history_c(abs(max_bin_var_idx_history_c(:, 5) - most_frequent_lags) <= threshold, :);
% reset the LIST and max_bin_var_idx
history_list = most_frequent_lags_bin;

range = mode(history_list(1:min(10, size(history_list, 1)), 3));
target_list = history_list(1:min(10, size(history_list, 1)), 3:4);
target_list = target_list(target_list(:,1) == range, :);
angle_s = min(target_list(:,2));
angle_l = max(target_list(:,2));
% To be at same size
angle_m = floor(angle_s + (angle_l - angle_s)/2);
found_bin = [range, angle_m]

% save wrist data & log
data_name = sprintf("wrist_%s", "example");
wrist_bin_threshold = 10;

data = imgs(cut_index(1, 1):cut_index(1, 2),range,angle_m-wrist_bin_threshold:angle_m+wrist_bin_threshold);
save(fullfile("output", data_name), 'data');

fig = figure(356);clf;
fig.WindowState = 'maximized';
max_bin_var_idx = found_bin;
subplot(3,1,1)
[dist, t] = plot_waveformOnBin(max_bin_var_idx, imgs, cut_index(1, 1), cut_index(1, 2), ...
    "flipped", 1, "phase", 1, "filter", 0, "emd_i", [4, 6]);
hold on
indices = cut_peak_time_ppg - cut_index(1, 1)/fs;
for k = 1:length(indices)
    line([indices(k) indices(k)], ylim, 'Color', 'red', 'LineStyle', '--', 'LineWidth', 0.5);
end
hold off;
xlim(xl)

subplot(3,1,2)
plot(cut_time_axis_ppg ,cut_ppg_data)
title("PPG ground truth")
hold on
indices = cut_peak_time_ppg - cut_index(1, 1)/fs;
for k = 1:length(indices)
    line([indices(k) indices(k)], ylim, 'Color', 'red', 'LineStyle', '--', 'LineWidth', 0.5);
end
hold off;
xlim(xl)

subplot(3,1,3)
d_dist = computeFirstDerivative(dist, 1/500);
plot(t, d_dist)
hold on
indices = cut_peak_time_ppg - cut_index(1, 1)/fs;
for k = 1:length(indices)
    line([indices(k) indices(k)], ylim, 'Color', 'red', 'LineStyle', '--', 'LineWidth', 0.5);
end
hold off;
title("First Derivative of Wrist Phase Waveform")
xlim(xl)


%% Head
head_distance = 0.70; % meter
target_angle = -15;
target_bpm = data_BPM;
angleThresholdDegree = 20;
fr = [0.7, 4];
secondDerivativeFlag = 0;
[max_bin_var_idx_history, max_bin_var_idx, fs] = find_bin(imgs, head_distance, cut_index(i, 1), cut_index(i, 2), ...
                                                       'thresholdOnBpm', thresholdOnBpm, 'targetBpm', target_bpm,  ...
                                                       'targetAngleDegree', target_angle, 'angleThresholdDegree', angleThresholdDegree, ...
                                                       "filter_range", fr, 'secondDerivative', secondDerivativeFlag);
max_bin_var_idx_history_c = max_bin_var_idx_history(max_bin_var_idx_history(:, 5) ~= 0, :);
if isempty(max_bin_var_idx_history_c)
    secondDerivativeFlag = 1;
    [max_bin_var_idx_history, max_bin_var_idx, fs] = find_bin(imgs, head_distance, cut_index(i, 1), cut_index(i, 2), ...
                                                       'thresholdOnBpm', thresholdOnBpm, 'targetBpm', target_bpm,  ...
                                                       'targetAngleDegree', target_angle, 'angleThresholdDegree', angleThresholdDegree, ...
                                                       "filter_range", fr, 'secondDerivative', secondDerivativeFlag);
end
% Clean the list with no zero value rows
max_bin_var_idx_history_c = max_bin_var_idx_history(max_bin_var_idx_history(:, 5) ~= 0, :);
% reset the LIST and max_bin_var_idx
history_list = max_bin_var_idx_history_c;
range_s = min(history_list(:, 3));
range_l = max(history_list(:, 3));
angle_s = min(history_list(:, 4));
angle_l = max(history_list(:, 4));
range_m = floor(range_s + (range_l - range_s)/2);
angle_m = floor(angle_s + (angle_l - angle_s)/2);
found_bin = [range_m, angle_m]


% save head data & log
data_name = sprintf("head_%s", "example");
head_angleBin_threshold = 10;
head_rangeBin_threshold = 2;

data = imgs(cut_index(1, 1):cut_index(1, 2), ...
    range_m-head_rangeBin_threshold:range_m+head_rangeBin_threshold,angle_m-head_angleBin_threshold:angle_m+head_angleBin_threshold);
save(fullfile("output", data_name), 'data');

fig = figure(357);clf;
fig.WindowState = 'maximized';
max_bin_var_idx = found_bin;
subplot(3,1,1)
[dist, t] = plot_waveformOnBin(max_bin_var_idx, imgs, cut_index(1, 1), cut_index(1, 2), ...
    "flipped", 1, "phase", 1, "filter", 1, "filter_range", fr);
hold on
indices = cut_peak_time_bcg - cut_index(1, 1)/fs;
for k = 1:length(indices)
    line([indices(k) indices(k)], ylim, 'Color', 'red', 'LineStyle', '--', 'LineWidth', 0.5);
end
hold off;
xlim(xl) 

subplot(3,1,2)
[dist, t] = plot_waveformOnBin(max_bin_var_idx, imgs, cut_index(1, 1), cut_index(1, 2), ...
    "flipped", 0, "phase", 1, "filter", 1, "filter_range", fr);
d_dist = computeSecondDerivative(dist, 1/500);
plot(t, d_dist)
title("Second Derivative of Head Phase Waveform")
% show the gt lines
hold on
indices = cut_peak_time_bcg - cut_index(1, 1)/fs;
for k = 1:length(indices)
    line([indices(k) indices(k)], ylim, 'Color', 'red', 'LineStyle', '--', 'LineWidth', 0.5);
end
hold off;
xlim(xl)

subplot(3,1,3)
plot(cut_time_axis_bcg ,cut_bcg_data)
hold on
indices = cut_peak_time_bcg - cut_index(1, 1)/fs;
for k = 1:length(indices)
    line([indices(k) indices(k)], ylim, 'Color', 'red', 'LineStyle', '--', 'LineWidth', 0.5);
end
hold off;
title("BCG ground truth")
xlim(xl)


%% neck
fr = [0.5, 10];
neck_distance = 0.60; % meter
target_angle = -20;
target_bpm = data_BPM;
thresholdOnBpm = 10;
angleThresholdDegree = 20;
secondDerivativeFlag = 1;
[max_bin_var_idx_history, max_bin_var_idx, fs] = find_bin(imgs, neck_distance, cut_index(1, 1), cut_index(1, 2), ...
    'thresholdOnBpm', thresholdOnBpm, 'targetBpm', target_bpm, 'secondDerivative', secondDerivativeFlag, "filter_range", fr, "targetAngleDegree", target_angle, "angleThresholdDegree", angleThresholdDegree);

max_bin_var_idx_history_c = max_bin_var_idx_history(max_bin_var_idx_history(:, 5) ~= 0, :);
history_list = max_bin_var_idx_history_c;
N = 40;
if N > size(history_list, 1)
    history_list = sortrows(history_list, 2, 'descend');
else
    history_list = sortrows(history_list, 2, 'descend');
    history_list = history_list(1:N, :);
end
max_bin_var_idx = [history_list(1,3),history_list(1,4)];

found_bin = [history_list(1, 3), history_list(1, 4)]

% save neck data & log & find bin middle plot
data_name = sprintf("neck_%s", "example");
neck_bin_threshold = 10;
% The first min number
closest_index = 1;
% if find bin find at bottom, get all 21 bottom bins
if min(128,history_list(closest_index, 4)+neck_bin_threshold) == 128
    data = imgs(cut_index(1, 1):cut_index(1, 2),history_list(closest_index, 3), ...
    end-neck_bin_threshold*2:end);
else
    data = imgs(cut_index(1, 1):cut_index(1, 2),history_list(closest_index, 3), ...
    history_list(closest_index, 4)-neck_bin_threshold:min(128,history_list(closest_index, 4)+neck_bin_threshold));
end
save(fullfile("output", data_name), 'data');

fig = figure(358);clf;
fig.WindowState = 'maximized';
max_bin_var_idx = found_bin;
subplot(3,1,1)
[dist, t] = plot_waveformOnBin(max_bin_var_idx, imgs, cut_index(1, 1), cut_index(1, 2), ...
    "flipped", 0, "phase", 1, "filter", 0, "filter_range", fr);
hold on
indices = cut_peak_time_neck - cut_index(1, 1)/fs;
for k = 1:length(indices)
    line([indices(k) indices(k)], ylim, 'Color', 'red', 'LineStyle', '--', 'LineWidth', 0.5);
end
hold off;
xlim(xl) 

subplot(3,1,2)
[dist, t] = plot_waveformOnBin(max_bin_var_idx, imgs, cut_index(1, 1), cut_index(1, 2), ...
    "flipped", 0, "phase", 1, "filter", 1, "filter_range", fr);
d_dist = computeSecondDerivative(dist, 1/500);
plot(t, d_dist)
title("Second Derivative of neck Phase Waveform")
% show the gt lines
hold on
indices = cut_peak_time_neck - cut_index(1, 1)/fs;
for k = 1:length(indices)
    line([indices(k) indices(k)], ylim, 'Color', 'red', 'LineStyle', '--', 'LineWidth', 0.5);
end
hold off;
xlim(xl)

subplot(3,1,3)
plot(cut_time_axis_neck ,cut_neck_data)
hold on
indices = cut_peak_time_neck - cut_index(1, 1)/fs;
for k = 1:length(indices)
    line([indices(k) indices(k)], ylim, 'Color', 'red', 'LineStyle', '--', 'LineWidth', 0.5);
end
hold off;
title("Neck ground truth")
xlim(xl)
