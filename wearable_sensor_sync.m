cali_time = 90; %second
end_cutting_time = 170;
fs = 500;
addpath("utils\")

%% 
load("data\wearable_data.mat")

imu1 = detrend(y);
imu2 = detrend(y2);
imu3 = detrend(y3);
imu4 = detrend(y4);
t_i1 = double(time_stamp)/1000;
t_i2 = double(time_stamp2)/1000;
t_i3 = double(time_stamp3)/1000;
t_i4 = double(time_stamp4)/1000;
% drop potential outlier
imu1 = imu1(5:end-1);
imu2 = imu2(5:end-1);
imu3 = imu3(5:end-1);
imu4 = imu4(5:end-1);
t_i1 = t_i1(5:end-1);
t_i2 = t_i2(5:end-1);
t_i3 = t_i3(5:end-1);
t_i4 = t_i4(5:end-1);
t_i1 = t_i1-t_i1(1);
t_i2 = t_i2-t_i2(1);
t_i3 = t_i3-t_i3(1);
t_i4 = t_i4-t_i4(1);

% eliminate sensor error reading
[t_i1,imu1, o1] = outlinerRemoval4DataPairs(t_i1,imu1);
[t_i2,imu2, o2] = outlinerRemoval4DataPairs(t_i2,imu2);
[t_i3,imu3, o3] = outlinerRemoval4DataPairs(t_i3,imu3);
[t_i4,imu4, o4] = outlinerRemoval4DataPairs(t_i4,imu4);

t_i1 = t_i1-t_i1(1);
t_i2 = t_i2-t_i2(1);
t_i3 = t_i3-t_i3(1);
t_i4 = t_i4-t_i4(1);

figure(111)
plot(t_i1, imu1)
hold on 
plot(t_i2, imu2)
plot(t_i3, imu3)
plot(t_i4, imu4)
legend()
hold off

% find shift_bins -- sliding window to find the smallest corr_shift
step = 0.2; % second - should smaller than 1s
w_len = 1.5; % second
p_shift_bins12 = zeros(cali_time/step+1, 3);
p_shift_bins13 = zeros(cali_time/step+1, 3);
p_shift_bins14 = zeros(cali_time/step+1, 3);
k = 1;
for w_start = 0:step:cali_time
    xl = [w_start, w_start+w_len];
    cut_imu1 = imu1(t_i1 > xl(1) & t_i1 < xl(2));
    cut_t_i1 = t_i1(t_i1 > xl(1) & t_i1 < xl(2));
    cut_imu2 = imu2(t_i2 > xl(1) & t_i2 < xl(2));
    cut_t_i2 = t_i2(t_i2 > xl(1) & t_i2 < xl(2));
    cut_imu3 = imu3(t_i3 > xl(1) & t_i3 < xl(2));
    cut_t_i3 = t_i3(t_i3 > xl(1) & t_i3 < xl(2));
    cut_imu4 = imu4(t_i4 > xl(1) & t_i4 < xl(2));
    cut_t_i4 = t_i4(t_i4 > xl(1) & t_i4 < xl(2));
    % normalzied & interpolate the data
    % only if there are weird timestamp jump
    % smooth the time stamp
    updatedArray = cut_t_i1;
    for i = 2:length(updatedArray) % Loop through the array starting from the second element
        diffValue = updatedArray(i) - updatedArray(i-1); % Compute difference
        if diffValue < 0 || diffValue > 0.01 % Check if out of range
            updatedArray(i) = updatedArray(i-1) + 0.002; % Adjust the value
        end
    end
    cut_t_i1 = updatedArray;
    [cut_n_imu_i1,~] = interpolationWithTimestamp(cut_imu1,cut_t_i1,fs);
    [cut_n_imu_i2,~] = interpolationWithTimestamp(cut_imu2,cut_t_i2,fs);
    [cut_n_imu_i3,~] = interpolationWithTimestamp(cut_imu3,cut_t_i3,fs);
    [cut_n_imu_i4,~] = interpolationWithTimestamp(cut_imu4,cut_t_i4,fs);
    cut_n_imu_i1 = detrend(cut_n_imu_i1);
    cut_n_imu_i2 = detrend(cut_n_imu_i2);
    cut_n_imu_i3 = detrend(cut_n_imu_i3);
    cut_n_imu_i4 = detrend(cut_n_imu_i4);
    % make them equal length
    [s_len, ~] = min([length(cut_n_imu_i1), length(cut_n_imu_i2), length(cut_n_imu_i3), length(cut_n_imu_i4)]);
    cut_n_imu_i1 = cut_n_imu_i1(1:s_len);
    cut_n_imu_i2 = cut_n_imu_i2(1:s_len);
    cut_n_imu_i3 = cut_n_imu_i3(1:s_len);
    cut_n_imu_i4 = cut_n_imu_i4(1:s_len);

    [acf12, lags12] = xcorr(cut_n_imu_i1, cut_n_imu_i2);
    [val12, index12] = findpeaks(acf12);
    if ~isempty(val12)
        [vm12, im12] = max(val12);
        p_shift_bins12(k,1) = lags12(index12(im12));
        p_shift_bins12(k,2) = max(acf12);
    else
        p_shift_bins12(k,2) = -inf; % discard this step with no peak
    end
    
    [acf13, lags13] = xcorr(cut_n_imu_i1, cut_n_imu_i3);
    [val13, index13] = findpeaks(acf13);
    if ~isempty(val13)
        [vm13, im13] = max(val13);
        p_shift_bins13(k,1) = lags13(index13(im13));
        p_shift_bins13(k,2) = max(acf13);
    else
        p_shift_bins13(k,2) = -inf; % discard this step with no peak
    end

    [acf14, lags14] = xcorr(cut_n_imu_i1, cut_n_imu_i4);
    [val14, index14] = findpeaks(acf14);
    if ~isempty(val14)
        [vm14, im14] = max(val14);
        p_shift_bins14(k,1) = lags14(index14(im14));
        p_shift_bins14(k,2) = max(acf14);
    else
        p_shift_bins14(k,2) = -inf; % discard this step with no peak
    end

    p_shift_bins12(k,3) = w_start;
    p_shift_bins13(k,3) = w_start;
    p_shift_bins14(k,3) = w_start;
    k = k + 1;
end

[~, shift12_index] = max(p_shift_bins12(:,2));
shift_bins12 = p_shift_bins12(shift12_index, 1);
shift_bins12_start = p_shift_bins12(shift12_index, 3)
[~, shift13_index] = max(p_shift_bins13(:,2));
shift_bins13 = p_shift_bins13(shift13_index, 1);
shift_bins13_start = p_shift_bins13(shift13_index, 3)
[~, shift14_index] = max(p_shift_bins14(:,2));
shift_bins14 = p_shift_bins14(shift14_index, 1);
shift_bins14_start = p_shift_bins14(shift14_index, 3)

% If locating one of the cali phase wrong - both sensor cali at same time
% find the next cali-time step_away second away from current point
step_away = 1;
if abs(shift_bins13_start - shift_bins12_start) < 1 % 1 second
    if p_shift_bins12(shift12_index, 2) > p_shift_bins13(shift13_index, 2) % 12 is correct
        [~, shift13_index] = max(p_shift_bins13(:,2));
        while abs(p_shift_bins13(shift13_index, 3)-shift_bins12_start) < step_away
            p_shift_bins13(shift13_index,2) = 0;
            [~, shift13_index] = max(p_shift_bins13(:,2));
        end
        shift_bins13 = p_shift_bins13(shift13_index, 1);
        shift_bins13_start = p_shift_bins13(shift13_index, 3)
    end
    if p_shift_bins12(shift12_index, 2) < p_shift_bins13(shift13_index, 2) % 13 is correct
        [~, shift12_index] = max(p_shift_bins12(:,2));
        while abs(p_shift_bins12(shift12_index, 3)-shift_bins13_start) < step_away
            p_shift_bins12(shift12_index,2) = 0;
            [~, shift12_index] = max(p_shift_bins12(:,2));
        end
        shift_bins12 = p_shift_bins12(shift12_index, 1);
        shift_bins12_start = p_shift_bins12(shift12_index, 3)
    end
end
if abs(shift_bins13_start - shift_bins14_start) < 1 % 1 second
    if p_shift_bins14(shift14_index, 2) > p_shift_bins13(shift13_index, 2) % 14 is correct
        [~, shift13_index] = max(p_shift_bins13(:,2));
        while abs(p_shift_bins13(shift13_index, 3)-shift_bins14_start) < step_away
            p_shift_bins13(shift13_index,2) = 0;
            [~, shift13_index] = max(p_shift_bins13(:,2));
        end
        shift_bins13 = p_shift_bins13(shift13_index, 1);
        shift_bins13_start = p_shift_bins13(shift13_index, 3)
    end
    if p_shift_bins14(shift14_index, 2) < p_shift_bins13(shift13_index, 2) % 13 is correct
        [~, shift14_index] = max(p_shift_bins14(:,2));
        while abs(p_shift_bins14(shift14_index, 3)-shift_bins13_start) < step_away
            p_shift_bins14(shift14_index,2) = 0;
            [~, shift14_index] = max(p_shift_bins14(:,2));
        end
        shift_bins14 = p_shift_bins14(shift14_index, 1);
        shift_bins14_start = p_shift_bins14(shift14_index, 3)
    end
end
if abs(shift_bins12_start - shift_bins14_start) < 1 % 1 second
    if p_shift_bins14(shift14_index, 2) > p_shift_bins12(shift12_index, 2) % 14 is correct
        [~, shift13_index] = max(p_shift_bins13(:,2));
        while abs(p_shift_bins12(shift12_index, 3)-shift_bins14_start) < step_away
            p_shift_bins12(shift12_index,2) = 0;
            [~, shift12_index] = max(p_shift_bins12(:,2));
        end
        shift_bins12 = p_shift_bins12(shift12_index, 1);
        shift_bins12_start = p_shift_bins12(shift12_index, 3)
    end
    if p_shift_bins14(shift14_index, 2) < p_shift_bins12(shift12_index, 2) % 12 is correct
        [~, shift14_index] = max(p_shift_bins14(:,2));
        while abs(p_shift_bins14(shift14_index, 3)-shift_bins12_start) < step_away
            p_shift_bins14(shift14_index,2) = 0;
            [~, shift14_index] = max(p_shift_bins14(:,2));
        end
        shift_bins14 = p_shift_bins14(shift14_index, 1);
        shift_bins14_start = p_shift_bins14(shift14_index, 3)
    end
end

fprintf("The sensor shift12 is %d\n", shift_bins12);
fprintf("The sensor shift13 is %d\n", shift_bins13);
fprintf("The sensor shift14 is %d\n", shift_bins14);

% apply the shift
t_i2 = t_i2 + shift_bins12 * 1/500;
t_i3 = t_i3 + shift_bins13 * 1/500;
t_i4 = t_i4 + shift_bins14 * 1/500;

% check the shift
figure(335)
xl = [shift_bins12_start, shift_bins12_start+w_len];
cut_imu1 = imu1(t_i1 > xl(1) & t_i1 < xl(2));
cut_t_i1 = t_i1(t_i1 > xl(1) & t_i1 < xl(2));
cut_imu2 = imu2(t_i2 > xl(1) & t_i2 < xl(2));
cut_t_i2 = t_i2(t_i2 > xl(1) & t_i2 < xl(2));
% normalzied & interpolate the data
[cut_n_imu_i1,cut_t_i_i1] = interpolationWithTimestamp(cut_imu1,cut_t_i1,fs);
[cut_n_imu_i2,cut_t_i_i2] = interpolationWithTimestamp(cut_imu2,cut_t_i2,fs);
cut_n_imu_i1 = detrend(normalize(cut_n_imu_i1, 'range'));
cut_n_imu_i2 = detrend(normalize(cut_n_imu_i2, 'range'));
subplot(3,1,1)
plot(cut_t_i_i1, cut_n_imu_i1)
hold on
plot(cut_t_i_i2, cut_n_imu_i2)
hold off
title("Sensor 1 & Sensor 2")
xl = [shift_bins13_start, shift_bins13_start+w_len];
cut_imu1 = imu1(t_i1 > xl(1) & t_i1 < xl(2));
cut_t_i1 = t_i1(t_i1 > xl(1) & t_i1 < xl(2));
cut_imu3 = imu3(t_i3 > xl(1) & t_i3 < xl(2));
cut_t_i3 = t_i3(t_i3 > xl(1) & t_i3 < xl(2));
% normalzied & interpolate the data
[cut_n_imu_i1,cut_t_i_i1] = interpolationWithTimestamp(cut_imu1,cut_t_i1,fs);
[cut_n_imu_i3,cut_t_i_i3] = interpolationWithTimestamp(cut_imu3,cut_t_i3,fs);
cut_n_imu_i1 = detrend(normalize(cut_n_imu_i1, 'range'));
cut_n_imu_i3 = detrend(normalize(cut_n_imu_i3, 'range'));
subplot(3,1,2)
plot(cut_t_i_i1, cut_n_imu_i1)
hold on
plot(cut_t_i_i3, cut_n_imu_i3)
hold off
title("Sensor 1 & Sensor 3")
xl = [shift_bins14_start, shift_bins14_start+w_len];
cut_imu1 = imu1(t_i1 > xl(1) & t_i1 < xl(2));
cut_t_i1 = t_i1(t_i1 > xl(1) & t_i1 < xl(2));
cut_imu4 = imu4(t_i4 > xl(1) & t_i4 < xl(2));
cut_t_i4 = t_i4(t_i4 > xl(1) & t_i4 < xl(2));
% normalzied & interpolate the data
[cut_n_imu_i1,cut_t_i_i1] = interpolationWithTimestamp(cut_imu1,cut_t_i1,fs);
[cut_n_imu_i4,cut_t_i_i4] = interpolationWithTimestamp(cut_imu4,cut_t_i4,fs);
cut_n_imu_i1 = detrend(normalize(cut_n_imu_i1, 'range'));
cut_n_imu_i4 = detrend(normalize(cut_n_imu_i4, 'range'));
subplot(3,1,3)
plot(cut_t_i_i1, cut_n_imu_i1)
hold on
plot(cut_t_i_i4, cut_n_imu_i4)
hold off
title("Sensor 1 & Sensor 4")

%% groundtruth labeling 
post_sensor_cali_time = floor(max([shift_bins12_start, shift_bins13_start, shift_bins14_start]))%% Get the groundtruth labels
% time data points
time_axis_scg = t_i1;
time_axis_ppg = t_i2;
time_axis_bcg = t_i3;
time_axis_neck = t_i4;

time_duration_gt = 180;
scg_data = z(5:end-1);
scg_data(o1) = [];
ppg_data = double(ADC_values2); ppg_data = ppg_data(5:end-1);
ppg_data(o2) = [];
bcg_data = x3(5:end-1);
bcg_data(o3) = [];
neck_data = z4(5:end-1);
neck_data(o4) = [];
[time_axis_scg,scg_data,~] = outlinerRemoval4DataPairs(time_axis_scg,scg_data);
[time_axis_ppg,ppg_data,~] = outlinerRemoval4DataPairs(time_axis_ppg,ppg_data);
[time_axis_bcg,bcg_data,~] = outlinerRemoval4DataPairs(time_axis_bcg,bcg_data);
[time_axis_neck,neck_data,~] = outlinerRemoval4DataPairs(time_axis_neck,neck_data);


fs_scg = length(time_stamp)/time_duration_gt;
fs_ppg = length(time_stamp2)/time_duration_gt;
fs_bcg = length(time_stamp3)/time_duration_gt;
fs_neck = length(time_stamp4)/time_duration_gt;

% PPG data processing
% Use the intersection tangent to find the peaks of ppg signal
fc = 10;
[b, a] = butter(4, fc/(fs_ppg/2), 'low'); 
smoothed_ppg = filtfilt(b, a, ppg_data);
first_derivative = computeFirstDerivative(smoothed_ppg, 1/fs_ppg);
[val_ppg, peak_index_ppg, width, prom] = findpeaks(normalize(first_derivative, 'range'), "MinPeakDistance", 0.4*fs_ppg, "MinPeakProminence", 0.001);
% get desire minPeakDistance & minPeakProminence
HBV_threshold = 0.2; %second
refined_MPD = mean(diff(peak_index_ppg))/fs - HBV_threshold;
refined_MPP = median(prom)*0.5;
[val_ppg, peak_index_ppg] = findpeaks(normalize(first_derivative, 'range'), "MinPeakDistance", refined_MPD*fs_ppg, "MinPeakProminence", refined_MPP);
% only check the data after cali_time second
peak_index_ppg = peak_index_ppg(peak_index_ppg > cali_time*fs_ppg);

valley_index_ppg = zeros(1, length(peak_index_ppg)-1);
for i = 1:length(peak_index_ppg)-1
    range = peak_index_ppg(i):peak_index_ppg(i+1);
    [~, min_idx] = min(smoothed_ppg(range));
    valley_index_ppg(i) = range(min_idx);
end
% find the intersection
intersection_points = zeros(1, length(valley_index_ppg));
for i = 1:length(valley_index_ppg)
    peak_idx = peak_index_ppg(i+1);
    valley_idx = valley_index_ppg(i);
    % Tangent at the peak
    m1 = (smoothed_ppg(peak_idx+1) - smoothed_ppg(peak_idx-1)) / (2); % Numerical derivative
    b1 = smoothed_ppg(peak_idx) - m1 * peak_idx; % y-intercept
    intersection_points(i) = round((smoothed_ppg(valley_idx)-b1)/m1);
end
% Reassign peak indices to intersection points
peak_index_ppg = intersection_points;
% delete all the NaN points
peak_index_ppg = peak_index_ppg(~isnan(peak_index_ppg));


% delete the outliers for ppg
% due to the tipical range of IBI is 100-300ms even for eilte athletes 
ppg_hbi = diff(peak_index_ppg);
outlier_index = ppg_hbi < median(ppg_hbi) - 0.3*fs_ppg | ppg_hbi > median(ppg_hbi) + 0.3*fs_ppg;
outlier_index = logical([outlier_index, 0]);
peak_index_ppg(outlier_index) = [];

% SCG data processing
fc = 0.7;
[b, a] = butter(1, fc/(fs_scg/2), 'high'); 
scg_data = filtfilt(b, a, scg_data);
fc = 30;
[b, a] = butter(1, fc/(fs_scg/2), 'low'); 
scg_data = filtfilt(b, a, scg_data);

% find the threshold with the SCG
refPPG_delay_range = [0.1, 0.3];
step_refPPG = 0.02;
sumation_RG_SCG = zeros(floor(diff(refPPG_delay_range)/step_refPPG), 1);
peak_index_scg = zeros(1, length(peak_index_ppg));
k = 1;
for SCG_RG = refPPG_delay_range(1):step_refPPG:refPPG_delay_range(2)
    % Loop through each peak in PPG and find the closest peak in SCG
    for i = 1:length(peak_index_ppg)
        % Corresponding timestamp for ppg peak
        ppg_peak_time = time_axis_ppg(peak_index_ppg(i));
        scg_peak_search_indices = find((ppg_peak_time - SCG_RG < time_axis_scg) & time_axis_scg < ppg_peak_time);
        if isempty(scg_peak_search_indices)
            last_nonzero = find(peak_index_scg ~= 0, 1, 'last');
            peak_index_scg = peak_index_scg(1:last_nonzero);
            fprintf('scg_peak_search_indices is empty. Stopping the loop and trim the ending zeros of peak finder at %d\n', last_nonzero);
            break; 
        end
        if length(scg_peak_search_indices) < 3
            peak_index_scg_section = 0;
        else
            [~, peak_index_scg_section] = findpeaks(scg_data(scg_peak_search_indices), "MinPeakDistance", length(scg_peak_search_indices) - 2, "MinPeakProminence", 0.05);
        end
        % Adjustment for initial data where findpeaks doesn't find anything
        if isempty(peak_index_scg_section)
            peak_index_scg_section = 0;
        end
        peak_index_scg(i) = scg_peak_search_indices(1) + peak_index_scg_section;
    end
    sumation_RG_SCG(k) = sum(scg_data(peak_index_scg));
    % sum(scg_data(peak_index_scg))
    % peak_index_scg
    % pause
    k = k + 1;
end
[~, SCG_RG_index] = max(sumation_RG_SCG);
SCG_RG = refPPG_delay_range(1) + step_refPPG * SCG_RG_index;

for i = 1:length(peak_index_ppg)
    % Corresponding timestamp for ppg peak
    ppg_peak_time = time_axis_ppg(peak_index_ppg(i));
    scg_peak_search_indices = find((ppg_peak_time - SCG_RG < time_axis_scg) & time_axis_scg < ppg_peak_time);
    if isempty(scg_peak_search_indices)
        last_nonzero = find(peak_index_scg ~= 0, 1, 'last');
        peak_index_scg = peak_index_scg(1:last_nonzero);
        fprintf('scg_peak_search_indices is empty. Stopping the loop and trim the ending zeros of peak finder at %d\n', last_nonzero);
        break; 
    end
    [~, peak_index_scg_section] = findpeaks(scg_data(scg_peak_search_indices), "MinPeakDistance", length(scg_peak_search_indices) - 2, "MinPeakProminence", 0.05);
    % Adjustment for initial data where findpeaks doesn't find anything
    if isempty(peak_index_scg_section)
        [~, max_index] = max(scg_data(scg_peak_search_indices));
        peak_index_scg_section = max_index;
    end
    peak_index_scg(i) = scg_peak_search_indices(1) + peak_index_scg_section;
end

% BCG data processing
fc = 0.5;
[b, a] = butter(1, fc/(fs_bcg/2), 'high'); 
bcg_data = filtfilt(b, a, bcg_data);
fc = 50;
[b, a] = butter(1, fc/(fs_bcg/2), 'low'); 
bcg_data = filtfilt(b, a, bcg_data);

% flipped it 
bcg_data = -bcg_data;

% find the threshold with the BCG
sumation_RG_BCG = zeros(floor(diff(refPPG_delay_range)/step_refPPG), 1);
peak_index_bcg = zeros(1, length(peak_index_ppg));
k = 1;
for BCG_RG = refPPG_delay_range(1):step_refPPG:refPPG_delay_range(2)
    for i = 1:length(peak_index_ppg)
        % Corresponding timestamp for ppg peak
        ppg_peak_time = time_axis_ppg(peak_index_ppg(i));
        bcg_peak_search_indices = find((ppg_peak_time - BCG_RG < time_axis_bcg) & time_axis_bcg < ppg_peak_time + 0.1); % plus 0.1 as bcg peak could be after ppg
        if length(bcg_peak_search_indices) < 3
            last_nonzero = find(peak_index_bcg ~= 0, 1, 'last');
            peak_index_bcg = peak_index_bcg(1:last_nonzero);
            fprintf('bcg_peak_search_indices is empty. Stopping the loop and trim the ending zeros of peak finder at %d\n', last_nonzero);
            break; 
        end
        [~, peak_index_bcg_section] = findpeaks(bcg_data(bcg_peak_search_indices), "MinPeakDistance", length(bcg_peak_search_indices) - 2, "MinPeakProminence", 0.05);
        % Adjustment for initial data where findpeaks doesn't find anything
        if isempty(peak_index_bcg_section)
            peak_index_bcg_section = 0;
        end
        peak_index_bcg(i) = bcg_peak_search_indices(1) + peak_index_bcg_section;
    end
    sumation_RG_BCG(k) = sum(bcg_data(peak_index_bcg));
    % sum(scg_data(peak_index_scg))
    % peak_index_scg
    % pause
    k = k + 1;
end
[~, BCG_RG_index] = max(sumation_RG_BCG);
BCG_RG = refPPG_delay_range(1) + step_refPPG * BCG_RG_index;

for i = 1:length(peak_index_ppg)
    % Corresponding timestamp for ppg peak
    ppg_peak_time = time_axis_ppg(peak_index_ppg(i));
    bcg_peak_search_indices = find((ppg_peak_time - BCG_RG < time_axis_bcg) & time_axis_bcg < ppg_peak_time + 0.1);
    if isempty(bcg_peak_search_indices)
        last_nonzero = find(peak_index_bcg ~= 0, 1, 'last');
        peak_index_bcg = peak_index_bcg(1:last_nonzero);
        fprintf('bcg_peak_search_indices is empty. Stopping the loop and trim the ending zeros of peak finder at %d\n', last_nonzero);
        break; 
    end
    [~, peak_index_bcg_section] = findpeaks(bcg_data(bcg_peak_search_indices), "MinPeakDistance", length(bcg_peak_search_indices) - 2, "MinPeakProminence", 0.05);
    % Adjustment for initial data where findpeaks doesn't find anything
    if isempty(peak_index_bcg_section)
        % peak_index_bcg_section = 0;
        [~, max_index] = max(bcg_data(bcg_peak_search_indices));
        peak_index_bcg_section = max_index;
    end
    peak_index_bcg(i) = bcg_peak_search_indices(1) + peak_index_bcg_section;
end

% Neck data processing
fc = 0.7;
[b, a] = butter(1, fc/(fs_neck/2), 'high'); 
neck_data = filtfilt(b, a, neck_data);
fc = 30;
[b, a] = butter(1, fc/(fs_neck/2), 'low'); 
neck_data = filtfilt(b, a, neck_data);

% find the threshold with the neck
sumation_RG_neck = zeros(floor(diff(refPPG_delay_range)/step_refPPG), 1);
peak_index_neck = zeros(1, length(peak_index_ppg));
k = 1;
for neck_RG = refPPG_delay_range(1):step_refPPG:refPPG_delay_range(2)
    % Loop through each peak in PPG and find the closest peak in neck
    for i = 1:length(peak_index_ppg)
        % Corresponding timestamp for ppg peak
        ppg_peak_time = time_axis_ppg(peak_index_ppg(i));
        neck_peak_search_indices = find((ppg_peak_time - neck_RG < time_axis_neck) & time_axis_neck < ppg_peak_time + 0.1); % plus 0.1 as neck peak could be after ppg
        if isempty(neck_peak_search_indices)
            last_nonzero = find(peak_index_neck ~= 0, 1, 'last');
            peak_index_neck = peak_index_neck(1:last_nonzero);
            fprintf('neck_peak_search_indices is empty. Stopping the loop and trim the ending zeros of peak finder at %d\n', last_nonzero);
            break; 
        end
        [~, peak_index_neck_section] = findpeaks(neck_data(neck_peak_search_indices), "MinPeakDistance", length(neck_peak_search_indices) - 2, "MinPeakProminence", 0.05);
        % Adjustment for initial data where findpeaks doesn't find anything
        if isempty(peak_index_neck_section)
            peak_index_neck_section = 0;
        end
        peak_index_neck(i) = neck_peak_search_indices(1) + peak_index_neck_section;
    end
    sumation_RG_neck(k) = sum(neck_data(peak_index_neck));
    % sum(neck_data(peak_index_neck))
    % peak_index_neck
    % pause
    k = k + 1;
end
[~, neck_RG_index] = max(sumation_RG_neck);
neck_RG = refPPG_delay_range(1) + step_refPPG * neck_RG_index;

for i = 1:length(peak_index_ppg)
    % Corresponding timestamp for ppg peak
    ppg_peak_time = time_axis_ppg(peak_index_ppg(i));
    neck_peak_search_indices = find((ppg_peak_time - neck_RG < time_axis_neck) & time_axis_neck < ppg_peak_time + 0.1);
    if isempty(neck_peak_search_indices)
        last_nonzero = find(peak_index_neck ~= 0, 1, 'last');
        peak_index_neck = peak_index_neck(1:last_nonzero);
        fprintf('neck_peak_search_indices is empty. Stopping the loop and trim the ending zeros of peak finder at %d\n', last_nonzero);
        break; 
    end
    [~, peak_index_neck_section] = findpeaks(neck_data(neck_peak_search_indices), "MinPeakDistance", length(neck_peak_search_indices) - 2, "MinPeakProminence", 0.05);
    % Adjustment for initial data where findpeaks doesn't find anything
    if isempty(peak_index_neck_section)
        [~, max_index] = max(neck_data(neck_peak_search_indices));
        peak_index_neck_section = max_index;
    end
    peak_index_neck(i) = neck_peak_search_indices(1) + peak_index_neck_section;
end

% Calculate the rough PTT
peak_time_scg = time_axis_scg(peak_index_scg);
peak_time_scg = peak_time_scg(peak_time_scg > cali_time);
peak_time_ppg = time_axis_ppg(peak_index_ppg);
peak_time_ppg = peak_time_ppg(peak_time_ppg > cali_time);
peak_time_bcg = time_axis_bcg(peak_index_bcg);
peak_time_bcg = peak_time_bcg(peak_time_bcg > cali_time);
peak_time_neck = time_axis_neck(peak_index_neck);
peak_time_neck = peak_time_neck(peak_time_neck > cali_time);

% In case peaks are too close that some over the ppg foot: 
% Especially after the exercise
ppgOverRange = 0.2;
peak_time_scg = peak_time_scg(peak_time_scg < max(peak_time_ppg)+ppgOverRange);
peak_time_bcg = peak_time_bcg(peak_time_bcg < max(peak_time_ppg)+ppgOverRange);
peak_time_neck = peak_time_neck(peak_time_neck < max(peak_time_ppg)+ppgOverRange);
n_peaks = min([length(peak_time_scg),length(peak_time_ppg),length(peak_time_bcg),length(peak_time_neck)]);
peak_time_scg = peak_time_scg(end-n_peaks+1:end);
peak_time_ppg = peak_time_ppg(end-n_peaks+1:end);
peak_time_bcg = peak_time_bcg(end-n_peaks+1:end);
peak_time_neck = peak_time_neck(end-n_peaks+1:end);

gt_ptt1 = peak_time_ppg - peak_time_scg;
gt_ptt2 = peak_time_bcg - peak_time_scg;
gt_ptt3 = peak_time_ppg - peak_time_bcg;
gt_ptt4 = peak_time_ppg - peak_time_neck;

% manual shift
% Refine the SCG, BCG Peaks from current draft result
% SCG
bin_width = 0.01;
[bin_counts, bin_edges] = histcounts(gt_ptt1, 'BinWidth', bin_width);
[max_count, max_index] = max(bin_counts);
highest_bin_start = bin_edges(max_index);          % Start of the range
highest_bin_end = bin_edges(max_index + 1);        % End of the range

SCG_RG = highest_bin_start + (highest_bin_end-highest_bin_start)/2 + 0.0;
% Display the results
fprintf('avg PTT-PPG&SCG is %.3f seconds.\n', SCG_RG);
SCG_peak_width_lim = 0.03;
for i = 1:length(peak_index_ppg)
    % Corresponding timestamp for ppg peak
    ppg_peak_time = time_axis_ppg(peak_index_ppg(i));
    scg_peak_search_indices = find((ppg_peak_time - SCG_RG - SCG_peak_width_lim < time_axis_scg) & time_axis_scg < ppg_peak_time - SCG_RG + SCG_peak_width_lim);
    if isempty(scg_peak_search_indices)
        last_nonzero = find(peak_index_scg ~= 0, 1, 'last');
        peak_index_scg = peak_index_scg(1:last_nonzero);
        fprintf('scg_peak_search_indices is empty. Stopping the loop and trim the ending zeros of peak finder at %d\n', last_nonzero);
        break; 
    end
    [~, peak_index_scg_section] = findpeaks(scg_data(scg_peak_search_indices), "MinPeakDistance", length(scg_peak_search_indices) - 2, "MinPeakProminence", 0.05);
    % Adjustment for initial data where findpeaks doesn't find anything
    if isempty(peak_index_scg_section)
        [~, max_index] = max(scg_data(scg_peak_search_indices));
        peak_index_scg_section = max_index;
    end
    peak_index_scg(i) = scg_peak_search_indices(1) + peak_index_scg_section;
end

% BCG
bin_width = 0.01;
[bin_counts, bin_edges] = histcounts(gt_ptt3, 'BinWidth', bin_width);
[max_count, max_index] = max(bin_counts);
highest_bin_start = bin_edges(max_index);          % Start of the range
highest_bin_end = bin_edges(max_index + 1);        % End of the range

BCG_RG = highest_bin_start + (highest_bin_end-highest_bin_start)/2;
% Display the results
fprintf('avg PTT-PPG&BCG is %.3f seconds.\n', BCG_RG);
BCG_peak_width_lim = 0.05;

for i = 1:length(peak_index_ppg)
    % Corresponding timestamp for ppg peak
    ppg_peak_time = time_axis_ppg(peak_index_ppg(i));
    bcg_peak_search_indices = find((ppg_peak_time - BCG_RG - BCG_peak_width_lim < time_axis_bcg) & time_axis_bcg < ppg_peak_time - BCG_RG + BCG_peak_width_lim);
    if isempty(bcg_peak_search_indices)
        last_nonzero = find(peak_index_bcg ~= 0, 1, 'last');
        peak_index_bcg = peak_index_bcg(1:last_nonzero);
        fprintf('bcg_peak_search_indices is empty. Stopping the loop and trim the ending zeros of peak finder at %d\n', last_nonzero);
        break; 
    end
    [~, peak_index_bcg_section] = findpeaks(bcg_data(bcg_peak_search_indices), "MinPeakDistance", length(bcg_peak_search_indices) - 2, "MinPeakProminence", 0.05);
    % Adjustment for initial data where findpeaks doesn't find anything
    if isempty(peak_index_bcg_section)
        [~, max_index] = max(bcg_data(bcg_peak_search_indices));
        peak_index_bcg_section = max_index;
    end
    peak_index_bcg(i) = bcg_peak_search_indices(1) + peak_index_bcg_section;
end

% neck
bin_width = 0.01;
[bin_counts, bin_edges] = histcounts(gt_ptt4, 'BinWidth', bin_width);
[max_count, max_index] = max(bin_counts);
highest_bin_start = bin_edges(max_index);          % Start of the range
highest_bin_end = bin_edges(max_index + 1);        % End of the range

neck_RG = highest_bin_start + (highest_bin_end-highest_bin_start)/2 + 0.0;
% Display the results
fprintf('avg PTT-PPG&neck is %.3f seconds.\n', neck_RG);
neck_peak_width_lim = 0.04;
for i = 1:length(peak_index_ppg)
    % Corresponding timestamp for ppg peak
    ppg_peak_time = time_axis_ppg(peak_index_ppg(i));
    neck_peak_search_indices = find((ppg_peak_time - neck_RG - neck_peak_width_lim < time_axis_neck) & time_axis_neck < ppg_peak_time - neck_RG + neck_peak_width_lim);
    if isempty(neck_peak_search_indices)
        last_nonzero = find(peak_index_neck ~= 0, 1, 'last');
        peak_index_neck = peak_index_neck(1:last_nonzero);
        fprintf('neck_peak_search_indices is empty. Stopping the loop and trim the ending zeros of peak finder at %d\n', last_nonzero);
        break; 
    end
    [~, peak_index_neck_section] = findpeaks(neck_data(neck_peak_search_indices), "MinPeakDistance", length(neck_peak_search_indices) - 2, "MinPeakProminence", 0.05);
    % Adjustment for initial data where findpeaks doesn't find anything
    if isempty(peak_index_neck_section)
        [~, max_index] = max(neck_data(neck_peak_search_indices));
        peak_index_neck_section = max_index;
    end
    peak_index_neck(i) = neck_peak_search_indices(1) + peak_index_neck_section;
end

xl = [cali_time, time_duration_gt];
xl = [cali_time+20, cali_time+30];
figure(102);clf;
subplot(6, 1, 2);
plot(time_axis_ppg, ppg_data);
hold on 
plot(time_axis_ppg(peak_index_ppg), ppg_data(peak_index_ppg), 'ro', 'MarkerFaceColor', 'r');
hold off
title('PPG Data');
xlabel('Time in seconds');
ylabel('Amplitude');
xlim(xl)
subplot(6, 1, 1);
plot(time_axis_scg, scg_data);
hold on 
plot(time_axis_scg(peak_index_scg), scg_data(peak_index_scg), 'ro', 'MarkerFaceColor', 'r');
hold off
title('SCG Data');
xlabel('Time in seconds');
ylabel('m/s^2');
xlim(xl)
subplot(6, 1, 3);
plot(time_axis_bcg, bcg_data);
hold on 
plot(time_axis_bcg(peak_index_bcg), bcg_data(peak_index_bcg), 'ro', 'MarkerFaceColor', 'r');
hold off
title('BCG flipped Data');
xlabel('Time in seconds');
ylabel('m/s^2');
xlim(xl)
subplot(6, 1, 4);
plot(time_axis_neck, neck_data);
hold on 
plot(time_axis_neck(peak_index_neck), neck_data(peak_index_neck), 'ro', 'MarkerFaceColor', 'r');
hold off
title('neck Data');
xlabel('Time in seconds');
ylabel('m/s^2');
xlim(xl)

% calculate PTT
peak_time_scg = time_axis_scg(peak_index_scg);
peak_time_scg_op = peak_time_scg;
peak_time_scg = peak_time_scg(peak_time_scg > cali_time);
peak_time_ppg = time_axis_ppg(peak_index_ppg);
peak_time_ppg_op = peak_time_ppg;
peak_time_ppg = peak_time_ppg(peak_time_ppg > cali_time);
peak_time_bcg = time_axis_bcg(peak_index_bcg);
peak_time_bcg_op = peak_time_bcg;
peak_time_bcg = peak_time_bcg(peak_time_bcg > cali_time);
peak_time_neck = time_axis_neck(peak_index_neck);
peak_time_neck_op = peak_time_neck;
peak_time_neck = peak_time_neck(peak_time_neck > cali_time);

% In case peaks are too close that some over the ppg foot: 
% Especially after the exercise
peak_time_scg = peak_time_scg(peak_time_scg < max(peak_time_ppg)+ppgOverRange);
peak_time_bcg = peak_time_bcg(peak_time_bcg < max(peak_time_ppg)+ppgOverRange);
peak_time_neck = peak_time_neck(peak_time_neck < max(peak_time_ppg)+ppgOverRange);
n_peaks = min([length(peak_time_scg),length(peak_time_ppg),length(peak_time_bcg),length(peak_time_neck)]);
peak_time_scg = peak_time_scg(end-n_peaks+1:end);
peak_time_ppg = peak_time_ppg(end-n_peaks+1:end);
peak_time_bcg = peak_time_bcg(end-n_peaks+1:end);
peak_time_neck = peak_time_neck(end-n_peaks+1:end);

gt_ptt1 = peak_time_ppg - peak_time_scg;
gt_ptt2 = peak_time_bcg - peak_time_scg;
gt_ptt3 = peak_time_ppg - peak_time_bcg;
gt_ptt4 = peak_time_ppg - peak_time_neck;
gt_ptt5 = peak_time_neck - peak_time_scg;
gt_ptt6 = peak_time_bcg - peak_time_neck;

median_gt_heart2Wrist_ptt = median(gt_ptt1)
median_gt_heart2Head_ptt = median(gt_ptt2)
median_gt_head2Wrist_ptt = median(gt_ptt3)
median_gt_neck2Wrist_ptt = median(gt_ptt4)
median_gt_heart2Neck_ptt = median(gt_ptt5)

subplot(6, 1, 5);
dscg = diff(peak_time_scg);
dppg = diff(peak_time_ppg);
dbcg = diff(peak_time_bcg);
dneck = diff(peak_time_neck);
hold on
plot(diff(peak_time_scg))
% to eliminate the shift if the cutoff is between first SCG and PPG peak
if peak_time_ppg(1) < peak_time_scg(1)
    plot(diff(peak_time_ppg(2:end)))
else
    plot(diff(peak_time_ppg))
end
plot(diff(peak_time_bcg))
plot(diff(peak_time_neck))
ylim([0.3, 1])
data_BPM = 1/median(dppg)*60;
title(sprintf("IBI GroundTruth with PPG bpm: %0.2f", data_BPM))
xlabel('data points');
ylabel('Time (seconds)');
legend("IBI-SCG", "IBI-PPG", "IBI-BCG", "IBI-neck")

subplot(6, 1, 6);
hold on
plot(gt_ptt1)
% plot(gt_ptt2)
% plot(gt_ptt3)
plot(gt_ptt4)
plot(gt_ptt5)
% plot(gt_ptt6)
ylim([-0.2, 0.4])
title(sprintf("PTT GroundTruth with PPG bpm: %0.2f\n heart2wrist: %0.4f  &  heart2neck: %0.4f", data_BPM, median_gt_heart2Wrist_ptt, median_gt_heart2Neck_ptt))
xlabel('data points');
ylabel('Time (seconds)');
legend("PTT-PPG&SCG", "PPT-PPG&neck", "PPT-SCG&neck")

numStd = 2;
% Delete outlier for SCG
fc = 0.6;
[b, a] = butter(4, fc/(fs_scg/2), 'high'); 
scg_data_l = filtfilt(b, a, scg_data);
pv_scg = scg_data_l(peak_index_scg);
outlier_index = pv_scg > median(pv_scg) + numStd*std(pv_scg) | pv_scg < median(pv_scg) - numStd*std(pv_scg);
outlier_index = logical([outlier_index, 0]);
peak_index_scg(outlier_index) = [];

subplot(6, 1, 1);
plot(time_axis_scg, scg_data);
hold on  
plot(time_axis_scg(peak_index_scg), scg_data(peak_index_scg), 'ro', 'MarkerFaceColor', 'r');
hold off
title('SCG Data');
xlabel('Time in seconds');
ylabel('m/s^2');
xlim(xl)

% Delete outlier for PPG 
fc = 0.6;
[b, a] = butter(4, fc/(fs_ppg/2), 'high'); 
ppg_data_l = filtfilt(b, a, ppg_data);
pv_ppg = ppg_data_l(peak_index_ppg);
outlier_index = pv_ppg > median(pv_ppg) + numStd*std(pv_ppg) | pv_ppg < median(pv_ppg) - numStd*std(pv_ppg);
outlier_index = logical([outlier_index, 0]);
peak_index_ppg(outlier_index) = [];

subplot(6, 1, 2);
plot(time_axis_ppg, ppg_data);
hold on 
plot(time_axis_ppg(peak_index_ppg), ppg_data(peak_index_ppg), 'ro', 'MarkerFaceColor', 'r');
hold off
title('PPG Data');
xlabel('Time in seconds');
ylabel('Amplitude');
xlim(xl)

% Delete outlier for BCG 
fc = 0.6;
[b, a] = butter(4, fc/(fs_bcg/2), 'high'); 
bcg_data_l = filtfilt(b, a, bcg_data);
pv_bcg = bcg_data_l(peak_index_bcg);
outlier_index = pv_bcg > median(pv_bcg) + numStd*std(pv_bcg) | pv_bcg < median(pv_bcg) - numStd*std(pv_bcg);
outlier_index = logical([outlier_index, 0]);
peak_index_bcg(outlier_index) = [];

subplot(6, 1, 3);
plot(time_axis_bcg, bcg_data);
hold on 
plot(time_axis_bcg(peak_index_bcg), bcg_data(peak_index_bcg), 'ro', 'MarkerFaceColor', 'r');
hold off
title('BCG flipped Data');
xlabel('Time in seconds');
ylabel('m/s^2');
xlim(xl)

% Delete outlier for Neck 
fc = 0.6;
[b, a] = butter(4, fc/(fs_neck/2), 'high'); 
neck_data_l = filtfilt(b, a, neck_data);
pv_neck = neck_data_l(peak_index_neck);
outlier_index = pv_neck > median(pv_neck) + numStd*std(pv_neck) | pv_neck < median(pv_neck) - numStd*std(pv_neck);
outlier_index = logical([outlier_index, 0]);
peak_index_neck(outlier_index) = [];

subplot(6, 1, 4);
plot(time_axis_neck, neck_data);
hold on 
plot(time_axis_neck(peak_index_neck), neck_data(peak_index_neck), 'ro', 'MarkerFaceColor', 'r');
hold off
title('neck Data');
xlabel('Time in seconds');
ylabel('m/s^2');
xlim(xl)
