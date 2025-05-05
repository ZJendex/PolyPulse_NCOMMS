% load radar example data
load("data\radar_data.mat")

% load radar config
a_radar_config;

% data process
radar_cali = z2(5:end-1);
radar_cali(o2) = [];
imu = detrend(radar_cali);
t_i = t_i2;
hard_code_radar_range = 7;

% focus range
% in case at radar sensor cali step, sensors had better sync waveform
sync_backward_time = 3; % second
xl_r = [post_sensor_cali_time-sync_backward_time,cali_time];
% xl_r = [10,cali_time];
angle_bin_range = [50, 60];
angles_shift_bins = zeros(diff(angle_bin_range), 3);
angle_index = 1;
for angle_bin = angle_bin_range(1):angle_bin_range(2)
    [dist, t] = plot_waveformOnBin([hard_code_radar_range, angle_bin], imgs, xl_r(1)*fs, xl_r(2)*fs, ...
                        'flipped', 0, "filter", 1, "filter_range", [0.7, 30], "phase", 1, "plot", 0);
    d_dist = computeSecondDerivative(dist, 1/fs);
    % filter for griping vertical movement
    fc = 10;
    [b, a] = butter(4, fc/(fs/2), 'low'); 
    d_dist = filtfilt(b, a, d_dist);

    k = 1;
    % find shift_bins -- sliding window to find the smallest corr_shift
    step = 0.5; % second
    w_len = 3; % second
    p_shift_bins = zeros((xl_r(2)-xl_r(1)-step)/step, 3);
    for w_start = 0:step:xl_r(2)-xl_r(1)-step
        xl = [w_start, w_start+w_len];

        % cut the data
        cut_d_dist = d_dist(t > xl(1) & t < xl(2));
        cut_t = t(t > xl(1) & t < xl(2));
        cut_imu = imu(t_i > xl(1) + xl_r(1) & t_i < xl(2) + xl_r(1));
        cut_t_i = t_i(t_i > xl(1) + xl_r(1) & t_i < xl(2) + xl_r(1));
        % normalzied the data
        [cut_n_imu_i,cut_t_i_i] = interpolationWithTimestamp(cut_imu,cut_t_i,fs); % could do extra modify to make sure each gt start at x.002 which align with 
        cut_n_d_dist = detrend(normalize(cut_d_dist, 'range'));
        cut_n_imu = detrend(normalize(cut_n_imu_i, 'range'));

        [s_len, ~] = min([length(cut_n_d_dist), length(cut_n_imu)]);
        cut_n_d_dist = cut_n_d_dist(1:s_len);
        cut_n_imu = cut_n_imu(1:s_len);

        % figure(3142);clf;
        % hold on 
        % plot(cut_n_d_dist)
        % plot(cut_n_imu)
        % hold off
        % 
        [acf, lags] = xcorr(cut_n_d_dist, cut_n_imu);

        [val, index] = findpeaks(acf);
        if ~isempty(val)
            [vm, im] = max(val);
            shift_bins = lags(index(im));
            p_shift_bins(k,1) = shift_bins;
            p_shift_bins(k,3) = w_start;
            % p_shift_bins(k,2) = sum((cut_n_d_dist - cut_n_imu).^2) / sum(cut_n_d_dist.^2);
            p_shift_bins(k,2) = max(acf);
        else
            p_shift_bins(k,2) = -inf; % discard this step with no peak
        end
        k = k + 1;
    end
    NMSEs = p_shift_bins(:, 2);
    % use the consecutively three values to ensure at the middle of the radar sensor calibration
    sliding_NMSEs_3sum = NMSEs(1:end-2) + NMSEs(2:end-1) + NMSEs(3:end); 
    [~, shift_index] = max(sliding_NMSEs_3sum);
    shift_index = shift_index + 1; % three sum has ignore the start&end bin value
    % [~, shift_index] = min(p_shift_bins(:,2));
    angles_shift_bins(angle_index,1) = p_shift_bins(shift_index, 1); % save the best shift
    angles_shift_bins(angle_index,2) = p_shift_bins(shift_index, 2); % save the best acf value
    angles_shift_bins(angle_index,3) = p_shift_bins(shift_index, 3); % save the window location
    angle_index = angle_index + 1;
end
% [shift_bins, selected_angle_index] = min(abs(angles_shift_bins(:,1)));
[shift_bins, selected_angle_index] = max(angles_shift_bins(:,2));
shift_bins = angles_shift_bins(selected_angle_index,1)
tar_win_start = angles_shift_bins(selected_angle_index,3);

angle_bin = angle_bin_range(1) + selected_angle_index - 1; % first index in matlab
fig = figure(334);clf;
fig.WindowState = 'maximized';
subplot(3,1,1)
[dist, t] = plot_waveformOnBin([hard_code_radar_range, angle_bin], imgs, xl_r(1)*fs, xl_r(2)*fs, ...
                        'flipped', 0, "filter", 1, "filter_range", [0.7,30], "phase", 1);
xlim(xl_r)

subplot(4,1,1)
d_dist = computeSecondDerivative(dist, 1/fs);
% filter for griping vertical movement
fc = 10;
[b, a] = butter(4, fc/(fs/2), 'low'); 
d_dist = filtfilt(b, a, d_dist);
plot(t+xl_r(1), d_dist)
title("Radar second derivative")
hold on 
xline(xl_r(1)+tar_win_start, 'r', 'LineWidth', 2);      % First red line at w
xline(xl_r(1)+tar_win_start + w_len, 'r', 'LineWidth', 2);
hold off
xlim(xl_r)

subplot(4,1,2)
plot(t_i2, detrend(radar_cali))
title("Y axis Wrist IMU")
xlim(xl_r)
hold on 
xline(xl_r(1)+tar_win_start, 'r', 'LineWidth', 2);      % First red line at w
xline(xl_r(1)+tar_win_start + w_len, 'r', 'LineWidth', 2);
hold off

% before shift
subplot(4,1,4)
xl = [tar_win_start, tar_win_start + w_len];
% cut the data
t_i = time_axis_ppg;
cut_d_dist = d_dist(t > xl(1) & t < xl(2));
cut_t = t(t > xl(1) & t < xl(2));
cut_imu = radar_cali(t_i > xl(1) + xl_r(1) & t_i < xl(2) + xl_r(1));
cut_t_i = t_i(t_i > xl(1) + xl_r(1) & t_i < xl(2) + xl_r(1));
% normalzied the data
[cut_n_imu_i,cut_t_i_i] = interpolationWithTimestamp(cut_imu,cut_t_i,fs);
cut_n_d_dist = detrend(normalize(cut_d_dist, 'range'));
cut_n_imu = detrend(normalize(cut_n_imu_i, 'range'));
[s_len, ~] = min([length(cut_n_d_dist), length(cut_n_imu)]);
cut_n_d_dist = cut_n_d_dist(1:s_len);
cut_t = cut_t(1:s_len);
cut_n_imu = cut_n_imu(1:s_len);
cut_t_i_i = cut_t_i_i(1:s_len);

plot(cut_t + xl_r(1), cut_n_d_dist)
hold on 
plot(cut_t_i_i, cut_n_imu)
hold off
legend("radar", "ground truth")

% plot the acf
subplot(4,1,3)
[s_len, ~] = min([length(cut_n_d_dist), length(cut_n_imu)]);
cut_n_d_dist = cut_n_d_dist(1:s_len);
cut_n_imu = cut_n_imu(1:s_len);
[acf, lags] = xcorr(cut_n_d_dist, cut_n_imu);
[val, index] = findpeaks(acf);
[vm, im] = max(val);
shift_bins_recalculate = lags(index(im))
if shift_bins_recalculate ~= shift_bins
    disp("Shift bin is different")
    shift_bins = shift_bins_recalculate
end
plot(lags, acf)
title("acf plot")

% apply the shift to gt peaks and gt time
time_axis_scg_s = time_axis_scg + shift_bins * 1/fs;
time_axis_ppg_s = time_axis_ppg + shift_bins * 1/fs;
time_axis_bcg_s = time_axis_bcg + shift_bins * 1/fs;
time_axis_neck_s = time_axis_neck + shift_bins * 1/fs;

peak_time_scg = peak_time_scg_op + shift_bins * 1/fs;
peak_time_ppg = peak_time_ppg_op + shift_bins * 1/fs;
peak_time_bcg = peak_time_bcg_op + shift_bins * 1/fs;
peak_time_neck = peak_time_neck_op + shift_bins * 1/fs;

% after shift
subplot(4,1,4)
hold on 
plot(cut_t_i_i + shift_bins * 1/fs, cut_n_imu, 'LineStyle','--')
hold off
legend("radar", "ground truth", "ground truth with shift")
xlim([xl_r(1)+tar_win_start, xl_r(1)+tar_win_start + w_len])