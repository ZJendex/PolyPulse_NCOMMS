function ...
    [max_bin_var_idx_history, max_bin_var_idx, fs] = find_bin(imgs, obj_distance, cutoff_start, cutoff_end, varargin)
a_radar_config;
% init value
max_bin_var = -inf;
max_bin_var_idx = -inf;

max_bin_var_idx_in_xyAxis = -inf;
max_bin_var_idx_history = zeros(num_of_frames,8);
m = 1;
%% varargin init
thresholdOnBpm = 30; %Hz
targetBpm = 60; %bpm
targetAngle = 0; %degree
r_t = 2; % #range bin threshold
a_t = 180; % angle degree threshold
secondDerivative = false;
filter_range = [0.5, 100];
flipped = false;
%% Check optional inputs
for i = 1:2:length(varargin)
    if strcmp(varargin{i}, 'thresholdOnBpm')
        thresholdOnBpm = double(varargin{i+1});  
    end
    if strcmp(varargin{i}, 'targetBpm')
        targetBpm = double(varargin{i+1});  
    end
    if strcmp(varargin{i}, 'targetAngleDegree')
        targetAngle = double(varargin{i+1});  
    end
    if strcmp(varargin{i}, 'rangeThresholdIndex')
        r_t = int8(varargin{i+1});  
    end
    if strcmp(varargin{i}, 'angleThresholdDegree')
        a_t = double(varargin{i+1});  
    end
    if strcmp(varargin{i}, 'secondDerivative')
        secondDerivative = boolean(varargin{i+1});  
    end
    if strcmp(varargin{i}, "filter_range")
        filter_range = varargin{i+1};
    end
    if strcmp(varargin{i}, "flipped")
        flipped = boolean(varargin{i+1});
    end
end

bpmFreq_upperbound = 60/(targetBpm+thresholdOnBpm);
bpmFreq_lowerbound = 60/(targetBpm-thresholdOnBpm);

%% Cut the data
if ~flipped
    imgs_ra = imgs(cutoff_start:cutoff_end,:,:);
else
    imgs_ra = -imgs(cutoff_start:cutoff_end,:,:);
end


%% Focus range in Cartesian Coordinates mapped range plot
calibration_bins = 1+2; % move one bin
target_distance = obj_distance; 

range_start_bin = max(2, round(target_distance/rangeResolution)-r_t+calibration_bins);
range_end_bin = min(RANGE_FFT, round(target_distance/rangeResolution)+r_t+calibration_bins);
for index = range_start_bin:range_end_bin
    for a_b = 1:1:angle_NFFT
        angle_bin_idx = a_b;
        range_bin_idx = index;
        vertical_pos_h = y_axis(index,a_b);
        horizontal_pos_h = x_axis(index,a_b);
        cur_bin_angle = rad2deg(angle(horizontal_pos_h + vertical_pos_h*1i));
        if cur_bin_angle > targetAngle + a_t || cur_bin_angle < targetAngle - a_t
            continue
        end

        % Get Bins
        data_on_bin = imgs_ra(:, range_bin_idx, angle_bin_idx);
        
        % Get Magnitude
        magnitude_on_bin = sum(abs(imgs(round(length(imgs)/2):round(length(imgs)/2 + 2*fs), range_bin_idx, angle_bin_idx))) ./ (2 * fs);
        % Pre-process signal
        [output_data_on_bin] = siganlProcessing_basic(data_on_bin, fs, "filter_range", filter_range, "flipped", 0);
        if secondDerivative
            output_data_on_bin = computeSecondDerivative(output_data_on_bin, 1/fs);
        end
        %% Calculate xcorr
        [acf, lags] = xcorr(output_data_on_bin,output_data_on_bin, "normalized");
        % Targeting period estimate
        startIdx = ceil(length(acf)/2) + round(bpmFreq_upperbound*fs);
        endIdx = ceil(length(acf)/2) + round(bpmFreq_lowerbound*fs);
        target_acf = (acf(startIdx:endIdx));
        % Find Peaks
        [acf_peaks, peak_lag_indices] = findpeaks(target_acf,'MinPeakWidth', 5);

        % check if peak exist
        if length(acf_peaks) < 1
            continue 
        end

        % find the highest correlation and the shift
        [max_acf, max_acf_index] = max(acf_peaks);
        peak_lags = lags(startIdx:endIdx);
        peak_lags = peak_lags(peak_lag_indices);
        peak_lags = peak_lags(max_acf_index);

        if max_acf > -10
            max_bin_var_idx_history(m,1) = max_acf;
            max_bin_var_idx_history(m,2) = magnitude_on_bin;
            max_bin_var_idx_history(m,3) = range_bin_idx;
            max_bin_var_idx_history(m,4) = angle_bin_idx;
            max_bin_var_idx_history(m,5) = peak_lags;
            max_bin_var_idx_history(m,6) = cur_bin_angle;
            m = m + 1;
        end
    end
end
max_bin_var_idx_history = sortrows(max_bin_var_idx_history, -1);

% init max_bin_var_index with the largest correlation value
max_bin_var_idx = [max_bin_var_idx_history(1,3),max_bin_var_idx_history(1,4)];

end

