function [dist, t] = plot_waveformOnBin(max_bin_var_idx, imgs, cutoff_start, cutoff_end, varargin)
%% radar configuration
a_radar_config;

%% varargin
rawX = false;  % By default, this feature is disabled
filter = true;
flipped = false;
phase = true;
filter_range = [0.5, 100];
emd_i = [-1, -1];
doPlot = true;
%% Check optional inputs
for i = 1:2:length(varargin)
    if strcmp(varargin{i}, 'rawX')
        rawX = varargin{i+1};  % Expect the next argument to be the value for 'rawX'
    end
    if strcmp(varargin{i}, "filter")
        filter = varargin{i+1};
    end
    if strcmp(varargin{i}, "flipped")
        flipped = boolean(varargin{i+1});
    end
    if strcmp(varargin{i}, "phase")
        phase = boolean(varargin{i+1});
    end
    if strcmp(varargin{i}, "filter_range")
        filter_range = varargin{i+1};
    end
    if strcmp(varargin{i}, "emd_i")
        emd_i = varargin{i+1};
    end
    if strcmp(varargin{i}, "plot")
        doPlot = boolean(varargin{i+1});
    end
end
%% retrieve the target bin
if nargin > 3 && cutoff_start > 0
    data_on_bin = imgs(cutoff_start:cutoff_end,max_bin_var_idx(1),max_bin_var_idx(2));
else
    data_on_bin = imgs(:,max_bin_var_idx(1),max_bin_var_idx(2));
end

% signal process
[output_data_on_bin] = siganlProcessing_basic(data_on_bin, fs, "filter", filter, "phase", phase, "filter_range", filter_range, "flipped", flipped, "emd_i", emd_i);

% convert to distance
dist = output_data_on_bin * c / (4*pi*Fc) * 1000; % mm
% x
t = (1:length(output_data_on_bin))/fs;


%% Plotting
% figure;
if doPlot
    if rawX
        plot(dist)
        xlabel("Samples")
    else
        plot(t, dist)
        xlabel("Time in Second")
    end

    ylabel("Distance (mm)")
    %% calculate angles and distance
    angle_bin_idx = max_bin_var_idx(2);
    range_bin_idx = max_bin_var_idx(1);
    vertical_pos_h = y_axis(range_bin_idx,angle_bin_idx);
    horizontal_pos_h = x_axis(range_bin_idx,angle_bin_idx);
    
    cur_bin_angle = rad2deg(angle(horizontal_pos_h + vertical_pos_h*1i));
    range_in_meter = range_bin_idx * rangeResolution;
    if phase
        title(sprintf('Angle\n Range on bin %d, Angle on bin %d\n fs %d\n DIST: %0.2fm ANGLE %0.2fdegree', max_bin_var_idx(1),max_bin_var_idx(2), fs, range_in_meter, cur_bin_angle));
    else
        title(sprintf('Magnitude\n Range on bin %d, Angle on bin %d\n fs %d\n DIST: %0.2fm ANGLE %0.2fdegree', max_bin_var_idx(1),max_bin_var_idx(2), fs, range_in_meter, cur_bin_angle));
    end
end
