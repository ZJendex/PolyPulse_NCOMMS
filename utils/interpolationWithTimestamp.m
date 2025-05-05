function [interpolated_data,uniform_time] = interpolationWithTimestamp(original_data,original_time,fs_target)
start_time = original_time(1);
end_time = original_time(end);

% Generate a regular time vector with 500 Hz frequency
uniform_time = (start_time:1/fs_target:end_time)';

% cube spline interpolation
interpolated_data = spline(original_time, original_data, uniform_time);

end

