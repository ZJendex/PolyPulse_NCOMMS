
function [output_data_on_bin] = siganlProcessing_basic(data_on_bin, fs, varargin)
    filter = true;
    phase = true;
    filter_range = [0.5, 60];
    flipped = false;
    emd_i = [-1, -1];

    for i = 1:2:length(varargin)
        if strcmp(varargin{i}, "filter")
            filter = varargin{i+1};
        end
        if strcmp(varargin{i}, "phase")
            phase = boolean(varargin{i+1});
        end
        if strcmp(varargin{i}, "filter_range")
            filter_range = varargin{i+1};
        end
        if strcmp(varargin{i}, "flipped")
            flipped = boolean(varargin{i+1});
        end
        if strcmp(varargin{i}, "emd_i")
            emd_i = varargin{i+1};
        end
    end

    %% Get phase data
    if phase
        single_fft_bin = angle(data_on_bin);
        single_fft_bin = squeeze(single_fft_bin);
        single_fft_bin = unwrap(single_fft_bin);
    else
        single_fft_bin = log10(abs(data_on_bin));
    end

    %% Emd wave decomposition
    if emd_i(1) ~= -1
        [imf,residual] = emd(single_fft_bin);
        single_fft_bin = sum(imf(:, emd_i(1):emd_i(2)), 2);
    end

    %% flipped the signal before the processing
    if flipped
        single_fft_bin = -single_fft_bin;
    end


    %% Apply butterworth filter
    if filter
        fc = filter_range(1);
        [b, a] = butter(2, fc/(fs/2), 'high'); 
        single_fft_bin = filtfilt(b, a, single_fft_bin);
        fc = filter_range(2);
        [b, a] = butter(2, fc/(fs/2), 'low'); 
        single_fft_bin = filtfilt(b, a, single_fft_bin);
    else
        fprintf("Without filter\n");
    end

    %% Detrend
    single_fft_bin = detrend(single_fft_bin, 4);

    %% Output Result
    output_data_on_bin = single_fft_bin;
end 


