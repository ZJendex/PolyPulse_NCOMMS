function [best_index, max_count] = most_appear_number(sorted_array, threshold)
% Given sorted array
array = sorted_array;

% Define the default threshold (window size)
if nargin == 1
    threshold = 1;
end

% Initialize variables to track the window with the most elements
max_count = 0;
best_index = 0;

% Two-pointer approach to slide the window
start_value = array(1);
end_value = array(end);
for i = start_value:end_value
    count = length(array(array > i - threshold-1 & array < i + threshold+1)); 
    
    % Update the max_count and best_window_start if this window has more elements
    if count > max_count
        max_count = count;
        best_index = i;
    end
end

end