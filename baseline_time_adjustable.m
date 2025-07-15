function [baseline_all] = baseline_time_adjustable(tau1, tau2, traces, fps)
% BASELINE time-dependent 
%   For each time_series data point (roi_t.traces), baseline is calculated as the minimum value of the
%   smoothed signal during the tau2 seconds before that point. Smooth filter is a mean
%   filter with a tau1 window (avoid edge effects)
%   tau2 = baseline window in seconds
%   tau1 = smoothing window in seconds

% Adjust tau1 and tau2 based on the fps
tau1 = floor(tau1 / (1/fps));
if mod(tau1, 2) == 0 % Ensure tau1 is odd to avoid edge effects in the filter
    tau1 = tau1 + 1;
end
tau2 = floor(tau2 / (1/fps)); % Convert tau2 to frame units

baseline = zeros(1, size(traces, 2)); % Initialize the baseline array

h = fspecial('average', [1 tau1]); % Create the smoothing filter
for n = 1:size(traces, 1)   % For each neuron
    data = traces(n, :);
    data_ext = horzcat(repmat(mean(data(1:ceil(tau1/2))), 1, ceil(tau1/2)), data, repmat(mean(data(end-ceil(tau1/2):end)), 1, ceil(tau1/2)));
    data_smooth_ext = filter2(h, data_ext); 
    data_smooth = data_smooth_ext(ceil(tau1/2)+1:end-ceil(tau1/2));   
    
    for t = 1:size(traces, 2) % For each time point
        if t < tau2
            F0 = data_smooth(1:t);     
        else
            F0 = data_smooth(t-tau2+1:t);
        end
        
        baseline(1, t) = mean(F0);
    end

    baseline_all(n, :) = baseline;

end

end
