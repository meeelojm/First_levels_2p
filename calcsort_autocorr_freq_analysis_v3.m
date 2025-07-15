function [acs, dominantFrequencies, Pxxs, Fss, finalSortedIndices,isoinds] = calcsort_autocorr_freq_analysis_v3(signalmat, numlags, fs, windowSize, maxpx, minPeakThreshold,freqRange)
    % calcsort_autocorr - Calculate and sort autocorrelations for each neuron
    %
    % Syntax: [acs, outperm, allIntervals, highAutocorrLowFreqProportion, cluster_indices, cluster_colors] = calcsort_autocorr_freq_analysis(signalmat)
    %
    % Inputs:
    %    signalmat - Matrix of neuronal activity signals (neurons x time)
    %    numlags - Number of lags for the autocorrelation calculation
    %    fs - Sampling frequency
    %    cthr - Threshold for clustering
    %    windowSize - Window size for PSD calculation
    %    maxpx - Maximum PSD bins to display
    %
    % Outputs:
    %    Visualizations of resorted autocorrelations and PSDs

    % Calculate autocorrelations for each sequence
%     figure;
%     make_it_tight = true;
%     subplot = @(m, n, p) subtightplot(m, n, p, [0.01 0.05], [0.1 0.1], [0.1 0.01]);
%     if ~make_it_tight, clear subplot; end
    
    % Initialize an empty array to store autocorrelations
    clear acs;
    
    % Loop through each neuron (row) in the signal matrix
    for i = 1:size(signalmat, 1)
        % Calculate the autocorrelation for each neuron with the specified lag
        ac = autocorr(signalmat(i, :), NumLags=numlags);
        acs(i, :) = ac;  % Store the autocorrelation in the 'acs' array
    end
    
    % Calculate peak counts and dominant frequencies for each sequence
    numNeurons = size(acs, 1);
    peakCounts = zeros(numNeurons, 1);
    dominantFrequencies = zeros(numNeurons, 1);
    Pxxs = [];
    Fss = [];
%     minPeakThreshold = 1; % Minimum number of peaks required
    
    for i = 1:numNeurons
        % Count peaks in the autocorrelation sequence
        [~, locations] = findpeaks(acs(i, :), 'MinPeakProminence', 0.1);
        peakCounts(i) = numel(locations);
        
        % Calculate the PSD and find the dominant frequency
        [Pxx, F] = pwelch(acs(i, :), windowSize, windowSize/2, windowSize, fs);
        %[Pxx, F] = pwelch(acs(i, :), windowSize, windowSize/2, [], fs);
        Pxxs(:, i) = Pxx;
        Fss(:, i) = F;
        [~, maxIdx] = max(Pxx);
        dominantFrequencies(i) = F(maxIdx);
    end
    % Separate sequences based on the peak count threshold
    sequencesWithPeaks = find(peakCounts >= minPeakThreshold);
    sequencesWithoutPeaks = find(peakCounts < minPeakThreshold);
    
    % Sort sequences with sufficient peaks by dominant frequencies and then peak count
    [~, orderWithPeaks] = sortrows([dominantFrequencies(sequencesWithPeaks), peakCounts(sequencesWithPeaks)]);
    sortedIndicesWithPeaks = sequencesWithPeaks(orderWithPeaks);
    
    % Sort sequences without sufficient peaks by dominant frequencies
    [~, orderWithoutPeaks] = sort(dominantFrequencies(sequencesWithoutPeaks));
    sortedIndicesWithoutPeaks = sequencesWithoutPeaks(orderWithoutPeaks);
    
    % Combine the sorted indices
    finalSortedIndices = [sortedIndicesWithPeaks; sortedIndicesWithoutPeaks];
    resortedAcs = acs(finalSortedIndices, :);
    resortedPsdMatrix = Pxxs(:, finalSortedIndices);
    resortedDominantFrequencies = dominantFrequencies(finalSortedIndices);

    figure('Renderer', 'painters', 'Position', [0 0 1352 1062]);
    
    % First subplot with 'gray' colormap (activity sorted by dominant frequencies)
    subplot(1,8,[1 3])
    %imagesc(dff_filtered(finalSortedIndices, start_frame:end_frame)); 
    imagesc(signalmat(finalSortedIndices,:))
    c = colorbar;  % Get the colorbar handle
    colormap(gca, flipud(bone));  % Apply 'gray' colormap only to this subplot
    caxis([-10 100]);
    title('Activity sorted by dominant frequencies');    
    xtick = get(gca, 'XTick');  % Get current x-tick values in frames
    xticks_in_seconds = xtick / fs;
    set(gca, 'XTickLabel', xticks_in_seconds);
    xlabel('Time (seconds)');
    ylabel('Neurons');
    xtickangle(90);
    c.Label.String = '%DF/F';  % Add label to the colorbar
    box off
    % Second subplot with 'turbo' colormap (autocorrelations sorted by dominant frequencies)
    subplot(1,8,[4 5])
    imagesc(resortedAcs); 
    c = colorbar;  % Get the colorbar handle
    colormap(gca, 'turbo');  % Apply 'turbo' colormap only to this subplot
    caxis([0 0.5]);
    title('Autocorrelations sorted by dominant frequencies');
    xtick = get(gca, 'XTick');  % Get current x-tick values in frames
    xticks_in_seconds = xtick / fs;
    xtickangle(90);
    set(gca, 'XTickLabel', xticks_in_seconds);
    xlabel('Time (seconds)');
    ylabel('Neurons');
    c.Label.String = 'Autocorrelation value';  % Add label to the colorbar
    box off
    % Third subplot with 'turbo' colormap (PSD sorted by dominant frequencies)
    subplot(1,8,[6 8])
    imagesc(resortedPsdMatrix(1:maxpx,:)');  % Plot only the first 10 frequency bins (adjust as necessary)
    ax = gca;
    temp = str2double(ax.XTickLabel)';  % Get the current x-tick labels
    xticks(temp);  % Set xticks to match the calculated frequencies
    set(gca, 'XTickLabel', round(F(temp), 3));  % Round frequency bins and apply
    caxis([0 0.5]); 
    c = colorbar;  % Get the colorbar handle
    colormap(gca, 'turbo');  % Apply 'turbo' colormap only to this subplot
    title('PSDs sorted by dominant frequencies');
    xtickangle(90);
    ylabel('Neurons');
    xlabel('Frequencies');
    c.Label.String = 'PSD';  % Add label to the colorbar
    box off
%%
    % Define the specific frequency range for isolating the neurons
%     freqRange = [0.03, 0.09];  % Adjust based on observation of circled frequencies
    
    % Step 3: Isolate cells based on the defined criteria
    isolated_indices = find(dominantFrequencies >= freqRange(1) & dominantFrequencies <= freqRange(2));
    
    % Output the isolated indices
    disp(['Number of cells isolated: ', num2str(length(isolated_indices))]);
    
    % Sort the isolated indices based on dominant frequencies (optional)
    [~, sorted_isolated_indices] = sort(dominantFrequencies(isolated_indices));
    

     isoinds = isolated_indices(sorted_isolated_indices);
end
