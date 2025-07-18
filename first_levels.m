%% first_level_analysis.m
% Basic first-level analysis pipeline for Yagnur neural data
% This script performs multiple analyses to characterize neural activity patterns:
% 1) Correlation vs. Distance: tests if neurons that are closer in space have more similar activity
% 2) Hopkins Statistic: quantifies clusterability of activity in high-dimensional space
% 3) Elbow Method: estimates optimal cluster number for k-means
% 4) k-means Clustering: groups neurons by activity patterns
% 5) PCA Sorting: assesses spatial organization along principal functional axes
% 6) Cross-Correlation Classification: examines temporal coordination between neuron subsets
% 7) Frequency Analysis: identifies dominant oscillatory components in activity

baseFolder = 'C:\Users\emiliajm\OneDrive - NTNU\Documents\test';
defaultFPS = 1.5;
[fishPaths, fpss] = discover_folders_ini(baseFolder, defaultFPS);
fprintf('Found %d experiments.\n',numel(fishPaths));

% --- Section 1: Load data and compute ΔF/F ---
% Rationale: Converting raw fluorescence to ΔF/F normalizes signals and removes baseline drift,
% ensuring comparisons across neurons and sessions are meaningful.
% Load tracing and position data
for f=1:numel(fishPaths)
    S              = load(fullfile(fishPaths{f}, 'results.mat'));  % contains .trace (raw fluorescence) and .position (XYZ)
    rawTraces     = double(S.trace);                              % matrix: neurons × time frames
    positionsOrig = S.position;                                 % matrix: neurons × 3 (X, Y, Z)
    [DeconvMat,FiltMat]=CaDeconvNew(rawTraces,3,30,0,100,0,200,1/31,0.1,0);
    % Baseline parameters for fluorescence
    tauFast       = 1;      % fast time constant (seconds)
    tauSlow       = 19.4;   % slow time constant adjusted for 30.9 Hz window
    % Compute baseline fluorescence with adjustable time constants
    baselineFluorescence = baseline_time_adjustable(tauFast, tauSlow, FiltMat, fpss(f));
    
    % Calculate ΔF/F (%)
    dffMovWindow = ((FiltMat - baselineFluorescence) ./ baselineFluorescence) * 100;
    
    % Remove neurons with NaNs or all-zero traces
    validRowIndices = find(all(~isnan(dffMovWindow), 2));
    filteredDff      = dffMovWindow(validRowIndices, :);
    dffMovWindow         = filteredDff(any(filteredDff, 2), :);  % final ΔF/F matrix
    positionsOrig = positionsOrig(validRowIndices, :);

    % Remove flyback plane
    temp = (unique(positionsOrig(:,5)));
    temp = max([temp]);
    idx_nonflyback = find(positionsOrig(:,5)~=temp);
    positionsOrig = positionsOrig(idx_nonflyback,:);
    dffMovWindow  = dffMovWindow(idx_nonflyback,:);
    % Define analysis window (frames)
    % make this 10 mins independent of the fps
    startFrame = 1;  % insert initial frame, usually skip the first 100 frames
    endFrame   = startFrame + round(20*60*fpss(f)) - 1;
end



%% Section 2: Correlation vs Distance (real vs shuffled) 150 microns
% Rationale: Spatial correlation analysis tests whether proximal neurons share more similar activity.
% By comparing real vs shuffled positions, we control for global correlations unrelated to space.
nFiles               = numel(fishPaths);
steps                = floor(linspace(0,200,51));
nBins                = length(steps) - 1;
subsampleSize        = size(dffMovWindow,1);                % or round(0.15 * nNeurons)
allCorrDistReal      = zeros(2*nFiles, nBins);
allCorrDistShuffled  = zeros(2*nFiles, nBins);

for iFile = 1:nFiles
    dffFull           = dffMovWindow;
    posFull           = positionsOrig;
    
    % --- pick a random subsample of neurons ---
    nNeurons          = size(dffFull,1);
    idxSample         = randperm(nNeurons, min(subsampleSize,nNeurons));
    dff               = dffFull(idxSample, startFrame:endFrame);
    pos               = posFull(idxSample, :);
    tic()
    % --- Real spatial correlation (2×nBins) ---
    [~, corrReal, ~, ~] = corr_vs_distance_AO( ...
        dff, pos, 1, 1, [], steps);
    allCorrDistReal(2*iFile-1 : 2*iFile, :) = corrReal;
    toc()
    tic()
    % --- Shuffled spatial correlation (2×nBins) ---
    posShuff      = pos(randperm(numel(idxSample)), :);
    dffShuff      = dff(:, randperm(size(dff,2)));
    [~, corrShuff, ~, ~] = corr_vs_distance_AO( ...
        dff, posShuff, 1, 1, [], steps);
    allCorrDistShuffled(2*iFile-1 : 2*iFile, :) = corrShuff;
    toc()
end

figure;plot(mean(allCorrDistReal(:,:)))
hold on;plot(mean(allCorrDistShuffled(:,:)))

%% Section 3: Hopkins statistic for real vs. uniform data
% Rationale: Hopkins statistic H > 0.5 indicates cluster tendency; H ≈ 0.5 is random.
hopkinsReal    = [];
hopkinsUniform = [];
for iFile = 1:nFiles
    m = round(0.25 * size(dffMovWindow,1));  % 10% of neurons
    % Real data Hopkins
    hopkinsReal(iFile) = hopkins(dffMovWindow, m);
    % Uniform random control data
    minVal = min(dffMovWindow(:));
    maxVal = max(dffMovWindow(:));
    uniformDff = minVal + (maxVal - minVal) * rand(size(dffMovWindow));
    hopkinsUniform(iFile) = hopkins(uniformDff, m);
end
display(hopkinsReal)
display(hopkinsUniform)

%% Section 4: Elbow method to estimate cluster count
% Rationale: Identify the k where adding clusters yields diminishing returns in within-cluster variance.
% Evaluate clustering up to 15 clusters on the clean ΔF/F data
dataWindow = dffMovWindow(:, startFrame:endFrame);
[maxClusters] = 25;  % maximum clusters to test
[sumDistAll, randSumDistAll, shuffCompSumDist] = evalclust_EY(dataWindow, maxClusters);

figure;
% 1) Sum of within-cluster distances vs. k
subplot(1,4,1);
plot(sumDistAll', 'k'); hold on;
plot(mean(sumDistAll,1), 'k', 'LineWidth',4);
plot(randSumDistAll', 'r');
plot(mean(randSumDistAll,1), 'r', 'LineWidth',4);
title('Within-cluster sumD');

% 2) Difference from random control
subplot(1,4,2);
plot(mean(randSumDistAll - sumDistAll,1), 'k', 'LineWidth',4);
title('Δ sumD vs random');

% 3) Slope of decay
subplot(1,4,3);
plot(mean(diff(sumDistAll,1,2),1), 'k', 'LineWidth',4);
title('Slope of sumD');

% 4) Difference in slopes
subplot(1,4,4);
plot(mean(diff(randSumDistAll,1,2) - diff(sumDistAll,1,2),1), 'k', 'LineWidth',4);
title('Δ slope vs random');


%% Section 5: k-means clustering and visualization
% Rationale: Partition neurons into k functional groups to reveal common response patterns.
% Define maximum clusters for analysis (user-defined)
maxClusters = 12%% INSERT number of clusters;

for fIdx = 1:numel(fishPaths)
    sliceData = double(dffMovWindow(:, startFrame:endFrame));
    
    % Generate colormap
    cmap = brewermap(maxClusters, 'Paired');
    % Run k-means with correlation distance
    defaultOpts = {'Distance','correlation','Replicates',50,'Start','sample'};
    [clusterIdx, clusterCenters, sumD, D] = kmeans(sliceData, maxClusters, defaultOpts{:});
    
    % 3D scatter: neurons colored by cluster assignment
    figure;
    scatter3(positionsOrig(:,1), positionsOrig(:,2), positionsOrig(:,3), 45, clusterIdx, 'filled');
    title('k-means Clusters'); colormap(cmap);
    view(-90, 90);
    c = colorbar('Direction','reverse');
    c.Ticks = 1:maxClusters;
    
    % Compute and plot cluster mean traces
    neuronTraces = cell(maxClusters,1);
    clusterMeans  = zeros(maxClusters, size(sliceData,2));
    for k = 1:maxClusters
        members = sliceData(clusterIdx==k, :);
        members = members(any(members,2), :);  % drop all-zero rows
        neuronTraces{k} = members;
        clusterMeans(k,:)  = mean(members,1);
    end
    
    % Plot each cluster mean
    figure;
    for k = 1:maxClusters
        subplot(maxClusters,1,k);
        plot(clusterMeans(k,:), 'LineWidth', 1.5, 'Color', cmap(k,:));
        ylabel(sprintf('Cluster %d', k));
        if k==1, title('Cluster Mean Traces'); end
        if k<maxClusters, set(gca, 'XTickLabel', []); end
        if k==maxClusters, xlabel('Frames'); end
    end

    figure;
    t = tiledlayout(ceil(sqrt(maxClusters)), ...
                   ceil(sqrt(maxClusters)), ...
                   'TileSpacing','compact','Padding','compact');
    sgtitle(sprintf('Fish %d: individual cluster views', fIdx));
    
    for k = 1:maxClusters
        nexttile;
        % Plot *all* neurons in light gray
        scatter3(pos(:,1), pos(:,2), pos(:,3), ...
                 20, 'MarkerEdgeColor',[.6 .6 .6], ...
                 'MarkerFaceColor',[.6 .6 .6], ...
                 'MarkerFaceAlpha',0.05);
        hold on;
        % Highlight the k-th cluster in its color
        idx_k = (clusterIdx == k);
        scatter3(pos(idx_k,1), pos(idx_k,2), pos(idx_k,3), ...
                 40, cmap(k,:), 'filled', ...
                 'MarkerFaceAlpha',0.8, ...
                 'MarkerEdgeColor','k');
        view(-90,90);
        axis tight off;
        title(sprintf('Cluster %d (%d cells)', k, sum(idx_k)));
    end
end


%% Section 6: Principal Component (PC) sorting
% Rationale: Tests if functional ordering (PC angle) aligns with anatomical A-P or L-R axes more than by chance.
nNeurons = size(dffMovWindow,1);
realCorrAP = zeros(numel(fishPaths),1);
realCorrLR = zeros(numel(fishPaths),1);
randCorrAP = zeros(numel(fishPaths),1);
randCorrLR = zeros(numel(fishPaths),1);

for fIdx = 1:numel(fishPaths)
    % Extract windowed data and positions
    dffWindow = dffMovWindow(:, startFrame:endFrame);
    pos = positionsOrig;
    
    % Compute PCs and real sorting
    coeff = pca(dffWindow', 'NumComponents', 2);
    angles = atan2(coeff(:,2), coeff(:,1));
    [~, sortOrder] = sort(angles);
    realCorrAP(fIdx) = corr((1:length(sortOrder))', pos(sortOrder,2), 'Type','Spearman');
    realCorrLR(fIdx) = corr((1:length(sortOrder))', pos(sortOrder,1), 'Type','Spearman');
    
    % Random control: random sort order
    randOrder = randperm(nNeurons);
    randCorrAP(fIdx) = corr((1:length(randOrder))', pos(randOrder,2), 'Type','Spearman');
    randCorrLR(fIdx) = corr((1:length(randOrder))', pos(randOrder,1), 'Type','Spearman');
end

% Plot sorted data
plot_mapTmap_sorted_neurons(sortOrder, dffMovWindow, pos, 0, 1, []);
view(-90, 90);
% Plot paired comparisons: Real vs Random for A-P and L-R
figure;
subplot(1,3,1);
pairedBoxScatter([abs(realCorrAP), abs(randCorrAP)]', {'Real','Random'});
title('A-P Spatial Organization');
subplot(1,3,2);
pairedBoxScatter([abs(realCorrLR), abs(randCorrLR)]', {'Real','Random'});
title('L-R Spatial Organization');
subplot(1,3,3);
pairedBoxScatter([abs(realCorrAP), abs(realCorrLR)]', {'A-P','L-R'});
title('A-P vs L-R');



%% Section 7: Cross-correlation classifier
% Rationale: Examines temporal coordination between early vs late firing neurons via cross-correlation peaks.
% Compute cross-correlation between first/last subsets of neurons
for f = 1:numel(fishPaths)
    fiveMinFrames = 300 * fpss(f);
    numLags      = floor(fiveMinFrames/4);
    windowWidth  = floor(fiveMinFrames/8);
    periodRange  = (numLags+1-windowWidth) : (numLags+1+windowWidth);

    nCells = size(dffMovWindow, 1);
    firstN  = round(nCells * 0.10); % First 10%
    lastN   = round(nCells * 0.10); % Last 10%

    % Indices for first and last 10%
    idxFirst = 1:firstN;
    idxLast  = (nCells-lastN+1):nCells;

    crossResults = crosscorss_classifier3(dffMovWindow(idxFirst, startFrame:endFrame), ...
                                          dffMovWindow(idxLast,  startFrame:endFrame), ...
                                          numLags, periodRange, 1);
end


%% Section 8: Frequency analysis of ΔF/F signals
% Rationale: Identifies dominant oscillatory frequencies and spatial distribution of low vs high–freq neurons.
% Convert original frame index to current fps index
for f = 1:numel(fishPaths)

    origStartFrame = 200; origFps = 2.5;
    startSec       = origStartFrame / origFps;
    startFrameFreq = round(startSec * fpss(f));
    endFrameFreq   = startFrameFreq + round(240 * fpss(f));  % 4-minute window
    numLagsFreq    = round(120 * fpss(f));                 % 2-minute window
    
    % Frequency analysis parameters
    windowSize         = 1024; %adjust according to the length of the segment 
    maxPeakCount       = 15;
    minPeakHeight      = 2;
    frequencyRange     = [0.01, 0.03];
    
    [autocorrs, domFreqs, powerSpecs, freqBins, isoIndices] = calcsort_autocorr_freq_analysis_v3(...
            dffMovWindow(:, startFrameFreq:endFrameFreq), numLagsFreq, fpss(f), ...
            windowSize, maxPeakCount, minPeakHeight, frequencyRange);
    
    % Categorize neurons by low vs. high frequency
    lowFreqMask  = domFreqs >= 0.01 & domFreqs <= 0.1;
    lowFreqValues= domFreqs(lowFreqMask);
    highFreqMask = ~lowFreqMask;
    highFreqValues= domFreqs(highFreqMask);
    
    % 3D scatter colored by low-frequency values, high-frequency in gray
    figure;
    % View 1\subplot(1,2,1);
    scatter3(neuronPositions(highFreqMask,1), neuronPositions(highFreqMask,2), neuronPositions(highFreqMask,3), ...
        43, 'ok', 'filled', 'MarkerFaceAlpha',0.1); hold on;
    scatter3(neuronPositions(lowFreqMask,1), neuronPositions(lowFreqMask,2), neuronPositions(lowFreqMask,3), ...
        43, lowFreqValues, 'o','filled');
    xlabel('X'); ylabel('Y'); zlabel('Z');
    colormap turbo; colorbar; caxis([0 0.1]);
    set(gca, 'ZDir','reverse'); view(90,90);
    % View 2\subplot(1,2,2);
    scatter3(neuronPositions(highFreqMask,1), neuronPositions(highFreqMask,2), neuronPositions(highFreqMask,3), ...
        43,'ok','filled','MarkerFaceAlpha',0.1); hold on;
    scatter3(neuronPositions(lowFreqMask,1), neuronPositions(lowFreqMask,2), neuronPositions(lowFreqMask,3), ...
        43, lowFreqValues, 'o','filled');
    xlabel('X'); ylabel('Y'); zlabel('Z');
    colormap turbo; colorbar; caxis([0 0.1]);
    set(gca, 'ZDir','normal','YDir','reverse','XDir','reverse'); view(198,42);
    view(-90, 90);
    title('Neurons according to frequencies')
end