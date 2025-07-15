function [crossres] = crosscorss_classifier3(signals1, signals2, NumLags, period, mode)
% CROSSCORSS_CLASSIFIER3 Computes and clusters cross-correlograms with asymmetry analysis.
%
% INPUTS:
%   signals1 - matrix (neurons x time)
%   signals2 - matrix (neurons x time)
%   NumLags  - number of lags to compute in cross-correlograms
%   period   - indices to extract from cross-correlograms for clustering
%   mode     - 1 = z-scored clustering, 2 = raw clustering
%
% OUTPUT:
%   crossres struct with fields:
%     .xcfs              - all raw cross-correlograms
%     .zs                - matrix used for clustering
%     .idxcfs            - [signal1_index, signal2_index, cluster]
%     .outperm           - sorting order from dendrogram
%     .clusters          - cell array of raw crosscorrelograms per cluster
%     .clusterz          - cell array of z-scored/used crosscorrelograms per cluster
%     .clustersids       - signal indices in each cluster
%     .asymmetry         - asymmetry index per correlogram
%     .global_asymmetry  - mean asymmetry across all pairs

    %------------------- Compute Cross-Correlograms -------------------%
    xcfs = [];
    idxcfs = [];
    
    % Demean signals
    signals1 = signals1 - mean(signals1, 2);
    signals2 = signals2 - mean(signals2, 2);
    
    for i = 1:size(signals1, 1)
        for j = 1:size(signals2, 1)
            [xcf, ~] = crosscorr(signals1(i,:), signals2(j,:), 'NumLags', NumLags);
            xcfs = [xcfs; xcf]; %#ok<AGROW>
            idxcfs = [idxcfs; i, j]; %#ok<AGROW>
        end
    end
    
    %------------------- Choose Feature Matrix -------------------%
    if mode == 1
        zs = zscore(xcfs(:, period), 0, 2);
        Z = linkage(zs, 'average', 'correlation');
    elseif mode == 2
        zs = xcfs(:, period);
        Z = linkage(zs, 'ward');
    else
        error('Mode must be 1 (z-scored) or 2 (raw).');
    end
    
    %------------------- Clustering and Visualization -------------------%
    figure('Renderer', 'painters', 'Position', [0 0 1100 300]);
    
    % Dendrogram
    subplot(1, 4, 1);
    [~, ~, outperm] = dendrogram(Z, 0, 'Orientation', 'left', 'ColorThreshold', 0.7);
    title('Hierarchical Clustering Dendrogram');
    xlabel('Distance'); ylabel('Cross-correlogram Index');
    
    % Heatmap of sorted correlograms
    sorted_crosscorrelograms = zs(outperm, :);
    subplot(1, 4, [2 3]);
    imagesc(sorted_crosscorrelograms); camroll(180); axis xy;
    if mode == 1
        title('Sorted Z-scored Cross-correlograms');
        caxis([-3 3]);
    else
        title('Sorted Raw Cross-correlograms');
        caxis([-0.1 0.25]);
    end 
    xlabel('Frames'); ylabel('Sorted Index');
    xline(find(period == NumLags + 1), 'k');  % vertical line at zero lag if included
    colormap(coolwarm); colorbar;
    
    % Similarity matrix
    subplot(1, 4, 4);
    temp = corrcoef(sorted_crosscorrelograms');
    imagesc(temp);  colormap(coolwarm); camroll(180); axis xy;
    title('Similarity (Correlation) Matrix');
    
    % Cluster assignment
    maxnumcl = 5;
    T = cluster(Z, 'MaxClust', maxnumcl);
    idxcfs(:, 3) = T;
    
    % Cluster organization
    clusters = cell(maxnumcl, 1);
    clusterz = cell(maxnumcl, 1);
    clustersids = cell(maxnumcl, 1);
    
    for k = 1:maxnumcl
        clusterIndices = (T == k);
        clusterz{k} = zs(clusterIndices, :);
        clusters{k} = xcfs(clusterIndices, :);
        clustersids{k} = idxcfs(clusterIndices, :);
    end
    
    %------------------- Asymmetry Analysis -------------------%
    allMid = ceil(size(xcfs, 2) / 2);
    neg_part = xcfs(:, 1:allMid-1);
    pos_part = xcfs(:, allMid+1:end);
    asymmetry = mean(pos_part, 2) - mean(neg_part, 2);
    global_asymmetry = mean(asymmetry);
    
    %------------------- Output Struct -------------------%
    crossres.xcfs = xcfs;
    crossres.zs = zs;
    crossres.idxcfs = idxcfs(outperm, :);
    crossres.outperm = outperm;
    crossres.clusters = clusters;
    crossres.clusterz = clusterz;
    crossres.clustersids = clustersids;
    crossres.asymmetry = asymmetry;
    crossres.global_asymmetry = global_asymmetry;
end
