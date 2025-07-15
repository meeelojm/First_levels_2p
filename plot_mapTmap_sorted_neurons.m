function plot_mapTmap_sorted_neurons(isort1, dff_data, positions, duration, colmap, stims_frames)
    % PLOT_MAPTMAP_SORTED_NEURONS Visualizes neurons sorted by mapTmap output
    %
    % Inputs:
    %   isort1       - sorted neuron indices from mapTmap
    %   dff_data     - dF/F matrix (neurons x time)
    %   positions    - 3D positions of neurons (neurons x 3)
    %   duration     - duration in frames (e.g., stimulus duration)
    %   colmap       - 1 for bone, 2 for turbo
    %   stims_frames - vector of stim frame indices (can be empty)
    
    % Normalize sorting rank into colormap values
    nNeurons = size(positions, 1);
    colors = zeros(nNeurons, 1);
    colors(isort1) = 1:nNeurons;
    colors = colors / max(colors);  % Normalize for color mapping
    
    % === Combined Raster Plot with Sort Rank Colorbar ===
    figure;
    
    % 1. Sorting index colorbar
    subplot(1, 12, 1);
    imagesc(colors(isort1));
    colormap(gca, jet);
    set(gca, 'YDir', 'normal', 'XTick', [], 'YTick', []);
    title('Sort Rank');
    freezeColors;
    
    % 2. Raster plot of sorted dF/F traces
    subplot(1, 12, 2:12);
    imagesc(dff_data);
    caxis([-20 100]);
    
    % Set colormap
    if colmap == 1
        colormap(flipud(bone)); colorbar
    elseif colmap == 2
        colormap turbo; colorbar
    end
    
    % Stimulus overlay logic
    if ~isempty(stims_frames) && isnumeric(stims_frames)
        disp('Plotting stimulus markers...');
        hold on
        % Try to infer offset from spacing
        if length(stims_frames) > 1
            increase = diff(stims_frames);
            offset_guess = increase(1);
        else
            offset_guess = duration;
        end
        frame_offset = stims_frames(1) - offset_guess;
    
        for s = 1:length(stims_frames)
            xline(stims_frames(s) - frame_offset + 1, 'r', 'LineWidth', 1.5);
        end
    else
        xline(duration*2, 'w', 'LineWidth', 3);  % fallback for single stim
    end
    
    title('All Responding Cells');
    xlabel('time (frames)');
    set(gca, 'YDir', 'normal', 'YTick', []);
    freezeColors;
    
    % === 3D scatter plot colored by sort rank ===
    figure;
    scatter3(positions(:,1), positions(:,2), positions(:,3)*10, 60, colors, 'filled');
    colormap(jet);
    colorbar;
    title('Neuron Positions Colored by mapTmap Sorting');
    xlabel('X'); ylabel('Y'); zlabel('Z');
    grid on;
    % view([0 90])
    set(gca, 'ZDir', 'reverse');
end
