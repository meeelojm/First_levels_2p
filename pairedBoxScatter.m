function pairedBoxScatter(temps, epochLabels)
%PAIREDBOXSCATTER Create boxplots with paired scatter lines
%   pairedBoxScatter(temps) plots boxplots and paired data
%   for a 2×N matrix temps, using default labels '1' and '2'.
%
%   pairedBoxScatter(temps, epochLabels) allows custom XTick labels:
%     epochLabels should be a 1×2 cell array of strings, e.g.
%       {'Epoch 1','Epoch 2'}
%
%   Example:
%     temps = randn(2,10) + [0;1];
%     pairedBoxScatter(temps, {'Before','After'});

    % Validate input
    if nargin < 2 || isempty(epochLabels)
        epochLabels = {'1','2'};
    elseif ~iscell(epochLabels) || numel(epochLabels)~=2
        error('epochLabels must be a 1×2 cell array of strings');
    end

    [rows, n] = size(temps);
    if rows~=2
        error('temps must be a 2×N matrix');
    end
    x = [1 2];

    % Create figure
    hold on;

    % 1) Boxplots
    data  = [temps(1,:) temps(2,:)];
    group = [ones(1,n), 2*ones(1,n)];
    boxplot(data, group, ...
        'Positions', x, ...    % force boxes at x=1 and x=2
        'Widths',    0.2, ...   % narrow boxes
        'Colors',    'k', ...   % black outlines
        'Symbol',    '');       % no outlier markers

    % 2) Paired lines and scatter
    for i = 1:n
        plot(x, temps(:,i), '-', 'Color', 0.85*[1 1 1], 'LineWidth', 0.2);
    end
    scatter(repmat(x(1),n,1), temps(1,:)',  50, 'b', 'filled');
    scatter(repmat(x(2),n,1), temps(2,:)', 50, 'r', 'filled');
    
    % Compute per-epoch means
    epochMeans = mean(temps,2);
    plot(x, epochMeans, '-ok', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'k');

    % 3) Format axes
    xlim([0.5 2.5]);
    set(gca, 'XTick', x, 'XTickLabel', epochLabels);
    ylabel('Value');
    box on; hold off;
end
