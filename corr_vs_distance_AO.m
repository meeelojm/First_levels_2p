function [all_corr_dist, all_corr_dist_pos, all_corr_dist_neg, all_dist_corr_pairs, steps] = corr_vs_distance_AO(dff_for_period, cells_positions, multiplane, please_split, varargin)
%corr_vs_distance_AO - Calculates correlation vs distance
%   Author: Anna Maria Ostenrath
%
%   Syntax:
%       [all_corr_dist, all_corr_dist_pos, all_corr_dist_neg] = corr_vs_distance_AO(dff_for_period, cells_positions, 1)
%       [all_corr_dist, all_corr_dist_pos, all_corr_dist_neg] = corr_vs_distance_AO(dff_for_period, cells_positions, 0, hemisphereindex)
%   Description:
%       function() - description
%    
%   Inputs:
%       dff_for_period - dff of the region/all the cells of interest in
%       form neuron x time
%       cells_positions - position of the cells in 2 or 3d (uses the second
%       dimension for splitting hemispheres)
%       multiplane - 1 or 0, so you can calculate in 3D (1) or 2D(0) 
%       please_split - 1 or 0, automatically split the brain in two
%       hemispheres using the y pos (kmeans)
%       varargin - extra arguements: 1) hemisphere index if you have
%       already split the brain 2) the steps you want to bin the distance
%       into default:  steps = floor(linspace(0,100,51));
%       input2 - Description
%       input3 - Description
%
%   Outputs:
%       all_corr_dist - list of correlations vs distance averaged over both
%       pos and neg correlations
%       all_corr_dist_pos - only the correlations >=0 
%       all_corr_dist_neg - only correlations <0 
%       steps - The steps that the distance was binned into
%
%   Examples: 
%       Line 1 of example
%       Line 2 of example
%       Line 3 of example
%
%   Other m-files required: none
%   Subfunctions: none
%   MAT-files required: none
%
%   See also: OTHER_FUNCTION_NAME1,  OTHER_FUNCTION_NAME2
%   Author: Anna Maria Ostenrath 
%   Date : Nov 2022

% First we need to split the brain/brain region into two sides
if please_split
    [pos_ind,~,~,~] = kmeans(cells_positions(:,2),2,'Replicates', 50,'Distance','sqeuclidean'); 

    [max_y, max_ind] = max(cells_positions(:,2)); 
    %     side = IDX_b(max_ind);
    upper_region = pos_ind(max_ind); 
    [min_y, min_ind] = min(cells_positions(:,2)); 
    lower_region = pos_ind(min_ind);
else 
    pos_ind = varargin{1}; % or use the hemisphere index provided 
end 
all_corr_dist = [];
all_corr_dist_pos  = [];
all_corr_dist_neg = [];

% now we make the edges for binning the distance or use the list provided
if size(varargin,2) > 1
    steps = varargin{2}; 
else
    steps = floor(linspace(0,100,26));
end

% now we loop over the two sides of the brain to calculated the corr vs
% distance
for side = 1:2
    % selecting only the neurons of the right size
    current_position = cells_positions(find(pos_ind == side), :);
    current_dff = dff_for_period(find(pos_ind == side),:);
    mat_dist=NaN(size(current_dff,1),size(current_dff,1));
    if multiplane
        for cell_one=1:size(current_position,1)
            for cell_two=1:size(current_position,1)
                if cell_one<cell_two
                    mat_dist(cell_one,cell_two)=sqrt((current_position(cell_one,1)-current_position(cell_two,1))^2 +(current_position(cell_one,2)-current_position(cell_two,2))^2+(current_position(cell_one,3)-current_position(cell_two,3))^2); % recalculate x and y beforehand
                end
            end
         end

    else
         for cell_one=1:size(current_position,1)
            for cell_two=1:size(current_position,1)
                if cell_one<cell_two
                    mat_dist(cell_one,cell_two)=sqrt((current_position(cell_one,1)-current_position(cell_two,1))^2 +(current_position(cell_one,2)-current_position(cell_two,2))^2); % recalculate x and y beforehand
                end
            end
         end
    end
    % now corr vs distance while I am at it 
     % here with hab
   [mat_act ,p_r]= corrcoef(current_dff(:,:)') ;% for studying right to right interactions 
    mat_act=triu(mat_act,1);

    new_dist = discretize(mat_dist,steps);

    mean_corr_values_general = zeros(size(steps,2),1); 
    mean_corr_only_pos = zeros(size(steps,2),1); 
    mean_corr_only_neg = zeros(size(steps,2),1); 
    for dist = 1:size(steps,2)
%         disp(steps(dist))
    % histogram(new_dist, mat_act_2); 
%         mean_corr_values_general(dist) = nanmean(mat_act_2(find(new_dist == steps(dist))));
%         all_pairs = mat_act_2(find(new_dist == steps(dist))); 
        mean_corr_values_general(dist) = nanmean(mat_act(find(new_dist == dist)));
        all_pairs = mat_act(find(new_dist == dist));
        posit_pair = [];
        neg_pair = []; 
        for pair = 1:length(all_pairs)
            if all_pairs(pair) < 0
                neg_pair = [neg_pair; all_pairs(pair)];
            elseif all_pairs(pair) >= 0
                posit_pair = [posit_pair; all_pairs(pair)];
            end
        end
        mean_corr_only_pos(dist) = nanmean(posit_pair);
        mean_corr_only_neg(dist) = nanmean(neg_pair);

    end
    all_corr_dist = [all_corr_dist; mean_corr_values_general(2:end)'];
    all_corr_dist_pos = [all_corr_dist_pos; mean_corr_only_pos(2:end)'];
    all_corr_dist_neg = [all_corr_dist_neg; mean_corr_only_neg(2:end)'];
    
    %% I just want a list for each cell pair now... 
    all_dist_corr_pairs = []; 
    for neuron1 = 1:size(current_position,1)
        for neuron2=1:size(current_position,1)
            if neuron1<neuron2
                all_dist_corr_pairs = [all_dist_corr_pairs; [mat_dist(neuron1, neuron2), mat_act(neuron1, neuron2), p_r(neuron1, neuron2)]];


            end

        end
    end
    
end

end