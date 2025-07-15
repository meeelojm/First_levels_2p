function H = hopkins(X, m)
    % X is your data (N x d), m is the number of samples (e.g., round(0.1*N))
    N = size(X,1);
    d = size(X,2);
    rand_idx = randperm(N, m);
    sample = X(rand_idx,:);
    
    % 1. Data-to-data (leave-one-out nearest neighbor)
    D_x = pdist2(sample, X);
    for i = 1:m
        D_x(i, rand_idx(i)) = inf; % exclude self
    end
    x = min(D_x, [], 2);

    % 2. Random-to-data
    minX = min(X);
    maxX = max(X);
    random_points = rand(m,d) .* (maxX - minX) + minX;
    D_y = pdist2(random_points, X);
    y = min(D_y, [], 2);

    % 3. Hopkins statistic
    H = sum(y) / (sum(x) + sum(y));
end
