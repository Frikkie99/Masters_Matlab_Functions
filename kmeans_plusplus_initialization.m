function centroids = kmeans_plusplus_initialization(X, K)
    % X: Input data matrix (each row represents a data point)
    % K: Number of centroids
    
    % Step 1: Randomly select the first centroid
    rng(1);
    centroids = zeros(K, size(X, 2));
    centroids(1, :) = X(randi(size(X, 1)), :);
    
    % Step 2: Calculate distances and choose next centroids
    for k = 2:K
        distances = pdist2(X, centroids(1:k-1, :), 'squaredeuclidean');
        minDistances = min(distances, [], 2);
        probabilities = minDistances / sum(minDistances);
        
        % Choose the next centroid based on probability
        nextCentroidIdx = randsample(size(X, 1), 1, true, probabilities);
        centroids(k, :) = X(nextCentroidIdx, :);
    end
end