function [C_idx_tc] = TCKMEANS2(data, K, tol, tau)
% Initialize the centroids randomly
centroids = kmeans_plusplus_initialization(data(:,2:end),K);
diff = 1;

while diff > tol
    a = centroids;
    % Assign each data point to the closest centroid
    distances = (pdist2(data(:,2:end), centroids).^2);
    [~, idx] = min(distances, [], 2);

    bf = idx;

    for i = 1 : K
        for j = 1 : length(data)
            alpha = (sum(exp(-tau^2*(pdist2(data(j,1),data(idx == i,1))).^2)) - 1)./(sum(exp(-tau^2*(pdist2(data(j,1),data(:,1))).^2)) - 1);
            distances(j,i) = distances(j,i)*(1 - alpha);
        end
    end
    
    [~, idx] = min(distances, [], 2);
    af = idx;

    change = sum(abs(bf - af));
    
    % Update the centroids
    for j = 1:K
        centroids(j, :) = mean(data(idx == j, 2:end));
    end
    diff = sum(sum(abs(centroids - a)));
end
[~, C_idx_tc] = min(distances, [], 2);
tau
end