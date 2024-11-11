function [C] = TCDBSCAN(X, eps, minPts, tau)
% N is the number of observations
[N, ~] = size(X);
% D is the distance matrix
D = zeros(N, N);
 
for i = 1:N
    for j = i+1:N
        % d is the euclidean distance between all data points minus the
        % exponential decay of their respective differences in time
        d = norm(X(i,2:3)-X(j,2:3)).*(1 - exp((-tau*norm(X(i,1) - X(j,1)))));
        % D is a symmetrical matrix
        D(i,j) = d;
        D(j,i) = d;
    end
end

% Initialise all points as unvisited
visited = zeros(N, 1);

% Find all core points (points that have minPts in their respective eps
% radii
% Calculate the number of neighbours for each data point, neglecting the
% data point itself
numNeighbors = sum(D <= eps & D ~= 0, 2);
% Logical boolean of core points
isCore = numNeighbors >= minPts;
% Indexes of core points
idx = find(isCore);

%Assign cluster labels to all points
C = zeros(N, 1);
clusterIdx = 1;
for i = idx'
    neighbors = find(D(i, :) <= eps & D(i,:) ~= 0);
    j = 1;
    while j <= length(neighbors) 
        k = neighbors(j);
        if ~visited(k)
           visited(k) = true;
           C(k) = clusterIdx;
           if isCore(k)
              new_neighbors = find(D(k,:) <= eps & D(k,:) ~= 0);
              neighbors = unique([neighbors new_neighbors], 'stable');
           end
        end
        j = j + 1;
   end
   clusterIdx = clusterIdx + 1;
end
end