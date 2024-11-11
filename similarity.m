function J = similarity(X, C)
% S is the Silhouette score of each data point. S can be value [-1,1] with
% 1 being the best and -1 being the worst.
[S,~] = silhouette(X, C);
J = sum(S)/length(X)*100;
end