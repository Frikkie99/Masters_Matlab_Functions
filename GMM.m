function [C, pi, mu, sigma] = GMM(X, K, maxIter)
% Generate sample data
N = length(X);    % number of data points
[~,D] = size(X);   % number of dimensions

% Initialize GMM parameters
[~, mu] = kmeans(X,K,"Start","plus"); % Means
% mu = [8;13];
sigma = repmat(eye(D), [1 1 K]);      % covanriaces matrices
pi = ones(1,K)/K;                     % mixing coefficients, adds to 1

% E-M algorithm
log_likelihood = zeros(maxIter, 1);
for iter = 1:maxIter
    % E-step: compute responsibilities
    r = zeros(N, K);
    for k = 1:K
%        if sum(sum(sigma(:,:,k))) == 0 
%             sigma(:,:,k) = eye(D);
%         end
        r(:,k) = pi(k) * mvnpdf(X, mu(k,:), (sigma(:,:,k) + eye(D)*1e-12)) + ones(N,1)*1E-15;
    end
    
    % Normalise responsibilities
    if sum(r(1,:), 2) == 0
        r(1,1) = 1;
    end

    r = r ./ sum(r, 2); % Sum rows
    
    % M-step: update parameters
    Nk = sum(r, 1); % Sum columns
    mu = (r'*X) ./ Nk';
    for k = 1:K
        Xc = X - mu(k,:);
        sigma(:,:,k) = (Xc' * diag(r(:,k)) * Xc) / Nk(k);
    end

    pi = Nk/ N;
    
    % Compute log-likelihood
    log_likelihood(iter) = sum(log(sum(r, 2)));
    
    % Check for convergence
%     if iter > 1 && abs(log_likelihood(iter) - log_likelihood(iter-1)) < 1e-6
%         break;
%     end
end
[~,C] = max(r,[],2);
end

